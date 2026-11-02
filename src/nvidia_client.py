#!/usr/bin/env python3
"""Thin HTTP clients for the NVIDIA BioNeMo NIM endpoints used by this project.

One function per endpoint, all hitting the hosted NIMs on
``https://health.api.nvidia.com``:

====================  ==================================================  ==========================
function              path                                                model
====================  ==================================================  ==========================
rfdiffusion_generate  ``/v1/biology/ipd/rfdiffusion/generate``            RFdiffusion 2.3 (IPD)
proteinmpnn_design    ``/v1/biology/ipd/proteinmpnn/predict``             ProteinMPNN (IPD)
boltz2_predict        ``/v1/biology/mit/boltz2/predict``                  Boltz-2 (MIT)
====================  ==================================================  ==========================

The DeepMind AlphaFold2 NIM (``/v1/biology/deepmind/alphafold2``) is retired and
returns HTTP 410, so Boltz-2 does the complex scoring of pipeline section 5.

Authentication
--------------
The API key is read from the ``NVIDIA_API_KEY`` environment variable **only**
(see ``.env.example``); it is never taken from an argument default, written to
disk, or included in any log line or exception message emitted here.

Asynchronous responses
----------------------
The hosted NIMs answer small requests synchronously (HTTP 200 with
``Nvcf-Status: fulfilled``). Larger ones may answer ``202 Accepted`` with an
``nvcf-reqid`` header; this module then polls
``/v2/nvcf/pexec/status/{req_id}`` until the result is ready, so callers always
get the finished payload back.

Retries
-------
429 and 5xx are retried with exponential backoff plus jitter, honouring
``Retry-After`` when the server sends it. 4xx other than 429 are raised
immediately as :class:`NvidiaApiError` with the server's validation message,
which is what the NIMs use to describe a malformed body.
"""
from __future__ import annotations

import contextlib
import os
import random
import threading
import time
from typing import Any, Callable, Dict, Iterator, List, Mapping, Optional, Sequence

import requests

__all__ = [
    "NvidiaApiError",
    "NvidiaAuthError",
    "NvidiaClient",
    "RateLimiter",
    "Throttle",
    "default_throttle",
    "configure_throttle",
    "rfdiffusion_generate",
    "proteinmpnn_design",
    "boltz2_predict",
    "BASE_URL",
    "RFDIFFUSION_PATH",
    "PROTEINMPNN_PATH",
    "BOLTZ2_PATH",
]

BASE_URL = "https://health.api.nvidia.com"
RFDIFFUSION_PATH = "/v1/biology/ipd/rfdiffusion/generate"
PROTEINMPNN_PATH = "/v1/biology/ipd/proteinmpnn/predict"
BOLTZ2_PATH = "/v1/biology/mit/boltz2/predict"
STATUS_PATH = "/v2/nvcf/pexec/status/{req_id}"

API_KEY_ENV = "NVIDIA_API_KEY"

DEFAULT_TIMEOUT = 300.0          # seconds, per HTTP request
DEFAULT_POLL_TIMEOUT = 1800.0    # seconds, total wall clock for a 202 poll loop
DEFAULT_POLL_INTERVAL = 5.0      # seconds between status polls
DEFAULT_MAX_RETRIES = 4          # retries after the first attempt
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})
DEFAULT_MAX_BACKOFF = 120.0      # seconds, cap on one backoff sleep

# --- client-side rate limiting ---------------------------------------------
# The hosted NIMs answer single calls happily but return HTTP 429 as soon as
# several are in flight: batch-01 logged 235 of them against Boltz-2 from a
# 6-worker pool. Retrying harder makes that worse, so the throttle below is the
# primary defence and the retry loop is only the backstop.
#
# Every request passes a *global* token bucket (all three endpoints share the
# account's quota) and then a per-endpoint bucket. Boltz-2 additionally takes a
# process-wide lock, so Boltz-2 calls are strictly serial no matter how many
# threads the caller runs. All of it is overridable from the environment, so a
# run can be re-tuned without a code change:
#
#   NVIDIA_GLOBAL_RPS         requests/second summed over all endpoints
#   NVIDIA_GLOBAL_BURST       bucket depth for the global bucket
#   NVIDIA_BOLTZ2_MIN_INTERVAL  seconds between the *starts* of Boltz-2 calls
#   NVIDIA_RFDIFFUSION_MIN_INTERVAL / NVIDIA_PROTEINMPNN_MIN_INTERVAL
#   NVIDIA_THROTTLE=0         disable throttling entirely (tests, dry runs)
#
# The Boltz-2 default is the rate measured by `scripts/calibrate_boltz2.py` on
# 2026-11-02 and recorded in docs/PIPELINE.md:
#   6 workers, no spacing  -> 31x 429 in 12 calls, 6 of 12 calls lost
#   serial, no spacing     -> 12x 429 in 20 calls, ~9.8 calls/min
#   serial, 8 s spacing    ->  0x 429 in 12 calls, 7.9 calls/min  <- the default
# i.e. the hosted Boltz-2 NIM sustains ~8 requests/minute from one key, and
# punishes concurrency far harder than volume. See that section of
# docs/PIPELINE.md before changing this number.
DEFAULT_GLOBAL_RPS = 1.0
DEFAULT_GLOBAL_BURST = 2.0
DEFAULT_BOLTZ2_MIN_INTERVAL = 8.0
DEFAULT_RFDIFFUSION_MIN_INTERVAL = 2.0
DEFAULT_PROTEINMPNN_MIN_INTERVAL = 2.0

# Endpoints that may never overlap with themselves, whatever the caller does.
SERIAL_PATHS = frozenset({BOLTZ2_PATH})


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        return float(raw)
    except ValueError:
        return default


class RateLimiter:
    """Token bucket with an optional minimum spacing between grants.

    Thread-safe and *deterministic*: grant times are tracked on an internal
    virtual clock that only ever moves forward, so the limiter behaves the same
    whether `sleep` really sleeps (production) or is a no-op (tests).

    `rate` is tokens per second, `burst` the bucket depth (how many requests may
    go out back to back after an idle period) and `min_interval` a floor on the
    gap between two consecutive grants.
    """

    def __init__(
        self,
        rate: float,
        burst: float = 1.0,
        min_interval: float = 0.0,
        monotonic: Callable[[], float] = time.monotonic,
        name: str = "",
    ) -> None:
        if rate <= 0:
            raise ValueError("rate must be > 0 requests/second")
        self.rate = float(rate)
        self.burst = max(1.0, float(burst))
        self.min_interval = max(0.0, float(min_interval))
        self.name = name
        self._monotonic = monotonic
        self._lock = threading.Lock()
        self._tokens = self.burst
        now = monotonic()
        self._updated = now
        self._virtual_now = now
        self._next_earliest = now
        self.grants = 0
        self.waited_seconds = 0.0

    def reserve(self, tokens: float = 1.0) -> float:
        """Book the next slot and return how long the caller must wait for it."""
        with self._lock:
            now = max(self._monotonic(), self._virtual_now)
            self._tokens = min(
                self.burst, self._tokens + (now - self._updated) * self.rate
            )
            self._updated = now
            shortfall = tokens - self._tokens
            token_time = now + (shortfall / self.rate if shortfall > 0 else 0.0)
            grant = max(token_time, self._next_earliest)
            self._tokens = min(
                self.burst, self._tokens + (grant - self._updated) * self.rate
            ) - tokens
            self._updated = grant
            self._virtual_now = grant
            self._next_earliest = grant + self.min_interval
            self.grants += 1
            wait = grant - self._monotonic()
            if wait > 0:
                self.waited_seconds += wait
            return max(0.0, wait)

    def acquire(self, sleep: Callable[[float], Any] = time.sleep,
                tokens: float = 1.0) -> float:
        """Wait (via `sleep`) until a slot is free; returns the seconds waited."""
        wait = self.reserve(tokens)
        if wait > 0:
            sleep(wait)
        return wait

    def stats(self) -> Dict[str, float]:
        return {
            "name": self.name,
            "rate_per_second": self.rate,
            "burst": self.burst,
            "min_interval": self.min_interval,
            "grants": self.grants,
            "waited_seconds": round(self.waited_seconds, 2),
        }


class Throttle:
    """Global + per-endpoint limiters, plus serial locks for some endpoints."""

    def __init__(
        self,
        global_limiter: Optional[RateLimiter] = None,
        endpoint_limiters: Optional[Mapping[str, RateLimiter]] = None,
        serial_paths: Sequence[str] = tuple(SERIAL_PATHS),
    ) -> None:
        self.global_limiter = global_limiter
        self.endpoint_limiters: Dict[str, RateLimiter] = dict(endpoint_limiters or {})
        self.serial_paths = frozenset(serial_paths)
        self._locks: Dict[str, threading.Lock] = {
            path: threading.Lock() for path in self.serial_paths
        }

    @contextlib.contextmanager
    def hold(self, path: str) -> Iterator[None]:
        """Serialise `path` against itself for the whole call, polling included."""
        lock = self._locks.get(path)
        if lock is None:
            yield
            return
        with lock:
            yield

    def acquire(self, path: str, sleep: Callable[[float], Any] = time.sleep) -> float:
        """Block until one request to `path` may be issued."""
        waited = 0.0
        if self.global_limiter is not None:
            waited += self.global_limiter.acquire(sleep)
        limiter = self.endpoint_limiters.get(path)
        if limiter is not None:
            waited += limiter.acquire(sleep)
        return waited

    def stats(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        if self.global_limiter is not None:
            out["global"] = self.global_limiter.stats()
        out["endpoints"] = {
            path: limiter.stats()
            for path, limiter in self.endpoint_limiters.items()
        }
        return out


def build_throttle() -> Throttle:
    """A `Throttle` configured from the environment (see the constants above)."""
    return Throttle(
        global_limiter=RateLimiter(
            rate=_env_float("NVIDIA_GLOBAL_RPS", DEFAULT_GLOBAL_RPS),
            burst=_env_float("NVIDIA_GLOBAL_BURST", DEFAULT_GLOBAL_BURST),
            name="global",
        ),
        endpoint_limiters={
            BOLTZ2_PATH: _interval_limiter(
                "NVIDIA_BOLTZ2_MIN_INTERVAL", DEFAULT_BOLTZ2_MIN_INTERVAL, "boltz2"
            ),
            RFDIFFUSION_PATH: _interval_limiter(
                "NVIDIA_RFDIFFUSION_MIN_INTERVAL",
                DEFAULT_RFDIFFUSION_MIN_INTERVAL,
                "rfdiffusion",
            ),
            PROTEINMPNN_PATH: _interval_limiter(
                "NVIDIA_PROTEINMPNN_MIN_INTERVAL",
                DEFAULT_PROTEINMPNN_MIN_INTERVAL,
                "proteinmpnn",
            ),
        },
    )


def _interval_limiter(env_name: str, default: float, name: str) -> RateLimiter:
    interval = max(1e-6, _env_float(env_name, default))
    return RateLimiter(rate=1.0 / interval, burst=1.0,
                       min_interval=interval, name=name)


_THROTTLE_LOCK = threading.Lock()
_THROTTLE: Optional[Throttle] = None


def default_throttle() -> Optional[Throttle]:
    """The process-wide throttle, built on first use; None when disabled."""
    global _THROTTLE
    if os.environ.get("NVIDIA_THROTTLE", "").strip() in ("0", "off", "false"):
        return None
    with _THROTTLE_LOCK:
        if _THROTTLE is None:
            _THROTTLE = build_throttle()
        return _THROTTLE


def configure_throttle(throttle: Optional[Throttle]) -> Optional[Throttle]:
    """Replace the process-wide throttle (used by the batch runner and tests)."""
    global _THROTTLE
    with _THROTTLE_LOCK:
        _THROTTLE = throttle
    return _THROTTLE


class NvidiaApiError(RuntimeError):
    """An NVIDIA NIM returned an error, or the result never arrived.

    `status_code` is the HTTP status (None for a client-side timeout) and
    `payload` the decoded body when there was one.
    """

    def __init__(
        self,
        message: str,
        status_code: Optional[int] = None,
        payload: Any = None,
        request_id: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.payload = payload
        self.request_id = request_id


class NvidiaAuthError(NvidiaApiError):
    """No API key in the environment, or the server rejected it (401/403)."""


def _api_key() -> str:
    key = os.environ.get(API_KEY_ENV, "").strip()
    if not key:
        raise NvidiaAuthError(
            f"{API_KEY_ENV} is not set in the environment. Copy .env.example to "
            f".env and put your build.nvidia.com key there, or export it in the "
            f"shell. The key value is never read from anywhere else."
        )
    return key


def _error_message(payload: Any, fallback: str) -> str:
    """Pull the human-readable message out of a NIM error body."""
    if isinstance(payload, Mapping):
        for field in ("error", "detail", "message", "title"):
            value = payload.get(field)
            if isinstance(value, str) and value:
                return value
            if value is not None:
                return str(value)
    if isinstance(payload, str) and payload.strip():
        return payload.strip()[:2000]
    return fallback


class NvidiaClient:
    """Session-backed client for the three BioNeMo NIMs.

    Parameters are all about transport, never about biology: see
    :meth:`rfdiffusion_generate`, :meth:`proteinmpnn_design` and
    :meth:`boltz2_predict` for the per-endpoint request schemas.
    """

    def __init__(
        self,
        base_url: str = BASE_URL,
        timeout: float = DEFAULT_TIMEOUT,
        max_retries: int = DEFAULT_MAX_RETRIES,
        poll_interval: float = DEFAULT_POLL_INTERVAL,
        poll_timeout: float = DEFAULT_POLL_TIMEOUT,
        session: Optional[requests.Session] = None,
        sleep=time.sleep,
        throttle: Any = "default",
        max_backoff: float = DEFAULT_MAX_BACKOFF,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.max_retries = max_retries
        self.poll_interval = poll_interval
        self.poll_timeout = poll_timeout
        self.session = session or requests.Session()
        self._sleep = sleep
        self.max_backoff = max_backoff
        # "default" means the process-wide throttle; None disables throttling
        # for this client only (both are useful in tests).
        self.throttle = default_throttle() if throttle == "default" else throttle
        # Validate the key eagerly: better to fail before a long pipeline run.
        self._key = _api_key()

    # -- transport ----------------------------------------------------------

    def _headers(self) -> Dict[str, str]:
        return {
            "Authorization": f"Bearer {self._key}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        }

    def _backoff(self, attempt: int, response: Optional[requests.Response]) -> float:
        """Seconds to wait before retry `attempt` (1-based).

        `Retry-After` wins when the server sends one (seconds or an HTTP date).
        Otherwise it is exponential backoff with *decorrelated jitter*: the base
        doubles per attempt and the actual wait is a random point in
        [base, base * 1.5], which spreads a thundering herd of retries instead
        of re-synchronising them. Capped at `max_backoff`.
        """
        if response is not None:
            retry_after = self._retry_after_seconds(response)
            if retry_after is not None:
                return min(self.max_backoff, retry_after)
        base = min(self.max_backoff, float(2 ** (attempt - 1)))
        return min(self.max_backoff, base * (1.0 + random.random() * 0.5))

    @staticmethod
    def _retry_after_seconds(response: requests.Response) -> Optional[float]:
        """`Retry-After` as seconds, accepting both the delta and date forms."""
        raw = response.headers.get("Retry-After") or response.headers.get(
            "retry-after"
        )
        if not raw:
            return None
        raw = str(raw).strip()
        try:
            return max(0.0, float(raw))
        except ValueError:
            pass
        try:
            from email.utils import parsedate_to_datetime
            import datetime as _dt

            when = parsedate_to_datetime(raw)
            if when is None:
                return None
            if when.tzinfo is None:
                when = when.replace(tzinfo=_dt.timezone.utc)
            delta = (when - _dt.datetime.now(_dt.timezone.utc)).total_seconds()
            return max(0.0, delta)
        except Exception:  # noqa: BLE001 - a malformed header is just "no hint"
            return None

    def _decode(self, response: requests.Response) -> Any:
        try:
            return response.json()
        except ValueError:
            return response.text

    def post(self, path: str, body: Mapping[str, Any]) -> Dict[str, Any]:
        """POST `body` to `path`, following 202/nvcf-reqid polling and retries.

        Rate limited: endpoints in `SERIAL_PATHS` (Boltz-2) are held under a
        process-wide lock for the whole call, so two of them never overlap, and
        every attempt - the first and each retry - must first be granted a slot
        by the global and per-endpoint token buckets.
        """
        if self.throttle is not None:
            with self.throttle.hold(path):
                return self._post_attempts(path, body)
        return self._post_attempts(path, body)

    def _post_attempts(self, path: str, body: Mapping[str, Any]) -> Dict[str, Any]:
        url = f"{self.base_url}{path}"
        last_error: Optional[NvidiaApiError] = None

        for attempt in range(1, self.max_retries + 2):
            retry_source: Optional[requests.Response] = None
            if self.throttle is not None:
                self.throttle.acquire(path, self._sleep)
            try:
                response = self.session.post(
                    url, headers=self._headers(), json=dict(body), timeout=self.timeout
                )
            except requests.Timeout as exc:
                last_error = NvidiaApiError(f"POST {path} timed out: {exc}")
            except requests.RequestException as exc:
                last_error = NvidiaApiError(f"POST {path} failed: {exc}")
            else:
                status = response.status_code
                req_id = response.headers.get("nvcf-reqid") or response.headers.get(
                    "Nvcf-Reqid"
                )
                if status == 200:
                    return self._as_mapping(self._decode(response), path, req_id)
                if status == 202:
                    if not req_id:
                        raise NvidiaApiError(
                            f"POST {path} returned 202 with no nvcf-reqid header, "
                            f"so the result cannot be polled for",
                            status_code=202,
                        )
                    return self._poll(req_id, path)
                if status in (401, 403):
                    raise NvidiaAuthError(
                        f"POST {path} rejected the API key (HTTP {status}): "
                        f"{_error_message(self._decode(response), 'unauthorized')}. "
                        f"Check {API_KEY_ENV} and that the key is entitled to "
                        f"this NIM.",
                        status_code=status,
                        request_id=req_id,
                    )
                payload = self._decode(response)
                error = NvidiaApiError(
                    f"POST {path} returned HTTP {status}: "
                    f"{_error_message(payload, f'HTTP {status}')}",
                    status_code=status,
                    payload=payload,
                    request_id=req_id,
                )
                if status not in RETRY_STATUS:
                    raise error
                last_error = error
                retry_source = response

            if attempt <= self.max_retries:
                self._sleep(self._backoff(attempt, retry_source))

        assert last_error is not None
        raise last_error

    def _poll(self, req_id: str, path: str) -> Dict[str, Any]:
        """Poll an accepted (202) request until it is fulfilled."""
        url = f"{self.base_url}{STATUS_PATH.format(req_id=req_id)}"
        deadline = time.monotonic() + self.poll_timeout
        while True:
            response = self.session.get(
                url, headers=self._headers(), timeout=self.timeout
            )
            if response.status_code == 200:
                return self._as_mapping(self._decode(response), path, req_id)
            if response.status_code != 202:
                payload = self._decode(response)
                raise NvidiaApiError(
                    f"polling {path} (request {req_id}) returned HTTP "
                    f"{response.status_code}: "
                    f"{_error_message(payload, 'polling failed')}",
                    status_code=response.status_code,
                    payload=payload,
                    request_id=req_id,
                )
            if time.monotonic() >= deadline:
                raise NvidiaApiError(
                    f"{path} (request {req_id}) was still running after "
                    f"{self.poll_timeout:.0f}s",
                    status_code=202,
                    request_id=req_id,
                )
            self._sleep(self.poll_interval)

    @staticmethod
    def _as_mapping(payload: Any, path: str, req_id: Optional[str]) -> Dict[str, Any]:
        if not isinstance(payload, Mapping):
            raise NvidiaApiError(
                f"{path} returned a {type(payload).__name__} body, expected a "
                f"JSON object",
                payload=payload,
                request_id=req_id,
            )
        return dict(payload)

    # -- endpoints ----------------------------------------------------------

    def rfdiffusion_generate(
        self,
        input_pdb: str,
        contigs: str,
        hotspot_res: Optional[Sequence[str]] = None,
        diffusion_steps: int = 15,
        random_seed: Optional[int] = None,
    ) -> str:
        """Generate a backbone with RFdiffusion; returns the output PDB text.

        Confirmed request schema (``/v1/biology/ipd/rfdiffusion/generate``):

        * ``contigs`` (string, **required**) - e.g. ``"A1-180/0 C1-10/0 70-70"``:
          keep chain A residues 1-180, chain break, keep chain C residues 1-10,
          chain break, then build a new 70-residue chain (the binder).
        * ``input_pdb`` (string) - PDB text. Either this or ``input_pdb_asset``
          is required; an empty body answers
          ``422 {"error":"Validation Errors: contigs: Field required"}`` and a
          body with only ``contigs`` answers ``422 {"error":"Either input_pdb or
          input_pdb_asset must be specified"}``.
        * ``hotspot_res`` (array of strings) - ``["C4","C5","C6"]`` style
          chain+residue labels the binder must contact.
        * ``diffusion_steps`` (int, 1-50, default 50), ``random_seed`` (int).

        Response: ``{"output_pdb": "<PDB text>", "elapsed_ms": <int>}``. In the
        output, the **designed binder is chain A** and the kept target segments
        are concatenated into chain B in contig order.
        """
        body: Dict[str, Any] = {
            "input_pdb": input_pdb,
            "contigs": contigs,
            "diffusion_steps": int(diffusion_steps),
        }
        if hotspot_res:
            body["hotspot_res"] = list(hotspot_res)
        if random_seed is not None:
            body["random_seed"] = int(random_seed)

        payload = self.post(RFDIFFUSION_PATH, body)
        output_pdb = payload.get("output_pdb")
        if not isinstance(output_pdb, str) or "ATOM" not in output_pdb:
            raise NvidiaApiError(
                "RFdiffusion response had no usable 'output_pdb'",
                payload={k: v for k, v in payload.items() if k != "output_pdb"},
            )
        return output_pdb

    def proteinmpnn_design(
        self,
        input_pdb: str,
        input_pdb_chains: Optional[Sequence[str]] = None,
        num_seq_per_target: int = 8,
        sampling_temp: Sequence[float] = (0.1,),
        ca_only: bool = False,
        use_soluble_model: bool = True,
        random_seed: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Design sequences for `input_pdb` with ProteinMPNN.

        Confirmed request schema (``/v1/biology/ipd/proteinmpnn/predict``):

        * ``input_pdb`` (string) - PDB text; required in practice (an empty body
          answers ``422 {"error":"No input PDB provided"}``, an unparseable one
          ``422 {"error":"Found no valid chains in input PDB"}``). May instead be
          a filename alongside ``input_pdb_asset``.
        * ``input_pdb_chains`` (array) - chains to **design**; every other chain
          is held fixed. Pass ``["A"]`` on RFdiffusion output to design the
          binder and freeze the target.
        * ``ca_only`` (bool, default false), ``use_soluble_model`` (bool,
          default false), ``random_seed`` (int),
          ``num_seq_per_target`` (int 1-100, default 1),
          ``sampling_temp`` (array of float 0-1, 0.1-0.3 recommended),
          plus the PSSM options (``pssm_jsonl``, ``pssm_multi``,
          ``pssm_threshold``).

        Response: ``{"mfasta": "<multi-FASTA>", "scores": [...], "probs": [...]}``.
        The first FASTA record is the *input* sequence (poly-G for an
        RFdiffusion backbone), not a design; this method strips it.

        Returns ``{"mfasta": str, "sequences": [{"header","sequence","score",
        "sample"}...], "scores": [...]}``.
        """
        body: Dict[str, Any] = {
            "input_pdb": input_pdb,
            "ca_only": bool(ca_only),
            "use_soluble_model": bool(use_soluble_model),
            "num_seq_per_target": int(num_seq_per_target),
            "sampling_temp": [float(t) for t in sampling_temp],
        }
        if input_pdb_chains:
            body["input_pdb_chains"] = list(input_pdb_chains)
        if random_seed is not None:
            body["random_seed"] = int(random_seed)

        payload = self.post(PROTEINMPNN_PATH, body)
        mfasta = payload.get("mfasta")
        if not isinstance(mfasta, str) or ">" not in mfasta:
            raise NvidiaApiError(
                "ProteinMPNN response had no usable 'mfasta'",
                payload={k: v for k, v in payload.items() if k != "probs"},
            )
        return {
            "mfasta": mfasta,
            "sequences": parse_mpnn_mfasta(mfasta),
            "scores": payload.get("scores", []),
        }

    def boltz2_predict(
        self,
        polymers: Sequence[Mapping[str, Any]],
        recycling_steps: int = 3,
        sampling_steps: int = 50,
        diffusion_samples: int = 1,
        step_scale: float = 1.638,
        without_potentials: bool = False,
    ) -> Dict[str, Any]:
        """Predict a complex with Boltz-2.

        Confirmed request schema (``/v1/biology/mit/boltz2/predict``):

        * ``polymers`` (array, **required**, 1-12 items). An empty body answers
          ``422 {"error":"body -> polymers: Field required"}`` and ``[]`` answers
          ``422 {"error":"body -> polymers: List should have at least 1 item
          after validation, not 0"}``. Each item:
          ``molecule_type`` (**required**, ``"protein"``/``"dna"``/``"rna"``),
          ``sequence`` (**required**, 1-4096 chars), ``id`` (optional, 1-4
          alphanumerics), ``cyclic`` (optional bool), ``msa`` (optional).
        * ``recycling_steps`` (1-10, default 3), ``sampling_steps``
          (10-1000, default 50), ``diffusion_samples`` (1-25, default 1),
          ``step_scale`` (0.5-5.0, default 1.638),
          ``without_potentials`` (bool), ``output_format`` (``"mmcif"`` only).

        Response keys: ``structures`` (list of
        ``{"structure","format","name","source"}``, mmCIF text),
        ``confidence_scores``, ``iptm_scores``, ``ptm_scores``,
        ``complex_plddt_scores``, ``complex_iplddt_scores``,
        ``complex_pde_scores``, ``complex_ipde_scores``,
        ``chains_ptm_scores``, ``pair_chains_iptm_scores``, ``metrics``.
        ``pae``/``pde`` come back **null** on the hosted endpoint even with
        ``write_full_pae``: full matrices are only written to disk by a
        self-hosted NIM, so per-residue PAE is not available here.
        """
        if not polymers:
            raise ValueError("boltz2_predict needs at least one polymer")
        if len(polymers) > 12:
            raise ValueError(
                f"Boltz-2 accepts at most 12 polymers, got {len(polymers)}"
            )
        cleaned: List[Dict[str, Any]] = []
        for polymer in polymers:
            sequence = str(polymer.get("sequence", "")).strip().upper()
            if not sequence:
                raise ValueError("every polymer needs a non-empty 'sequence'")
            if len(sequence) > 4096:
                raise ValueError(
                    f"polymer {polymer.get('id', '?')} is {len(sequence)} residues; "
                    f"Boltz-2 accepts at most 4096"
                )
            item: Dict[str, Any] = {
                "molecule_type": str(polymer.get("molecule_type", "protein")),
                "sequence": sequence,
            }
            if polymer.get("id"):
                item["id"] = str(polymer["id"])
            if polymer.get("cyclic"):
                item["cyclic"] = True
            if polymer.get("msa"):
                item["msa"] = polymer["msa"]
            cleaned.append(item)

        body = {
            "polymers": cleaned,
            "recycling_steps": int(recycling_steps),
            "sampling_steps": int(sampling_steps),
            "diffusion_samples": int(diffusion_samples),
            "step_scale": float(step_scale),
            "without_potentials": bool(without_potentials),
        }
        payload = self.post(BOLTZ2_PATH, body)
        if not payload.get("structures"):
            raise NvidiaApiError(
                "Boltz-2 response contained no structures",
                payload={k: v for k, v in payload.items() if k != "structures"},
            )
        return payload


def parse_mpnn_mfasta(mfasta: str) -> List[Dict[str, Any]]:
    """Parse a ProteinMPNN multi-FASTA into design records.

    Drops the leading ``>input, ...`` record, which echoes the input sequence
    rather than a design. Each returned record has ``header``, ``sequence``,
    ``score`` (float or None) and ``sample`` (int or None).
    """
    records: List[Dict[str, Any]] = []
    header: Optional[str] = None
    chunks: List[str] = []

    def flush() -> None:
        if header is None:
            return
        sequence = "".join(chunks).strip()
        if not sequence:
            return
        if header.lstrip().startswith("input"):
            return  # the echoed input backbone sequence, not a design
        fields: Dict[str, str] = {}
        for part in header.split(","):
            if "=" in part:
                key, _, value = part.partition("=")
                fields[key.strip()] = value.strip()
        def _num(key, cast):
            try:
                return cast(fields[key])
            except (KeyError, ValueError, TypeError):
                return None
        records.append(
            {
                "header": header.strip(),
                "sequence": sequence,
                "score": _num("score", float),
                "sample": _num("sample", int),
            }
        )

    for line in mfasta.splitlines():
        if line.startswith(">"):
            flush()
            header = line[1:]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    flush()
    return records


# -- module-level conveniences (one call, one client) ------------------------

def rfdiffusion_generate(input_pdb: str, contigs: str, **kwargs) -> str:
    """One-shot :meth:`NvidiaClient.rfdiffusion_generate`."""
    client_kwargs = {k: kwargs.pop(k) for k in ("timeout", "max_retries") if k in kwargs}
    return NvidiaClient(**client_kwargs).rfdiffusion_generate(
        input_pdb, contigs, **kwargs
    )


def proteinmpnn_design(input_pdb: str, **kwargs) -> Dict[str, Any]:
    """One-shot :meth:`NvidiaClient.proteinmpnn_design`."""
    client_kwargs = {k: kwargs.pop(k) for k in ("timeout", "max_retries") if k in kwargs}
    return NvidiaClient(**client_kwargs).proteinmpnn_design(input_pdb, **kwargs)


def boltz2_predict(polymers: Sequence[Mapping[str, Any]], **kwargs) -> Dict[str, Any]:
    """One-shot :meth:`NvidiaClient.boltz2_predict`."""
    client_kwargs = {k: kwargs.pop(k) for k in ("timeout", "max_retries") if k in kwargs}
    return NvidiaClient(**client_kwargs).boltz2_predict(polymers, **kwargs)
