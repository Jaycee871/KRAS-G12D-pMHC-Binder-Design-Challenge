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

import os
import random
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence

import requests

__all__ = [
    "NvidiaApiError",
    "NvidiaAuthError",
    "NvidiaClient",
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
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.max_retries = max_retries
        self.poll_interval = poll_interval
        self.poll_timeout = poll_timeout
        self.session = session or requests.Session()
        self._sleep = sleep
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
        """Seconds to wait before retry `attempt` (1-based)."""
        if response is not None:
            retry_after = response.headers.get("Retry-After")
            if retry_after:
                try:
                    return max(0.0, float(retry_after))
                except ValueError:
                    pass
        return min(60.0, (2 ** (attempt - 1))) * (1.0 + random.random() * 0.25)

    def _decode(self, response: requests.Response) -> Any:
        try:
            return response.json()
        except ValueError:
            return response.text

    def post(self, path: str, body: Mapping[str, Any]) -> Dict[str, Any]:
        """POST `body` to `path`, following 202/nvcf-reqid polling and retries."""
        url = f"{self.base_url}{path}"
        last_error: Optional[NvidiaApiError] = None

        for attempt in range(1, self.max_retries + 2):
            retry_source: Optional[requests.Response] = None
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
