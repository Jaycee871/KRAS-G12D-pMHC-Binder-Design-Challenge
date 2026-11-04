#!/usr/bin/env python3
"""Measure a sustainable Boltz-2 request rate against the hosted NIM.

Problem 1 of the batch-01 post-mortem: a 6-worker pool drew 235 HTTP 429s from
`/v1/biology/mit/boltz2/predict` and only finished 50 of 56 candidates. Single
calls had always been fine, so the limit bites on *concurrency*, not on volume.

This script issues strictly serial Boltz-2 calls at a configurable minimum
interval and reports, per interval: wall time per call, the realised request
rate, and how many 429s the client had to back off on. Run it with a couple of
intervals and take the fastest one that logs zero 429s; that number is the
`NVIDIA_BOLTZ2_MIN_INTERVAL` default in `src/nvidia_client.py` and the figure
recorded in `docs/PIPELINE.md`.

The payload defaults to the real four-chain G12D complex plus a binder, so the
numbers describe the calls the batch actually makes rather than a toy job.

    python3 scripts/calibrate_boltz2.py --calls 6 --intervals 0,6
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import nvidia_client  # noqa: E402
import run_pipeline  # noqa: E402
from nvidia_client import NvidiaApiError, NvidiaClient, RateLimiter, Throttle  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]

# A short, cheap stand-in binder: calibration is about request *rate*, but the
# default payload keeps the real target so latency is representative.
PROBE_BINDER = "ADFKSDITKAVLASAPNPTAAFASIAEQEGDVATVELRSGIRYTVRRGVWTGACEAEGPCVQRGE"


class CountingClient(NvidiaClient):
    """Counts every status the retry loop had to back off on."""

    def __init__(self, *args, **kwargs):
        self.retry_statuses: dict = {}
        super().__init__(*args, **kwargs)

    def _backoff(self, attempt, response=None):
        if response is not None:
            code = response.status_code
            self.retry_statuses[code] = self.retry_statuses.get(code, 0) + 1
            print(f"    [retry] HTTP {code} on attempt {attempt}", file=sys.stderr)
        return super()._backoff(attempt, response)


def probe_polymers(target_dir: Path, state: str, binder: str):
    states = run_pipeline.load_target_states(target_dir)
    return list(states[state]) + [
        {"id": run_pipeline.BINDER_CHAIN, "molecule_type": "protein",
         "sequence": binder}
    ]


def measure(interval: float, calls: int, polymers, boltz_kwargs, timeout: float,
            workers: int = 1):
    """`calls` Boltz-2 calls spaced at least `interval` seconds apart.

    `workers` > 1 deliberately removes the serial lock, which is how the
    batch-01 429 storm is reproduced for comparison.
    """
    nvidia_client.configure_throttle(
        Throttle(
            global_limiter=RateLimiter(rate=10.0, burst=1.0, name="global"),
            endpoint_limiters={
                nvidia_client.BOLTZ2_PATH: RateLimiter(
                    rate=1.0 / max(interval, 1e-6), burst=1.0,
                    min_interval=interval, name="boltz2",
                )
            },
            serial_paths=(nvidia_client.BOLTZ2_PATH,) if workers == 1 else (),
        )
    )
    client = CountingClient(timeout=timeout, max_retries=4)
    latencies, failures = [], []

    def one(index):
        call_started = time.monotonic()
        try:
            client.boltz2_predict(polymers, **boltz_kwargs)
        except NvidiaApiError as exc:
            failures.append(str(exc)[:200])
            print(f"  call {index + 1}/{calls} FAILED: {str(exc)[:160]}",
                  file=sys.stderr)
            return
        latency = time.monotonic() - call_started
        latencies.append(latency)
        print(f"  call {index + 1}/{calls}: {latency:.1f}s", flush=True)

    started = time.monotonic()
    if workers == 1:
        for index in range(calls):
            one(index)
    else:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(one, range(calls)))
    elapsed = time.monotonic() - started
    ok = len(latencies)
    return {
        "min_interval_s": interval,
        "workers": workers,
        "calls": calls,
        "succeeded": ok,
        "failed": len(failures),
        "elapsed_s": round(elapsed, 1),
        "mean_latency_s": round(statistics.fmean(latencies), 1) if ok else None,
        "median_latency_s": round(statistics.median(latencies), 1) if ok else None,
        "max_latency_s": round(max(latencies), 1) if ok else None,
        "realised_rate_per_min": round(60.0 * ok / elapsed, 2) if elapsed else None,
        "retry_statuses": dict(client.retry_statuses),
        "http_429": client.retry_statuses.get(429, 0),
        "failures": failures,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--calls", type=int, default=6,
                        help="serial calls per interval")
    parser.add_argument("--intervals", default="0,6",
                        help="comma-separated minimum intervals in seconds")
    parser.add_argument("--workers", type=int, default=1,
                        help="concurrent callers (1 = the serial default; >1 "
                             "reproduces the batch-01 429 storm)")
    parser.add_argument("--state", default="g12d", choices=run_pipeline.STATES)
    parser.add_argument("--binder", default=PROBE_BINDER)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--sampling-steps", type=int, default=50)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--target-dir", default=str(run_pipeline.DEFAULT_TARGET_DIR))
    parser.add_argument("--out", default=str(REPO_ROOT / "results" / "calibration"
                                             / "boltz2_rate.json"))
    args = parser.parse_args(argv)

    polymers = probe_polymers(Path(args.target_dir), args.state, args.binder)
    boltz_kwargs = {
        "recycling_steps": args.recycling_steps,
        "sampling_steps": args.sampling_steps,
        "diffusion_samples": 1,
        "without_potentials": True,
    }
    results = []
    for raw in args.intervals.split(","):
        interval = float(raw.strip())
        print(f"[calibrate] {args.calls} serial Boltz-2 calls, "
              f"min interval {interval:g}s, {args.workers} worker(s)")
        results.append(measure(interval, args.calls, polymers, boltz_kwargs,
                               args.timeout, workers=args.workers))
        print(f"[calibrate] {json.dumps(results[-1], indent=2)}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "measured_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "state": args.state,
        "boltz_kwargs": boltz_kwargs,
        "serial": args.workers == 1,
        "workers": args.workers,
        "runs": results,
    }, indent=2))
    print(f"[calibrate] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
