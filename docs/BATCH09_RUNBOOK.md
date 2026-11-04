# Batch-09 execution and recovery

The 2026-11-04 model run completed 12 calls and 36 samples using the repository
secret. Results are in `results/batch-09/MATCHED_ANCHOR_RESULTS.md`.
The complete checkpoint, analysis, structures and compressed raw response
archive are preserved in this review branch. Historical batches 06–08 remain
unchanged. Reanalysis of the completed panel needs no model requests.

## Frozen experiment

| Item | Setting |
| --- | --- |
| Candidates | `bb50-g1`, `bb54-g2`, `bb25` — historical sequences |
| Target states | G12D, WT, apo HLA+B2M, unrelated-peptide decoy |
| Shared regional anchor | HLA A69 and A155, maximum distance 10 Å |
| `bb50-g1`, `bb25` | 8 recycling steps, 50 sampling steps, potentials enabled |
| `bb54-g2` | 3 recycling steps, 50 sampling steps, without potentials |
| Per candidate/state | 1 request containing 3 distinct structures |
| Full panel | 12 model requests, 36 structures |

These are the settings implemented and tested in the merged batch-09 code.
The old PR body describing five candidates, Q70/R114 anchors and one common
sampling setting does not describe that final implementation.

## Run on GitHub

The initial push of this workflow to `run/batch09-matched-anchor` starts the
actual 12-call model panel. The push trigger is restricted to that branch and
changes to this workflow file, so committing result files does not run it again.

For subsequent manual preflight runs after merging, open **Actions → Batch 09 matched-anchor
counter-screen → Run workflow**. The default `dry_run=true` checks all inputs
and tests without requiring a key or making a model request.

For the actual panel, use `dry_run=false`. The workflow reads the repository's
`NVIDIA_API_KEY` secret and fails before requesting a model if it is absent.
It never prints the key. No Proteinbase submission is performed.

Download the `batch09-matched-anchor-<run-id>-<attempt>` artifact after the run.
It includes the complete raw model responses, coordinates, checkpoint and run
metadata, even if scoring or a later request failed. The full analysis JSON
is produced only when all 36 samples pass data-integrity checks. These checks
are not biological selectivity thresholds.

**Do not use “Re-run jobs” to resume a paid model run.** Each fresh Actions job
starts a fresh filesystem. Instead, recover the downloaded artifact locally
as below. The model workflow does not commit or submit results. A separate, narrowly
scoped preservation workflow imports the reviewed, completed artifact into the
existing PR branch without making model calls.

## Local preflight and execution

```bash
python -m pip install -r requirements.txt
python scripts/batch09_matched_anchor.py --dry-run
# With NVIDIA_API_KEY already set in the execution environment:
python scripts/batch09_matched_anchor.py --max-boltz-calls 12
# New runs have versioned full readout flags. For the archived panel:
python scripts/audit_carboxylate_readout.py
```

The dry run prints the exact sequences, target polymers, sampling parameters
and regional constraint. It creates no result files and requires no key.
The cap applies to new model requests in this invocation; cached response
recovery does not use that cap. HTTP retries follow the existing NVIDIA client.

## Recovery after interruption

Keep `results/batch-09/responses/`, `structures/` and
`matched_anchor_checkpoint.jsonl` together. For an Actions artifact, restore its
`results/batch-09/` directory into the same repository-relative path in a local
checkout. Run the same command and selection again.

The runner saves each full response before scoring its first draw, then skips
only draws already checkpointed from that exact response. A call is complete
only after all requested draws are present. Sampling parameters, sequence,
constraint and coordinate hashes must match before resuming. A partial legacy
checkpoint with no original response stops with a recovery message; it is not
silently skipped or combined with a fresh model call. A malformed trailing
JSONL record also stops for explicit recovery; preserve the original file.

If model scoring fails, the response and structure remain available, and the
analyzer refuses a final report until the scoring issue is resolved. Empty,
short and duplicate-structure responses are preserved for inspection and are
not counted as a complete panel.

## Evidence to review before committing results

Review `matched_anchor_checkpoint.jsonl`, `matched_anchor_run.json`,
`matched_anchor_analysis.json` and the matching coordinates together. Raw
responses are also archived in `matched_anchor_responses.tar.gz` with a SHA-256
checksum; extract that archive in `results/batch-09/` for response-level recovery. Preserve the artifact or an equivalent
response backup to allow recovery and reanalysis.

The matched comparison remains diagnostic: no apo-area cutoff is introduced,
WT footprint loss remains descriptive, and ipTM is not an affinity or
selectivity measurement. Consult `results/batch-09/SUMMARY.md` for the control
audit and candidate priorities before deciding which sequences to advance.

The original unversioned checkpoint is historical. Current readout and pose
measurements are in `results/readout-correction/`; see its README before using
any archived gate count.
