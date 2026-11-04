# Batch-10 plan — counter-screen first, generation second

**Status:** protocol only. No batch-10 model call has been made by this commit.

## Why this batch exists

The corrected submission-review pool has 11 retained direct-readout sequences.
Ranks 1–3 already have a matched G12D / WT / apo / decoy panel. Rank 4 has one
older WT diagnostic, but ranks 4–11 do **not** have a matched apo panel. That is
a major evidence gap because peptide-free HLA binding is an explicit competition
objective.

Before generating another backbone, batch-10 asks a cheaper and more useful
question:

> Is one of the eight lower-ranked retained sequences already less
> apo-prone than the three headline leads?

The experiment reuses batch-09's HLA regional anchor (A69 + A155, 10 Å) on
**both G12D and apo** so the docking prior is matched across the comparison.

## Phase 1 — screen

Eight frozen sequences: review ranks **4–11**.

For each design:

- one G12D call, 3 diffusion samples;
- one apo call, 3 diffusion samples;
- corrected oxygen-specific Asp6 readout from `selectivity_panel.panel_metrics`;
- no RFdiffusion;
- no ProteinMPNN;
- no gate change;
- no Track-1 design slot spent.

Exactly **16 Boltz-2 calls / 48 structures** if run to completion.

The predeclared advance heuristic is:

1. corrected direct Asp6 readout on G12D in **at least 2 of 3** draws; and
2. median apo total buried area / median G12D total buried area **< 1.0**.

This is deliberately a triage heuristic, **not** a claim that buried area is a
calibrated affinity predictor. It only asks for the correct direction before
spending more calls.

## Phase 2 — confirm

Only phase-1 advancers get:

- WT, 3 diffusion samples;
- decoy, 3 diffusion samples.

Maximum another 16 calls if all eight advance. In practice the phase is
adaptive: zero advancers means zero confirm calls.

## Frozen per-design sampling

The sampling regime comes from the best previously observed G12D condition for
each retained sequence where such a sweep exists. It is then applied identically
to every arm for that sequence.

The code refuses to run if review ranks 4–11 or their sequence hashes change.

## Run

Inspect without an API call:

```bash
python scripts/batch10_counter_screen.py --phase screen --dry-run
```

Run phase 1:

```bash
python scripts/batch10_counter_screen.py --phase screen --max-boltz-calls 16
```

After phase 1 is complete, inspect the fixed advance decision:

```bash
python scripts/batch10_counter_screen.py --phase confirm --dry-run
```

Then, only if there are advancers:

```bash
python scripts/batch10_counter_screen.py --phase confirm --max-boltz-calls 16
```

All responses are cached before scoring, checkpoints are fsync'd, structure
hashes are verified on resume, and duplicate diffusion structures are rejected.
