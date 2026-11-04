# KRAS G12D pMHC Binder Design — Challenge 5

De novo binder design for the **KRAS G12D peptide–MHC complex**
(HLA-A\*11:01, PDB [7OW6](https://www.rcsb.org/structure/7OW6)) —
**Anthropic × Adaptyv Protein Design Competition, Challenge 5**, hosted on
Proteinbase.

[![tests](https://github.com/Jaycee871/KRAS-G12D-pMHC-Binder-Design-Challenge/actions/workflows/tests.yml/badge.svg)](https://github.com/Jaycee871/KRAS-G12D-pMHC-Binder-Design-Challenge/actions/workflows/tests.yml)

> **Status: design batches 01–08 complete; batch-09 control audit complete,
> matched-anchor model panel pending. No design has been submitted.**
> The control audit retired WT geometry-loss as a hard selectivity gate.
> `bb50-g1`, `bb54-g2` and `bb25` are the fixed subjects for the next matched
> G12D/WT/apo/decoy diagnostic. See `results/batch-09/SUMMARY.md` and
> [`docs/BATCH09_RUNBOOK.md`](docs/BATCH09_RUNBOOK.md). No batch-09 model
> outputs are currently committed.

---

## The target

| | |
| --- | --- |
| Peptide | **VVVGADGVGK** — KRAS G12D 10-mer |
| Presented on | **HLA-A\*11:01** |
| Structure | **PDB 7OW6** |
| Full complex | HLA-A\*11:01 heavy chain + β2-microglobulin (B2M) + the 10-residue peptide |

Chains extracted from 7OW6 by `src/fetch_target.py` (verified by an actual run):
heavy chain `A` (276 aa), B2M `B` (99 aa), peptide `C` (10 aa, `VVVGADGVGK`).
The TCR chains present in the crystal structure are stripped — the design
target is the pMHC surface.

## The three objectives

A design must:

1. **bind** the G12D complex — HLA + B2M + `VVVGADGVGK`
2. **not bind** the wild-type complex — HLA + B2M + `VVVGAGGVGK`
3. **not bind** peptide-free HLA + B2M

The WT and G12D peptides differ **only at peptide position 6** (`D` vs `G`).
That single side chain is the entire specificity problem: a binder that grips
the α1/α2 groove rim will bind all three targets and score nothing. The
paratope has to sit over position 6. See
[`docs/PIPELINE.md`](docs/PIPELINE.md#0-the-central-difficulty).

```
G12D   V V V G A D G V G K
WT     V V V G A G G V G K
                 ^
            position 6
```

## Ranking

Designs are ranked by, in order:

1. **mutation specificity** — binds G12D, not WT
2. **peptide-dependent binding** — binds the complex, not peptide-free HLA
3. **lowest K_D** — affinity is the last tie-break, not the goal

The repository's internal ranking mirrors this with a single
`specificity_gap = score(G12D) − max(score(WT), score(apo))`, which is the
column `src/make_submission.py` sorts by.

## Submission rules

| rule | value |
| --- | --- |
| Deadline | **Sun 8 Nov 2026, 23:59 Anywhere on Earth** (= 9 Nov 11:59 UTC = 9 Nov 19:59 Asia/Taipei) |
| Format | **CSV**, ordered by preferred rank |
| Columns | `name`, `sequence`, `molecule_class` required |
| Chain length | **10–250 aa** per chain |
| `molecule_class` | `protein`, `nanobody`, `scfv`, `fab_kappa`, or `fab_lambda`; Fab uses `{VH}:{VL}` |
| Originality | designs must be **unique, de novo and zero-shot**; known binders may calibrate filters but cannot be the starting binder |
| Per-participant caps | Track 1: submit **20–40**; Tracks 2/3: up to **20** |
| Track 1 screening | the **top 20** ranked designs that pass filtering are screened |

Submit from the live Challenge 5 page:
<https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/pmhc>.
The current deadline plan is in
[`docs/SUBMISSION_PLAN.md`](docs/SUBMISSION_PLAN.md).

## Pipeline overview

The project is now an **active, measured design campaign**, not a planned
pipeline.

1. **Target preparation** — G12D, WT, apo HLA+B2M, and an unrelated
   HLA-A*11:01 decoy peptide are represented under `data/targets/`.
2. **Backbone generation** — RFdiffusion has been run across multiple design
   hypotheses through batch-08.
3. **Sequence design** — ProteinMPNN plus targeted donor-biased variants has
   been run across the surviving backbones.
4. **Complex prediction** — hosted Boltz-2 is the current workhorse, with
   repeated constrained G12D/WT and diagnostic apo/decoy arms.
5. **Interface analysis** — `src/complex_analysis.py` and
   `src/groove_geometry.py` measure chain-aware buried area, position-6
   contacts, direct Asp6 read-out, interface charge, pose spread, and groove
   partitioning.
6. **Decision gates** — later batches use reproducibility and measured
   geometry/chemistry. The old independent-docking ipTM specificity gap and
   paired MM-GBSA ddG are **retired as decision metrics**; historical code is
   retained for provenance.
7. **Submission tooling** — `src/make_submission.py` and
   `src/validate_submission.py` build/validate the current Proteinbase CSV
   schema.

Latest decision audit: `results/batch-09/SUMMARY.md`; latest model results:
`results/batch-08/SUMMARY.md`. Batch-09 matched-anchor results are pending.
No sequence has been submitted, and selectivity has not been established.

## Repository layout

```
.
├── README.md
├── LICENSE                      MIT — Pack Kwan Low (Jaycee)
├── NOTICE                       third-party data/software attribution
├── requirements.txt             CPU deps; GPU tools installed separately
├── pytest.ini
├── .env.example                 key names only — real keys are never committed
├── .gitignore                   Python + structures, weights, outputs, .env
├── .github/workflows/tests.yml  pytest on push and PR
├── src/
│   ├── fetch_target.py          7OW6 → three target states
│   ├── validate_submission.py   enforce every submission rule, exit non-zero
│   └── make_submission.py       scored table → compliant top-N CSV
├── tests/                       pytest suite
├── docs/
│   ├── PIPELINE.md              the design pipeline, concretely
│   └── SUBMISSION_PLAN.md       day-by-day plan to Nov 8
├── examples/                    format-only examples, not designs
└── data/targets/                generated FASTAs + TARGETS.md (PDBs git-ignored)
```

## How to run

```bash
git clone https://github.com/Jaycee871/KRAS-G12D-pMHC-Binder-Design-Challenge.git
cd KRAS-G12D-pMHC-Binder-Design-Challenge
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

**Build the three targets** (downloads 7OW6 from RCSB; needs network):

```bash
python src/fetch_target.py
# writes data/targets/{g12d_complex,wt_complex_template,apo_hla_b2m}.pdb
#        data/targets/{g12d_complex,wt_complex,apo_hla_b2m}.fasta
#        data/targets/TARGETS.md
```

The small FASTA files are committed; the PDB coordinates are git-ignored and
regenerated by this command.

**Validate a submission file** (exits non-zero with per-row messages):

```bash
python src/validate_submission.py submissions/batch1.csv
```

**Build a submission from a scored designs table:**

```bash
python src/make_submission.py scored_designs.csv -o submissions/batch1.csv -n 20
# ranks by specificity_gap (then score), truncates to N, and refuses to write
# anything that would fail validation
```

**Run the tests:**

```bash
python -m pytest -q
```

## Secrets

`.env` is git-ignored. `.env.example` lists only the *names* of the variables
used (`NVIDIA_API_KEY`, `HF_TOKEN`, `NATURE_API_KEY`). Real keys live in the
local environment or in GitHub Actions secrets and are **never** committed,
printed or logged.

## Status

| component | state |
| --- | --- |
| Target preparation (`fetch_target.py`) | ✅ working — run against RCSB, 3 targets produced, peptide `VVVGADGVGK` extracted |
| Submission validator (`validate_submission.py`) | ✅ working, tested |
| Submission builder (`make_submission.py`) | ✅ working, tested |
| Test suite | ✅ automated CI; run `python -m pytest -q` for current count |
| CI (pytest on push/PR) | ✅ configured |
| Pipeline documentation | ✅ written (`docs/PIPELINE.md`) |
| Submission plan | ✅ written (`docs/SUBMISSION_PLAN.md`) |
| Backbone generation (RFdiffusion) | ✅ run through batch-08 |
| Sequence design (ProteinMPNN + donor-bias variants) | ✅ run through batch-08 |
| Complex prediction / counter-screening (hosted Boltz-2) | ✅ run across G12D / WT / apo / decoy diagnostics |
| Latest model result / audit | ⚠️ batch-08 model results; batch-09 control audit retires the WT geometry-loss hard gate |
| Matched-anchor counter-screen | ⏳ 12-call / 36-sample panel prepared; model results pending |
| **Designs produced** | **✅ multiple batches generated and scored** |
| **Designs submitted** | **⬜ none** |

## Licence

MIT — see [`LICENSE`](LICENSE). Author: **Pack Kwan Low (Jaycee)**.
Third-party attribution in [`NOTICE`](NOTICE). Not affiliated with Anthropic,
Adaptyv Bio, Proteinbase or the RCSB PDB.
