# KRAS G12D pMHC Binder Design — Challenge 5

De novo binder design for the **KRAS G12D peptide–MHC complex**
(HLA-A\*11:01, PDB [7OW6](https://www.rcsb.org/structure/7OW6)) —
**Anthropic × Adaptyv Protein Design Competition, Challenge 5**, hosted on
Proteinbase.

[![tests](https://github.com/Jaycee871/KRAS-G12D-pMHC-Binder-Design-Challenge/actions/workflows/tests.yml/badge.svg)](https://github.com/Jaycee871/KRAS-G12D-pMHC-Binder-Design-Challenge/actions/workflows/tests.yml)

> **Status: scaffold only. No binder has been designed, folded, scored or
> submitted.** What works today is target preparation, submission building and
> submission validation (all CPU-only, all tested). The design and
> structure-prediction stages are written out concretely in
> [`docs/PIPELINE.md`](docs/PIPELINE.md) but have not been run.

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
| Format | CSV or FASTA, **≤ 5 MB** |
| Columns | `name`, `sequence` required; `molecule_class` optional |
| Designs per submission | **≤ 20** |
| Rate limit | **one submission per 24 hours** |
| Chain length | **10–250 aa** per chain |
| Molecule formats | single chain, nanobody, scFv, Fab — scFv/Fab written as `{VH}:{VL}` |
| Originality | designs must be **de novo**; novelty score **≥ 3/4**, checked automatically |
| Per-participant caps | up to **40** designs Track 1; up to **20** each for Tracks 2 and 3 |

Submit at
<https://proteinbase.com/competitions/anthropic-adaptyv-2026/rounds/5/submit>.
The day-by-day batching plan is in
[`docs/SUBMISSION_PLAN.md`](docs/SUBMISSION_PLAN.md).

## Pipeline overview

Implemented here (CPU, tested) in **bold**; everything else is planned.

1. **Target preparation** — `src/fetch_target.py` downloads 7OW6 from RCSB,
   identifies the heavy chain / B2M / peptide, and writes the three target
   states (G12D complex, WT complex, peptide-free HLA+B2M).
2. Backbone generation — RFdiffusion and BindCraft against a peptide-centred
   epitope (peptide residues 4–8, weighted on residue 6).
3. Sequence design — ProteinMPNN, 8–16 sequences per backbone.
4. Monomer refolding filter — AlphaFold2, keep pLDDT ≥ 85 and RMSD < 2 Å to
   the design backbone.
5. Complex prediction — AlphaFold-Multimer or Boltz, run **three times per
   design** (G12D, WT, peptide-free), scoring ipTM and PAE-interface.
6. Ranking — by `specificity_gap`, behind hard gates on ipTM, interface
   coverage of peptide position 6, and PAE.
7. **Submission** — `src/make_submission.py` ranks and truncates,
   `src/validate_submission.py` enforces every documented rule.

Full detail, including thresholds and filters: [`docs/PIPELINE.md`](docs/PIPELINE.md).

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
├── tests/                       pytest suite (57 tests)
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
| Test suite | ✅ 57 passed |
| CI (pytest on push/PR) | ✅ configured |
| Pipeline documentation | ✅ written (`docs/PIPELINE.md`) |
| Submission plan | ✅ written (`docs/SUBMISSION_PLAN.md`) |
| Backbone generation (RFdiffusion / BindCraft) | ⬜ not started — needs GPU |
| Sequence design (ProteinMPNN) | ⬜ not started |
| AF2 monomer filter | ⬜ not started |
| Complex prediction / three-target scoring | ⬜ not started |
| **Designs produced** | **⬜ none** |
| **Designs submitted** | **⬜ none** |

## Licence

MIT — see [`LICENSE`](LICENSE). Author: **Pack Kwan Low (Jaycee)**.
Third-party attribution in [`NOTICE`](NOTICE). Not affiliated with Anthropic,
Adaptyv Bio, Proteinbase or the RCSB PDB.
