# Batch-04 — donor-biased, pos6-centred design under the Asp6 read-out gate

**Date:** 2026-11-03 → 2026-11-04 (Asia/Taipei).
**Scripts:** `scripts/run_batch04.py`, `scripts/rank_batch04.py`,
`scripts/replicate_batch04.py`, `scripts/replicate_rmsd_batch04.py`.
**Raw data:** `checkpoint.jsonl` (authoritative, 54 rows), `candidates.json`,
`backbone_poses.json`, `run_stats.json`, `gate_results.json`, `designs.csv`,
`replicates.json`, `replicate_rmsd.json`, `scored_structures/`,
`replicate_structures/`.
**Boltz-2 calls:** 54 scoring (one constrained G12D call per candidate) + 6
replicate calls = **60**. Track-1 design slots consumed: **0** (nothing
submitted).

## Note on the decision this file replaces

The previous version of this file recorded the decision **not to run
batch-04**. That decision was correct on its own terms: the pocket-constrained
Boltz-2 protocol had just failed the 7OW6 positive control on the half that
mattered for ranking (`specificity_gap` = −0.006 for a known selective TCR,
and ±0.3 swings between byte-identical calls), and batch-04 as then designed
would have ranked ~80 new designs by that discredited quantity. See
[`results/validation/CONTROL_7OW6.md`](../validation/CONTROL_7OW6.md).

**Why it was reversed:** the gate audit
([`results/gate-audit/ASP6_GATE_AUDIT.md`](../gate-audit/ASP6_GATE_AUDIT.md))
then established that **0 of the 88 designs** ever scored with the Asp6
read-out form any hydrogen bond or salt bridge to the Asp6 carboxylate, and
that no existing design reaches the peptide-burial floor. Re-ranking the
existing library therefore cannot produce a gate-passing submission — the
shortfall is in the *designs*, not in the scoring. A gate-passing submission
requires new generation. The user approved the run on that basis.

Crucially, this batch is **not** the cancelled batch-04. It changes the
*designing*, not the scoring: nothing here is ranked on `specificity_gap` or on
the paired MM-GBSA `ddG_sel`, both of which remain retired. ipTM appears only
as a pose-quality read-out, never as a gap and never as a tie-breaker.

## What was run

1. **Hotspots collapsed onto peptide position 6** (`C6`, `C6-C7`, `C5-C7`)
   across four strategies, with the design-stage pose filter tightened on
   position 6 so only backbones that actually reach the mutated side chain are
   paid for.
2. **Donor bias applied locally.** The hosted ProteinMPNN NIM ignores
   `bias_AA_jsonl` and returns HTTP 500 for `pssm_jsonl` / `bias_by_res_jsonl`,
   so API-side biasing was unavailable. Instead: over-sample ProteinMPNN,
   *select* on donor content at the positions pointing at Asp6, and derive one
   deterministic **donor-anchored** graft per backbone (nearest shell position →
   Arg, second → Tyr, Asp/Glu in the pos6 shell → their neutral amides).
3. **Constrained Boltz-2**, strictly serial at ≥ 8 s intervals, once per
   candidate, with the validated flat pocket term
   `{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`.
4. **Gate:** a binder donor (Arg/Lys/His/Ser/Thr/Tyr) heavy atom within 3.5 Å
   of Asp6 OD1/OD2 **and** `asp6_buried_frac > 0`.

The run was interrupted twice by the platform and resumed from
`checkpoint.jsonl` both times; `run_stats_resumed.json` is the record of the
final segment only (17 calls). **54 of 54 candidates are scored, with 0
analysis errors.** Re-running the resume path confirms nothing is outstanding.

## Result 1 — the donor bias worked, chemically

For the first time in this project, designs read the Asp6 carboxylate
directly:

| quantity | batches 02+03 (88 designs) | **batch-04 (54 designs)** |
| --- | --- | --- |
| Asp6 hydrogen bond | **0** | **2** |
| Asp6 salt bridge | **0** | **3** |
| any donor within 3.5 Å of OD1/OD2 | **0** | **2** |
| `asp6_buried_frac > 0` | 5 / 88 | **23 / 54** |
| max `asp6_buried_frac` | 0.538 | **0.888** |
| designs with ≥ 1 Asp6 **side-chain** contact | — | 24 / 54 |

**2 of 54 pass the Asp6 donor gate**, at
`asp6_min_donor_distance` = **2.94 Å** and **2.70 Å** — both genuine
salt-bridge-length Arg NH2···OD2 contacts.

Both survivors are **donor-anchored grafts**. None of the 38 MPNN-selected
sequences passed. The deterministic Arg graft is what produced the result; the
donor-content *selection* over MPNN samples did not.

## Result 2 — the interface charge moved, as intended

| | median net interface charge |
| --- | --- |
| batch-02 | −1.0 |
| batch-03 | −1.0 |
| **batch-04, all 54** | **0.0** (mean −0.5) |
| batch-04, donor-anchored (16) | **0.0** (mean **+0.08**) |
| batch-04, mpnn-selected (38) | −0.45 (mean −0.74) |

31 of 54 designs are at or above neutral, against a −1.0 median in both prior
batches. Distribution: min −4.0, max +3.0, 23 negative / 31 zero-or-positive.
**Yes — the donor bias did what it was built to do**, and the effect is
concentrated in exactly the subset it was applied to.

## Result 3 — the ranked survivors

Ranked on measured geometry only (peptide burial, Asp6 **side-chain**
contacts, charge, topological diversity). Both survivors are in separate
70 %-identity clusters.

| # | design | strategy | len | `peptide_buried_frac` | Asp6 side-chain contacts | `asp6_min_donor_distance` | partner | `asp6_buried_frac` | net charge | peptide area | ipTM* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `B4-c6-62_bb25` (donor-anchored) | B4-c6-62 | 62 | **0.229** | **13** | **2.94 Å** | ARG29:NH2 | **0.853** | **+1.0** | 205 Å² | 0.419 |
| 2 | `B4-c67-65_bb26` (donor-anchored) | B4-c67-65 | 65 | 0.173 | 6 | **2.70 Å** | ARG58:NH2 | 0.643 | 0.0 | 176 Å² | 0.245 |

\* ipTM is a pose-quality read-out only. Neither reaches 0.75; the batch max is
0.708 and the median 0.348.

**Neither survivor meets the `peptide_buried_frac` ≥ 0.40 ranking floor, and
neither does any of the other 52** (batch max 0.353). This floor has not been
loosened. It is worth recording — without acting on it here — that the 7OW6
control measures **0.224 predicted / 0.207 crystal**, i.e. the real,
experimentally validated selective TCR also fails this floor, and
`CONTROL_7OW6.md` proposes 0.18 on that evidence. Rank 1's 0.229 is slightly
*above* the real TCR's value. The recalibration remains unapplied and is the
user's call.

## Result 4 — replicate pose stability: this is what kills them

The one ranking axis the 7OW6 post-mortem left standing is whether a pose
reproduces across identical constrained calls. Three byte-identical re-runs per
survivor, plus binder C-alpha RMSD measured in the pMHC frame (HLA + peptide
superposed, binder **not** re-fitted; frame RMSD 0.06–0.36 Å throughout, so the
numbers below are real).

| | `B4-c6-62_bb25` | `B4-c67-65_bb26` |
| --- | --- | --- |
| **gate passes / draws** | **3 of 4** | **1 of 4** |
| `asp6_min_donor_distance` over draws | 2.94, 2.72, —, 3.25 | 2.70, —, —, — |
| `asp6_buried_frac` | 0.853, 0.876, 0.792, 0.904 (SD 0.048) | 0.643, 0.201, 0.045, 0.000 (SD 0.086) |
| Asp6 side-chain contacts | 13, 15, 10, 15 (SD 2.4) | 6, **0, 0, 0** |
| `peptide_buried_frac` | 0.229, 0.235, 0.198, 0.246 (SD 0.021) | 0.173, 0.222, 0.187, 0.228 |
| net interface charge | +1.0, +1.0, −2.0, +1.0 | 0.0, −2.0, −3.0, −2.0 |
| ipTM | 0.419, 0.402, 0.438, 0.440 (SD 0.018) | 0.245, 0.430, 0.432, 0.320 |
| **binder C-alpha RMSD, pairwise mean / max** | **5.67 Å / 10.38 Å** | **16.92 Å / 18.13 Å** |

Pairwise binder RMSDs:

* `bb25` — 0.42, 1.40, 1.30 Å among three of the four poses, and
  **10.15–10.38 Å** from the fourth. Three draws land in one tight, reproducible
  docked mode that reads the carboxylate (Asp6 burial 0.79–0.90 in *every* draw,
  including the failing one); one draw flips to a completely different mode.
* `bb26` — **15.16–18.13 Å between every single pair.** There is no
  reproducible pose at all. Its gate pass was one sampling draw, and the three
  replicates give **zero** Asp6 side-chain contacts.

This is precisely the mechanism `CONTROL_7OW6.md` predicted: with
`constraints` and `without_potentials: false`, a weak binder's landscape has
several docked modes that the sampling flips between, while a real TCR's funnel
never flips.

**`B4-c67-65_bb26` should be dropped.** Its gate pass does not survive
replication.

## Per-strategy breakdown

| strategy | n scored | median charge | gate survivors | `asp6_buried_frac` > 0 | max `peptide_buried_frac` |
| --- | --- | --- | --- | --- | --- |
| `B4-c6-62` (C6 only, 62 aa) | 19 | 0.0 | **1** | 8 | 0.319 |
| `B4-c67-65` (C6-C7, 65 aa) | 10 | −0.5 | 1 (fails replication) | 5 | 0.269 |
| `B4-pos6-68` (C5-C7, 68 aa) | 12 | −0.5 | 0 | 4 | 0.228 |
| `B4-pos6-72` (C5-C7, 72 aa) | 13 | +1.0 | 0 | 6 | 0.353 |

The single-hotspot `C6` strategy at 62 aa is the one that produced the only
design whose carboxylate read-out reproduces. `B4-pos6-68` — the setting that
produced `bb13` in batch-03 — produced no survivor here.

## Verdict

* The batch did what it was designed to do: **the Asp6 read-out problem is
  solved at the chemistry level.** Two designs form Arg···Asp6 salt bridges at
  2.7–2.9 Å, 23 of 54 bury the aspartate at all (vs 5 of 88 before), and the
  interface charge moved off −1.0 to 0.0. That is real, measured progress and it
  validates the local donor-graft approach over API-side biasing.
* **Exactly one design — `B4-c6-62_bb25` — survives both the gate and
  replication**, and even it flips to an unrelated docked mode on one draw in
  four, and sits at `peptide_buried_frac` 0.229 against a 0.40 floor it cannot
  meet (nor can the real TCR, at 0.224).
* Nothing here reaches the standard of a confident submission.

## Recommendation

**Do not submit at 07:30.** Reasons, in order of weight:

1. The Track-1 cap is 40 designs, **0 are spent**, and the deadline is
   2026-11-08 (AoE) — five days away. There is no scarcity argument for
   spending a slot on a design with a 1-in-4 pose flip.
2. `bb26` is disqualified by replication. That leaves a one-design submission,
   which wastes the ≤ 20-designs-per-submission slot and the 24-hour cadence.
3. `bb25`'s honest status is "the first design in this project to read the
   G12D aspartate reproducibly", not "a good binder". ipTM 0.42 and
   `peptide_buried_frac` 0.229 are a weak complex.

**Do instead, in priority order:**

1. **Re-run the `B4-c6-62` strategy at volume.** It is the only strategy that
   produced a replicating survivor, and the donor-anchored graft is the only
   variant type that ever passes the gate (2/16 vs 0/38). Generate
   donor-anchored variants *for every backbone*, not one per backbone.
2. **Fold replication into the gate itself.** A design should pass the Asp6
   gate in ≥ 3 of 3 draws, not 1 of 1. On this batch that costs 3× the calls on
   survivors only, which is affordable — the whole 54-candidate run was 54
   calls.
3. **Settle the `peptide_buried_frac` floor.** 0.40 is not reachable by the
   real 7OW6 TCR. Either apply the control-calibrated 0.18 or accept that this
   column cannot be a hard gate. This is a user decision and has not been made
   here.
4. Keep `bb25` as the held reference design to beat.

No submission has been made, no PR merged, nothing pushed to `main`.
