# Reproducibility inventory of `results/`

A factual, batch-by-batch record of what each result directory tested, what it
concluded, which artefacts are committed, what is missing, and whether someone
with only this repository could rerun it. Written 2026-11-04 (Asia/Taipei) on
branch `docs/methodology-record`; no existing result file was modified.

**How "rerunnable" is used below.** Three levels:

* **Re-analysable** — the committed artefacts are enough to recompute the stated
  numbers without any model call.
* **Re-runnable with API budget** — the script and its settings are committed, so
  the experiment can be repeated against the hosted NIMs, but the original
  structures cannot be recovered exactly (the endpoints are stochastic and
  unversioned here).
* **Not reproducible as run** — something needed is absent and no command in the
  repo regenerates it.

**Repository-wide facts that affect every row.**

* `*.pdb` / `*.cif` are gitignored except `results/*/structures/*`. The bulk dumps
  `all_structures/`, `scored_structures/`, `replicate_structures/` and
  `results/batch-01/structures/` are ignored by policy; they exist on this machine
  (about 230 MB in total) but are **not in git**.
* No hosted-model version, digest or function id is recorded in any artefact. Only
  request bodies are kept.
* `requirements.txt` has lower-bound pins only; OpenMM (needed for the MM-GBSA
  work) is not listed at all.
* `logs/` is gitignored by policy; checkpoints are the authoritative record.
* Only batch-09 preserves raw API responses.

---

## batch-01 — first real design batch

* **Tested:** four strategies (`S1-core-short` 65 aa, `S2-core-long` 95,
  `S3-pos6-tight` 80, `S4-full-groove` 85) scored against all four states, under
  the original hard-gate set including the ipTM `specificity_gap`.
* **Concluded:** 0 of 56 scored candidates cleared every gate; 11 of 56 buried any
  peptide surface, 11 touched nothing; mean buried area 568 Å² total vs 21 Å² on
  the peptide; `peptide_buried_frac ≥ 0.40` passed by 0 of 56. The 7OW6 control
  row (groove 0.90 Asp6 burial, TCR 0.00 with 7 position-6 contacts) already
  showed the `asp6_buried_frac ≥ 0.60` gate was unreachable by construction.
* **Committed:** `SUMMARY.md`, `checkpoint.jsonl`, `designs.csv`,
  `designs.details.json`, `candidates.json`, `run_stats.json` (6 tracked files).
* **Missing / ignored:** `structures/` is explicitly gitignored (it holds
  overlapping top-5 dumps of two re-runs, not a curated set);
  `all_structures/` (18 MB) and `scored_structures/` (20 MB) ignored. No raw API
  responses. No structure hashes.
* **Sampling settings recorded?** Yes — `run_stats.json` carries `boltz_kwargs`,
  `diffusion_steps`, `strategies`, throttle and retry statistics, and
  `src/run_batch.py --seed` defaults to 1, though the seed actually used is not
  written into `run_stats.json`.
* **Rerunnable?** Re-analysable for all tabulated metrics (checkpoint is complete).
  Not reproducible as run at the structure level: no coordinates are in git and
  the per-run seed is not recorded.

## batch-02 — pose-filtered, Asp6 read-out gated

* **Tested:** five S1 variants (58/65/72 aa core, 65/68 aa pos6) over 120
  backbones with the two new pose filters; ipTM demoted from hard gate to ranking
  target; `asp6_buried_frac ≥ 0.60` replaced by the Asp6 H-bond/salt-bridge
  read-out.
* **Concluded:** 0 of 60 cleared all hard gates; 11 of 60 survived the predicted
  pose filter; 0 of 60 formed an Asp6 hydrogen bond or salt bridge; 0 of 60 reached
  `peptide_buried_frac` 0.40 (best 0.27).
* **Committed:** `SUMMARY.md`, `checkpoint.jsonl`, `designs.csv`,
  `designs.details.json`, `candidates.json`, `backbone_poses.json`,
  `run_stats.json`, and 5 curated mmCIF in `structures/` (12 tracked files).
* **Missing:** `all_structures/` (8.1 MB) and `scored_structures/` (21 MB) ignored;
  no raw responses; no structure hashes; 5 structures retained of 93 predictions.
* **Sampling settings recorded?** Yes, in `run_stats.json` (`boltz_kwargs`,
  `diffusion_steps`, `pose_filter`, throttle). Seed not recorded in the artefact.
* **Rerunnable?** Re-analysable; re-runnable with API budget. Structure-level
  reproduction of the 88 discarded predictions is not possible from the repo.

## batch-03 — pos6-centred, and the source of the rescore cohort

* **Tested:** `S1-pos6-65` / `S1-pos6-68` over 120 backbones.
* **Concluded:** 0 of 28 cleared all hard gates; 8 of 28 survived the predicted
  pose filter; 0 of 28 formed an Asp6 bond; best `peptide_buried_frac` 0.31.
* **Committed:** same file set as batch-02 (12 tracked files, 5 curated mmCIF).
* **Missing:** bulk dumps ignored (8.0 MB + 9.5 MB); no raw responses; no hashes.
* **Sampling settings recorded?** Yes, as batch-02.
* **Rerunnable?** Re-analysable; re-runnable with API budget.

## calibration — Boltz-2 throttle study

* **Tested:** sustainable hosted-Boltz-2 request rate: serial 5, serial 20,
  8 s-interval, and 6 parallel workers.
* **Concluded:** ~8 requests/minute is sustainable from one key; 6 workers gave
  31 × HTTP 429 and 6 lost calls out of 12; serial-20 with no interval realised
  9.76 calls/min with 12 × 429. This is why every later batch ran strictly serial
  at ≥ 8 s.
* **Committed:** four JSON files with timestamps, `boltz_kwargs`, latencies,
  retry statuses and failures. `scripts/calibrate_boltz2.py` is committed.
* **Missing:** nothing material; no structures were retained (none needed).
* **Rerunnable?** Re-runnable with API budget; the measurement is of a hosted
  service at one point in time (2026-11-02) and will not reproduce exactly.

## validation — the 7OW6 positive control

* **Tested:** does the pocket-constrained Boltz-2 protocol (a) place a known
  binder where it really sits and (b) give the ipTM specificity gap a clearly
  positive value for a known-selective TCR?
* **Concluded:** geometry **passes** (peptide buried fraction 0.224 predicted vs
  0.207 crystal; 825 vs 813 Å² total; variable-domain Cα RMSD 1.29 / 1.51 Å in a
  0.85 Å pMHC frame); discrimination **fails** (gap = −0.006; design-regime gaps
  swing ±0.3 between byte-identical calls). The control also shows the real TCR
  failing two of six hard gates. The gate recalibration proposed here (0.40 → 0.18,
  −2.0 → −3.5) was recorded as evidence and **not applied**.
* **Committed:** `CONTROL_7OW6.md`, `7ow6_constrained_control.json`,
  `7ow6_constrained_checkpoint.jsonl`, `7ow6_replicates.jsonl`,
  `7ow6_replicates_summary.json`, `7ow6_control_rmsd.json`,
  `7ow6_buried_area.json` (7 tracked). Scripts `control_7ow6_constrained.py`,
  `control_7ow6_replicates.py`, `control_7ow6_rmsd.py`, `validate_7ow6.py`.
* **Missing:** **the five control mmCIFs are untracked** —
  `7ow6_control_{g12d,wt,decoy}_constrained.cif`,
  `7ow6_control_{g12d,apo}_unconstrained.cif` sit in `results/validation/` on disk
  but are caught by the global `*.cif` rule, which whitelists only
  `results/*/structures/`. Also missing: the crystal `data/raw/7OW6.pdb`
  (gitignored, re-fetchable), any raw API response, and any hash of either.
* **Sampling settings recorded?** Yes — constraint form and per-call settings are
  in the checkpoint and described in the markdown.
* **Rerunnable?** Numbers are re-analysable from JSON, but the **RMSD check cannot
  be re-run** from a clone: `control_7ow6_rmsd.py` needs both the predicted mmCIFs
  and the crystal. Fix is cheap (~3 MB): move the five files under
  `results/validation/structures/`.

## ddg-scoring — paired same-pose MM-GBSA, and its retirement

* **Tested:** `ddG_sel = dG_bind(G12D pose) − dG_bind(same pose, pos6→Gly)` on the
  7OW6 crystal and prediction, an off-site p3 negative, a decoy negative, an
  off-peptide design, plus a groove control and a 2×2 parameter sweep.
* **Concluded:** wrong sign at all four corners (+1.27, +2.32, +7.13, +9.60
  kcal/mol); replicate SD 0.8–8.7 kcal/mol; the decoy gives the strongest
  wrong-direction value (−3.53 ± 0.12); the groove's apolar Val3→Gly (−7.32)
  matches the real Arg114–Asp6 salt bridge (−7.42). Diagnosis: the metric measures
  side-chain burial rather than recognition, and 7OW6 is ill-posed for it (no TCR
  atom within 5 Å of the Asp6 side chain). Retired without a fix;
  `scripts/ddg_rescore_designs.py` was never run on a design.
* **Committed:** `FINDINGS.md`, `control_7ow6_ddg.json`,
  `control_diagnostics.jsonl`, `groove_control.jsonl`; `src/ddg_scoring.py` and
  `tests/test_ddg_scoring.py` (10 tests) in CI.
* **Missing:** OpenMM is not in `requirements.txt` and no OpenMM / force-field
  version is recorded in the outputs; the input `data/raw/7OW6.pdb` is gitignored
  with no hash; no minimised structures were kept.
* **Sampling settings recorded?** The physics knobs (eps_in, L-BFGS iterations,
  jitter σ, 16 Å truncation, restraint constant) are recorded in `FINDINGS.md` and
  in the JSON; the software stack is not.
* **Rerunnable?** CPU-only and in principle re-runnable with no API budget, but
  **not reproducible as run** until the environment is pinned and the input
  structure is hash-identified — the findings themselves show the result is
  sensitive to numerical detail.

## gate-audit — Asp6 read-out audit of the existing library

* **Tested:** whether anything already scored (batches 01–03) qualifies under the
  only surviving filter.
* **Concluded:** **0 of the 88 designs scored with the Asp6 read-out (batches 02
  and 03) form any hydrogen bond or salt bridge to Asp6**;
  `asp6_min_donor_distance` is unpopulated in every row; batch-01 was never scored
  for the read-out; best `peptide_buried_frac` anywhere is 0.361. Re-ranking the
  library cannot produce a gate-passing submission.
* **Committed:** `ASP6_GATE_AUDIT.md` only (1 file).
* **Missing:** **no script**. The audit has no committed generator, so its counts
  can be recomputed by hand from `results/batch-0{1,2,3}/designs.csv` but not by
  command. No intermediate table was kept.
* **Rerunnable?** Re-analysable by hand from committed CSVs; not reproducible as a
  command.

## constrained-rescore — the pocket-constraint discovery

* **Tested:** rescoring all 28 batch-03 sequences under the newly discovered flat
  pocket constraint; also the schema itself (nested forms and `[chain, index]`
  pairs return HTTP 422).
* **Concluded:** the binder lands on the peptide in 26 of 28 under the constraint
  (vs 8 of 28), Asp6 bonds 3 of 28 (vs 0), but `peptide_buried_frac ≥ 0.40` still
  0 of 28 (max 0.32); and the top gaps did not reproduce (`bb18` +0.298 → −0.300),
  which is why batch-04 was not launched at the time.
* **Committed:** `FINDINGS.md`, `checkpoint.jsonl`, `replicate_top3.jsonl`,
  `g12d_repeat_spread.json` (4 tracked files). Scripts
  `constrained_rescore.py`, `replicate_top3_gaps.py`, `g12d_repeat_spread.py`.
* **Missing:** **28 per-design constrained mmCIFs are present on disk but
  untracked** (same `*.cif` rule as `validation/`). `FINDINGS.md` also points at a
  driver path `/workspace/constrained_rescore.py`, i.e. an out-of-repo location;
  the in-repo equivalent is `scripts/constrained_rescore.py`. No raw responses, no
  hashes.
* **Sampling settings recorded?** The constraint is recorded; per-call
  `boltz_kwargs` are in the checkpoint rows.
* **Rerunnable?** Re-analysable; re-runnable with API budget. Coordinates not
  recoverable from a clone.

## constraint_test — the four-backbone constraint probe

* **Tested:** before/after effect of the pocket constraint on four batch-03
  backbones.
* **Concluded:** three backbones that were off-peptide landed on it
  (e.g. `bb51` 0.0 → 0.187 peptide buried fraction, 0 → 9 position-6 contacts);
  `bb13` gained an Asp6 H-bond and ipTM 0.192 → 0.801.
* **Committed:** `report.json` only (1 tracked file).
* **Missing:** the four `*_pocket.cif` files on disk are untracked; no script is
  named in the directory and none is obviously dedicated to it; no raw responses.
* **Rerunnable?** The JSON is self-contained as a record; the experiment is not
  reproducible as a command from the repo alone.

## batch-04 — donor-biased design, Asp6 read-out gate

* **Tested:** hotspots collapsed onto position 6; local donor bias (the hosted
  ProteinMPNN ignores `bias_AA_jsonl` and 500s on `pssm_jsonl`), one deterministic
  donor-anchored graft per backbone; 54 constrained G12D calls + 6 replicates.
* **Concluded:** the chemistry worked (2 Asp6 H-bonds, 3 salt bridges, 23 of 54
  with `asp6_buried_frac > 0` vs 5 of 88 before; median interface charge moved
  −1.0 → 0.0), but only `B4-c6-62_bb25` survived gate **and** replication (3 of 4
  draws, pairwise binder RMSD 5.67 Å mean / 10.38 Å max); `bb26` had no
  reproducible pose (15.2–18.1 Å between every pair) and was dropped. 0 of 54
  reached `peptide_buried_frac` 0.40.
* **Committed:** `SUMMARY.md`, `checkpoint.jsonl` (54 rows, authoritative),
  `candidates.json`, `backbone_poses.json`, `designs.csv`, `gate_results.json`,
  `replicates.json`, `replicate_rmsd.json`, `run_stats.json`,
  `run_stats_resumed.json`, 8 curated mmCIF (18 tracked files).
* **Missing / defective:** `scored_structures/` (ignored, 11 MB) holds **36 files
  for 54 rows** — the name `{backbone}_{variant}` collides for the two
  MPNN-selected candidates per backbone and the second overwrites the first;
  recorded in `run_stats.json:structure_file_caveat`. Checkpoint metrics are
  unaffected; both gate survivors are donor-anchored and intact.
  `replicate_structures/` ignored. No raw responses, no structure hashes.
  `run_stats_resumed.json` covers only the final 17 of 54 scoring calls (the run
  was interrupted twice).
* **Sampling settings recorded?** Yes at batch level (`run_stats.json`:
  constraint, `boltz_kwargs`, gate definition, pose filter, strategies). **Not**
  per row: `checkpoint.jsonl` rows carry only `row`, `scored_at`,
  `structure_file` — no `boltz_kwargs`, no `structure_sha1`, no condition label.
  `scripts/run_batch04.py --seed` defaults to 4042 but the value used is not
  written into the artefacts.
* **Rerunnable?** Re-analysable for all reported metrics; re-runnable with API
  budget. Per-draw settings must be inferred from the batch-level record.

## batch-05 — `B5-c6-62` at volume, replication inside the gate

* **Tested:** batch-04's only replicating strategy at volume (200 backbones, 60
  grafts), with a 3-of-3 replication requirement folded into the gate.
* **Concluded:** 4 of 60 cleared the Asp6 gate on the screening draw; **0 of 4
  cleared it 3 of 3**. Chemistry scaled (39 of 60 with `asp6_buried_frac > 0`, max
  0.966, two 2.03 Å donor contacts); charge held at median 0.0 / mean +0.39. Every
  gate failure coincided with a flip between docked modes
  (`bb50-g1` 14.26 Å mean / 20.37 Å max over all four poses). 0 of 60 reached 0.40;
  42 of 60 clear the unapplied 0.18.
* **Committed:** `SUMMARY.md`, `checkpoint.jsonl` (60 rows), `candidates.json`,
  `backbone_poses.json`, `designs.csv`, `gate_results.json`, `replicates.json`,
  `replicate_rmsd.json`, `final_ranking.json`, `survivors.csv`, `run_stats.json`,
  16 curated mmCIF (27 tracked files).
* **Missing:** `all_structures/` (14 MB), `scored_structures/` (19 MB),
  `replicate_structures/` (3.7 MB) ignored; no raw responses; no structure hashes.
* **Sampling settings recorded?** Yes at batch level — `run_stats.json` records the
  constraint, `boltz_kwargs` (recycling 3 / sampling 50 / `diffusion_samples` 1 /
  `without_potentials` false), the pose filter and the gate definition verbatim.
  Again not per checkpoint row. Seed default 5050, not written to artefacts.
* **Rerunnable?** Re-analysable; re-runnable with API budget.

## sampling-study — are the pose flips a property of the designs or of the sampler?

* **Tested:** five existing designs × four sampling cells (baseline, `nopot`,
  `recyc8`, `combo`), 3 draws each, 39 new Boltz-2 calls, no new designs.
* **Concluded:** the spread collapses and **the best cell is design-specific**:
  `bb50-g1` 14.44 / 20.37 Å → **0.70 / 1.05 Å** at `recycling_steps` 8; `bb25`
  7.27 → 1.96 Å; `bb135-g1` 4.54 → 2.04 Å under combo; `bb173-g1` improves in no
  cell. Three designs reached the Asp6 gate 3 of 3 where batch-05 had none. The
  study's own caveats are recorded: n = 3 per cell, G12D only, one axis per cell,
  and "a collapsed spread is self-consistency, not accuracy".
* **Committed:** `SUMMARY.md`, `checkpoint.jsonl` (39 rows), `analysis.json`,
  `structures/` (39 mmCIF, all tracked) — 42 tracked files. This directory keeps a
  structure for every call.
* **Missing:** no raw responses; **no `structure_sha1` and no sequence string in
  the checkpoint rows**, so the byte-identical-response artefact that the later
  selectivity panel documents cannot be audited retrospectively here; the baseline
  cell was not re-run but re-used from batch-04/05 replicate draws, so the baseline
  rows live in a different directory.
* **Sampling settings recorded?** Yes, per row (`boltz_kwargs`, `condition`,
  `constraint`) — the first directory in the project to do this.
* **Rerunnable?** Fully re-analysable including RMSDs, because every structure is
  committed. Re-runnable with API budget.

## selectivity-panel — replication at n = 10 and the first off-target panel

* **Tested:** does the collapse hold at n = 10, and does the Asp6 read-out survive
  WT / apo / decoy arms? 57 Boltz-2 calls.
* **Concluded:** `bb50-g1` 10/10 (RMSD 1.58 Å mean / 3.78 Å max), `bb135-g1` 10/10
  (2.17 / 4.79), `bb25` 9/10 with one ~14 Å mode flip. **No design is selective:**
  all three bury more on peptide-free HLA (1101–1454 Å²) than on G12D
  (627–1076 Å²); `bb25` buries 196 Å² of WT peptide vs 190 Å² of G12D;
  `bb135-g1` salt-bridges the decoy's Asp8 in 3 of 4 draws and scores higher ipTM
  on WT than G12D in every draw. `bb50-g1`'s decoy arm is the single clean
  rejection (238 → 62 Å²).
* **Committed:** `SUMMARY.md`, `checkpoint.jsonl` (57 rows), `analysis.json`,
  `structures/` (57 mmCIF, all tracked) — 60 tracked files.
* **Missing:** no raw responses; **no `structure_sha1`** in the rows, which matters
  here because the SUMMARY itself reports that repeat calls returned byte-identical
  structures (6–8 distinct of 10 on G12D; all 4 identical on apo for two designs) —
  the distinct-count claim cannot be re-verified from the checkpoint alone, only by
  re-hashing the committed mmCIFs.
* **Sampling settings recorded?** Yes, per row (`boltz_kwargs`, `condition`,
  `constrained`, `constraint`, `arm`, `sequence`).
* **Rerunnable?** Fully re-analysable (structures committed); the distinct-structure
  audit is recoverable by hashing the tracked files. Re-runnable with API budget.

## batch-06 — shrink the binder, WT arm inside the gate

* **Tested:** 34/40/46 aa binders, a peptide-centred design filter, the ratio
  `buried_area_peptide / buried_area_hla` as the optimised quantity, and
  `WT_CRITERION` (≥ 30 % loss on both median peptide area and median position-6
  contacts) placed inside the gate.
* **Concluded:** null. Peptide share 0.209 → **0.213** while every component shrank
  ~21 %; 28 of 264 backbones survived the filter; 0 of 41 designs cleared 3-of-3.
  The WT arm was run as an explicit diagnostic on the two 2-of-3 near-misses:
  `bb54-g2` lost **59.2 %** of peptide area and **73.1 %** of position-6 contacts
  on WT (the first WT discrimination in the project) and rejected the decoy
  (−87 %), but its apo burial (1116–1220 Å²) still exceeded on-target (532–793 Å²).
* **Committed:** `SUMMARY.md`, `checkpoint.jsonl` (41 rows),
  `panel_checkpoint.jsonl` (30), `analysis.json`, `panel_results.json`,
  `wt_diagnostic.json`, `candidates.json` + `candidates_w2.json`,
  `backbone_poses.json` + `_w2`, `run_stats.json` + `_w2`, 30 mmCIF (42 tracked).
* **Missing:** `all_structures/` (17 MB) and `scored_structures/` (12 MB) ignored;
  no raw responses.
* **Sampling settings recorded?** Yes, and best-in-class for the design batches:
  per-row `boltz_kwargs`, `condition`, `constraint`, `structure_sha1`,
  `distinct_structures_in_call`, plus batch-level `run_stats.json` with the filter
  and gate definitions. Seed default 6060 is **not** written into `run_stats.json`
  (batch-08 does record its seed; batch-06 and -07 do not).
* **Rerunnable?** Re-analysable, including the distinct-structure claim (hashes are
  in the checkpoint). Re-runnable with API budget.

## batch-07 — sampling sweep on bb54-g2 and the rigidified graft

* **Tested:** 17 sampling cells on `bb54-g2`; a rigidifying redesign protecting the
  six Asp6-shell positions; 21 new candidates in paired graft / rigid twins;
  72 Boltz-2 calls, 216 structures.
* **Concluded:** `bb54-g2` reaches 3-of-3 reproducibly only at `r3np`
  (recycling 3 + `without_potentials`), pooled 8 of 9 draws — i.e. the earlier
  2-of-3 plateau was a sampling artefact. **Pose tightness and gate reproducibility
  are anti-correlated** (tightest cell `r10np` 2.25/2.81 Å gates 2/3; the 3-of-3
  cell sits at 13.93/16.97 Å). Rigidification failed: 1 of 21 hit 3-of-3 once and
  gave 1/3 on repeat; self-rigidifying the parent degraded it; Track B 0 of 6. No
  staple was formable, so no `rigid-gps` variant exists. Pooled WT margin
  −33.4 % / −52.4 % with no strict separation; apo unchanged (1102–1239 Å² vs
  379–758 Å²). 0 of 183 draws reach `peptide_buried_frac` 0.40.
* **Committed:** `SUMMARY.md`, `panel_checkpoint.jsonl`, `analysis.json`,
  `screen_summary.jsonl`, `screen_stats.json`, `cross_condition.json`,
  `depth_results.json`, `pooled_wt.json`, `candidates*.json`, `backbone_poses*.json`
  (including the preserved wave-1 `*_w1_ds30.json` artefacts), 18 curated mmCIF
  (34 tracked files).
* **Missing:** 198 of 216 structures are in the ignored `all_structures/` (65 MB),
  while `panel_checkpoint.jsonl` still names them under `structures/` — the paths
  in the checkpoint do not resolve in a clone (the SUMMARY says so explicitly). No
  raw responses. No seed in `screen_stats.json`.
* **Sampling settings recorded?** Yes, per row, including `structure_sha1` and
  `distinct_structures_in_call`; the duplicate-call pair (`r3np_rep`/`r3np_rep3`)
  is documented and excluded from pooled statistics.
* **Rerunnable?** Re-analysable for every tabulated metric; **pose spreads for the
  198 untracked structures are not recomputable in a clone**. Re-runnable with API
  budget.

## batch-08 — shape steering, de-wetting, and the pooled WT decay

* **Tested:** four shape-steered strategies (192 RFdiffusion calls) against five
  pre-fixed shape thresholds; helix-apolar penalty in MPNN selection; 13 paired
  de-wetted twins; a controlled de-wetting of `bb54-g2` itself; a pooled WT
  recomputation; an apo/decoy probe on the two most peptide-centred designs.
* **Concluded:** null on all counts. **0 of 192 backbones cleared all five shape
  thresholds** (median approach height **−20.4 Å**; all scored designs came from a
  reserve list); peptide share median **0.195**, helix share **0.744**; 0 of 36
  designs reached 3-of-3. De-wetting moved the peptide share by **+0.002** across
  13 pairs and, applied to `bb54-g2`, **destroyed the gate 3/3 → 0/3** while moving
  the helix share +0.007. Pooled over n = 12 G12D vs 15 WT draws, `bb54-g2`'s WT
  area loss is **29.4 %**, below the 30 % leg — the criterion applied verbatim to
  more data. Apo burial ~1.9× on-target (1206.7 vs 627.3 Å² medians).
* **Committed:** `SUMMARY.md`, `panel_checkpoint.jsonl`, `analysis.json`,
  `design_stats.json` + `_w2` (including `seed: 8080` and the shape thresholds),
  `dewet_results.json`, `apo_probe.json`, `pooled_wt.json`, `panel_results.json`,
  `candidates*.json`, `backbone_poses*.json`, 28 curated mmCIF (41 tracked).
* **Missing:** `results/batch-08/structures/` is itself gitignored with a
  whitelist of the curated set, so 101 of the 129 structures on disk are not in
  git; the SUMMARY states the rest are reproducible from
  `panel_checkpoint.jsonl`, which records every structure's SHA-1. No raw
  responses.
* **Sampling settings recorded?** Yes, per row, plus the RFdiffusion seed and
  `diffusion_steps` at batch level. This is the most completely parameterised of
  the generation batches.
* **Rerunnable?** Re-analysable for all reported medians; re-runnable with API
  budget (an RFdiffusion repeat at seed 8080 is the closest thing in the project to
  a deterministic regeneration, though the hosted endpoint is unversioned).

## batch-09 — control-calibrated audit and the matched-anchor counter-screen

* **Tested (audit, 0 model calls):** whether the WT geometry-loss criterion is
  control-calibrated. **Tested (panel, 12 calls / 36 structures):** whether apo
  over-burial survives when all four arms get the same HLA regional anchor (binder
  within 10 Å of A69 and A155).
* **Concluded (audit):** the criterion is retired as a hard gate, because the
  known-selective JDIa41b1 TCR (>4000-fold, 743 pM vs 3.0 µM) binds WT in a nearly
  identical pose (0.29 Å) via indirect electrostatics. Historical batches are left
  unchanged and continue to report the criterion they used.
  **Concluded (panel):** apo over-burial persists under a matched prior —
  apo/G12D total-area ratios **1.337** (`bb50-g1`), **1.654** (`bb54-g2`),
  **1.850** (`bb25`); WT peptide-area margins **−23.7 % / −17.1 % / −3.1 %**;
  G12D peptide shares **0.218 / 0.276 / 0.268**; Asp6 read-out 3/3, 1/3, 2/3. The
  document states plainly that buried area is not an affinity measurement.
* **Committed:** `SUMMARY.md`, `MATCHED_ANCHOR_RESULTS.md`,
  `matched_anchor_checkpoint.jsonl` (36 rows), `matched_anchor_analysis.json`,
  `matched_anchor_run.json`, **`matched_anchor_responses.tar.gz` (3.0 MB) with a
  `.sha256` manifest**, and all 36 mmCIF in `structures/` (43 tracked files), plus
  `docs/BATCH09_RUNBOOK.md`, `scripts/batch09_matched_anchor.py`,
  `scripts/batch09_analyze.py` and two test modules. The source commit
  (`75176b2`), the GitHub Actions run id and the artifact SHA-256 are recorded.
* **Missing / incomplete:**
  * **No RMSD of any kind.** Verified: zero `rmsd` fields across the 36 records,
    and `batch09_analyze.py` computes none. The "1.58 Å / 3.78 Å" figure quoted for
    `bb50-g1` in `SUMMARY.md` is recycled from the selectivity panel (different
    constraint, different settings) and is not a measurement of these poses.
  * **No `asp6_buried_frac`.** Verified: zero occurrences across the 36 records.
    The surviving gate has two legs (donor ≤ 3.5 Å **and** `asp6_buried_frac > 0`);
    only the donor leg was evaluated, under the field `asp6_gate_pass`.
  * `results/batch-09/responses/` and `recovered-evidence/` are gitignored; the
    archive is the committed form.
* **Sampling settings recorded?** Yes, completely: per-row `sampling`,
  `constraint`, `constraint_kind`, `structure_sha1`,
  `distinct_structures_in_call`, plus `matched_anchor_run.json` with per-design
  settings and the call count.
* **Rerunnable?** **Fully re-analysable with no model calls** — the only directory
  in the project of which this is true end to end (structures, hashes and raw
  responses are all committed). Re-runnable with API budget via the Actions
  workflow or `scripts/batch09_matched_anchor.py`.

## submission-review — the candidate pool

* **Tested:** an inventory of 213 distinct sequences recorded downstream in batches
  04–09, the sampling study and the selectivity panel, under an explicit curation
  policy; no model calls.
* **Concluded:** **11 candidates retained**, 9 short of the Track-1 minimum of 20.
  147 had no recorded successful direct Asp6 read-out; 37 retired with the
  batch-08 shape/de-wetting family; 16 with the rigidified family; 2 excluded as a
  conservative bb135 family curation. `format_preview_only.csv` passes the CSV
  schema validator, which does not enforce the Track-1 count.
* **Committed:** `README.md`, `review_shortlist.csv`, `candidate_audit.csv`,
  `evidence_pool.json` (462 KB, with source file SHA-256s and per-record
  locations), `format_preview_only.csv` (5 tracked files);
  `scripts/build_evidence_pool.py` and `tests/test_evidence_pool.py`.
* **Missing:** batches 01–03 are deliberately out of scope, so this is not a census
  of every sequence generated; the README says so. Rows lacking retrievable
  coordinates remain explicitly unhashed.
* **Rerunnable?** Yes — `python scripts/build_evidence_pool.py` regenerates it from
  committed artefacts only, and the CSV exports were re-generated for consistent
  line endings in commit `9caefe7`.

---

## The gaps that matter most

1. **The positive control's structures are not in git.** Five mmCIFs in
   `results/validation/` (plus 28 in `constrained-rescore/` and 4 in
   `constraint_test/`) are caught by the global `*.cif` ignore, which whitelists
   only `results/*/structures/`. The control is the pivot of the whole record and
   its RMSD check cannot be re-run from a clone. Cheapest high-value fix in the
   repo (~15 MB total, or ~3 MB for the control alone).
2. **Raw model responses exist for batch-09 only.** Batches 01–08, the sampling
   study, the selectivity panel and both metric post-mortems kept derived metrics
   and (sometimes) coordinates, never the responses. Any question that needs a
   field the analysis code did not extract at the time is unanswerable.
3. **No model versioning anywhere.** Not one artefact records which hosted
   RFdiffusion / ProteinMPNN / Boltz-2 build produced it. Re-running any batch
   today compares against an unknown.
4. **Seeds are partially recorded.** `--seed` defaults exist (1, 4042, 5050, 6060,
   7070, 8080) but only batch-08 writes the seed into its artefacts; and the
   Boltz-2 calls themselves are stochastic with no exposed seed, so a byte-exact
   rerun is impossible in principle for the scoring half.
5. **Environment is unpinned.** Lower-bound requirements only; buried-area numbers
   depend on Biopython's Shrake–Rupley, and the MM-GBSA work needs OpenMM, which is
   not declared at all.
6. **batch-09 lacks half its own gate and all pose data** (no `asp6_buried_frac`,
   no RMSD), and its SUMMARY quotes an RMSD recycled from an earlier protocol.
7. **Structure coverage is uneven.** Fully covered: sampling study (39/39),
   selectivity panel (57/57), batch-09 (36/36). Partially covered: batch-06
   (30 panel structures), batch-07 (18 of 216), batch-08 (28 of 129). Effectively
   uncovered: batches 01–05, where only 0–16 curated files are tracked and the
   bulk dumps are ignored — and in batch-04 36 files stand for 54 rows because of
   a filename collision.
8. **`results/gate-audit/` has no generator script**, so one of the audit results
   the record leans on cannot be reproduced by command.
9. **Early checkpoints carry no per-row settings or hashes.** Batches 01–05 record
   `boltz_kwargs` only at batch level and no `structure_sha1`; the
   byte-identical-response artefact discovered later cannot be audited for them.
10. **A retired metric is still the code default.**
    `src/make_submission.py:32` keeps `DEFAULT_RANK_COLUMNS = ("specificity_gap",
    "score")`, and `src/design_gates.py` still enforces `SPECIFICITY_GAP_MIN` and
    the ipTM thresholds. Documented here as an open defect; not fixed on this
    branch.
