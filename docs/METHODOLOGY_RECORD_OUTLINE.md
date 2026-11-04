# Methodology record — proposed section structure

**Status:** proposal only. Nothing in `results/` is altered by this document, no
number or conclusion in any existing file is restated with a different value, and
`src/design_gates.py` is untouched.

**What the record is for.** The campaign produced no submitted design (0 of 40
Track-1 slots used; deadline 2026-11-08 23:59 AoE) and no candidate with credible
evidence of G12D selectivity. The contribution being written up is therefore not a
binder but a decision record: what we proposed to rank designs by, how a positive
control was used to decide whether a ranking quantity was allowed to make
decisions, and what happened when hypotheses about the interface were tested.

**Scope limit, to be held in every section.** Every negative statement in this
record is a statement about computational diagnostics — Boltz-2 co-folding
predictions, Shrake–Rupley buried-area partitions, contact censuses, and a
single-structure MM-GBSA calculation. The defensible form is "we did not establish
credible evidence of selectivity". The indefensible form, which must not appear,
is "these sequences have been shown not to bind". No sequence in this project has
been tested experimentally.

**Precision about criteria changes.** The record must not claim that no threshold
was ever loosened. One hard criterion was retired: the WT geometry-loss gate
(≥ 30 % loss of both median peptide buried area and median position-6 contacts on
WT), retired in `results/batch-09/SUMMARY.md`. The principle to state instead is:

> Historical results are reported under the criteria in force at the time, and
> every later revision of a criterion carries an explicit reason and a timestamp.

That, and not an unbroken threshold, is what separates a legitimate scientific
correction from moving the goalposts after seeing the results.

---

## Section 1 — The target and the central difficulty

**Content.** KRAS G12D decamer `VVVGADGVGK` on HLA-A\*11:01 (PDB 7OW6); comparators
WT `VVVGAGGVGK` (differs at peptide position 6 only), peptide-free HLA+β2m ("apo")
and the unrelated HLA-A\*11:01 decoy `AVFDRKSDAK`. Why the problem is hard: the
discriminating side chain is one aspartate lying inside a groove whose α1/α2
helices present most of the available surface.

**Evidence in the repo.** `README.md` (target table, three objectives);
`docs/PIPELINE.md` §0 "The central difficulty" and §1 "Target preparation";
`src/fetch_target.py`; `data/targets/` (FASTAs committed, `TARGETS.md`);
`docs/LITERATURE.md` (Poole et al. 2022, Xie et al. 2025).

**Not reproducible from what is committed.** `data/raw/` and all `*.pdb`/`*.cif`
target coordinates are gitignored; `data/targets/*.pdb` must be regenerated with
`python src/fetch_target.py`, which needs network access to RCSB. No checksum of
the downloaded 7OW6 entry is recorded anywhere, so a future RCSB revision would
change the inputs silently. **To fix:** commit a SHA-256 of the fetched PDB/mmCIF
and the fetch date next to `data/targets/TARGETS.md`.

---

## Section 2 — What the pipeline actually was

**Content.** RFdiffusion → ProteinMPNN (plus deterministic donor-anchored grafts)
→ hosted Boltz-2 co-folding, optionally with a flat pocket constraint →
interface analysis. State that all model calls went to hosted NVIDIA BioNeMo NIMs,
not to local installs.

**Evidence.** `docs/PIPELINE.md` §§2–6 and "Confirmed endpoint schemas";
`src/nvidia_client.py` (202 + `nvcf-reqid` polling, retry handling);
`src/run_batch.py`, `src/run_pipeline.py`, `src/complex_analysis.py`,
`src/groove_geometry.py`, `src/pose_filter.py`, `src/design_gates.py`;
`scripts/run_batch04.py` … `scripts/run_batch08.py`,
`scripts/batch09_matched_anchor.py`. Throttle calibration:
`results/calibration/boltz2_rate_serial{5,20}.json`,
`boltz2_rate_interval8.json`, `boltz2_rate_parallel6.json` — e.g. 20 serial calls
with no minimum interval gave 9.76 calls/min with 12 × HTTP 429
(`boltz2_rate_serial20.json`), which is why every later batch ran strictly serial
at ≥ 8 s spacing. The validated flat pocket-constraint schema (and the three
rejected forms returning HTTP 422) is in
`results/constrained-rescore/FINDINGS.md`.

**Not reproducible.** No version or digest of the hosted RFdiffusion, ProteinMPNN
or Boltz-2 NIM is recorded in any run artefact — only the request bodies. A
hosted-endpoint update would be invisible in this repo. `requirements.txt` pins
nothing beyond lower bounds (`numpy>=1.24`, `pandas>=2.0`, `biopython>=1.81`,
`requests>=2.31`), and the buried-area numbers depend on the Biopython
Shrake–Rupley implementation. **To fix:** record the NIM response metadata /
function id per call, and add a lock file or exact pins for the analysis stack.

---

## Section 3 — Ranking metric 1: the ipTM specificity gap, and its retirement

**Arc.** Proposed as the ranking quantity → checked against the 7OW6 positive
control → the known-selective Poole TCR scores a *negative* gap → retired.

**Evidence.** `results/validation/CONTROL_7OW6.md`:
`specificity_gap = ipTM(G12D) − max(WT, apo, decoy) = 0.9268 − 0.9325 = −0.006`;
replicate spreads ipTM(G12D) 0.9268 ± 0.0010, WT 0.9325 ± 0.0001, decoy
0.9136 ± 0.0005, apo 0.7844; restricting to the binder↔peptide pair does not
change the sign (0.8966 vs 0.9001). Design-regime irreproducibility, same file:
`S1-pos6-68_bb18` +0.298 → −0.300, `bb19` +0.221 → −0.041, `bb13` +0.183 → +0.151;
repeat ipTM(G12D) 0.486/0.619/0.524 (SD 0.056) and 0.277/0.758/0.758 (SD 0.227).
Raw data `results/validation/7ow6_constrained_control.json`,
`7ow6_constrained_checkpoint.jsonl`, `7ow6_replicates.jsonl`,
`7ow6_replicates_summary.json`, `results/constrained-rescore/replicate_top3.jsonl`,
`g12d_repeat_spread.json`. The geometry half of the same control *passed*:
peptide buried fraction 0.224 predicted vs 0.207 crystal, total buried area 825 vs
813 Å², position-6 contacts 7 vs 7, variable-domain Cα RMSD 1.29 Å (α) / 1.51 Å (β)
in a pMHC frame of 0.85 Å (`results/validation/7ow6_control_rmsd.json`).

**Not reproducible.** The five control mmCIFs
(`results/validation/7ow6_control_*.cif`) are **present on disk but untracked** —
`*.cif` is gitignored and the whitelist only covers `results/*/structures/`. A
clone therefore has the control's numbers but not its coordinates, and
`scripts/control_7ow6_rmsd.py` cannot be re-run without new API calls. No raw API
responses were kept for this control. **To fix:** move the five control structures
under `results/validation/structures/` (whitelisted) or add an explicit
`.gitignore` exception; they are ~3 MB in total.

---

## Section 4 — Ranking metric 2: the paired MM-GBSA ddG, and its retirement

**Arc.** Built specifically to cancel the variance that killed metric 1 (one pose,
two states, Asp6→Gly as a pure deletion) → same positive control → wrong sign at
all four parameter corners → retired without a fix.

**Evidence.** `results/ddg-scoring/FINDINGS.md`: crystal 7OW6 `ddG_sel` **+1.27 ±
0.81** where negative is required; predicted 7OW6 +0.81 ± 1.90; off-site p3 control
−0.83 ± 1.58; **decoy −3.53 ± 0.12 — the strongest and most confident
"G12D-favouring" value in the panel, on a peptide that cannot have a G12D
preference**; wrong sign at all four corners (+1.27, +2.32, +7.13, +9.60 kcal/mol).
Diagnosis to carry over verbatim in substance: (1) the quantity measures
side-chain **burial**, not recognition — the groove's apolar Val3→Gly gives −7.32
against the real Arg114–Asp6 salt bridge's −7.42
(`results/ddg-scoring/groove_control.jsonl`); (2) 7OW6 is **ill-posed** for any
single-pose direct-contact energy, because no TCR atom comes within 5 Å of the
Asp6 side chain (nearest 5.73 Å, 0 atoms < 4 Å) while HLA Arg114 NH2 sits at
2.95 Å; (3) the numerical floor (SD 0.8–8.7, parameter sensitivity 8.3 kcal/mol)
sits above the ~1 kcal/mol signal. `scripts/ddg_rescore_designs.py` is committed
and **was never run against a design**. Raw data `control_7ow6_ddg.json`,
`control_diagnostics.jsonl`, `groove_control.jsonl`; code `src/ddg_scoring.py`,
tests `tests/test_ddg_scoring.py`.

**Not reproducible.** The MM-GBSA path needs OpenMM, which is not in
`requirements.txt` at all, and no force-field/OpenMM version is recorded in the
result JSONs. The contact census was computed from `data/raw/7OW6.pdb`, which is
gitignored. **To fix:** pin the OpenMM/ff14SB/GBn2 versions in the output JSON and
add the 7OW6 raw-file hash.

---

## Section 5 — The positive control as a decision instrument

**Content.** This is the methodological core and should be its own section, not an
aside. The control did not merely annotate the metrics; it changed research
decisions: batch-04 as originally designed (~80 designs ranked by
`specificity_gap`) was **cancelled** on the control's evidence, then
**reinstated with a different purpose** (change the designing, not the scoring)
only after the gate audit showed the shortfall was in the designs.

**Evidence.** `results/validation/CONTROL_7OW6.md` ("no new design batch was
launched"); `results/batch-04/SUMMARY.md` §"Note on the decision this file
replaces"; `results/gate-audit/ASP6_GATE_AUDIT.md`: **0 of the 88 designs scored
with the Asp6 read-out (batches 02 and 03) form any hydrogen bond or salt bridge
to Asp6**, `asp6_min_donor_distance` unpopulated in every row, best
`peptide_buried_frac` 0.361, and batch-01 was never scored for the read-out at
all.

**Also record here:** what the control costs. The same control shows the real,
experimentally validated TCR **fails two of six hard gates and three of four
ranking targets** (`peptide_buried_frac` 0.224 vs a 0.40 floor; net interface
charge −2.9 against a −2.0 floor).

**Not reproducible.** As §3: the control's coordinates are untracked. The audit
itself (`ASP6_GATE_AUDIT.md`) has **no script committed** — it is a prose file with
no generator, so its per-batch counts cannot be regenerated by command, only
recomputed by hand from `results/batch-0{1,2,3}/designs.csv`. **To fix:** add the
audit script, or a short recomputation snippet, next to the markdown.

---

## Section 6 — The one unresolved criterion: the 0.40 peptide-share floor

**Content.** `src/design_gates.PEPTIDE_BURIED_FRACTION_MIN` is **0.40** and has
never been changed. `results/validation/CONTROL_7OW6.md` proposed **0.18** on the
evidence that the real TCR measures 0.224 predicted / 0.207 crystal, i.e. the floor
rejects the positive control itself. The recalibration was deliberately left
unapplied and recorded as the user's decision in batch-04, batch-05, the sampling
study, batch-06, batch-07 and batch-08. Batch-09 further demoted it to
"diagnostic / ranking evidence, not a hard gate" in prose while the constant in
code remains 0.40. This should be written up as an **open item with a visible
tension**, not as a resolved question: a threshold that is simultaneously (a) in
force in code, (b) known to reject the positive control, and (c) described in the
latest batch document as unsuitable as an absolute gate.

**Evidence.** `src/design_gates.py` lines 27–50 (`PEPTIDE_BURIED_FRACTION_MIN =
0.40`, `NET_INTERFACE_CHARGE_MIN = -2.0`, `SPECIFICITY_GAP_MIN = 0.25` still
present); `results/validation/CONTROL_7OW6.md` "Proposed gate recalibration
(evidence, not applied)"; per-batch restatements in
`results/batch-04/SUMMARY.md`, `results/batch-05/SUMMARY.md` (0 of 60 clear 0.40,
42 clear 0.18), `results/sampling-study/SUMMARY.md`,
`results/batch-06/SUMMARY.md` (0 of 41 clear 0.40, 26 clear 0.18),
`results/batch-07/SUMMARY.md` (0 of 183 draws clear 0.40, 137 clear 0.18),
`results/batch-08/SUMMARY.md` (1 of 114 draws clears 0.40 — max 0.48 — and that
draw is a 386 Å² near-collapsed interface with no donor within 3.5 Å).

**Reproducible?** Yes, from the committed checkpoints; `src/design_gates.evaluate`
and `tests/test_design_gates.py` are committed and run in CI.

---

## Section 7 — The retirement of the WT geometry-loss gate (batch-09)

**Content.** The criterion: `scripts/batch06_panel.WT_CRITERION` — median peptide
buried area **and** median position-6 heavy-atom contacts both ≥ 30 % lower on WT,
over 3 draws per state in the same cell, with the Gly-by-construction Asp6
side-chain loss excluded. It was pre-registered before batch-06's results and
applied verbatim in batches 06, 07 and 08. It was **retired** in
`results/batch-09/SUMMARY.md` with this reason: the known-selective JDIa41b1 TCR
(Poole et al., Nat Commun 2022, DOI 10.1038/s41467-022-32811-1) has 743 ± 18 pM
affinity for G12D vs 3.0 µM for WT (>4000-fold) while binding WT in a nearly
identical pose — **0.29 Å** Cα RMSD over the variable domains — with selectivity
attributed to indirect electrostatics. A criterion that demands a selective binder
lose 30 % of its WT footprint is therefore not control-calibrated.

**The timing must be stated plainly in the same paragraph.** The retirement
happened in the batch in which all three leads were failing that criterion:
bb54-g2's pooled WT area loss had decayed to **−29.4 %** against the 30 % leg
(`results/batch-08/SUMMARY.md`, `results/batch-08/pooled_wt.json`), and in the
batch-09 matched-anchor panel the WT peptide-area margins are **−23.7 %
(bb50-g1), −17.1 % (bb54-g2), −3.1 % (bb25)**
(`results/batch-09/matched_anchor_analysis.json`, ratios 0.763 / 0.829 / 0.969).
A reader is entitled to weigh the argument and the timing together, and the
document should give them both rather than only the argument.

**Date and authorship of the change:** `results/batch-09/SUMMARY.md`, dated
2026-11-04 (Asia/Taipei), commit `a26f784`.

**Evidence.** `scripts/batch06_panel.py` (criterion text), `results/batch-06`,
`-07`, `-08` SUMMARY.md WT sections, `results/batch-08/pooled_wt.json`,
`results/batch-09/SUMMARY.md`, `docs/LITERATURE.md`.

**Not reproducible.** The Poole affinity and 0.29 Å figures are literature values
cited from the paper, not computed here; the record should mark them as such.

---

## Section 8 — Three improvement hypotheses, each tested, each kept

Each sub-section states the hypothesis, the pre-stated test, the measured outcome
and the verdict. All three are negative and are reported as such.

**8.1 Shrink the binder (batch-06).** 34/40/46 aa against 62–72 aa; a
peptide-centred design-stage pose filter (≥ 0.25 of target contact on the peptide,
≤ 90 heavy atoms on the groove). Result: every part of the interface shrank ~21 %
(total 960.1 → 760.8 Å², HLA 772.9 → 591.3, peptide 201.2 → 158.2) and the peptide
share of buried area moved **0.209 → 0.213**. 28 of 264 backbones survived the
filter; 0 of 41 designs cleared the Asp6 gate 3 of 3. The most peptide-centred
design (`bb38-g1`, share 0.32) lost the aspartate in 6 of 6 panel draws. Evidence:
`results/batch-06/SUMMARY.md`, `analysis.json`, `checkpoint.jsonl` (41 rows),
`panel_checkpoint.jsonl` (30), `wt_diagnostic.json`.

**8.2 Rigidify the graft (batch-07).** bb54-g2's donor anchor was already on a
helix (`shell_ss` HHHH, donor-loop length 0), so the edits targeted connectors and
termini with the six Asp6-shell positions protected; no Cys–Cys staple was formable
on the backbone, so no `rigid-gps` variant was ever built. Result: 1 of 21 designs
hit 3-of-3 once (`bb54_rigid-gp-p2`) and gave 1/3 on an identical repeat;
self-rigidifying the parent made it worse (1/3 vs 3/3); Track B fresh rigid
scaffolds were 0 of 6. The batch's useful positive finding is separate: cell
`r3np` (recycling 3, `without_potentials` true) gives bb54-g2 8 of 9 gated draws,
and **pose tightness and gate reproducibility are anti-correlated** — the tightest
cell `r10np` (2.25/2.81 Å) gates only 2/3 while the 3-of-3 cell sits at
13.93/16.97 Å. Evidence: `results/batch-07/SUMMARY.md`, `analysis.json`,
`panel_checkpoint.jsonl`, `screen_summary.jsonl`, `screen_stats.json`,
`cross_condition.json`, `depth_results.json`, `pooled_wt.json`,
`scripts/rigidify.py`.

**8.3 Steer the interface by shape (batch-08).** Five shape thresholds fixed before
the run and calibrated on batch-07's 97 backbones. Result: **0 of 192 backbones
cleared all five**, median approach height **−20.4 Å** (binder centroid below the
peptide); every scored design came from a ranked reserve list. Peptide share over
114 on-target draws **0.195** (below batch-06's 0.213), helix share **0.744**; 0 of
36 designs reached 3-of-3. The paired de-wetting twins moved the peptide share by
**+0.002** (helix +0.017) across 13 pairs, and the controlled de-wetting of
bb54-g2 itself destroyed the gate **3/3 → 0/3** while moving the helix share by
+0.007. Evidence: `results/batch-08/SUMMARY.md`, `analysis.json`,
`design_stats.json`, `dewet_results.json`, `apo_probe.json`,
`panel_checkpoint.jsonl`, `pooled_wt.json`, `src/groove_geometry.py`.

**Framing note for 8.1–8.3.** The comparator that makes these three results
interpretable is that the real 7OW6 TCR itself sits at peptide share 0.207 and
helix share 0.793. A groove-dominated interface is what a working binder on this
target looks like.

---

## Section 9 — Replication as a gate, and the sampling-settings discovery

**Content.** Batch-05 folded replication into the gate (3 of 3 independent
constrained draws) and produced **0 of 60** survivors, with every gate failure
coinciding with a flip between docked modes. The sampling-settings study then
showed the instability is a property of the **sampling settings**, not of the
designs: `bb50-g1` goes from 14.44/20.37 Å (baseline) to **0.70/1.05 Å** under
`recycling_steps` 8, `bb25` 7.27 → 1.96 Å, `bb135-g1` 4.54 → 2.04 Å under combo;
and the best cell is design-specific and not additive (combo rescued `bb135-g1`
but re-opened `bb25`). Three designs reached 3-of-3 that had reached it zero times
in batch-05. Batch-07 then added the anti-correlation between pose tightness and
gate reproducibility. The honest caveat to carry: **a collapsed spread is
self-consistency, not accuracy.**

**Evidence.** `results/batch-05/SUMMARY.md`, `gate_results.json`,
`replicate_rmsd.json`, `replicates.json`, `final_ranking.json`, `survivors.csv`;
`results/sampling-study/SUMMARY.md`, `checkpoint.jsonl` (39 rows),
`analysis.json`, `structures/` (39 mmCIF, tracked);
`results/batch-07/SUMMARY.md` sweep table.

**Not reproducible.** Sampling-study rows carry `boltz_kwargs` and `constraint`
but **no `structure_sha1` and no sequence string** (unlike batches 06–09), so the
duplicate-response artefact cannot be audited retrospectively for this study; the
sequences must be joined back from batch-04/05 `designs.csv` by design key.

---

## Section 10 — The off-target panels and the apo result

**Content.** The selectivity panel is where the project's central negative result
first appears: reading the aspartate is not the same as needing it. All three
designs bury more surface on peptide-free HLA (1101–1454 Å²) than on the G12D
complex (627–1076 Å²); `bb135-g1` salt-bridges the decoy's Asp8 in 3 of 4 draws;
`bb25` buries 196 Å² of WT peptide against 190 Å² of G12D. The one genuinely good
sign is `bb50-g1`'s decoy behaviour (peptide footprint 238 → 62 Å², position-6
contacts 9–17 → 0–2), and it should be reported at exactly that weight: one arm of
one design, not evidence of selectivity.

**Evidence.** `results/selectivity-panel/SUMMARY.md`, `analysis.json`,
`checkpoint.jsonl` (57 rows), `structures/` (57 mmCIF, tracked);
`results/batch-06/SUMMARY.md` §4; `results/batch-08/SUMMARY.md` §6
(apo 1110–1287 Å² vs 542–758 Å² on target).

**Methodological caveat to carry through.** Until batch-09, the apo arm was
**unconstrained** (no chain C to anchor the pocket term to), so it was never
like-for-like with the three constrained arms. Both the selectivity panel and
batches 06–08 state this.

---

## Section 11 — The matched-anchor counter-screen (batch-09)

**Content.** The one experiment that removes the unequal-prior objection: all four
arms given the same HLA regional anchor (binder within 10 Å of HLA A69 and A155).
Result: apo/G12D total buried-area ratios **1.337 (bb50-g1), 1.654 (bb54-g2),
1.850 (bb25)**; WT peptide-area margins **−23.7 % / −17.1 % / −3.1 %**; G12D
peptide shares **0.218 / 0.276 / 0.268**; Asp6 read-out 3/3, 1/3, 2/3. 12 Boltz-2
calls, 36 structures, all SHA-1 distinct, ~150 s in a GitHub Actions runner.

**Two defects of this batch that the record must flag.**

1. **No RMSD was computed in batch-09.** `matched_anchor_checkpoint.jsonl`
   contains no RMSD field of any kind (verified: zero matches for `rmsd` across
   all 36 records), and `scripts/batch09_analyze.py` produces none. The
   "1.58 Å mean / 3.78 Å max" figure quoted for `bb50-g1` in
   `results/batch-09/SUMMARY.md` is a **recycled** number from the n = 10
   selectivity panel under a different constraint and different settings; it is
   not a measurement of the matched-anchor poses. The record should say so, and
   `MATCHED_ANCHOR_RESULTS.md` is correct in warning that these draws must not be
   pooled with the peptide-C6-constrained draws.
2. **`asp6_buried_frac` is absent from every batch-09 record** (verified: zero
   matches across the checkpoint). The surviving Asp6 gate has two legs — a donor
   within 3.5 Å of OD1/OD2 **and** `asp6_buried_frac > 0` — so only the donor leg
   was evaluated in that batch. The reported 3/3, 1/3 and 2/3 counts are
   donor-distance counts, under the field name `asp6_gate_pass`.

**Evidence.** `results/batch-09/MATCHED_ANCHOR_RESULTS.md`, `SUMMARY.md`,
`matched_anchor_analysis.json`, `matched_anchor_checkpoint.jsonl` (36 rows),
`matched_anchor_run.json`, `matched_anchor_responses.tar.gz` + `.sha256`,
`structures/` (36 mmCIF, tracked), `docs/BATCH09_RUNBOOK.md`,
`scripts/batch09_matched_anchor.py`, `scripts/batch09_analyze.py`,
`tests/test_batch09_*.py`.

**Reproducibility.** This is the best-documented batch in the project: frozen
protocol, source commit `75176b2`, Actions run id, raw responses archived with a
SHA-256 manifest, every structure hashed, and re-analysis needs no model calls.

---

## Section 12 — Decay of a margin with sample size

**Content.** A short, self-contained methodological section. `bb54-g2`'s WT
peptide-area margin: **−59.2 % at n = 3** (batch-06, under `recyc8`), **−33.4 % at
n = 9/12 distinct draws** (batch-07, pooled under `r3np`), **−29.4 % at n = 12/15**
(batch-08 pooled, below the 30 % criterion), **−17.1 %** in the batch-09
matched-anchor panel (a different protocol, so not a continuation of the same
series but the same direction). The point is general: a margin first seen at n = 3
is an excursion until more draws say otherwise, and the decay direction here was
monotone.

**Evidence.** `results/batch-06/SUMMARY.md` §3, `results/batch-07/SUMMARY.md`
pooled WT table, `results/batch-08/pooled_wt.json` and SUMMARY §5,
`results/batch-09/matched_anchor_analysis.json`.

**Caveat to state.** The three pooled values come from different numbers of draws
*and* the last from a different constraint, so the series is suggestive of a
small-n excursion, not a controlled sample-size experiment.

---

## Section 13 — Hypotheses we hold, and must not state as findings

**Content.** One short section, to keep claims and observations apart.

* **Observed:** RFdiffusion backbones for this target place the binder below the
  peptide (median approach height −20.4 Å over 192 backbones) and the resulting
  interfaces are helix-dominated (median helix share 0.744 over 114 draws).
* **Hypothesis, awaiting test:** that this happens *because* the α1/α2 helices
  present the largest contiguous surface, i.e. a generator-prior explanation. The
  repo contains the contact bias, not its cause. Nothing here tests an
  alternative, such as hotspot specification, the diffusion-step setting (30 steps
  gave 0 of 96 usable backbones in batch-07, 15 steps is what every other batch
  used), or the target representation.
* **Also a hypothesis:** that the apo over-burial reflects genuine peptide-
  independent groove affinity. What is measured is predicted buried area under a
  matched anchor, which is not a free energy and not a K_D.

**Evidence.** `results/batch-08/SUMMARY.md` §1 and "What this batch establishes";
`results/batch-07/SUMMARY.md` wave-1 note; `results/batch-09/MATCHED_ANCHOR_RESULTS.md`.

---

## Section 14 — Known defects, left open

To be listed as defects, not fixed inside the methodology-record branch.

1. **`src/make_submission.py` line 32:** `DEFAULT_RANK_COLUMNS = ("specificity_gap",
   "score")`. The `--rank-column` default therefore still falls back to the
   **retired** ipTM metric when no column is given; only `README.md` was patched to
   tell the user to pass `--rank-column evidence_rank --ascending`. A scored table
   containing a `specificity_gap` column will silently be ranked on it.
2. **batch-04 structure-file collision:** `scored_structures/` holds 36 files for
   54 rows because the name is `{backbone}_{variant}` and the two MPNN-selected
   candidates per backbone collide; the second overwrites the first. Recorded in
   `results/batch-04/run_stats.json` (`structure_file_caveat`). Checkpoint metrics
   are unaffected and both gate survivors are donor-anchored, so their files are
   intact. The directory is gitignored in any case.
3. **`SPECIFICITY_GAP_MIN = 0.25` and the ipTM gates remain live in
   `src/design_gates.py`** even though the quantity is retired in prose. Anyone
   calling `evaluate()` with a `specificity_gap` key still gets it enforced.
4. **batch-09 gaps:** no RMSD, no `asp6_buried_frac` (see §11).

---

## Section 15 — What was not done, and what it would take

A closing section listing the things a reader might reasonably expect and will not
find: no experimental validation of any kind; no positive control that a
*direct-contact* metric could pass (identified as the single most valuable missing
control in `results/ddg-scoring/FINDINGS.md`); no ensemble or alchemical free
energy; no alternative backbone generator; no submission. Plus the reproducibility
gaps collected in `docs/REPRODUCIBILITY_INVENTORY.md`.

---

## Evidence index for the record

| Claim class | Primary artefacts |
| --- | --- |
| Positive control, geometry + ipTM | `results/validation/` (7 tracked files; 5 mmCIF untracked) |
| MM-GBSA retirement | `results/ddg-scoring/` (4 files) |
| Library audit | `results/gate-audit/ASP6_GATE_AUDIT.md` (no script) |
| Constraint discovery and rescore | `results/constrained-rescore/`, `results/constraint_test/` |
| Throttle calibration | `results/calibration/` (4 files) |
| Design batches | `results/batch-01` … `results/batch-09` |
| Replication and settings | `results/batch-05`, `results/sampling-study` |
| Selectivity | `results/selectivity-panel`, `results/batch-06` … `-09` |
| Candidate pool | `results/submission-review/` (11 of 213 sequences retained) |

Per-batch reproducibility is assessed in
[`REPRODUCIBILITY_INVENTORY.md`](REPRODUCIBILITY_INVENTORY.md).
