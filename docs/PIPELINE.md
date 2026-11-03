# Design pipeline

Pipeline for Challenge 5. Sections 1-6 are implemented and have been run
end to end against the hosted NVIDIA BioNeMo NIMs — see
[Running it](#running-it) at the end of this document for the commands, the
confirmed API schemas and what is still missing. Sections 2-5 are *not* run
locally on a GPU: `src/nvidia_client.py` calls the hosted models instead.

## 0. The central difficulty

The three objectives are not three separate problems, they are one problem:

| state | peptide | required outcome |
| --- | --- | --- |
| G12D pMHC | `VVVGADGVGK` | **bind** |
| WT pMHC | `VVVGAGGVGK` | **do not bind** |
| peptide-free HLA + B2M | — | **do not bind** |
| decoy pMHC | `AVFDRKSDAK` | **do not bind** |

The G12D and WT complexes differ by **one side chain**: peptide position 6,
aspartate vs glycine. Everything else — the entire HLA-A\*11:01 α1/α2 groove
rim, B2M, the peptide backbone, the other nine side chains — is identical.

**Therefore: a binder whose interface sits on the groove rim or the helices
will bind all four target states and score zero.** The paratope has to be centred on
the solvent-exposed face of the peptide, directly over position 6, and it has
to read the Asp carboxylate specifically — a hydrogen-bond donor / salt-bridge
partner (Arg, Lys, His, Ser/Thr/Tyr pairs) buried against it, so that the
glycine version loses contact *and* pays a desolvation penalty rather than
simply losing one weak contact. Peptide-contact area should dominate the
interface; HLA contact should be the minimum needed to position the binder
(roughly: aim for ≥ 40 % of buried surface on the peptide chain).

The apo (peptide-free) requirement follows from the same design rule — a binder
that needs the peptide's backbone and position-6 side chain to complete its
interface has nothing to grip on empty HLA. The decoy row is the fourth state
added after the literature review (sources and rationale in section 1): an unrelated, well
characterised HLA-A\*11:01 10-mer. It catches the failure the first three states
miss — a binder that reads the groove rim plus generic peptide backbone and
clears the WT gate only by luck.

Two cautions from the literature review ([`LITERATURE.md`](LITERATURE.md)) bear
directly on how the interface is judged:

* **Contact counting alone can mislead.** In Poole et al. 2022 — the one paper
  on exactly this pair — the selective TCR shows *few discernible differences in
  direct contacts* between the WT and G12D complexes; the selectivity is
  attributed to **indirect electrostatic** effects of the Asp side chain. So a
  design that looks non-selective by contact count may still be selective
  energetically, and vice versa. Treat the contact count as a necessary filter,
  not as evidence of specificity, and let the buried-Asp fraction (section 6)
  and the four-way ipTM comparison carry the judgement.
* **Excess negative interface charge hurts.** Xie et al. 2025 found that the
  refinement round which raised their experimental success rate also **reduced
  excessive negative charge** at the interface. That matters here more than
  usual, because the residue we are discriminating on *is itself an aspartate*:
  the temptation is to decorate the paratope with carboxylates, which both
  repels the target and costs solubility. Present complementary **positive
  charge locally** — Arg/Lys/His, or Ser/Thr/Tyr donors — directly over position
  6, while keeping the binder's **net interface charge** from running strongly
  negative. Both are checked as a gate in section 6.

## 1. Target preparation (implemented — `src/fetch_target.py`)

Downloads RCSB **7OW6**, strips the TCR chains, waters and heteroatoms, and
writes three target states to `data/targets/`:

* `g12d_complex.pdb` / `.fasta` — HLA-A\*11:01 heavy chain (chain A, 276 aa) +
  B2M (chain B, 99 aa) + peptide `VVVGADGVGK` (chain C, 10 aa)
* `wt_complex_template.pdb` + `wt_complex.fasta` — same complex with the
  peptide sequence changed at position 6 (D→G). The PDB is a *template*
  carrying G12D coordinates; the WT side chain must be rebuilt before scoring,
  and the FASTA is the authoritative input for complex prediction.
* `apo_hla_b2m.pdb` / `.fasta` — peptide-free HLA + B2M

The counter-screen decoy is kept alongside them as
`decoy_a1101_complex.fasta` — the same HLA-A\*11:01 heavy chain and B2M with
the peptide replaced by `AVFDRKSDAK`, the EBV EBNA-3B (EBNA-4) 399–408 "AVF"
decamer, one of the two immunodominant HLA-A\*11:01 epitopes in A11-positive
donors and a standard A\*11:01 reference binder
([de Campos-Lima et al. 1994, *J Exp Med* 179:1297,
10.1084/jem.179.4.1297](https://doi.org/10.1084/jem.179.4.1297);
[Midgley et al. 2003, *J Virol* 77:11517,
10.1128/JVI.77.21.11517-11530.2003](https://doi.org/10.1128/JVI.77.21.11517-11530.2003)).
It was chosen because it is a 10-mer — same length and therefore the same
groove register and bulge geometry as `VVVGADGVGK` — is unrelated in sequence,
and carries its own acidic side chains (Asp4, Asp8), so it also catches a
paratope that has merely learned to like a carboxylate on the peptide face
rather than *this* carboxylate at *this* position. There is no PDB entry for it
on A\*11:01, so the complex is predicted rather than templated; where a
template is wanted instead, `QVPLRPMTYK` (HIV-1 Nef 73–82, PDB 1QVO) is a
crystallographically characterised A\*11:01 decamer and is a drop-in
alternative.

Hotspot definition for the design stage: peptide chain C residues **4–8**,
weighted on residue 6 (Asp), plus the small number of α1/α2 rim residues within
8 Å of that patch. Define the epitope from the peptide outward, never from the
helices inward.

## 2. Backbone generation — RFdiffusion / BindCraft

* RFdiffusion binder design (`scaffoldguided` or free hallucination) against
  the `g12d_complex` target with the hotspot residues above.
* Length sweep 55–120 aa — small enough to sit over the peptide without
  draping onto the helices, large enough to fold. Keep to the 10–250 aa
  submission window regardless.
* 3,000–6,000 backbones per length band; discard any whose predicted contact
  map does not place ≥ 3 residues within 5 Å of peptide position 6.
* BindCraft run in parallel on the same hotspot as an independent prior — its
  AF2-in-the-loop hallucination tends to produce different topologies from
  RFdiffusion, which matters for the novelty requirement.
* Topology preference: compact mini-binders (three-helix bundle, ferredoxin-like
  αβ) presenting a concave or loop-rich face. For the antibody-format tracks,
  nanobody scaffolds with an extended CDR3 that can insert into the groove.

## 3. Sequence design — ProteinMPNN

* **8–16 sequences per backbone**, sampling temperature 0.1–0.2.
* Fixed-chain mode with the target chains held fixed; design only the binder.
* Bias against free cysteine (no unpaired Cys) and against long hydrophobic
  surface patches; soft bias toward Arg/Lys/His/Tyr at positions whose Cβ points
  at peptide position 6 — complementary donors for the Asp carboxylate, and
  deliberately *not* more carboxylates (Xie et al. 2025; section 0).
* Filters: net charge −6…+6, net *interface* charge −2…+4 (section 6), no
  N-glycosylation sequon (N-X-S/T), no more than 4 consecutive identical
  residues, no Met at position 1 issues for expression.

## 4. Monomer refolding filter — AlphaFold2 (single sequence / monomer)

> **Not wired up.** The hosted DeepMind AlphaFold2 NIM
> (`/v1/biology/deepmind/alphafold2`) was retired and now answers HTTP 410, and
> Boltz-2 has no monomer-refold mode that returns a design-vs-prediction RMSD.
> `src/run_pipeline.py` therefore skips straight from section 3 to section 5 and
> relies on the four-state panel to do the filtering, which costs four complex
> predictions on sequences a cheap refold would have dropped. Closing this gap
> needs a local ColabFold/Boltz monomer run or an equivalent NIM.

Refold each designed sequence on its own and keep only:

* **pLDDT ≥ 85** (mean over the binder), and
* **RMSD < 2 Å** (Cα) to the design backbone it came from.

This removes the bulk of ProteinMPNN output cheaply and is the single most
predictive cheap filter for expression and folding. Expect ~5–15 % survival.

## 5. Complex prediction — run four times per design

For every surviving design, run complex prediction **four times**, against:

1. the **G12D** complex (HLA + B2M + `VVVGADGVGK`),
2. the **WT** complex (HLA + B2M + `VVVGAGGVGK`),
3. **peptide-free** HLA + B2M,
4. the **decoy** complex (HLA + B2M + `AVFDRKSDAK`).

Implemented with the hosted **Boltz-2** NIM
(`/v1/biology/mit/boltz2/predict`), one call per design per state, since the
AlphaFold2 NIM is retired (section 4). Record per run:

* `ipTM` (and pTM), the primary interface confidence
* **PAE-interface**: mean predicted aligned error across the binder↔target
  residue pairs (lower is better; < 10 Å is the usual cut)
* buried surface area, split into peptide-buried vs HLA-buried
* number of binder heavy atoms within 5 Å of peptide position 6
* **buried fraction of the mutated Asp at peptide position 6** (definition in
  section 6)
* net charge of the binder residues forming the interface (section 6)

The fourth run is one extra prediction per design — the cheapest available
insurance against the dominant failure mode, an allele binder that passes a
three-state panel. Run it on every design that survives section 4; there is no
point ranking a design before it is known not to be a generic A\*11:01 binder.

## 6. Scoring and ranking — the specificity gap

Primary ranking quantity, mirroring the competition's own ordering
(mutation specificity → peptide-dependent binding → affinity):

```
specificity_gap = score(G12D) - max( score(WT), score(apo), score(decoy) )
```

where `score` is ipTM (or a normalised ipTM / PAE-interface composite). Rank
descending on `specificity_gap`; this column is what `src/make_submission.py`
sorts by out of the box.

### Buried fraction of the mutated Asp (peptide position 6)

A second design objective, alongside the specificity gap. Define

```
asp6_buried_frac = ( SASA_free - SASA_complex ) / SASA_free
```

where `SASA_free` is the solvent-accessible surface area of the **Asp6 side
chain** (peptide chain C residue 6, side-chain atoms CB/CG/OD1/OD2) in the
**free pMHC** — HLA + B2M + peptide, binder absent — and `SASA_complex` is the
same atoms' SASA in the **binder-bound complex**. Compute both with the same
probe radius (1.4 Å) and the same algorithm (FreeSASA / Shrake–Rupley) on the
predicted G12D complex, so the two numbers differ only by the binder.

* **Ranking term only (since batch-02): preferred `asp6_buried_frac ≥ 0.75`.**
  It is reported for every design and used behind the specificity gap, but it
  no longer vetoes anything — see the measured caveat below.

> **Measured caveat (2026-11-02), acted on in batch-02.** Run on PDB 7OW6
> itself (`scripts/validate_7ow6.py`), the HLA groove buries 0.90 of the Asp6
> side chain and 0.44 of the whole peptide, which is the sanity check this
> definition needs to pass. But with the *selective TCR* of Poole et al.
> (chains D+E) as the binder, `asp6_buried_frac = 0.00`: the aspartate is
> already ~90 % buried by the groove, leaving only 8.8 Å² of side-chain SASA,
> and the TCR buries none of the remainder while still making 7–11 heavy-atom
> contacts with position 6. **A gate of 0.60 measured this way is unreachable
> even by a known selective binder**, so batch-02 deleted
> `ASP6_BURIED_FRACTION_MIN` and replaced it with the read-out gate below.

### Asp6 read-out gate (replaces `asp6_buried_frac ≥ 0.60`)

What a selective binder has to do is *read* the carboxylate, not bury it. The
gate is therefore chemical plus geometric, computed on the predicted G12D
complex by `src/complex_analysis.py`:

* **hydrogen bond** — a donor nitrogen, or a hydroxyl/thiol oxygen-equivalent
  (Ser OG, Thr OG1, Tyr OH, Cys SG), on the **binder** within **3.5 Å** of
  Asp6 OD1 or OD2;
* **salt bridge** — a charged donor (Arg NE/NH1/NH2, Lys NZ, His ND1/NE2)
  within **4.0 Å** of OD1 or OD2.

A design **passes** when it makes at least one hydrogen bond *or* salt bridge
to the carboxylate, **or** it makes **≥ 7 heavy-atom contacts with peptide
position 6** — the *Poole escape clause*, calibrated on the 7OW6 TCR, which
reads the aspartate indirectly because HLA Arg114/Gln70 already cap it. In both
branches ≥ 3 position-6 contacts are still required. Columns:
`asp6_hbond`, `asp6_salt_bridge`, `asp6_min_donor_distance`,
`asp6_hbond_partner`.

### Pose filters (new in batch-02)

Batch-01's finding was that the problem is the **pose**, not the sequence:
Boltz-2 re-docks most binders onto the HLA α3 face or β2m (mean 24 Å² buried on
the peptide). Two filters in `src/pose_filter.py` enforce the pose and, more
importantly, stop the budget being spent on designs that have already failed:

1. **Design-stage filter** — on the RFdiffusion complex itself (binder chain A,
   concatenated target chain B), before ProteinMPNN runs: ≥ 10 binder heavy
   atoms within 5 Å of the peptide, ≥ 3 distinct peptide residues contacted,
   ≥ 1 contact with position 6, and closest approach to position 6 ≤ 6.0 Å.
   An off-peptide backbone costs zero ProteinMPNN and zero Boltz-2 calls.
2. **Predicted-pose filter** — after the **G12D** Boltz-2 state and before the
   other three: ≥ 25 Å² buried on the peptide and ≥ 1 contact with position 6.
   A binder that Boltz-2 has docked away from the peptide costs one call
   instead of four. This is the main budget saver: a batch's cost is
   `candidates + 3 × (candidates surviving the filter)`.

The reasoning is energetic, not geometric. A binder that genuinely *buries* the
carboxylate makes the glycine variant pay a real desolvation and cavity penalty;
a binder that merely touches position 6 loses one weak contact and keeps
binding. Where compute allows, confirm with Xie et al.'s reverse-mutation score
— rescore the top designs with the D→G substitution in place and keep those with
the largest predicted ΔΔG.

Hard gates applied before ranking (geometry and chemistry only, since
batch-02):

* ≥ 40 % of buried interface area on the peptide chain
* the **Asp6 read-out gate** above (H-bond or salt bridge, or ≥ 7 position-6
  contacts)
* ≥ 3 binder heavy atoms within 5 Å of peptide position 6 (necessary, not
  sufficient — see the Poole caution in section 0)
* net interface charge in **−2 … +4**
* PAE-interface(G12D) < 10 Å
* the predicted-pose filter above

**Ranking targets, not gates (demoted in batch-02).** Boltz-2 reports ipTM for
the binder against its best-scoring chain pair, so a high ipTM is compatible
with the binder sitting on β2m; until the pose filter is proven these are
reported and ranked on, never used to veto:

* `ipTM(G12D) ≥ 0.75`
* `max( ipTM(WT), ipTM(apo), ipTM(decoy) ) ≤ 0.50`
* `specificity_gap ≥ 0.25`
* `asp6_buried_frac ≥ 0.75`

These thresholds are mirrored in code in [`src/design_gates.py`](../src/design_gates.py)
so the document and the tooling cannot drift apart; `evaluate()` returns the
list of hard gates a scored design fails and `ranking_shortfalls()` the ranking
targets it misses.

Tie-breakers, in order: larger `specificity_gap`, then higher
`asp6_buried_frac`, then lower PAE-interface on G12D, then higher monomer
pLDDT, then topological diversity (do not submit 20 near-identical sequences —
cluster at 70 % identity and take cluster representatives).

## 7. Novelty

Submissions need a novelty score **≥ 3/4**, checked automatically. Designs must
be de novo. Before each batch: BLAST / foldseek every candidate against the PDB
and against known TCR/antibody germlines, and drop anything that looks like a
retrieved natural binder or a lightly mutated known TCR from the 7OW6 entry
itself. De novo hallucinated scaffolds normally clear this easily; grafted CDR
loops from a natural TCR do not.

## 8. Submission

`src/make_submission.py` takes the scored table, ranks it, truncates to ≤ 20,
and validates the result against every documented rule before writing. Then
upload at the submit page and log the batch in `docs/SUBMISSION_PLAN.md`.

## Running it

Everything in sections 1-6 runs from this repository against the hosted NVIDIA
BioNeMo NIMs. No GPU, no model weights, no container.

```bash
pip install -r requirements.txt
export NVIDIA_API_KEY=...          # from build.nvidia.com; never commit it
python3 src/fetch_target.py        # section 1: writes data/targets/*.pdb|fasta

# section 2-6, no network and no quota: stub responses, same code paths
python3 src/run_pipeline.py --dry-run --out designs/dryrun.csv

# the real thing: 2 backbones x 4 sequences, scored against all four states
python3 src/run_pipeline.py \
    --backbones 2 --binder-lengths 70,85 --seqs-per-backbone 8 \
    --out designs/scored_designs.csv

python3 src/make_submission.py --scored designs/scored_designs.csv   # section 8
```

Batches run through `src/run_batch.py`, which wraps the same code with the
strategy table, the pose filters, the throttle and checkpointing:

```bash
# batch-02: five S1-core-short variants, pose filter on, budget capped
python3 src/run_batch.py --strategy-set batch02 \
    --backbones-per-strategy 24 --seqs-per-backbone 8 \
    --seqs-scored-per-backbone 2 --seed 2026 \
    --out-dir results/batch-02 --max-boltz-calls 430 --time-budget 10200
python3 src/summarize_batch.py --dir results/batch-02 --title "Batch 02"
```

`--strategy-set` picks `batch01`, `batch02` or `all`; `--no-pose-filter`
disables both pose filters (for a controlled comparison); `--max-boltz-calls`
and `--time-budget` stop the run cleanly at a budget, leaving the checkpoint
resumable with `--resume`.

`run_pipeline.py` writes a CSV with exactly the columns of
[`examples/scored_designs_example.csv`](../examples/scored_designs_example.csv),
plus a `*.details.json` sidecar holding the raw per-state Boltz-2 confidences
(chain-pair ipTM matrices, pLDDT, ipLDDT, interface PDE) that the CSV columns
are derived from. Designs are sorted by `specificity_gap`, and every gate
failure from [`src/design_gates.py`](../src/design_gates.py) is printed.

Cost per design is **four Boltz-2 calls** (G12D, WT, apo, decoy) plus a share of
one RFdiffusion and one ProteinMPNN call. On 2026-11-02 a one-backbone,
one-candidate live run took 22 s wall clock in total; RFdiffusion on the trimmed
target answered in ~6 s, ProteinMPNN in ~2 s, and each Boltz-2 call on the
455-residue complex in ~4 s. Use `--max-candidates` to cap spend, and
`--recycling-steps 1 --sampling-steps 30` for cheap smoke runs.

### Rate limits and the measured sustainable request rate

The hosted NIMs are shared and rate limited per key. Batch-01 was run with six
worker threads and no pacing and paid for it: **235 HTTP 429 retries**, and only
50 of 56 candidates finished a four-state score in 1506 s. Everything is now
funnelled through the token-bucket throttle in
[`src/nvidia_client.py`](../src/nvidia_client.py), which enforces a global rate
plus a per-endpoint minimum interval, holds a process-wide lock so Boltz-2 calls
are **strictly serial**, honours `Retry-After` (seconds or HTTP-date), and backs
off exponentially with full jitter up to 120 s.

The sustainable rate was measured with
[`scripts/calibrate_boltz2.py`](../scripts/calibrate_boltz2.py) against the live
endpoint on 2026-11-02 (raw JSON in `results/calibration/`):

| configuration | calls | HTTP 429s | lost calls | achieved rate |
| --- | --- | --- | --- | --- |
| 6 workers, no pacing | 12 | **31** | 6 | — (unusable) |
| serial, no pacing | 20 | 12 | 0 | 9.76 calls/min |
| **serial, 8 s minimum interval** | 12 | **0** | 0 | **7.87 calls/min** |

**The sustainable Boltz-2 rate from one API key is ~8 requests per minute
(one request every 8 s, serial).** That is the default:
`DEFAULT_BOLTZ2_MIN_INTERVAL = 8.0`, overridable with `NVIDIA_BOLTZ2_MIN_INTERVAL`
or `--boltz-min-interval`. Even at 8 s a small residue of 429s survives (the
limit is shared, so neighbours matter): the 56-candidate rescore saw roughly one
retry per seven requests, every one of them absorbed by a single backoff with no
lost candidate, against 235 retries and six lost candidates before.
RFdiffusion and ProteinMPNN are paced at 2 s and are
not serialised; the global bucket is 1 request/s with a burst of 2. Setting
`NVIDIA_THROTTLE=0` disables all of it (tests only).

Practical consequence for planning: a four-state screen is four Boltz-2 calls,
so **~30-35 s per candidate** end to end. A 56-candidate batch is ~33 minutes,
and that is the number to size batches against - not the worker count.
`src/run_batch.py` therefore checkpoints every scored candidate to
`checkpoint.jsonl` (plus its mmCIF) with an `fsync` before continuing, and
`--resume` picks up where an interrupted run stopped without re-spending calls.

### Confirmed endpoint schemas

Verified live on 2026-11-02 by reading the 422 validation bodies and NVIDIA's
NIM reference docs; the full field lists live in the docstrings of
[`src/nvidia_client.py`](../src/nvidia_client.py). Base URL
`https://health.api.nvidia.com`, `Authorization: Bearer $NVIDIA_API_KEY`.

| NIM | path | required | key optional fields | response |
| --- | --- | --- | --- | --- |
| RFdiffusion 2.3 | `/v1/biology/ipd/rfdiffusion/generate` | `contigs`, and one of `input_pdb` / `input_pdb_asset` | `hotspot_res` (`["C4",...]`), `diffusion_steps` (1-50), `random_seed` | `{"output_pdb", "elapsed_ms"}` |
| ProteinMPNN | `/v1/biology/ipd/proteinmpnn/predict` | `input_pdb` (or asset) | `input_pdb_chains` (chains to **design**), `num_seq_per_target` (1-100), `sampling_temp` (array), `ca_only`, `use_soluble_model`, `random_seed`, `pssm_*` | `{"mfasta", "scores", "probs"}` |
| Boltz-2 | `/v1/biology/mit/boltz2/predict` | `polymers` (1-12, each `molecule_type` + `sequence` ≤ 4096 aa) | `id`, `cyclic`, `msa`, `recycling_steps` (1-10), `sampling_steps` (10-1000), `diffusion_samples` (1-25), `step_scale`, `without_potentials` | `{"structures" (mmCIF), "iptm_scores", "ptm_scores", "pair_chains_iptm_scores", "complex_plddt_scores", "complex_iplddt_scores", "complex_ipde_scores", ...}` |

Three things about these endpoints shape the code:

* **RFdiffusion relabels chains.** In the output PDB the designed binder is
  chain **A** and every kept target segment is concatenated into chain **B**, in
  contig order — so the contig `A1-180/0 C1-10/0 70-70` comes back as a
  70-residue chain A plus a 190-residue chain B whose residues 181-190 are the
  peptide. ProteinMPNN is therefore called with `input_pdb_chains: ["A"]`, which
  designs the binder and reports `fixed_chains=['B']`.
* **ProteinMPNN echoes the input.** The first record of `mfasta` is the input
  backbone sequence (poly-glycine for an RFdiffusion backbone), not a design;
  `parse_mpnn_mfasta` drops it.
* **The global ipTM is useless here, use the chain-pair matrix.** For a
  pMHC + binder prediction `iptm_scores` is dominated by the three *native*
  interfaces (HLA/B2M/peptide) and sits near 0.79 even when the binder misses
  completely. The per-state score in the CSV is instead the largest
  binder↔target-chain entry of `pair_chains_iptm_scores`.

### Known gaps before a real design run

* **No per-residue PAE.** The hosted Boltz-2 returns `pae: null` even with
  `write_full_pae: true` (full matrices are only written to disk by a
  self-hosted NIM), so the section 6 `PAE-interface < 10 Å` gate cannot be
  evaluated as documented. `pae_interface_g12d` in the CSV is currently the
  **interface PDE** (`complex_ipde_scores`, Å) as a stand-in; it is a related
  but not identical quantity, and should be replaced by a real interface PAE
  before the gate is trusted.
* **No MSA.** Boltz-2 is called in single-sequence mode. Accuracy on the pMHC
  side would improve with an MSA for the HLA and B2M chains (the endpoint
  accepts up to three `msa` database entries per protein polymer).
* **No monomer refold filter** (section 4), and therefore no design-backbone
  RMSD or monomer pLDDT column.
* **No AlphaFold-Multimer cross-check** — the agreement requirement in
  section 5 is unmet while only one predictor is available.
* **`asp6_buried_frac` is measured against the binder-deleted complex**, not
  against a separately predicted free pMHC. That isolates the binder's
  contribution exactly, which is what the definition wants, but it does not
  capture any peptide relaxation on binder release.
* **Interface analysis now fails loudly.** All SASA, contact and charge
  measurement lives in [`src/complex_analysis.py`](../src/complex_analysis.py),
  which identifies the binder, peptide, HLA heavy chain and B2M **by content**
  (sequence and length) rather than by hardcoded chain IDs, strips waters,
  ligands, hydrogens and alternate locations, checks that peptide position 6 is
  actually an aspartate, and raises `StructureAnalysisError` when a chain or
  side-chain atom it needs is absent. Previously every one of those paths
  returned `0.0`, which is indistinguishable from a binder that genuinely buries
  nothing - and in batch-01 the zeros were genuine, so the silence cost a whole
  batch of diagnosis. The CSV now also carries `buried_area_total`,
  `buried_area_peptide`, `buried_area_hla`, `contacts_peptide`, `contacts_hla`
  and `analysis_error` so *where* a binder landed is visible without reopening
  the structures.
