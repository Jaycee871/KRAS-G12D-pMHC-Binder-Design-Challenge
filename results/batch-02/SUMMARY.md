# Batch 02 - pose-filtered, Asp6-readout gated

Run with `src/run_batch.py` against the hosted NVIDIA BioNeMo NIMs (RFdiffusion 2.3, ProteinMPNN, Boltz-2), 4630 s wall clock, 260 API calls (boltz2 93, proteinmpnn 47, rfdiffusion 120).

## What was produced

| quantity | count |
| --- | --- |
| RFdiffusion backbones | 120 |
| candidates offered for scoring | 60 |
| RFdiffusion backbones passing the design pose filter | 47 |
| candidates scored (>= the G12D state) | 60 |
| candidates surviving the Boltz-2 pose filter | 11 |
| candidates with a complete four-state score | 11 |
| Boltz-2 predictions | 93 |
| designs clearing **every** section 6 hard gate | 0 |
| designs whose interface analysis raised an error | 0 |

Retryable HTTP statuses hit: 429 x41.

Boltz-2 throttle: strictly serial, 8.0 s minimum interval between requests, 137 slots granted, 514.59 s spent waiting on the throttle (the alternative, measured in `results/calibration/`, is 31 x HTTP 429 and 6 lost calls out of 12 at 6 workers).

## Strategies

| strategy | binder length | hotspots | MPNN T | rationale |
| --- | --- | --- | --- | --- |
| `S1-core-65` | 65 | C4 C5 C6 C7 C8 | 0.1 | the batch-01 S1-core-short setting exactly (65 aa, C4-C8, T=0.1): the only strategy that produced a peptide-centred design, re-run over many more backbones |
| `S1-core-58` | 58 | C4 C5 C6 C7 C8 | 0.1 | shorter S1 variant (58 aa): less scaffold to drape onto the a1/a2 helices, so a larger share of the paratope has to land on the peptide |
| `S1-core-72` | 72 | C4 C5 C6 C7 C8 | 0.1 | longer S1 variant (72 aa): enough scaffold for a stable fold while staying in the batch-01 length band that worked |
| `S1-pos6-65` | 65 | C5 C6 C7 | 0.1 | S1 length with hotspots centred on C5-C7, i.e. on the mutated Asp6 and its two neighbours, to pull the paratope onto the mutation itself |
| `S1-pos6-68` | 68 | C5 C6 C7 | 0.15 | C5-C7 hotspots at 68 aa and a slightly higher MPNN temperature, for sequence diversity over the same position-6-centred pose |

## Gate survival (of the scored candidates)

| gate | passing | % |
| --- | --- | --- |
| peptide_buried_frac >= 0.40 | 0/60 | 0 % |
| Asp6 H-bond / salt bridge (or >= 7 pos6 contacts) | 2/60 | 3 % |
| pos6 contacts >= 3 | 8/60 | 13 % |
| net interface charge in -2..+4 | 41/60 | 68 % |
| PAE/PDE-interface(G12D) < 10 A | 60/60 | 100 % |
| pose filter (binder on the peptide) | 11/60 | 18 % |
| **all hard gates** | 0/60 | 0 % |

Ranking targets (ipTM; demoted from hard gates in batch-02 because the reported ipTM is the strongest binder-to-*any*-chain interface and so says nothing about where the binder sits):

| ranking target | meeting it | % |
| --- | --- | --- |
| ipTM(G12D) >= 0.75 | 1/60 | 2 % |
| max(WT, apo, decoy) <= 0.50 | 54/60 | 90 % |
| specificity_gap >= 0.25 | 0/60 | 0 % |
| asp6_buried_frac >= 0.75 (preferred) | 0/60 | 0 % |

## Per-strategy results

| strategy | scored | mean ipTM(G12D) | best ipTM(G12D) | mean gap | best gap | mean peptide buried frac | mean buried area on peptide | designs touching the peptide | passing all gates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `S1-core-58` | 11 | 0.213 | 0.506 | -0.405 | -0.405 | 0.03 | 21 A^2 | 2/11 | 0 |
| `S1-core-65` | 12 | 0.354 | 0.823 | -0.046 | -0.027 | 0.02 | 23 A^2 | 3/12 | 0 |
| `S1-core-72` | 13 | 0.314 | 0.732 | -0.210 | -0.126 | 0.02 | 24 A^2 | 2/13 | 0 |
| `S1-pos6-65` | 13 | 0.437 | 0.734 | -0.169 | +0.214 | 0.03 | 34 A^2 | 3/13 | 0 |
| `S1-pos6-68` | 11 | 0.267 | 0.552 | -0.371 | -0.133 | 0.05 | 51 A^2 | 4/11 | 0 |

## Where the binders actually landed

Buried area is the change in solvent-accessible surface on the *target* side: the same predicted coordinates with the binder chain present and with it deleted. A design can score a respectable ipTM by gripping the HLA alpha-3 face or beta-2m and never touch the peptide, which is what most of this batch does, so the split matters far more than the total.

| quantity | value |
| --- | --- |
| designs burying **any** peptide surface | 14/60 |
| designs with >= 3 heavy-atom contacts to peptide position 6 | 8/60 |
| designs that touch the target nowhere at all | 7/60 |
| mean buried area, total | 978 A^2 |
| mean buried area on the peptide | 30 A^2 |
| mean buried area on the HLA heavy chain | 823 A^2 |
| best peptide-buried fraction seen | 0.27 |

### Positive control: PDB 7OW6

The same analysis code run on the real HLA-A*11:01 / beta-2m / VVVGADGVGK complex (`scripts/validate_7ow6.py`), treating a known binding partner as the binder:

| reference binder | peptide buried frac | Asp6 side-chain buried frac | pos6 contacts |
| --- | --- | --- | --- |
| HLA groove (chain A) | 0.44 | **0.90** (85 -> 9 A^2) | 19 |
| TCR (chains D+E) | 0.21 | 0.00 (8.8 -> 8.8 A^2) | 7 |

The second row is the uncomfortable one: in the real complex the G12D aspartate is **already ~90 % buried by the groove itself** (only 8.8 A^2 of side chain left exposed), and the selective TCR of Poole et al. 2022 buries none of what remains while still making 7 heavy-atom contacts with position 6. Measured this way the `asp6_buried_frac >= 0.60` gate was unreachable by construction, and batch-02 replaced it with the Asp6 hydrogen-bond / salt-bridge read-out criterion.

## Pose filter

Two filters, both new in batch-02. The first reads the *RFdiffusion* complex and rejects a backbone whose binder is not on the peptide before any sequence is designed or any Boltz-2 call is spent. The second reads the Boltz-2 **G12D** prediction and stops there when the binder has been re-docked away from the peptide, saving the WT, apo and decoy calls.

| quantity | value |
| --- | --- |
| enabled | True |
| design filter: peptide contacts / residues / pos6 contacts / pos6 distance | >= 10 / >= 3 / >= 1 / <= 6.0 A |
| backbones generated | 120 |
| backbones passing the design filter | 47 |
| predicted-pose filter: peptide area / pos6 contacts | >= 25.0 A^2 / >= 1 |
| candidates surviving the predicted-pose filter | 11/60 |
| Boltz-2 calls saved by it | 147 |

## Leaderboard - top 10 by specificity gap

| # | name | strategy | ipTM G12D | WT | apo | decoy | gap | peptide buried | buried on peptide | pos6 contacts | Asp6 H-bond | salt bridge | asp6 buried | iface charge | gates failed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `KRASG12D_NIM_001_4ea249` | `S1-pos6-65` | 0.6686 | 0.2469 | 0.4548 | 0.2891 | 0.214 | 0.163 | 119.9 | 5 | False | False | 0.0 | 2.0 | 2 |
| 2 | `KRASG12D_NIM_002_57f728` | `S1-core-65` | 0.3972 | 0.3295 | 0.4245 | 0.3097 | -0.027 | 0.097 | 146.3 | 5 | False | False | 0.0 | -3.8 | 3 |
| 3 | `KRASG12D_NIM_003_4f4b4b` | `S1-core-65` | 0.3371 | 0.1884 | 0.4027 | 0.2809 | -0.066 | 0.046 | 61.1 | 4 | False | False | 0.0 | -4.9 | 3 |
| 4 | `KRASG12D_NIM_004_d455dc` | `S1-core-72` | 0.3071 | 0.4333 | 0.3539 | 0.2082 | -0.126 | 0.215 | 179.3 | 1 | False | False | 0.0 | -3.0 | 4 |
| 5 | `KRASG12D_NIM_005_e390d7` | `S1-pos6-68` | 0.3471 | 0.4009 | 0.4799 | 0.1044 | -0.133 | 0.274 | 228.0 | 9 | False | False | 0.538 | -3.0 | 2 |
| 6 | `KRASG12D_NIM_006_d32344` | `S1-pos6-65` | 0.3005 | 0.54 | 0.4117 | 0.5794 | -0.279 | 0.098 | 176.9 | 5 | False | False | 0.524 | -1.9 | 2 |
| 7 | `KRASG12D_NIM_007_00ecf4` | `S1-core-72` | 0.3897 | 0.6324 | 0.3123 | 0.6839 | -0.294 | 0.062 | 130.5 | 2 | False | False | 0.0 | -5.0 | 4 |
| 8 | `KRASG12D_NIM_008_1fc587` | `S1-core-58` | 0.226 | 0.2071 | 0.6315 | 0.4532 | -0.405 | 0.229 | 177.8 | 8 | False | False | 0.193 | 1.0 | 1 |
| 9 | `KRASG12D_NIM_009_05993b` | `S1-pos6-65` | 0.2502 | 0.492 | 0.6911 | 0.3121 | -0.441 | 0.098 | 150.7 | 4 | False | False | 0.081 | -1.8 | 2 |
| 10 | `KRASG12D_NIM_010_243fb9` | `S1-pos6-68` | 0.3055 | 0.4122 | 0.7607 | 0.1026 | -0.455 | 0.168 | 153.6 | 3 | False | False | 0.0 | -2.0 | 2 |

### Most common gate failures

* `only 0` - 95 designs
* `peptide-buried fraction` - 60 designs
* `no hydrogen` - 58 designs
* `predicted pose:` - 49 designs
* `net interface` - 19 designs
* `only 2` - 2 designs
* `only 1` - 1 designs

## What to change next

1. **Pose, again.** 14/60 scored designs bury any peptide surface and 11/60 survived the predicted-pose filter; 11 were worth a full four-state panel. The filter is doing its job as a budget saver, but it cannot create a pose Boltz-2 does not predict: the next step is to score the *designed* complex directly (a predictor that accepts a template or a docked pose), rather than asking a co-folder to rediscover it from sequence.
2. **Asp6 read-out.** 0/60 designs hydrogen bond or salt bridge to the Asp6 carboxylate. Note what the 7OW6 control says: HLA Arg114 already caps that carboxylate, so a binder has to compete with the groove for it - which is why the Poole TCR reads the aspartate indirectly and why a dense position-6 contact footprint substitutes for the bond in the gate.
3. **ipTM stays a ranking term.** Best ipTM(G12D) this batch was 0.823 and the best specificity gap +0.214. Until a design clears the geometry gates there is nothing for an ipTM threshold to filter.
4. **Best strategy so far: `S1-pos6-65`.** Concentrate backbones there, and keep spending the budget on *more backbones*, not more sequences per backbone - same-backbone sequences give correlated poses.
5. **Keep the serial Boltz-2 throttle.** The calibrated sustainable rate is ~8 requests/minute from one key. With the pose filter an off-peptide candidate now costs one call instead of four, so batch size should be planned as `candidates + 3 x (candidates surviving the pose filter)`.

