# Batch-03

Run with `src/run_batch.py` against the hosted NVIDIA BioNeMo NIMs (RFdiffusion 2.3, ProteinMPNN, Boltz-2), 13629 s wall clock, 212 API calls (boltz2 52, proteinmpnn 40, rfdiffusion 120).

## What was produced

| quantity | count |
| --- | --- |
| RFdiffusion backbones | 120 |
| candidates offered for scoring | 28 |
| RFdiffusion backbones passing the design pose filter | 40 |
| candidates scored (>= the G12D state) | 28 |
| candidates surviving the Boltz-2 pose filter | 8 |
| candidates with a complete four-state score | 8 |
| Boltz-2 predictions | 52 |
| designs clearing **every** section 6 hard gate | 0 |
| designs whose interface analysis raised an error | 0 |

Retryable HTTP statuses hit: 429 x17.

Boltz-2 throttle: strictly serial, 8.0 s minimum interval between requests, 70 slots granted, 255.69 s spent waiting on the throttle (the alternative, measured in `results/calibration/`, is 31 x HTTP 429 and 6 lost calls out of 12 at 6 workers).

## Strategies

| strategy | binder length | hotspots | MPNN T | rationale |
| --- | --- | --- | --- | --- |
| `S1-pos6-65` | 65 | C5 C6 C7 | 0.1 | S1 length with hotspots centred on C5-C7, i.e. on the mutated Asp6 and its two neighbours, to pull the paratope onto the mutation itself |
| `S1-pos6-68` | 68 | C5 C6 C7 | 0.15 | C5-C7 hotspots at 68 aa and a slightly higher MPNN temperature, for sequence diversity over the same position-6-centred pose |

## Gate survival (of the scored candidates)

| gate | passing | % |
| --- | --- | --- |
| peptide_buried_frac >= 0.40 | 0/28 | 0 % |
| Asp6 H-bond / salt bridge (or >= 7 pos6 contacts) | 4/28 | 14 % |
| pos6 contacts >= 3 | 7/28 | 25 % |
| net interface charge in -2..+4 | 24/28 | 86 % |
| PAE/PDE-interface(G12D) < 10 A | 28/28 | 100 % |
| pose filter (binder on the peptide) | 8/28 | 29 % |
| **all hard gates** | 0/28 | 0 % |

Ranking targets (ipTM; demoted from hard gates in batch-02 because the reported ipTM is the strongest binder-to-*any*-chain interface and so says nothing about where the binder sits):

| ranking target | meeting it | % |
| --- | --- | --- |
| ipTM(G12D) >= 0.75 | 0/28 | 0 % |
| max(WT, apo, decoy) <= 0.50 | 25/28 | 89 % |
| specificity_gap >= 0.25 | 0/28 | 0 % |
| asp6_buried_frac >= 0.75 (preferred) | 0/28 | 0 % |

## Per-strategy results

| strategy | scored | mean ipTM(G12D) | best ipTM(G12D) | mean gap | best gap | mean peptide buried frac | mean buried area on peptide | designs touching the peptide | passing all gates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `S1-pos6-65` | 13 | 0.297 | 0.597 | -0.133 | -0.063 | 0.04 | 23 A^2 | 3/13 | 0 |
| `S1-pos6-68` | 15 | 0.375 | 0.641 | -0.189 | +0.076 | 0.09 | 60 A^2 | 6/15 | 0 |

## Where the binders actually landed

Buried area is the change in solvent-accessible surface on the *target* side: the same predicted coordinates with the binder chain present and with it deleted. A design can score a respectable ipTM by gripping the HLA alpha-3 face or beta-2m and never touch the peptide, which is what most of this batch does, so the split matters far more than the total.

| quantity | value |
| --- | --- |
| designs burying **any** peptide surface | 9/28 |
| designs with >= 3 heavy-atom contacts to peptide position 6 | 7/28 |
| designs that touch the target nowhere at all | 0/28 |
| mean buried area, total | 935 A^2 |
| mean buried area on the peptide | 43 A^2 |
| mean buried area on the HLA heavy chain | 687 A^2 |
| best peptide-buried fraction seen | 0.31 |

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
| backbones passing the design filter | 40 |
| predicted-pose filter: peptide area / pos6 contacts | >= 25.0 A^2 / >= 1 |
| candidates surviving the predicted-pose filter | 8/28 |
| Boltz-2 calls saved by it | 60 |

## Leaderboard - top 10 by specificity gap

| # | name | strategy | ipTM G12D | WT | apo | decoy | gap | peptide buried | buried on peptide | pos6 contacts | Asp6 H-bond | salt bridge | asp6 buried | iface charge | gates failed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `KRASG12D_NIM_001_db847d` | `S1-pos6-68` | 0.522 | 0.3555 | 0.4461 | 0.4112 | 0.076 | 0.128 | 197.3 | 4 | False | False | 0.0 | -2.0 | 2 |
| 2 | `KRASG12D_NIM_002_f1abe5` | `S1-pos6-65` | 0.2751 | 0.186 | 0.3377 | 0.1996 | -0.063 | 0.241 | 165.6 | 9 | False | False | 0.0 | 1.0 | 1 |
| 3 | `KRASG12D_NIM_003_a11e67` | `S1-pos6-65` | 0.1486 | 0.1794 | 0.3511 | 0.145 | -0.203 | 0.183 | 101.6 | 2 | False | False | 0.0 | 1.1 | 3 |
| 4 | `KRASG12D_NIM_004_b87d01` | `S1-pos6-68` | 0.3268 | 0.4995 | 0.541 | 0.2056 | -0.214 | 0.198 | 151.8 | 4 | False | False | 0.0 | -2.0 | 2 |
| 5 | `KRASG12D_NIM_005_6b375a` | `S1-pos6-68` | 0.1652 | 0.2673 | 0.3869 | 0.2491 | -0.222 | 0.277 | 125.3 | 9 | False | False | 0.0 | 0.0 | 1 |
| 6 | `KRASG12D_NIM_006_968f68` | `S1-pos6-68` | 0.1922 | 0.4209 | 0.2588 | 0.2079 | -0.229 | 0.311 | 81.3 | 8 | False | False | 0.0 | -2.0 | 1 |
| 7 | `KRASG12D_NIM_007_e6e356` | `S1-pos6-68` | 0.447 | 0.6987 | 0.6086 | 0.5784 | -0.252 | 0.205 | 105.2 | 3 | False | False | 0.0 | -2.0 | 2 |
| 8 | `KRASG12D_NIM_008_2cc093` | `S1-pos6-68` | 0.2823 | 0.5758 | 0.5592 | 0.2026 | -0.293 | 0.284 | 235.2 | 9 | False | False | 0.2 | -2.0 | 1 |
| 9 | `KRASG12D_NIM_009_64d756` | `S1-pos6-65` | 0.1421 |  |  |  |  | 0.055 | 31.2 | 0 | False | False | 0.0 | -0.9 | 4 |
| 10 | `KRASG12D_NIM_010_29e3dc` | `S1-pos6-65` | 0.5933 |  |  |  |  | 0.0 | 0.0 | 0 | False | False | 0.0 | -2.0 | 4 |

### Most common gate failures

* `only 0` - 39 designs
* `peptide-buried fraction` - 28 designs
* `no hydrogen` - 24 designs
* `predicted pose:` - 20 designs
* `net interface` - 4 designs
* `only 2` - 1 designs

## What to change next

1. **Pose, again.** 9/28 scored designs bury any peptide surface and 8/28 survived the predicted-pose filter; 8 were worth a full four-state panel. The filter is doing its job as a budget saver, but it cannot create a pose Boltz-2 does not predict: the next step is to score the *designed* complex directly (a predictor that accepts a template or a docked pose), rather than asking a co-folder to rediscover it from sequence.
2. **Asp6 read-out.** 0/28 designs hydrogen bond or salt bridge to the Asp6 carboxylate. Note what the 7OW6 control says: HLA Arg114 already caps that carboxylate, so a binder has to compete with the groove for it - which is why the Poole TCR reads the aspartate indirectly and why a dense position-6 contact footprint substitutes for the bond in the gate.
3. **ipTM stays a ranking term.** Best ipTM(G12D) this batch was 0.641 and the best specificity gap +0.076. Until a design clears the geometry gates there is nothing for an ipTM threshold to filter.
4. **Best strategy so far: `S1-pos6-68`.** Concentrate backbones there, and keep spending the budget on *more backbones*, not more sequences per backbone - same-backbone sequences give correlated poses.
5. **Keep the serial Boltz-2 throttle.** The calibrated sustainable rate is ~8 requests/minute from one key. With the pose filter an off-peptide candidate now costs one call instead of four, so batch size should be planned as `candidates + 3 x (candidates surviving the pose filter)`.

