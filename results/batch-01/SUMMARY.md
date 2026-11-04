# Batch 01 - first real design batch

Run with `src/run_batch.py` against the hosted NVIDIA BioNeMo NIMs (RFdiffusion 2.3, ProteinMPNN, Boltz-2), 2123 s wall clock, 100 API calls (boltz2 100).

## What was produced

| quantity | count |
| --- | --- |
| RFdiffusion backbones | 18 |
| candidates offered for scoring | 56 |
| candidates with a complete four-state score | 56 |
| Boltz-2 predictions | 100 |
| designs clearing **every** section 6 hard gate | 0 |
| designs whose interface analysis raised an error | 0 |

Retryable HTTP statuses hit: 429 x70.

Boltz-2 throttle: strictly serial, 8.0 s minimum interval between requests, 172 slots granted, 717.66 s spent waiting on the throttle (the alternative, measured in `results/calibration/`, is 31 x HTTP 429 and 6 lost calls out of 12 at 6 workers).

## Strategies

| strategy | binder length | hotspots | MPNN T | rationale |
| --- | --- | --- | --- | --- |
| `S1-core-short` | 65 | C4 C5 C6 C7 C8 | 0.1 | short minibinder arcing over the peptide; hotspots on the solvent-exposed peptide core C4-C8, low MPNN temperature |
| `S2-core-long` | 95 | C4 C5 C6 C7 C8 | 0.2 | longer scaffold on the same peptide-core hotspots, higher MPNN temperature for more sequence diversity |
| `S3-pos6-tight` | 80 | C5 C6 C7 | 0.15 | hotspots collapsed onto C5-C7, i.e. the mutated Asp6 and its immediate neighbours, to force a position-6-centred paratope (Xie et al. 2025 mutation-site hotspot) |
| `S4-full-groove` | 85 | C1 C3 C5 C6 C8 C10 A62 A65 A72 A146 A150 A155 A158 | 0.15 | control: whole peptide plus a1/a2 rim, the docking mode PIPELINE.md section 0 predicts will bind all four states |

## Gate survival (of the scored candidates)

| gate | passing | % |
| --- | --- | --- |
| ipTM(G12D) >= 0.75 | 1/56 | 2 % |
| max(WT, apo, decoy) <= 0.50 | 42/56 | 75 % |
| specificity_gap >= 0.25 | 1/56 | 2 % |
| peptide_buried_frac >= 0.40 | 0/56 | 0 % |
| asp6_buried_frac >= 0.60 | 2/56 | 4 % |
| pos6 contacts >= 3 | 4/56 | 7 % |
| net interface charge in -2..+4 | 52/56 | 93 % |
| PAE/PDE-interface(G12D) < 10 A | 56/56 | 100 % |
| **all gates** | 0/56 | 0 % |

## Per-strategy results

| strategy | scored | mean ipTM(G12D) | best ipTM(G12D) | mean gap | best gap | mean peptide buried frac | mean buried area on peptide | designs touching the peptide | passing all gates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `S1-core-short` | 25 | 0.286 | 0.870 | -0.202 | +0.254 | 0.04 | 23 A^2 | 6/25 | 0 |
| `S2-core-long` | 12 | 0.176 | 0.376 | -0.178 | -0.032 | 0.01 | 6 A^2 | 1/12 | 0 |
| `S3-pos6-tight` | 12 | 0.199 | 0.340 | -0.206 | +0.034 | 0.03 | 36 A^2 | 3/12 | 0 |
| `S4-full-groove` | 7 | 0.249 | 0.491 | -0.220 | +0.161 | 0.02 | 13 A^2 | 1/7 | 0 |

## Where the binders actually landed

Buried area is the change in solvent-accessible surface on the *target* side: the same predicted coordinates with the binder chain present and with it deleted. A design can score a respectable ipTM by gripping the HLA alpha-3 face or beta-2m and never touch the peptide, which is what most of this batch does, so the split matters far more than the total.

| quantity | value |
| --- | --- |
| designs burying **any** peptide surface | 11/56 |
| designs with >= 3 heavy-atom contacts to peptide position 6 | 4/56 |
| designs that touch the target nowhere at all | 11/56 |
| mean buried area, total | 568 A^2 |
| mean buried area on the peptide | 21 A^2 |
| mean buried area on the HLA heavy chain | 508 A^2 |
| best peptide-buried fraction seen | 0.36 |

### Positive control: PDB 7OW6

The same analysis code run on the real HLA-A*11:01 / beta-2m / VVVGADGVGK complex (`scripts/validate_7ow6.py`), treating a known binding partner as the binder:

| reference binder | peptide buried frac | Asp6 side-chain buried frac | pos6 contacts |
| --- | --- | --- | --- |
| HLA groove (chain A) | 0.44 | **0.90** (85 -> 9 A^2) | 19 |
| TCR (chains D+E) | 0.21 | 0.00 (8.8 -> 8.8 A^2) | 7 |

The second row is the uncomfortable one: in the real complex the G12D aspartate is **already ~90 % buried by the groove itself** (only 8.8 A^2 of side chain left exposed), and the selective TCR of Poole et al. 2022 buries none of what remains while still making 7 heavy-atom contacts with position 6. Measured this way the `asp6_buried_frac >= 0.60` gate is unreachable by construction - see *What to change next*.

## Leaderboard - top 10 by specificity gap

| # | name | strategy | ipTM G12D | WT | apo | decoy | gap | peptide buried | buried on peptide | pos6 contacts | asp6 buried | iface charge | gates failed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `KRASG12D_NIM_001_dd0baa` | `S1-core-short` | 0.6923 | 0.4379 | 0.3645 | 0.1263 | 0.254 | 0.34 | 215.9 | 9 | 0.627 | 1.0 | 2 |
| 2 | `KRASG12D_NIM_002_070e8e` | `S4-full-groove` | 0.491 | 0.2018 | 0.3297 | 0.2951 | 0.161 | 0.0 | 0.0 | 0 | 0.0 | 3.1 | 5 |
| 3 | `KRASG12D_NIM_003_b93fec` | `S1-core-short` | 0.3636 | 0.2229 | 0.2902 | 0.2398 | 0.073 | 0.0 | 0.0 | 0 | 0.0 | -1.0 | 5 |
| 4 | `KRASG12D_NIM_004_02d59d` | `S3-pos6-tight` | 0.3335 | 0.2246 | 0.2995 | 0.261 | 0.034 | 0.158 | 188.3 | 9 | 0.0 | -1.9 | 4 |
| 5 | `KRASG12D_NIM_005_c07ddd` | `S1-core-short` | 0.8699 | 0.853 | 0.3098 | 0.3821 | 0.017 | 0.0 | 0.0 | 0 | 0.0 | 2.0 | 5 |
| 6 | `KRASG12D_NIM_006_89717e` | `S1-core-short` | 0.3554 | 0.3496 | 0.3368 | 0.2087 | 0.006 | 0.0 | 0.0 | 0 | 0.0 | 0.0 | 5 |
| 7 | `KRASG12D_NIM_007_699719` | `S4-full-groove` | 0.3782 | 0.4005 | 0.365 | 0.1482 | -0.022 | 0.0 | 0.0 | 0 | 0.0 | -2.0 | 5 |
| 8 | `KRASG12D_NIM_008_b0d71c` | `S4-full-groove` | 0.2553 | 0.1955 | 0.2842 | 0.2183 | -0.029 | 0.113 | 91.2 | 1 | 0.0 | -1.0 | 5 |
| 9 | `KRASG12D_NIM_009_2aec78` | `S2-core-long` | 0.2808 | 0.2485 | 0.3129 | 0.2091 | -0.032 | 0.0 | 0.0 | 0 | 0.0 | 0.1 | 5 |
| 10 | `KRASG12D_NIM_010_7c1cc1` | `S2-core-long` | 0.2105 | 0.1622 | 0.28 | 0.1871 | -0.07 | 0.0 | 0.0 | 0 | 0.0 | 3.0 | 5 |

### Most common gate failures

* `peptide-buried fraction` - 56 designs
* `Asp6 buried` - 54 designs
* `only 0` - 50 designs
* `net interface` - 4 designs
* `ipTM(G12D) 0.210` - 2 designs
* `ipTM(G12D) 0.150` - 2 designs
* `ipTM(G12D) 0.152` - 2 designs
* `specificity_gap -0.188` - 2 designs
* `ipTM(G12D) 0.115` - 2 designs
* `ipTM(G12D) 0.136` - 2 designs
* `ipTM(G12D) 0.129` - 2 designs
* `specificity_gap -0.258` - 2 designs
* `specificity_gap -0.419` - 2 designs
* `ipTM(G12D) 0.692` - 1 designs
* `ipTM(G12D) 0.491` - 1 designs
* `specificity_gap 0.161` - 1 designs
* `ipTM(G12D) 0.364` - 1 designs
* `specificity_gap 0.073` - 1 designs
* `ipTM(G12D) 0.334` - 1 designs
* `specificity_gap 0.034` - 1 designs
* `ipTM(wt) 0.853` - 1 designs
* `specificity_gap 0.017` - 1 designs
* `ipTM(G12D) 0.355` - 1 designs
* `specificity_gap 0.006` - 1 designs
* `ipTM(G12D) 0.378` - 1 designs
* `specificity_gap -0.022` - 1 designs
* `ipTM(G12D) 0.255` - 1 designs
* `specificity_gap -0.029` - 1 designs
* `only 1` - 1 designs
* `ipTM(G12D) 0.281` - 1 designs
* `specificity_gap -0.032` - 1 designs
* `specificity_gap -0.070` - 1 designs
* `ipTM(G12D) 0.155` - 1 designs
* `specificity_gap -0.081` - 1 designs
* `ipTM(G12D) 0.292` - 1 designs
* `specificity_gap -0.083` - 1 designs
* `ipTM(G12D) 0.178` - 1 designs
* `specificity_gap -0.086` - 1 designs
* `ipTM(G12D) 0.591` - 1 designs
* `ipTM(wt) 0.505` - 1 designs
* `ipTM(apo) 0.712` - 1 designs
* `specificity_gap -0.121` - 1 designs
* `ipTM(G12D) 0.245` - 1 designs
* `specificity_gap -0.135` - 1 designs
* `specificity_gap -0.151` - 1 designs
* `ipTM(G12D) 0.205` - 1 designs
* `specificity_gap -0.165` - 1 designs
* `specificity_gap -0.177` - 1 designs
* `specificity_gap -0.178` - 1 designs
* `ipTM(G12D) 0.137` - 1 designs
* `specificity_gap -0.183` - 1 designs
* `ipTM(G12D) 0.190` - 1 designs
* `ipTM(G12D) 0.128` - 1 designs
* `specificity_gap -0.189` - 1 designs
* `ipTM(G12D) 0.376` - 1 designs
* `ipTM(apo) 0.567` - 1 designs
* `specificity_gap -0.192` - 1 designs
* `only 2` - 1 designs
* `ipTM(G12D) 0.123` - 1 designs
* `specificity_gap -0.197` - 1 designs
* `ipTM(G12D) 0.084` - 1 designs
* `specificity_gap -0.205` - 1 designs
* `ipTM(G12D) 0.107` - 1 designs
* `specificity_gap -0.213` - 1 designs
* `specificity_gap -0.216` - 1 designs
* `specificity_gap -0.231` - 1 designs
* `specificity_gap -0.234` - 1 designs
* `specificity_gap -0.238` - 1 designs
* `ipTM(G12D) 0.438` - 1 designs
* `ipTM(decoy) 0.681` - 1 designs
* `specificity_gap -0.243` - 1 designs
* `ipTM(G12D) 0.192` - 1 designs
* `specificity_gap -0.248` - 1 designs
* `ipTM(G12D) 0.124` - 1 designs
* `specificity_gap -0.252` - 1 designs
* `ipTM(G12D) 0.216` - 1 designs
* `specificity_gap -0.256` - 1 designs
* `ipTM(G12D) 0.106` - 1 designs
* `ipTM(G12D) 0.156` - 1 designs
* `specificity_gap -0.262` - 1 designs
* `ipTM(G12D) 0.085` - 1 designs
* `specificity_gap -0.269` - 1 designs
* `ipTM(G12D) 0.116` - 1 designs
* `specificity_gap -0.270` - 1 designs
* `ipTM(G12D) 0.174` - 1 designs
* `specificity_gap -0.275` - 1 designs
* `specificity_gap -0.277` - 1 designs
* `ipTM(G12D) 0.093` - 1 designs
* `specificity_gap -0.284` - 1 designs
* `ipTM(G12D) 0.340` - 1 designs
* `ipTM(apo) 0.636` - 1 designs
* `specificity_gap -0.295` - 1 designs
* `specificity_gap -0.329` - 1 designs
* `ipTM(G12D) 0.431` - 1 designs
* `ipTM(wt) 0.556` - 1 designs
* `ipTM(apo) 0.503` - 1 designs
* `ipTM(decoy) 0.765` - 1 designs
* `specificity_gap -0.334` - 1 designs
* `ipTM(G12D) 0.535` - 1 designs
* `ipTM(wt) 0.575` - 1 designs
* `ipTM(decoy) 0.872` - 1 designs
* `specificity_gap -0.337` - 1 designs
* `ipTM(G12D) 0.103` - 1 designs
* `specificity_gap -0.350` - 1 designs
* `ipTM(G12D) 0.126` - 1 designs
* `ipTM(apo) 0.520` - 1 designs
* `specificity_gap -0.394` - 1 designs
* `ipTM(G12D) 0.092` - 1 designs
* `ipTM(apo) 0.509` - 1 designs
* `specificity_gap -0.418` - 1 designs
* `ipTM(G12D) 0.208` - 1 designs
* `ipTM(apo) 0.627` - 1 designs
* `ipTM(G12D) 0.321` - 1 designs
* `ipTM(wt) 0.740` - 1 designs
* `ipTM(G12D) 0.242` - 1 designs
* `ipTM(apo) 0.662` - 1 designs
* `specificity_gap -0.421` - 1 designs
* `ipTM(G12D) 0.091` - 1 designs
* `ipTM(apo) 0.531` - 1 designs
* `specificity_gap -0.440` - 1 designs
* `ipTM(G12D) 0.125` - 1 designs
* `ipTM(apo) 0.749` - 1 designs
* `specificity_gap -0.624` - 1 designs
* `ipTM(G12D) 0.475` - 1 designs
* `specificity_gap 0.000` - 1 designs

## What to change next

1. **Score the designed pose, do not ask a co-folder to rediscover it.** 45/56 designs bury no peptide surface at all. RFdiffusion puts the backbone over the peptide, but Boltz-2 is handed sequences only, re-docks from scratch and usually parks the binder on the HLA alpha-3 / beta-2m face. Either score the RFdiffusion complex in place, or keep only predictions whose binder lands where it was designed to.
2. **Redefine the Asp6 gate so it is physically reachable.** Position 6 of VVVGADGVGK points into the groove; the real selective TCR scores 0.00 on `asp6_buried_frac` as currently defined. Gate instead on a hydrogen bond or salt bridge to the Asp6 carboxylate within 3.5 A, plus the position-6 contact count, and keep the buried fraction as a ranking term only.
3. **Treat ipTM as a ranking term, not a hard gate, until the pose problem is fixed.** Best ipTM(G12D) this batch was 0.870 and the best specificity gap +0.254; nothing came near the 0.75 / 0.25 pair, and the ipTM that is reported is the strongest binder-to-any-chain interface, so an off-peptide docking mode can still score 0.4-0.5.
4. **Spend the next batch's budget on `S1-core-short`-style backbones and on more backbones per strategy**, not on more sequences per backbone: sequences from the same backbone give strongly correlated poses, and the pose is the thing that is failing.
5. **Keep the serial Boltz-2 throttle.** The calibrated sustainable rate is ~8 requests/minute from one key; a four-state screen is therefore ~30 s per candidate. Plan batch size against that number rather than against concurrency.

