# Asp6 read-out gate audit of the existing design library

Generated 2026-11-03, after the paired MM-GBSA ddG metric was retired (PR #7).

Question: under the only filter that survived the 7OW6 post-mortem - a direct binder donor within 3.5 A of Asp6 OD1/OD2 with asp6_buried_frac > 0 - does anything already scored qualify for the next submission slot?

- `batch-01/designs.csv`: 56 scored designs, no Asp6 H-bond or salt-bridge columns were recorded in this batch (pre-dates the read-out gate), max peptide_buried_frac 0.361, max asp6_buried_frac 0.629, designs with asp6_buried_frac > 0: 2, median net_interface_charge +0.0
- `batch-02/designs.csv`: 60 scored designs, Asp6 H-bond 0, Asp6 salt bridge 0, max peptide_buried_frac 0.274, max asp6_buried_frac 0.538, designs with asp6_buried_frac > 0: 4, median net_interface_charge -1.0
- `batch-03/designs.csv`: 28 scored designs, Asp6 H-bond 0, Asp6 salt bridge 0, max peptide_buried_frac 0.311, max asp6_buried_frac 0.200, designs with asp6_buried_frac > 0: 1, median net_interface_charge -1.0

**Result: 0 of the 88 designs scored with the Asp6 read-out (batches 02 and 03) form any hydrogen bond or salt bridge to Asp6, and batch-01 was never scored for it, and `asp6_min_donor_distance` is unpopulated in every row (no donor was ever found in range). No design reaches the peptide_buried_frac >= 0.40 ranking floor either (best 0.361).**

Consequences:

1. Re-ranking the existing library cannot produce a gate-passing submission. The shortfall is in the designs, not in the scoring.
2. The interfaces are net acidic (median net interface charge negative in batches 02 and 03), which is the likely reason no donor ever approaches Asp6 - exactly the failure the cancelled batch-04 ProteinMPNN donor bias was meant to fix.
3. Any gate-passing submission before the 2026-11-08 deadline requires a new generation run (RFdiffusion backbones plus donor-biased ProteinMPNN plus constrained Boltz-2), not further analysis.
