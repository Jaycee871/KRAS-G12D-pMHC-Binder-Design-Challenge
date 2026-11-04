# Batch-11 — peptide-only generation with HLA re-insertion

**Date:** 2026-11-05 (Asia/Taipei)  
**Final status:** complete after infrastructure and measurement corrections.  
**Scientific RFdiffusion outputs:** 24 / 24 requested.  
**ProteinMPNN calls:** 0 (no backbone passed the full design gate).  
**Boltz-2 calls:** 0.  
**Track-1 slots spent:** 0.  
**Historical thresholds changed:** no.

## Why this batch was different

Batches 06–10 generated or evaluated binders with HLA present and repeatedly
found HLA-dominant / apo-prone interfaces. Batch-08 had already tried multiple
peptide hotspots plus shape steering, so another hotspot retune would not have
been a new mechanism.

Batch-11 therefore removed HLA and B2M from the generative target. RFdiffusion
saw only the crystallographic KRAS G12D peptide conformation. HLA was restored
after generation solely to ask whether a peptide-reading backbone could fit the
actual pMHC geometry.

The predeclared design gate required:

1. Asp6 side-chain engagement;
2. Val8 side-chain engagement;
3. at least two distinct binder residues across those contacts; and
4. at most two binder heavy atoms within 2.2 A of restored HLA A1–180.

## Infrastructure correction

The first workflow attempt was not scientific evidence. A workflow-escaping
bug wrote GitHub expressions as literal `${{ ... }}` strings, causing all
RFdiffusion requests to receive an invalid bearer token and return HTTP 401.
That attempt is preserved in `INVALID_RUN.md` and must not be interpreted as
a 0/24 design result.

The workflow escaping was fixed and the second execution returned all 24
RFdiffusion backbones with zero generation-service failures.

## Side-chain measurement correction

RFdiffusion's returned PDB keeps the target peptide sequence but emits only
N/CA/C/O atoms for those kept residues. The first design-stage scorer therefore
reported Asp6 and Val8 side-chain distances as unavailable.

No threshold was changed. The 24 preserved backbones were reanalysed offline,
without model calls, by rigidly aligning their peptide CA atoms to 7OW6 and
restoring the crystallographic peptide side-chain coordinates.

Use `design_reanalysis.json` and `backbone_reanalysis.json` as the
authoritative batch-11 geometry result. The original `design_run.json` is
retained as an execution record.

## Corrected result

- **9 / 24** backbones satisfy the Asp6 + Val8 bidentate contact requirements
  before considering HLA.
- **0 / 24** satisfy the HLA steric-fit requirement.
- **0 / 24** satisfy the full frozen design gate.
- severe HLA clash atoms: **41 minimum, 154 maximum**, versus the frozen
  allowance of <=2.
- peptide CA alignment RMSD across the generated set is approximately
  **0.11–0.17 A**, so the clash result is not explained by a bad peptide
  superposition.

The strongest peptide-reading examples make the intended contacts but occupy
space that belongs to the HLA groove once the real pMHC is restored. For
example, `B11-peptide40_bb03` reaches Asp6 at 3.452 A and Val8 at 2.230 A,
yet still has 42 severe HLA clash atoms.

## Interpretation

The peptide-only intervention answered the intended mechanistic question.

Removing HLA from generation can make RFdiffusion focus on the peptide:
nine backbones satisfy the predeclared two-residue peptide-contact geometry.
But removing HLA also removes the physical boundary that makes those approach
directions accessible in a real pMHC. The resulting backbones commonly wrap
through the HLA-occupied volume.

This is therefore not evidence that peptide-focused generation is impossible.
It is evidence that **peptide focus without simultaneous HLA exclusion geometry
is physically underconstrained for this target**.

## Decision

Generation is frozen.

The predeclared stop rule said that a peptide-only batch producing no full-gate
advancer would end backbone generation rather than trigger another
result-driven threshold retune. We therefore do not loosen the 2.2 A clash
criterion, do not rerun the same mechanism with more seeds, and do not create
another generic rescue batch before the deadline.

The project now moves to a transparent, evidence-tiered fallback Track-1
package from already generated sequences.

## Evidence files

- `design_backbones.tar.gz` — the 24 RFdiffusion PDBs from the valid run.
- `backbone_poses.json` — original execution-time geometry record; peptide
  side-chain fields are not authoritative because RFdiffusion omitted them.
- `backbone_reanalysis.json` — corrected per-backbone geometry using restored
  7OW6 side chains.
- `design_reanalysis.json` — corrected aggregate result; 0 model calls.
- `INVALID_RUN.md` — audit record of the earlier workflow-authentication
  failure.
