import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import audit_carboxylate_readout as audit


def test_audit_preserves_original_and_distinguishes_missing_coordinates(tmp_path,monkeypatch):
    for name in ('scripts/selectivity_panel.py','scripts/audit_carboxylate_readout.py','src/complex_analysis.py'):
        target=tmp_path/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    source=tmp_path/'results/batch-09/matched_anchor_checkpoint.jsonl'
    source.parent.mkdir(parents=True)
    original=[json.loads(l) for l in (ROOT/'results/batch-09/matched_anchor_checkpoint.jsonl').read_text().splitlines()]
    row=next(r for r in original if r['design']=='bb50-g1' and r['arm']=='g12d' and r['draw']==1)
    dest=source.parent/row['structure_file']
    dest.parent.mkdir(parents=True)
    shutil.copy2(ROOT/'results/batch-09'/row['structure_file'],dest)
    missing={**row,'draw':99,'structure_file':'structures/missing.cif'}
    source.write_text('\n'.join(json.dumps(x) for x in (row,missing))+'\n')
    before=source.read_bytes()
    monkeypatch.setattr(audit,'sources',lambda root:[source])
    report=audit.audit(tmp_path)
    assert report['coverage']=={'checkpoint_records':2,'remeasured':1,'unavailable':1,'changed_flags':1}
    assert report['remeasured'][0]['corrected_metrics']['asp6_gate_pass'] is False
    assert source.read_bytes()==before
    assert report['source_sha256'][str(source.relative_to(tmp_path))]==hashlib.sha256(before).hexdigest()
    dest.write_text(dest.read_text()+'\n')
    with pytest.raises(ValueError,match='hash mismatch'):
        audit.audit(tmp_path)
