import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('evidence_pool', ROOT / 'scripts/build_evidence_pool.py')
pool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pool)


@pytest.fixture(scope='module')
def report():
    return pool.build(ROOT)


def shortlist(report):
    by_sha = {c['sequence_sha256']: c for c in report['candidates']}
    return [by_sha[s] for s in report['shortlist_sha256']]


def test_review_keeps_exact_recorded_sequences_and_flags_shortfall(report):
    rows = shortlist(report)
    assert [c['lead_name'] for c in rows[:3]] == ['bb50-g1','bb25','bb54-g2']
    assert len({c['sequence'] for c in rows}) == len(rows)
    assert [c['evidence_rank'] for c in rows] == list(range(1,len(rows)+1))
    assert all(c['any_direct_readout'] and not c['exclusion_reason'] for c in rows)
    assert all(c['origins'] for c in rows)
    assert report['submission_ready'] is False
    assert report['shortfall_to_20'] == max(0,20-len(rows))


def test_duplicates_do_not_become_independent_samples(report):
    lead = shortlist(report)[0]
    p = next(p for p in lead['protocols'] if p['source'] == 'results/selectivity-panel/checkpoint.jsonl' and p['arm']=='g12d')
    assert p['records'] == 7
    assert p['unique_hashed_structures'] == 4
    assert p['unique_hashed_successes'] == 4
    assert len({json.dumps(p,sort_keys=True) for p in lead['protocols']}) == len(lead['protocols'])
    assert any(p['condition']=='matched_regional_anchor' for p in lead['protocols'])
    assert any(p['condition']=='recyc8' for p in lead['protocols'])


def test_retired_families_are_excluded_not_silently_reintroduced(report):
    for c in shortlist(report):
        assert not any('rigid' in n.lower() for n in c['aliases'])
        assert not any(n.startswith('B8-') and n!='B8-bb54-g2-parent' for n in c['aliases'])
        assert not any('bb135' in n for n in c['aliases'])


def test_confidence_changes_cannot_reorder_evidence(monkeypatch, report):
    original = list(pool.source_records(ROOT))
    modified = copy.deepcopy(original)
    for _,_,r in modified:
        r['row']['iptm'] = -1000
        r['row']['iptm_g12d'] = 1000
        r['row']['specificity_gap'] = 9999
    monkeypatch.setattr(pool,'source_records',lambda root: iter(modified))
    assert pool.build(ROOT)['shortlist_sha256'] == report['shortlist_sha256']


def test_reproducible_exports_and_edit_similarity(report,tmp_path):
    pool.write_report(report,tmp_path/'a')
    pool.write_report(report,tmp_path/'b')
    for p in (tmp_path/'a').iterdir():
        assert p.read_bytes() == (tmp_path/'b'/p.name).read_bytes()
    assert pool.edit_similarity('ACDE','ACDE') == 1
    assert pool.edit_similarity('ACDE','ACD') == .75
    assert pool.edit_similarity('ACD','ACDE') == .75


def test_corrected_inventory_uses_new_counts_and_preserves_original(report):
    audit_path=ROOT/'results/readout-correction/audit.json'
    matched_path=ROOT/'results/readout-correction/batch09_corrected_analysis.json'
    before=(ROOT/'results/submission-review/evidence_pool.json').read_bytes()
    corrected=pool.build(ROOT,audit_path,matched_path)
    lead=shortlist(corrected)[0]
    assert lead['matched']['arms']['g12d']['asp6_gate_passes']==2
    assert shortlist(report)[0]['matched']['arms']['g12d']['asp6_gate_passes']==3
    assert corrected['shortlist_sha256']==report['shortlist_sha256']
    assert corrected['readout_correction_applied'] is True
    assert corrected['unverified_records_remain'] is True
    assert any(p['readout_status']=='coordinate-remeasured' for p in lead['protocols'])
    assert (ROOT/'results/submission-review/evidence_pool.json').read_bytes()==before
