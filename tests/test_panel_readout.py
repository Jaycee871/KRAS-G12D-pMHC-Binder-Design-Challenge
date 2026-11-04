"""Chemical counterexample and real-coordinate regressions for panel readout."""
import sys
from pathlib import Path
import numpy as np
import pytest
from Bio.PDB import Atom, Residue, NeighborSearch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import selectivity_panel as panel
import complex_analysis as ca
import run_batch04 as original


def atom(name, xyz, element):
    return Atom.Atom(name, np.array(xyz,dtype=float), 0, 1, ' ', name, 1, element)


def synthetic(acceptor='OD1', distance=3.0, resname='ASP'):
    residue=Residue.Residue((' ',6,' '),resname,' ')
    residue.add(atom('CB',(0,0,0),'C'))
    residue.add(atom(acceptor,(10,0,0),'O'))
    donor=Residue.Residue((' ',1,' '),'ARG',' ')
    donor.add(atom('NH1',(10,distance,0),'N'))
    return residue, NeighborSearch(list(donor))


def test_carbon_proximity_does_not_pass_oxygen_gate():
    residue, _ = synthetic()
    donor=Residue.Residue((' ',1,' '),'ARG',' ')
    donor.add(atom('NH1',(0,3,0),'N'))
    assert panel.nearest_carboxylate_donor(residue,NeighborSearch(list(donor))) is None


@pytest.mark.parametrize('resname,name',[('ASP','OD1'),('ASP','OD2'),('GLU','OE1'),('GLU','OE2')])
def test_only_actual_carboxylate_oxygens_are_accepted(resname,name):
    residue, search=synthetic(acceptor=name,resname=resname)
    match=panel.nearest_carboxylate_donor(residue,search)
    assert match['acceptor_atom']==name
    assert match['distance']==3.0


def test_rounding_does_not_turn_3504_into_a_pass():
    residue, search=synthetic(distance=3.504)
    result=panel.nearest_carboxylate_donor(residue,search)
    assert round(result['distance'],2)==3.5
    assert result['distance']>original.GATE_DONOR_DISTANCE


@pytest.mark.parametrize('design,expected',[('bb50-g1',[True,False,True]),
    ('bb54-g2',[True,False,False]),('bb25',[True,False,True])])
def test_real_batch09_matches_original_two_leg_gate(design,expected):
    import json
    rows=[json.loads(l) for l in (ROOT/'results/batch-09/matched_anchor_checkpoint.jsonl').read_text().splitlines()]
    selected=[r for r in rows if r['design']==design and r['arm']=='g12d']
    actual=[]
    for row in selected:
        text=(ROOT/'results/batch-09'/row['structure_file']).read_text()
        corrected=panel.panel_metrics(text,row['sequence'],'g12d')
        reference=ca.interface_metrics(text,binder_sequence=row['sequence'],binder_chain_hint='D')
        gate=bool(original.gate_donor_distance(reference)['donor'] and reference['asp6_buried_frac']>0)
        assert corrected['asp6_gate_pass']==gate
        assert corrected['asp6_buried_frac']==pytest.approx(reference['asp6_buried_frac'])
        actual.append(corrected['asp6_gate_pass'])
    assert actual==expected


def test_burial_leg_is_required(monkeypatch):
    import json
    rows=[json.loads(l) for l in (ROOT/'results/batch-09/matched_anchor_checkpoint.jsonl').read_text().splitlines()]
    row=next(r for r in rows if r['design']=='bb50-g1' and r['arm']=='g12d' and r['draw']==0)
    text=(ROOT/'results/batch-09'/row['structure_file']).read_text()
    # Equal bound/free SASA retains the oxygen contact but removes burial.
    monkeypatch.setattr(ca,'_atom_sasa', lambda model: {(c.id,r.id[1],a.id):1.0
        for c in model for r in c for a in r})
    result=panel.panel_metrics(text,row['sequence'],'g12d')
    assert result['asp6_donor_distance_pass'] is True
    assert result['asp6_buried_frac']==0
    assert result['asp6_gate_pass'] is False
