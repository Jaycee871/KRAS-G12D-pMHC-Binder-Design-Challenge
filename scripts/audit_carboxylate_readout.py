#!/usr/bin/env python3
"""Offline, additive remeasurement of checkpoints using the corrected chemistry.

No historical checkpoint or model response is overwritten. Missing structures
remain unverified. Different sources, conditions and arms are never pooled.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import itertools
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'),str(ROOT/'scripts')]
import selectivity_panel as panel
import replicate_rmsd_batch04 as rmsd
import selectivity_panel_analyze as frame_analysis


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources(root):
    return [root/'results/selectivity-panel/checkpoint.jsonl'] + [
        root/f'results/batch-{n:02}/panel_checkpoint.jsonl' for n in (6,7,8)] + [
        root/'results/batch-09/matched_anchor_checkpoint.jsonl']


def audit(root):
    report={'readout_version':panel.READOUT_VERSION,'model_calls':0,
            'historical_inputs_unchanged':True,'source_sha256':{},
            'groups':{},'remeasured':[],'unavailable':[]}
    for path in sources(root):
        relative=str(path.relative_to(root))
        report['source_sha256'][relative]=sha(path)
        for line_no,line in enumerate(path.read_text().splitlines(),1):
            if not line.strip(): continue
            r=json.loads(line)
            entry={'source':relative,'line':line_no,'design':r['design'],
                   'arm':r['arm'],'draw':r['draw'],'condition':r.get('condition','matched_regional_anchor')}
            coord=(path.parent/r['structure_file']).resolve()
            if not coord.is_relative_to(root.resolve()):
                raise ValueError('coordinate path escapes repository')
            group_key=json.dumps([relative,r['design'],r['arm'],entry['condition']])
            group=report['groups'].setdefault(group_key,{'source':relative,'design':r['design'],
                'arm':r['arm'],'condition':entry['condition'],'records':0,'remeasured':0,
                'unavailable':0,'changed_flags':0})
            group['records']+=1
            if not coord.is_file():
                report['unavailable'].append({**entry,'structure_file':r['structure_file'],
                    'recorded_structure_sha1':r.get('structure_sha1'), 'reason':'coordinate absent from clone'})
                group['unavailable']+=1
                continue
            recorded_hash=r.get('structure_sha1')
            if recorded_hash and hashlib.sha1(coord.read_bytes()).hexdigest()!=recorded_hash:
                raise ValueError(f'coordinate hash mismatch: {coord}')
            sequence=r.get('sequence') or r['row']['sequence']
            measured=panel.panel_metrics(coord.read_text(),sequence,r['arm'])
            keys=('readout_version','pos6_min_donor_distance','pos6_donor_partner',
                  'asp6_buried_frac','asp6_donor_distance_pass','asp6_gate_pass','donor_contacts')
            original=r['row'].get('asp6_gate_pass')
            corrected=measured['asp6_gate_pass']
            changed=original is not None and str(original).lower()!=str(corrected).lower()
            report['remeasured'].append({**entry,'coordinate':str(coord.relative_to(root.resolve())),
                'coordinate_sha256':sha(coord),'original_hash_verified':bool(recorded_hash),
                'original_asp6_gate_pass':original,'gate_changed':changed,
                'corrected_metrics':{k:measured.get(k) for k in keys}})
            group['remeasured']+=1
            group['changed_flags']+=int(changed)
    report['coverage']={'checkpoint_records':sum(g['records'] for g in report['groups'].values()),
        'remeasured':len(report['remeasured']),'unavailable':len(report['unavailable']),
        'changed_flags':sum(x['gate_changed'] for x in report['remeasured'])}
    # Include analyzer and chemistry hashes, not just input-data hashes.
    report['analysis_sha256']={str(p.relative_to(root)):sha(p) for p in (
        root/'scripts/selectivity_panel.py',root/'scripts/audit_carboxylate_readout.py',
        root/'src/complex_analysis.py')}
    return report


def corrected_batch09(root, report):
    path=root/'results/batch-09/matched_anchor_analysis.json'
    result=copy.deepcopy(json.loads(path.read_text()))
    result['readout_version']=panel.READOUT_VERSION
    result['supersedes_recorded_readout_counts']=True
    result['historical_analysis_sha256']=sha(path)
    rows=[json.loads(l) for l in (root/'results/batch-09/matched_anchor_checkpoint.jsonl').read_text().splitlines()]
    corrected=[r for r in report['remeasured'] if r['source']=='results/batch-09/matched_anchor_checkpoint.jsonl']
    if len(corrected)!=36: raise ValueError('full batch-09 correction needs all 36 structures')
    for design, info in result['designs'].items():
        for arm, summary in info['arms'].items():
            subset=[r for r in corrected if r['design']==design and r['arm']==arm]
            summary['asp6_gate_passes']=sum(r['corrected_metrics']['asp6_gate_pass'] for r in subset)
            values=[r['corrected_metrics']['asp6_buried_frac'] for r in subset
                    if r['corrected_metrics']['asp6_buried_frac'] is not None]
            summary['median_asp6_buried_fraction']=statistics.median(values) if values else None
            originals=[r for r in rows if r['design']==design and r['arm']==arm]
            models=[rmsd.model_of(root/'results/batch-09'/r['structure_file']) for r in originals]
            binder=frame_analysis.binder_chain_of(models[0],len(originals[0]['sequence']))
            frame=['A'] if arm=='apo' else ['A','C']
            pairs=[]
            for (i,a),(j,b) in itertools.combinations(enumerate(models),2):
                for chain in frame+[binder]:
                    if len(rmsd.ca_list(a,chain))!=len(rmsd.ca_list(b,chain)):
                        raise ValueError('unequal coordinate coverage in RMSD comparison')
                value,floor=frame_analysis.binder_rmsd_frame(a,b,frame,binder,binder)
                pairs.append({'draw_pair':[originals[i]['draw'],originals[j]['draw']],
                    'binder_ca_rmsd':value,'frame_ca_rmsd':floor})
            summary['pose']={'fit_chains':frame,'binder_chain':binder,'n':len(models),
                'pairwise':pairs,'mean':statistics.mean(p['binder_ca_rmsd'] for p in pairs),
                'max':max(p['binder_ca_rmsd'] for p in pairs)}
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'results/readout-correction')
    args=parser.parse_args(argv)
    args.out.mkdir(parents=True,exist_ok=True)
    report=audit(ROOT)
    corrected=corrected_batch09(ROOT,report)
    (args.out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    (args.out/'batch09_corrected_analysis.json').write_text(json.dumps(corrected,indent=2)+'\n')
    print(json.dumps(report['coverage']))
    for name,info in corrected['designs'].items():
        g=info['arms']['g12d']
        print(name,'full gate',g['asp6_gate_passes'],'/ 3; RMSD mean/max',
              round(g['pose']['mean'],2),round(g['pose']['max'],2))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
