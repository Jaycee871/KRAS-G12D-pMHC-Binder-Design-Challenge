#!/usr/bin/env python3
"""Build a review pool from recorded results; never generate or submit sequences."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def truth(value):
    return str(value).lower() == 'true'


def edit_similarity(a, b):
    """1 - Levenshtein distance / longer length; not structural similarity."""
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (x != y)))
        previous = current
    return 1 - previous[-1] / max(len(a), len(b), 1)


def source_records(root):
    paths = []
    for batch in range(4, 10):
        paths.extend(sorted((root / 'results' / f'batch-{batch:02}').glob('*checkpoint.jsonl')))
    for name in ('sampling-study', 'selectivity-panel'):
        paths.append(root / 'results' / name / 'checkpoint.jsonl')
    for path in paths:
        for line_no, line in enumerate(path.read_text().splitlines(), 1):
            if line.strip():
                yield path, line_no, json.loads(line)
    # These original repeat runs are not present in the JSONL checkpoints.
    for batch in (4, 5):
        path = root / 'results' / f'batch-{batch:02}' / 'replicates.json'
        for i, group in enumerate(json.loads(path.read_text())):
            for j, row in enumerate(group['runs']):
                yield path, f'[{i}].runs[{j}]', {
                    'row': row, 'sequence': group['sequence'],
                    'design': group['backbone'] + '_' + group['variant'],
                    'condition': 'original_replicates_settings_in_script',
                }


def exclusion(candidate):
    names = candidate['aliases']
    if any('bb135' in n and ('B5-' in n or n == 'bb135-g1') for n in names):
        return 'bb135 family deprioritized: measured decoy-Asp8 promiscuity in g1; other variants unvalidated'
    if any('rigid' in n.lower() for n in names):
        return 'batch-09 audit: rigidified family retired after replication failure'
    if any(n.startswith('B8-') and n != 'B8-bb54-g2-parent' for n in names):
        return 'batch-09 audit: batch-08 shape/de-wetting family retired'
    if not candidate['any_direct_readout']:
        return 'no recorded successful direct Asp6 readout'
    return ''


def build(root):
    candidates, sources = {}, {}
    for path, location, record in source_records(root):
        row = record.get('row') or {}
        sequence = record.get('sequence') or row.get('sequence')
        if not sequence:
            raise ValueError(f'missing sequence: {path}:{location}')
        if not 10 <= len(sequence) <= 250 or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError(f'invalid sequence: {path}:{location}')
        sha = hashlib.sha256(sequence.encode()).hexdigest()
        relative = str(path.relative_to(root))
        sources[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
        name = record.get('design') or row.get('backbone', '') + '_' + row.get('variant', '')
        c = candidates.setdefault(sha, {'sequence': sequence, 'sequence_sha256': sha,
            'aliases': set(), 'backbones': set(), 'protocols': {}, 'origins': [],
            'any_direct_readout': False, 'matched': None})
        c['aliases'].add(name)
        backbone = record.get('backbone') or row.get('backbone')
        if backbone:
            c['backbones'].add(backbone)
        arm = record.get('arm', 'g12d')
        protocol = {'source': relative, 'arm': arm,
                    'condition': record.get('condition', 'matched_regional_anchor' if 'batch-09' in relative else 'initial_screen'),
                    'sampling': record.get('sampling', record.get('boltz_kwargs')),
                    'constraint': record.get('constraint')}
        protocol_id = json.dumps(protocol, sort_keys=True)
        g = c['protocols'].setdefault(protocol_id, {**protocol, 'records': 0,
            'scoring_errors': 0, 'success_records': 0, 'structures': {},
            'unhashed_records': 0})
        g['records'] += 1
        c['origins'].append({'source': relative, 'location': location})
        if row.get('analysis_error'):
            g['scoring_errors'] += 1
            continue
        passed = arm == 'g12d' and truth(row.get('asp6_gate_pass', row.get('gate_pass')))
        g['success_records'] += int(passed)
        c['any_direct_readout'] |= passed
        digest = record.get('structure_sha1')
        coordinate = record.get('structure_file')
        if not digest and coordinate and (path.parent / coordinate).is_file():
            digest = hashlib.sha1((path.parent / coordinate).read_bytes()).hexdigest()
        if digest:
            previous = g['structures'].get(digest)
            if previous is not None and previous != passed:
                raise ValueError(f'conflicting readout for identical structure: {digest}')
            g['structures'][digest] = passed
        else:
            g['unhashed_records'] += 1
    matched_path = root / 'results/batch-09/matched_anchor_analysis.json'
    sources[str(matched_path.relative_to(root))] = hashlib.sha256(matched_path.read_bytes()).hexdigest()
    matched = json.loads(matched_path.read_text())['designs']
    for c in candidates.values():
        c['aliases'] = sorted(c['aliases'])
        c['backbones'] = sorted(c['backbones'])
        protocols = []
        for p in c['protocols'].values():
            structures = p.pop('structures')
            p['unique_hashed_structures'] = len(structures)
            p['unique_hashed_successes'] = sum(structures.values())
            protocols.append(p)
        c['protocols'] = protocols
        c['exclusion_reason'] = exclusion(c)
        c['family'] = next(iter(c['backbones']), c['aliases'][0])
        if any('bb54' in n for n in c['aliases']):
            c['family'] = 'bb54_parent_and_derivatives'
        for name in ('bb50-g1', 'bb25', 'bb54-g2'):
            if name in c['aliases']:
                c['matched'] = matched[name]
                c['lead_name'] = name
        c['review_warnings'] = ['selectivity_unvalidated', 'competition_novelty_unverified']
        arms = {p['arm'] for p in protocols}
        for arm in ('wt','apo','decoy'):
            if arm not in arms:
                c['review_warnings'].append(arm + '_untested_in_inventory')
        if c['matched']:
            c['review_warnings'].append('matched_apo_overburial_unresolved')
        else:
            c['review_warnings'].append('matched_counter_screen_not_run')
        repeats = [p for p in protocols if p['arm']=='g12d' and p['records']>=3]
        if repeats and all(p['success_records']==0 for p in repeats):
            c['review_warnings'].append('initial_hit_not_reproduced_in_recorded_repeats')
        c['tier'] = 'excluded'
        if not c['exclusion_reason']:
            c['tier'] = 'multi_sample_mixed' if any(
                p['arm'] == 'g12d' and p['unique_hashed_successes'] >= 2
                for p in protocols) else 'limited_readout'
            if c.get('lead_name'):
                c['tier'] = 'matched_lead'
    chosen = []
    for name in ('bb50-g1', 'bb25', 'bb54-g2'):
        chosen.append(next(c for c in candidates.values() if c.get('lead_name') == name))
    for tier in ('multi_sample_mixed', 'limited_readout'):
        pending = [c for c in candidates.values() if c['tier'] == tier]
        while pending:
            families = {c['family'] for c in chosen}
            def key(c):
                similarity = max((edit_similarity(c['sequence'], s['sequence']) for s in chosen), default=0)
                return c['family'] in families, similarity, c['sequence_sha256']
            pick = min(pending, key=key)
            chosen.append(pick)
            pending.remove(pick)
    for rank, c in enumerate(chosen, 1):
        c['evidence_rank'] = rank
        c['name'] = f'KRAS_REVIEW_{rank:02}_{c["sequence_sha256"][:8]}'
        c['max_edit_similarity_to_earlier'] = round(max((edit_similarity(c['sequence'], s['sequence'])
            for s in chosen[:rank-1]), default=0), 4)
    return {'submission_ready': False, 'validated_selective_candidates': 0,
        'track1_minimum_in_project_plan': 20, 'shortlist_count': len(chosen),
        'shortfall_to_20': max(0, 20-len(chosen)), 'model_calls': 0,
        'scope': 'recorded downstream batches 04-09, sampling study and selectivity panel; not batches 01-03',
        'ranking_policy': 'Fixed matched leads bb50-g1, bb25, bb54-g2; then multi-sample mixed and limited-readout tiers. Within a tier, prefer an unused generation-backbone family, then lower maximum edit similarity to earlier entries, then SHA256. This is a review heuristic, not a calibrated biological score.',
        'source_sha256': sources, 'candidates': sorted(candidates.values(), key=lambda c: c['sequence_sha256']),
        'shortlist_sha256': [c['sequence_sha256'] for c in chosen]}


def write_report(report, out):
    out.mkdir(parents=True, exist_ok=True)
    (out / 'evidence_pool.json').write_text(json.dumps(report, indent=2) + '\n')
    indexed = {c['sequence_sha256']: c for c in report['candidates']}
    chosen = [indexed[s] for s in report['shortlist_sha256']]
    fields = ['evidence_rank', 'name', 'sequence', 'molecule_class', 'tier', 'family',
              'max_edit_similarity_to_earlier', 'sequence_sha256','aliases','review_warnings']
    with (out / 'review_shortlist.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n'); w.writeheader()
        for c in chosen:
            w.writerow({**{k: ';'.join(c[k]) if isinstance(c[k],list) else c[k] for k in fields if k != 'molecule_class'}, 'molecule_class': 'protein'})
    with (out / 'candidate_audit.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['sequence_sha256','sequence','tier','aliases','exclusion_reason'], lineterminator='\n')
        w.writeheader()
        for c in report['candidates']:
            w.writerow({k: ';'.join(c[k]) if k == 'aliases' else c[k] for k in w.fieldnames})
    with (out / 'format_preview_only.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['name','sequence','molecule_class'], lineterminator='\n')
        w.writeheader()
        for c in chosen:
            w.writerow({'name': c['name'], 'sequence': c['sequence'], 'molecule_class': 'protein'})
    return chosen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/submission-review')
    args = parser.parse_args()
    report = build(ROOT)
    chosen = write_report(report, args.out)
    print(json.dumps({k: report[k] for k in ('submission_ready','shortlist_count','shortfall_to_20','model_calls')}))
    for c in chosen:
        print(c['evidence_rank'], c.get('lead_name', c['aliases'][0]), c['tier'])


if __name__ == '__main__':
    main()
