"""Summarize already-completed fixed comparisons; no molecular calculation."""
import csv
import json
from pathlib import Path
import statistics
import sys
from collections import Counter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import read_json, record, write_new

W = ROOT / 'workspaces/nikasha_recovery_20260922'


def outcome_counts(rows, method):
    c = Counter(r['methods'][method]['outcome'] for r in rows)
    return {'total': len(rows), **{k: c[k] for k in ('correct', 'wrong', 'inconclusive', 'unavailable')}}


def compare(rows, a, b):
    common = [r for r in rows if all(r['methods'][m]['R'] is not None for m in (a, b))]
    changes = [{k: r[k] for k in ('case_id', 'expected_class', 'source_conditioning_metal')} |
               {m: r['methods'][m] for m in (a, b)} for r in common
               if r['methods'][a]['decision'] != r['methods'][b]['decision']]
    return {'declared': len(rows), 'common': len(common), a: outcome_counts(common, a), b: outcome_counts(common, b), 'transitions': changes}


def main():
    data = read_json(W / 'proposal_comparison.json'); selection = read_json(W / 'final_selection.json')
    triple = read_json(W / 'final_DFT_triples/result.json')
    account = list(csv.DictReader((W / 'scheduler_accounting.txt').open(), delimiter='|'))
    jobs = {r['JobID']: r for r in account if '.' not in r['JobID']}
    batches = {'proposal': ['1204185', '1204188', '1204189', '1204190', '1204191'],
               'DFT': ['1203976', '1203977', '1203978', '1203979'], 'CC': ['1202429']}
    costs = {name: {'jobs': [jobs[j] for j in ids], 'allocated_core_seconds': sum(int(jobs[j]['AllocCPUS']) * int(jobs[j]['ElapsedRaw']) for j in ids),
                    'allocated_GPU_seconds': int(jobs['1204185']['ElapsedRaw']) if name == 'proposal' else 0} for name, ids in batches.items()}
    paired = {}
    for z in ('all225', 'La', 'Ca'):
        rows = [r for r in data['rows'] if z == 'all225' or r['source_conditioning_metal'] == z]
        paired[z] = {a + '_vs_' + b: compare(rows, a, b) for a, b in [('DFT', 'context_composite'),
            ('core_native', 'context_composite'), ('context_composite', 'proposal_old'), ('context_composite', 'proposal_new')]}
    spread = {}
    for arm in ('La4', 'Ca5'):
        eligible = [p for p in data['pools'] if p['pool'] == arm and all(p['methods'][m]['R'] is not None for m in ('context_composite', 'proposal_new'))]
        pairs = [(p['root_case_id'], p['methods']['context_composite']['within_protein_range'], p['methods']['proposal_new']['within_protein_range']) for p in eligible]
        spread[arm] = {'common_protein_groups': len(pairs), 'old_median_range': statistics.median(p[1] for p in pairs),
                       'proposal_median_range': statistics.median(p[2] for p in pairs),
                       'groups_range_decreased': sum(p[2] < p[1] for p in pairs), 'groups_range_increased': sum(p[2] > p[1] for p in pairs),
                       'group_ranges': pairs}
    cc_path = ROOT / 'workspaces/electronic_accuracy_20260919/autocollect_v1/final/collection.json'
    cc = read_json(cc_path)
    summary = {'comparison': record(W / 'proposal_comparison.json'), 'final_selection': record(W / 'final_selection.json'),
        'DFT_triples': record(W / 'final_DFT_triples/result.json'), 'scheduler_accounting': record(W / 'scheduler_accounting.txt'),
        'counts': data['counts'], 'paired': paired, 'spread': spread, 'costs': costs,
        'independent_CC': {'collection': record(cc_path), 'complete_endpoints': cc['complete_endpoints'], 'total_endpoints': cc['total_endpoints'],
                          'complete_pairs': 0, 'interpretation': 'All La endpoints completed, all Ca endpoints error-terminated in ORCA MDCI; no correlated La/Ca classification comparison available.'},
        'scientific_calls_during_recovery': 0, 'thresholds_changed': False, 'production_changed': False}
    if (W / 'summary.json').exists():
        if read_json(W / 'summary.json') != summary: raise ValueError('saved scientific summary differs; use a new version')
    else:
        write_new(W / 'summary.json', summary)
    table = []
    names = {'DFT': 'Preserved DFT', 'core_native': 'Native MACE core', 'context_native': 'Native MACE context',
             'core_composite': 'Core composite', 'context_composite': 'Original context composite',
             'proposal_old': 'Proposal, old bands', 'proposal_new': 'Proposal, new fixed bands'}
    order = list(names)
    def counts_string(t): return '/'.join(str(t[k]) for k in ('correct', 'wrong', 'inconclusive', 'unavailable'))
    for method in order:
        table.append('| ' + names[method] + ' | ' + ' | '.join(counts_string(data['counts']['single'][z][method]) for z in ('all225', 'La', 'Ca')) + ' | ' + counts_string(data['counts']['all100_triples'][method]) + ' |')
    report = ['# Keep the original composite: donor proposals did not improve discrimination', '',
        '**The original context composite gives the most correct single-fold calls on this consumed reference challenge.** It gives 203 correct calls, two wrong calls and two inconclusives among 207 available sources. Native MACE core gives fewer decisive calls but zero wrong calls. Preserved DFT, native core and the original context composite tie on all 94 available three-fold summaries. The new accommodation proposal method offers no overall gain, and its separate canonical-only bands worsen transfer.', '',
        'This is a structural robustness test on 25 already-consumed protein groups, not 225 independent biological controls or prospective validation. The 25 canonical geometries remain calibration-only. All 225 noncanonical sources remain in denominators; no thresholds, samples or production settings changed during recovery.', '',
        '## Completed fixed comparisons', '',
        'Cells are **correct / wrong / inconclusive / unavailable**. The La-conditioned arm is the main scanner use case.', '',
        '| Method | All225 | La100 | Ca125 | All100 La triples |',
        '|---|---:|---:|---:|---:|', *table, '',
        'The four triples per protein are correlated subsets, with every member required. DFT, native core and original context composite each have 23/25 proteins with all four triples correct; two groups contain unavailable preparations.', '',
        '| Method | Strict La4 (25) | Strict Ca5 (25) | Equal arm medians (25) |',
        '|---|---:|---:|---:|']
    for method in order:
        report.append('| ' + names[method] + ' | ' + ' | '.join(counts_string(data['counts']['pools'][p][method]) for p in ('La4', 'Ca5', 'balanced')) + ' |')
    report += ['', '## Common coverage and actual transitions', '']
    for arm in ('all225', 'La', 'Ca'):
        p = paired[arm]['DFT_vs_context_composite']
        report.append(f"- {arm}: on the same {p['common']} available sources, DFT {counts_string(p['DFT'])}; original context composite {counts_string(p['context_composite'])}.")
    report += ['', 'On La-conditioned sources, the composite makes one additional correct call and removes the one DFT wrong call. It does not establish broader La/Ca affinity accuracy. Raw per-case transitions and all matched subsets are in summary.json and proposal_comparison.json.', '',
        '| Case | Condition | DFT → original context composite |', '|---|---|---|']
    for c in paired['all225']['DFT_vs_context_composite']['transitions']:
        report.append(f"| {c['case_id']} | {c['source_conditioning_metal']} | {c['DFT']['decision']} → {c['context_composite']['decision']} |")
    report += ['', '## Accommodation utility', '',
        'All 415 eligible optimizers succeeded; 2,046 native MACE energy/force calls completed. There were 822 fresh proposal GFN2 attempts: 820 complete and two failed. Eight identical proposal=q0 endpoints reused the exact prior result. Selection chose 345 proposals and 68 origins; three metal endpoints remained unavailable. No failed proposal fell back to a baseline success.', '',
        'The proposal protocol has 205 available site pairs. On those same205 sources, the original composite is201 correct/2 wrong/2 inconclusive. Old-band proposals retain201 correct, defer one wrong call to inconclusive (201/1/3), and add two failures outside common coverage. The independently frozen new bands instead yield197/2/6: exactly four formerly correct calls become inconclusive, with no wrong call recovered. Both proposal interpretations retain20 unavailable cases. Changing calibration is therefore material and was not a solution to transfer accuracy.', '']
    for arm, s in spread.items():
        report.append(f"- {arm}, same {s['common_protein_groups']} complete protein groups: median within-protein raw-score range {s['old_median_range']:.6f} → {s['proposal_median_range']:.6f} model kcal/mol; range decreased in {s['groups_range_decreased']} groups and increased in {s['groups_range_increased']}.")
    report += ['', 'The exact missing proposal pairs are the pre-existing A8R3S4 Ca-conditioned sample3 solvent failure, plus new failures for P12293 Ca-conditioned sample4 (La vacuum) and Q60AR6 Ca-conditioned sample1 (La ALPB). The original 17 preparation failures remain separate. Lower energy or narrower geometric spread does not by itself establish better classification.', '',
        '## Cost and independent electronic reference', '',
        f"- Preserved DFT: 416/416 new endpoints completed; **{costs['DFT']['allocated_core_seconds']:,} allocated core-seconds**, no GPU. Archived25 pairs reused.",
        f"- Proposal continuation: **{costs['proposal']['allocated_core_seconds']:,} allocated core-seconds and {costs['proposal']['allocated_GPU_seconds']:,} GPU-seconds**, including the two failed solvent attempts and allocation/collection overhead. This is incremental to archived q0/context preparation and not full source-to-score cost or a matched-hardware speed comparison.",
        f"- Independent CC job1202429: **{costs['CC']['allocated_core_seconds']:,} allocated core-seconds**, no GPU; only3/6 endpoints complete. All three Ca outputs show MDCI error termination and MPI segmentation faults. No complete Ca/La pair exists, so it cannot adjudicate whether the electronic treatment improves classification. It is too costly for routine use. The failure cause was not isolated by this recovery.",
        '- Original preparation and parser/collection CPU were not separately metered. Recovery itself made zero molecular calls; scheduler costs above are actual Sept20–22 receipts, not estimates.', '',
        '## Recommendation and reusable artifacts', '',
        '**Retain the original context composite for this development comparison; retain the production baseline/default.** Reject promotion of the current proposal-selection method. Three-fold aggregation is robust on the available consumed references, with no gain over preserved DFT. The next experiment must provide a transferable accuracy benefit rather than a favorable energy shift or recalibrated canonical separation.', '',
        'Final reusable proposal selection: `workspaces/nikasha_recovery_20260922/final_selection.json`. Full fixed join: `proposal_comparison.json`; full DFT/MACE join: `final_DFT_MACE_comparison.json`; exact100-triple DFT comparison: `final_DFT_triples/result.json`; compact machine-readable summary and transitions: `summary.json`. All are under that recovery workspace and pin the original immutable runs. Parent separately owns the shared-candidate-pool experiment; no new chemistry belongs to this recovery.', '']
    (ROOT / 'diagnostics/nikasha_recovery_20260922/REPORT.md').write_text('\n'.join(report))
    print(json.dumps({'counts': data['counts'], 'spread': spread, 'costs': {k: {a: v[a] for a in ('allocated_core_seconds','allocated_GPU_seconds')} for k,v in costs.items()}}))


if __name__ == '__main__': main()
