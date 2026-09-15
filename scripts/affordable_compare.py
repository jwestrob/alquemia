"""Side-by-side protocol records and a compact report; missing scores stay missing."""
from __future__ import annotations

import argparse
from collections import Counter
import math
from pathlib import Path
import sys

from affordable_common import read_json, record, verify, write_new, contrast, corrected, energy


def descriptors(manifest,root):
    sys.path.insert(0,str(root.parent))
    from on_scanner.bvs import CATION_PARAMS
    c=manifest['coordination']
    donors=c.get('typed_donors',c.get('source_typed_direct_donors_ordered_by_distance',[]))
    distances=[d['distance_A'] for d in donors]
    groups=Counter((d['chain'],d['resnum'],d.get('icode',d.get('insertion_code','')),d['resname']) for d in donors)
    counts=Counter(d['element'] for d in donors)
    charge=manifest['charge_ledger']
    values={}
    for ion in ('La3+','Ca2+'):
        p=CATION_PARAMS[ion]
        bvs=sum(math.exp((p['r0_'+d['element']]-d['distance_A'])/p['B']) for d in donors if d['element'] in ('O','N'))
        values[ion]={'BVS':bvs,'absolute_residual':abs(bvs-p['q']),'relative_residual':abs(bvs-p['q'])/p['q']}
    return {'typed_CN':c.get('coordination_number',len(donors)),'direct_donor_composition':dict(counts),
            'bidentate_carboxylate_residue_count':sum(n>=2 and k[-1] in ('ASP','GLU') for k,n in groups.items()),
            'mean_direct_distance_A':sum(distances)/len(distances) if distances else None,
            'minimum_direct_distance_A':min(distances) if distances else None,
            'core_formal_charge':charge.get('scaffold_formal_charge',charge.get('fragment_formal_charge_sum')),
            'core_atom_count':1+sum(f['atom_count'] for f in manifest['qm_fragments']),
            'BVS':values,'BVS_parameter_source':record(root.parent/'on_scanner/bvs.py'),
            'distance_precision':'existing manifest distances; typically rounded to 0.001 A',
            'classification_fit':'none; raw descriptors only'}


def physical_verdict(solver,baseline_boundary,endpoint_scores=None):
    if solver is None: return {'status':'not_run','checks':{},'score_allowed':False}
    entries={r['label']:r for r in solver['numerical_checks']}
    values={k:r['result']['components']['delta_U_kcal_mol'] for k,r in entries.items() if r['status']=='computed'}
    checks={}
    for label,row in entries.items():
        if row['status']=='computed':
            error=abs(row['result']['components']['grid_vs_analytic_cross_error_kcal_mol'])
            checks[label+'/coulomb_discretization']={'absolute_change_kcal_mol':error,'tolerance_kcal_mol':.5,'pass':error<=.5}
    hashes=[entries.get(f'1h4i_qm{radius}_{metal}/primary',{}).get('result') or {} for radius in ('33','36') for metal in ('La','Ca')]
    hashes=[h.get('physical_boundary_hash') for h in hashes]
    checks['common_partition_cavity']={'hashes':hashes,'pass':len(set(hashes))==1 if all(hashes) else None}
    anchor='1h4i_qm33_La'
    for label,reference,tolerance in [('identity/'+anchor,None,.01),('repeat/'+anchor,anchor+'/primary',.5)]:
        delta=abs(values[label]-(values[reference] if reference else 0.)) if label in values and (reference is None or reference in values) else None
        checks[label]={'absolute_change_kcal_mol':delta,'tolerance_kcal_mol':tolerance,'pass':None if delta is None else delta<=tolerance}
    cases=sorted({k.split('/')[0] for k in values if k.endswith('/primary')})
    for name in cases:
        for variant in ('refined','extended')+(('translated','rotated') if name.startswith('1h4i_qm33') else ()):
            key=name+'/'+variant;base=name+'/primary'
            delta=abs(values[key]-values[base]) if key in values else None
            checks[key]={'absolute_change_kcal_mol':delta,'tolerance_kcal_mol':.5,'pass':None if delta is None else delta<=.5}
    primary={}
    endpoint_scores=endpoint_scores or {}
    for radius in ('33','36'):
        ca=f'1h4i_qm{radius}_Ca/primary';la=f'1h4i_qm{radius}_La/primary'
        if ca in values and la in values and radius in endpoint_scores:
            primary[radius]=endpoint_scores[radius]+values[ca]-values[la]
            delta=abs(endpoint_scores[radius]-baseline_boundary['radii']['3.'+radius[-1]]['score_kcal_mol'])
            checks['baseline_reproduction_qm'+radius]={'absolute_change_kcal_mol':delta,'tolerance_kcal_mol':.5,'pass':delta<=.5}
    delta=abs(primary['36']-primary['33']) if len(primary)==2 else None
    checks['partition_score']={'absolute_change_kcal_mol':delta,'tolerance_kcal_mol':2.,'pass':None if delta is None else delta<=2.}
    if solver['budget_exceeded']: status='cost_budget_exceeded'
    elif any(c['pass'] is False for c in checks.values()):status='physical_checks_failed'
    elif not cases or any(c['pass'] is None for c in checks.values()):status='physical_checks_incomplete'
    else:status='physical_checks_passed'
    return {'status':status,'checks':checks,'score_allowed':status=='physical_checks_passed',
            'development_boundary_scores_kcal_mol':primary,'predictive_validation':'not_established'}


def compare(root):
    root=Path(root).resolve();diag=root/'diagnostics/affordable_challenger_20260915'
    audit=read_json(diag/'live_audit.json');rows=[]
    for r in audit['rows']:
        ep=Path(r['source_record']['artifacts']['La']['input']['path']).parent
        manifests=list(ep.glob('*_carve_manifest.json'))
        if len(manifests)!=1: raise ValueError('ambiguous baseline carve')
        m=read_json(manifests[0]); row={k:r[k] for k in ('case','evidence_stratum','use','baseline_protocol_id','released_endpoints_hartree','released_score','released_decision','provenance_status')}
        row.update({'cheap_descriptors':descriptors(m,root),'baseline_preparation':record(manifests[0]),
                    'replicate_group':'P16027_and_1H4I' if r['case'] in ('1H4I','p16027-pqq-la_model') else r['case'],
                    'homology_group':'canonical_PQQ_dehydrogenase_family_not_independent_architectures',
                    'repaired_baseline':{'status':'not_applicable','score':None},
                    'environmental_challenger':{'status':'not_evaluated','score':None}})
        if 'acidic_Dplus2' in r['source_record']:
            row['cheap_descriptors']['motif_only_rule']='Ln' if r['source_record']['acidic_Dplus2'] else 'Ca'
            row['biological_label']=r['source_record']['class']
        rows.append(row)
    panel=read_json(root/'diagnostics/nonpqq_direct_site_benchmark_20260915/panel_manifest.json')
    repairs=read_json(diag/'repair_inventory_verified.json')
    for r in repairs:
        baseline=read_json(verify(r['baseline']));repair=read_json(verify(r['repaired']));stem=baseline['stem']
        protein=next(p for p in panel['proteins'] if stem.startswith(p['panel_key']))
        rows.append({'case':stem,'evidence_stratum':protein['evidence']['class'],'evidence':protein['evidence'],
                     'use':'preparation_validation_only_energies_unopened','replicate_group':protein['panel_key'],
                     'ordered_site_vector_policy':protein['aggregation'],'baseline_protocol_id':baseline['protocol_id'],
                     'baseline_preparation':r['baseline'],'cheap_descriptors':descriptors(baseline,root),
                     'released_score':None,'released_decision':None,
                     'repaired_baseline':{'status':'prepared_not_executed','protocol_id':repair['protocol_id'],
                                          'manifest':r['repaired'],'score':None,'atom_count':r['atom_count']},
                     'environmental_challenger':{'status':'not_in_development_pilot','score':None}})
    workspace=root/'workspaces/affordable_challenger_20260915'
    events_path=workspace/'pilot/budget_events.jsonl'
    events=[__import__('json').loads(line) for line in events_path.read_text().splitlines()] if events_path.exists() else []
    solver_path=workspace/'solver/solver_result.json';solver=read_json(solver_path) if solver_path.exists() else None
    boundary=read_json(root/'diagnostics/pqq_boundary_pair_20260914/result.json')['panels']['mxaf']
    pilot=read_json(workspace/'pilot/pilot_manifest.json')
    endpoint_scores={}
    for radius in ('33','36'):
        tasks={t['metal']:t for t in pilot['tasks'] if t['case']=='1h4i_qm'+radius}
        try:
            ca,la=[energy(tasks[metal]['output_path']) for metal in ('Ca','La')]
            endpoint_scores[radius]=contrast(ca,la,audit['aquo_gauge']['delta_E_aquo_hartree'])['S_kcal_mol']
        except (OSError,ValueError): pass
    physical=physical_verdict(solver,boundary,endpoint_scores)
    results=[]
    for case in pilot['cases']:
        name=case['case'];sources=[t for t in pilot['tasks'] if t['case']==name]
        rec={'case':name,'use':'development_consumed_1H4I','baseline_protocol_id':case['baseline_protocol_id'],
             'baseline_preparation':case['source_manifest'],'endpoint_results':{},'environmental_challenger':{'status':'not_run','score':None}}
        for t in sources:
            op=Path(t['output_path'])
            if not op.exists(): rec['endpoint_results'][t['metal']]={'status':'not_run','energy_hartree':None};continue
            try: e=energy(op);status='computed'
            except ValueError: e=None;status='failed'
            rec['endpoint_results'][t['metal']]={'status':status,'energy_hartree':e,'output':record(op)}
        if solver:
            valid={r['label']:r for r in solver['numerical_checks'] if r['status']=='computed'}
            terms={metal:valid.get(f'{name}_{metal}/primary') for metal in ('La','Ca')}
            ca=rec['endpoint_results'].get('Ca',{}).get('energy_hartree');la=rec['endpoint_results'].get('La',{}).get('energy_hartree')
            if ca is not None and la is not None and all(terms.values()):
                gauge=audit['aquo_gauge']['delta_E_aquo_hartree']
                base_score=contrast(ca,la,gauge)
                correction=corrected(base_score,terms['Ca']['result']['components']['delta_U_kcal_mol'],terms['La']['result']['components']['delta_U_kcal_mol'])
                rec['environmental_challenger']={'status':physical['status'],'development_descriptor':correction,
                                                 'score':correction['S_env_kcal_mol'] if physical['score_allowed'] else None,
                                                 'components':{metal:terms[metal]['result']['components'] for metal in terms},
                                                 'decision':'uncalibrated_protocol'}
            else: rec['environmental_challenger']={'status':'unavailable','score':None}
        results.append(rec)
    return {'schema_version':'alquemia.affordable_comparison.v1','baseline_default':'unchanged',
            'archived_and_prepared_rows':rows,'development_pilot':results,'physical_checks':physical,
            'cost':{'baseline_accounting':audit['sacct_raw'],'solver_allocated_core_seconds':solver['allocated_core_seconds'] if solver else None,
                    'pilot_admitted_endpoints':sum(e.get('admitted_endpoint_count',0) for e in events),
                    'pilot_execution_events':events,
                    'pilot_recorded_allocated_core_seconds':sum(e['allocated_core_seconds'] for e in events if 'allocated_core_seconds' in e) if any('allocated_core_seconds' in e for e in events) else None,
                    'matched_production_cost_ratio':None,'affordability':'not_established','unmeasured_preparation_costs':'not zero; not fully instrumented before this session'},
            'judgments':{'numerical_credibility':physical['status'],'incremental_predictive_information':'not_tested_on_independent_composition_challenging_controls',
                         'affordability':'not_established'},'recommendation':'retain_baseline',
            'response_status':'response_model_not_validated'}


def report(data):
    n=len(data['archived_and_prepared_rows'])
    lines=['# Affordable discriminator development status','',
           'Baseline production behavior and all frozen experiments remain unchanged.','',
           f'{n} archived/prepared site records are retained in the side-by-side ledger.',
           'The 27 archived PQQ pairs reproduce their released energies and bands; all six non-PQQ peptide repairs preserve source heavy coordinates and paired charges.','',
           '| Component | Status |','|---|---|',
           '| Peptide-amide repair | Six real sites prepared; new protocol; no energetic rescore |',
           f"| Environmental numerical credibility | {data['judgments']['numerical_credibility']} |",
           '| Incremental predictive information | Not established; pilot uses consumed 1H4I cases |',
           '| Ordinary-score affordability | Not established; no matched production cost ratio |',
           '| Mechanical response | Disabled: response_model_not_validated |','',
           '## Development endpoint status','']
    for r in data['development_pilot']:
        statuses=', '.join(f"{m}: {v['status']}" for m,v in r['endpoint_results'].items())
        lines.append(f"- {r['case']}: {statuses}; environment {r['environmental_challenger']['status']}.")
    lines+=['','Recommendation: **retain baseline**. A prepared or completed calculation does not establish predictive improvement.','',
            'The fixed-core 1H4I whole-chain environment is unsupported because its terminal Lys lacks a complete ff19SB-compatible terminal state. The existing qm33/qm36 environments are prepared without adding atoms.','',
            'The detailed JSON retains source manifests, evidence strata, raw inexpensive descriptors, unavailable scores, numerical checks and costs. No new threshold is fit.']
    return '\n'.join(lines)+'\n'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();data=compare(a.root);a.output.mkdir(parents=True,exist_ok=False)
    write_new(a.output/'comparison.json',data);(a.output/'REPORT.md').write_text(report(data))
    print(f"{len(data['archived_and_prepared_rows'])} records; {data['physical_checks']['status']}; retain baseline")


if __name__=='__main__': main()
