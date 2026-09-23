"""Read-only two-source geometry-by-rank inventory; no molecular execution."""
from pathlib import Path
import argparse,hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,xyz


def geometry_key(task):
    return (tuple(tuple(a) for a in xyz(verify(task['xyz']))),task['charge'],task['multiplicity'])


def inventory(comparison,output):
    result=read_json(comparison)
    q4=next(r for r in result['rows'] if r['case_id']=='q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-2')
    eligible=[r for r in result['rows'] if r['precision_numerical_comparison'] and r['precision_numerical_comparison']['delta_R']['operational'] is not None]
    maximum=max(eligible,key=lambda r:abs(r['precision_numerical_comparison']['delta_R']['operational']))
    if maximum['case_id']!='p38539-pqq-la_model__conditioned_La__seed-1_sample-2':raise ValueError('frozen largest-shift source differs')
    targets=[];lookups={}
    for row,candidate in [(q4,'adaptive_Ca'),(maximum,'adaptive_La')]:
        target={'case_id':row['case_id'],'metal':'La','candidate':candidate,'medium':'vacuum','versions':{},
                'selection_rule':'new_classification_abstention' if row is q4 else 'largest_absolute_pooled_R_change',
                'observed_pooled_delta_R':row['precision_numerical_comparison']['delta_R']['operational']}
        for version,pin in [('old',row['source_collection']),('new',row['precision_collection'])]:
            c=next(c for c in read_json(verify(pin))['cases'] if c['case_id']==row['case_id'])
            cell=c['matrix']['La'][c['aliases'][candidate]['representative']];low=cell['low']['vacuum']
            m=read_json(verify(low['manifest']));task=next(t for t in m['tasks'] if t['task_id']==low['task_id'])
            receipt=read_json(verify(low['receipt']));text=verify(low['output']).read_text()
            if not(receipt['normal_termination'] and receipt['scf_converged'] and receipt['returncode']==0):raise ValueError('source scalar incomplete')
            expected=8 if version=='old' else 1
            if receipt['parallelism']['nprocs']!=expected:raise ValueError('actual known rank differs')
            r={'task':task,'manifest':low['manifest'],'receipt':low['receipt'],'output':low['output'],
               'energy_hartree':low['energy_hartree'],'actual_ranks':expected,
               'SCF_cycles':int(re.search(r'SCF CONVERGED AFTER\s+(\d+)\s+CYCLES',text)[1]),
               'parameter_export':low['audit']['parameter_export'],'native_valence_electrons':low['audit']['native_valence_electrons'],
               'original_recipe':verify(task['input']).read_text(),'orca':m['orca'],'additional_exact_matches':[]}
            target['versions'][version]=r;lookups[(row['case_id'],version)]=geometry_key(task)
        a,b=target['versions'].values()
        if a['original_recipe']!=b['original_recipe'] or a['orca']!=b['orca'] or a['parameter_export']['sha256']!=b['parameter_export']['sha256']:
            raise ValueError('source Hamiltonian/state recipes differ')
        targets.append(target)
    prior_paths=[ROOT/'workspaces/native_gfn2_rank_panel_20260923/run_v2/manifest.json']
    prior_paths += [ROOT/f'workspaces/native_gfn2_ranks_20260923/run_v1/rank_{n}/manifest.json' for n in (1,4,8)]
    prior_paths += [ROOT/f'workspaces/native_pool_continuation_20260923/run_v1/stage{s}/manifest.json' for s in (1,2)]
    prior_paths += [ROOT/f'workspaces/native_solvent_force_20260923/run_v1/stage{s}/manifest.json' for s in (0,1,2)]
    scopes=[]
    for path in prior_paths:
        m=read_json(path);tasks=m.get('all_tasks',m['tasks']);matches=[];cases=set()
        for t in tasks:
            cases.add(t.get('case_id',t.get('case')))
            key=geometry_key(t)
            for target in targets:
                for version,row in target['versions'].items():
                    if key!=lookups[(target['case_id'],version)]:continue
                    ep=Path(t['output_path']+'.execution.json')
                    match={'target':target['case_id'],'geometry':version,'manifest':record(path),'task_id':t['task_id'],
                           'recipe_matches':verify(t['input']).read_text()==row['original_recipe'],
                           'receipt':record(ep) if ep.exists() else None}
                    row['additional_exact_matches'].append(match);matches.append(match)
        scopes.append({'manifest':record(path),'task_count':len(tasks),'cases':sorted(c for c in cases if c is not None),'exact_geometry_state_matches':matches})
    missing=[]
    for target in targets:
        for version,ranks in [('old',1),('new',8)]:
            row=target['versions'][version]
            if row['additional_exact_matches']:raise ValueError('exact prior match found; inspect reuse before declaring missing')
            missing.append({'case_id':target['case_id'],'candidate':target['candidate'],'metal':'La','medium':'vacuum',
                            'geometry':version,'requested_ranks':ranks,'source_task':row['task'],'source_manifest':row['manifest']})
    out={'status':'read_only_inventory_no_submission','comparison':record(comparison),'targets':targets,'searched_scopes':scopes,
         'existing_exact_cells':4,'missing_cells':missing,'maximum_new_scalar_calls_if_authorized':4,
         'new_MACE_calls':0,'new_optimization_calls':0,'new_DFT_calls':0,'new_molecular_calls':0,
         'matching_policy':'Exact parsed elemental coordinates and charge/multiplicity; recipe/ORCA checked separately. No approximate-coordinate reuse.',
         'proposed_policy':'Unchanged nativeGFN2/300K/MaxIter500/mixer/tolerances/NoAutostart; old and new geometry at1and8ranks. No continuation/seed variation or preferred-value selection.',
         'implementation':record(__file__)}
    write_new(output,out);print(json.dumps({'existing_cells':4,'missing_calls':len(missing),'scanned_prior_tasks':sum(s['task_count'] for s in scopes),'exact_prior_matches':sum(len(s['exact_geometry_state_matches']) for s in scopes)},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('comparison','output'):p.add_argument('--'+k,required=True,type=Path)
    inventory(**vars(p.parse_args()))
