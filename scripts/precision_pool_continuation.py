"""Inventory fixed precision pools for a uniform two-pass native qualification."""
from __future__ import annotations
import argparse,json,re
import numpy as np
from functools import lru_cache
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,energy,xyz
from compact_solvation import input_text,diagnostics
from run_orca_task_manifest import load_manifest_tasks,_completed_attempt_is_valid
from structure_informed_starts import scf_details
from native_pool_continuation import recipe

PROTOCOL='precision_three_geometry_pool_two_native_continuations_v1'
CANDIDATES=('origin','adaptive_Ca','adaptive_La')
FOLDS=('q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-2',
       'p38539-pqq-la_model__conditioned_La__seed-1_sample-2',
       'a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1',
       'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4')
CRYSTALS=('1H4I','4MAE','1KB0')


@lru_cache(None)
def pinned_json(path,sha):
    return read_json(verify({'path':path,'sha256':sha}))


def data(pin):return pinned_json(pin['path'],pin['sha256'])


@lru_cache(None)
def normalized(path,sha):
    verify({'path':path,'sha256':sha});m,ts=load_manifest_tasks(Path(path))
    return m,{t['task_id']:t for t in ts}


def fingerprint(r):
    return (r['xyz']['sha256'],r['charge'],r['multiplicity'],r['medium'],
            r['gbw']['sha256'],r['xtbw']['sha256'],r['parameter_export']['sha256'])


def prior_pairs(stage_pairs):
    index={};audited=[]
    for a,b in stage_pairs:
        ca,cb=read_json(a),read_json(b);ma=data(ca['manifest']);mb=data(cb['manifest'])
        if ca['stage']!=1 or cb['stage']!=2 or mb['previous']!=record(a):raise InvalidArtifact('prior two-pass linkage differs')
        first={r['task_id']:r for r in ca['rows']}
        for r in cb['rows']:
            p=first[r['task_id']]
            if p['status']!='confirmed_restart' or r['status']!='confirmed_restart':continue
            source=p['source']
            for stage,row,manifest in [(1,p,ma),(2,r,mb)]:
                t=next(t for t in manifest['tasks'] if t['task_id']==row['task_id'])
                body=verify(t['input']).read_text()
                if 'EnGrad' in body:continue
                if body!=recipe(t['charge'],t['multiplicity'],t['medium']):raise InvalidArtifact('prior continuation recipe differs')
                for pin in (row['actual']['receipt'],row['actual']['output'],row['seed_after']['preserved_after'],row['seed_after']['preserved_gbw_after']):verify(pin)
            # Scalar-only reuse; preserve gradient-bearing existing outputs separately.
            if any(next(t for t in m['tasks'] if t['task_id']==r['task_id'])['gradient_requested'] for m in (ma,mb)):continue
            fp=fingerprint(source);entry={'stage1_collection':record(a),'stage2_collection':record(b),'task_id':r['task_id'],
                'stage1':p,'stage2':r,'actual_rank':read_json(verify(r['actual']['receipt']))['parallelism']['nprocs']}
            index.setdefault(fp,[]).append(entry)
        audited.append({'stage1':record(a),'stage2':record(b)})
    return index,audited


def inventory(precision34,precision225,agreement,output):
    p34=read_json(precision34);p225=read_json(precision225)
    selected=[{'case_id':r['case_id'],'role':r['role'],'expected_class':r['known_class'],'collection':r['collection']}
              for r in p34['rows'] if r['role']=='calibration' or r['case_id'] in CRYSTALS]
    if len(selected)!=28 or sum(r['role']=='calibration' for r in selected)!=25:raise InvalidArtifact('designated25 plus3 crystals required')
    for cid in FOLDS:
        r=next(r for r in p225['rows'] if r['case_id']==cid)
        selected.append({'case_id':cid,'role':'consumed_noncanonical_development','expected_class':r['expected_class'],'collection':r['precision_collection']})
    if len({r['case_id'] for r in selected})!=32:raise InvalidArtifact('population duplicates or missing source')
    root=Path(__file__).resolve().parents[1]
    pairs=[(root/'workspaces/native_pool_continuation_20260923/run_v1/stage1/collection.json',root/'workspaces/native_pool_continuation_20260923/run_v1/stage2/collection.json'),
           (root/'workspaces/precision_geometry_continuation_20260923/run_v2/stage1/collection.json',root/'workspaces/precision_geometry_continuation_20260923/run_v2/stage2/collection.json')]
    reuse,prior=prior_pairs(pairs);rows=[];cases=[]
    for group in selected:
        case=next(c for c in data(group['collection'])['cases'] if c['case_id']==group['case_id'])
        if case['pool']['status']!='available':raise InvalidArtifact('declared precision pool unavailable: '+group['case_id'])
        case_rows=[]
        for metal in ('Ca','La'):
            for candidate in CANDIDATES:
                cell=case['matrix'][metal][candidate]
                for medium in ('vacuum','alpb'):
                    low=cell['low'][medium];row={**{k:group[k] for k in ('case_id','role','expected_class')},
                        'candidate':candidate,'metal':metal,'medium':medium,'status':'unavailable','reason':None,
                        'pool_collection':group['collection'],'source_cell':{k:low[k] for k in ('manifest','task_id','receipt','output','energy_hartree')},
                        'MACE_energy_eV':cell['components']['MACE_eV'],'cell_xyz':cell['xyz'],'continued_reuse':None}
                    try:
                        manifest=low['manifest'];m,ts=normalized(manifest['path'],manifest['sha256']);t=next(t for t in m['tasks'] if t['task_id']==low['task_id'])
                        receipt=data(low['receipt']);output_path=verify(low['output']);text=output_path.read_text()
                        if not _completed_attempt_is_valid(verify(low['receipt']),output_path,manifest_sha256=manifest['sha256'],task=ts[t['task_id']],runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
                            raise InvalidArtifact('actual source receipt incomplete/incompatible')
                        if energy(output_path)!=low['energy_hartree']:raise InvalidArtifact('source energy differs')
                        details=scf_details(text);source_input=verify(t['input']).read_text()
                        stripped=re.sub(r'^\s*MaxIter 500\s*\n','',source_input,flags=re.M)
                        if stripped!=input_text(t['charge'],t['multiplicity'],medium,'native'):raise InvalidArtifact('source Hamiltonian differs')
                        if (t['metal'],t['medium'],t['multiplicity'])!=(metal,medium,1):raise InvalidArtifact('source state differs')
                        audit=low.get('audit') or diagnostics(low,t)
                        for pin in (t['xyz'],audit['parameter_export']):verify(pin)
                        expected=xyz(verify(cell['xyz']));actual=xyz(verify(t['xyz']))
                        if [a[0] for a in expected]!=[a[0] for a in actual]:raise InvalidArtifact('pool/source atom order differs')
                        delta=float(np.max(np.abs(np.asarray([a[1:] for a in expected])-np.asarray([a[1:] for a in actual]))))
                        if delta>1e-12:raise InvalidArtifact('pool/source coordinates exceed existing reconciliation tolerance')
                        runtime=verify(receipt['artifacts']['runtime_input']);base=runtime.with_suffix('')
                        gbw=Path(str(base)+'.gbw');xtbw=Path(str(base)+'.xtbw')
                        row.update(xyz=t['xyz'],input=t['input'],output=low['output'],receipt=low['receipt'],
                            charge=t['charge'],multiplicity=t['multiplicity'],electron_count=details['electrons'],
                            parameter_export=audit['parameter_export'],source_energy_hartree=low['energy_hartree'],
                            pool_source_max_coordinate_component_difference_A=delta,
                            source_SCF_cycles=details['cycles'],source_ranks=receipt['parallelism']['nprocs'],
                            source_manifest=manifest,source_task_id=t['task_id'],orca=receipt['orca_executable'],
                            gbw=record(gbw) if gbw.exists() else None,xtbw=record(xtbw) if xtbw.exists() else None)
                        if not row['gbw'] or not row['xtbw']:raise InvalidArtifact('matching_cold_seed_pair_missing')
                        found=reuse.get(fingerprint(row),[])
                        if found:
                            # Prefer lowest path deterministically, never an energy-based choice.
                            found.sort(key=lambda p:(p['stage2_collection']['path'],p['task_id']))
                            row['continued_reuse']=found[0]
                        row.update(status='compatible_seed_pair_available',gbw_bytes=gbw.stat().st_size,xtbw_bytes=xtbw.stat().st_size)
                    except (InvalidArtifact,KeyError,OSError) as exc:row['reason']=str(exc)
                    rows.append(row);case_rows.append(row)
        cases.append({**group,'denominator':12,'available':sum(r['status']=='compatible_seed_pair_available' for r in case_rows),
            'two_pass_reuse_cells':sum(r['continued_reuse'] is not None for r in case_rows)})
        print(group['case_id'],cases[-1]['available'],'/12, reuse',cases[-1]['two_pass_reuse_cells'],flush=True)
    available=sum(r['status']=='compatible_seed_pair_available' for r in rows);reuse_count=sum(r['continued_reuse'] is not None for r in rows)
    result={'protocol_id':PROTOCOL,'precision34':record(precision34),'precision225':record(precision225),'agreement':record(agreement),
        'cases':cases,'rows':rows,'case_denominator':32,'cell_denominator':384,'available_seed_pairs':available,
        'missing_seed_or_unsupported_cells':384-available,'two_pass_reuse_cells':reuse_count,'reused_molecular_calls':2*reuse_count,
        'maximum_new_scalar_calls':2*(available-reuse_count),'logical_calls_before_reuse':768,'existing_continuation_inventories':prior,
        'implementation':record(__file__),'new_molecular_calls_in_inventory':0,'production_changed':False}
    write_new(output,result)
    return {k:v for k,v in result.items() if k not in ('rows','cases')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('precision34','precision225','agreement','output'):p.add_argument('--'+k,required=True,type=Path)
    print(json.dumps(inventory(**vars(p.parse_args())),indent=2))
