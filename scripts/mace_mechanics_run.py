"""Finite cheap-core/scaffold mechanics tasks through the existing MACE runner."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import check_atoms,accepted_attempt
from mace_global_benchmark import snapshot,numerical_parent_gate
from mace_curvature import actual_collection
from mace_gb import MODEL as GB_MODEL
from mace_short_engine import COMPONENT,ADAPTER

CORE='alquemia.mace_mechanics_core.v1'
GB='alquemia.mace_mechanics_gb.v1'
SHORT='alquemia.mace_mechanics_short.v1'
POLICY='source_original_H_coupled_metal_peptide_scaffold_response_v1'
TOL={'total_charge_e':1e-5,'energy_kcal_mol':.01,'force_max_eV_A':.001}
UNAVAILABLE={'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,
             'entropy_correction_kcal_mol':None,'calibrated_class':None}


def preparation(path):
    p=read_json(path)
    if p['policy_id']!=POLICY or set(p['cases'])!={'GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9'}:
        raise InvalidArtifact('declared physical preparation required')
    for ref in [p['agreement'],*p['implementation'].values(),p['DFT_tasks'],p['DFT_reference']]:verify(ref)
    for name,ref in p['cases'].items():
        c=read_json(verify(ref))
        for field in ('source_preparation','physical_atoms','core_jacobians','full_jacobians'):verify(c[field])
        if len(c['grids'])!=17:raise InvalidArtifact('physical grid count changed')
        for grid in c['grids'].values():
            la,ca=(grid['endpoints'][metal] for metal in ('La','Ca'))
            a,b=(xyz(verify(e['xyz'])) for e in (la,ca))
            if la['charge']-ca['charge']!=1 or a[0][1:]!=b[0][1:] or a[1:]!=b[1:]:
                raise InvalidArtifact('paired physical core state changed')
    if set(p['physical_cases'])!={'GGR_1GLG','ALPHA_1F6S','ALPHA_6IP9'}:
        raise InvalidArtifact('physical protein set changed')
    for full in p['physical_cases'].values():
        for grid in full['grids'].values():
            la,ca=(grid['endpoints'][metal] for metal in ('La','Ca'))
            a,b=(xyz(verify(e['xyz'])) for e in (la,ca))
            if la['charge']-ca['charge']!=1 or len(a)!=len(b):raise InvalidArtifact('full paired charge/count differs')
            for x,y in zip(a,b):
                if x[1:]!=y[1:] or (x[0]!=y[0] and (x[0],y[0])!=('La','Ca')):
                    raise InvalidArtifact('full paired atom identity/coordinates differ')
    return p


def expected_tasks(p,short=False):
    tasks=[]
    if short:
        for name,full in p['physical_cases'].items():
            for point,grid in full['grids'].items():
                for metal,e in grid['endpoints'].items():
                    tasks.append({'task_id':f'full__{name}_{point}_{metal}','kind':'full','case_id':name,'point':point,
                                  'metal':metal,**e,'energy_component':COMPONENT})
        selected={k:v for k,v in p['cases'].items() if k.startswith('ALPHA')}
    else:selected=p['cases']
    for name,ref in selected.items():
        c=read_json(verify(ref))
        for point,grid in c['grids'].items():
            if short and point!='center':continue
            for metal,e in grid['endpoints'].items():
                task_id=f'{name}_{point}_{metal}'
                if task_id in p['reused_GGR']:continue
                t={'task_id':('core__' if short else '')+task_id,'kind':'core','case_id':name,
                   'point':point,'metal':metal,**e}
                if short:t['energy_component']=COMPONENT
                tasks.append(t)
    for t in tasks:t['variant']='primary'
    if len(tasks)!=(106 if short else 116):raise InvalidArtifact('expected finite cheap inventory changed')
    return tasks


def prepare_core(prepared,reference_mace,reference_gb,short_gate,agreement,output):
    p=preparation(prepared);c,cm=actual_collection(reference_mace);g,gm=actual_collection(reference_gb)
    gate,sm=actual_collection(short_gate)
    if (not gate['numerical_checks_pass'] or cm['checkpoint_label']!='medium' or gm['source_mace_collection']!=record(reference_mace)
            or sm['model']['checkpoint']!=cm['model']['checkpoint']):raise InvalidArtifact('matching completed curvature/short gates required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,(*cm['implementation'],*p['implementation'],'mace_mechanics_run.py','mace_short_engine.py'))
    model=copy.deepcopy(cm['model']);model['preparation_policy']=POLICY
    m={'schema_version':CORE,'protocol_id':'mace_polar_1m_coupled_core_curvature_v1',
       'agreement':record(agreement),'preparation':record(prepared),'reference_mace':record(reference_mace),
       'reference_gb':record(reference_gb),'short_gate':record(short_gate),'numerical_reference':cm['numerical_reference'],
       'software':cm['software'],'model':model,'implementation':pins,'tasks':expected_tasks(p),'tolerances':TOL,
       'run_inventory':{'new_MACE_calls':116,'reused_MACE_calls':20,'new_DFT_calls':0},**UNAVAILABLE}
    seal(out,m);return validate(out/'manifest.json')


def prepare_short(prepared,short_gate,agreement,output):
    p=preparation(prepared);gate,gm=actual_collection(short_gate)
    if not gate['numerical_checks_pass']:raise InvalidArtifact('short energy/force gate failed')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,(*gm['implementation'],*p['implementation'],'mace_mechanics_run.py'))
    model=copy.deepcopy(gm['model']);model['preparation_policy']=POLICY
    m={'schema_version':SHORT,'protocol_id':'mace_polar_1m_coupled_scaffold_short_v1',
       'agreement':record(agreement),'preparation':record(prepared),'short_gate':record(short_gate),
       'numerical_reference':gm['numerical_reference'],'software':gm['software'],'model':model,
       'implementation':pins,'tasks':expected_tasks(p,True),'tolerances':TOL,
       'run_inventory':{'new_short_calls':106,'reused_GGR_short_center_gradients':4,'new_DFT_calls':0},**UNAVAILABLE}
    seal(out,m);return validate(out/'manifest.json')


def seal(out,m):
    for t in m['tasks']:
        t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
    write_new(out/'manifest.json',m)


def prepare_gb(collection,solver_validation,agreement,output):
    c,m=actual_collection(collection);validate(verify(c['manifest']))
    if m['schema_version']!=CORE:raise InvalidArtifact('actual core mechanics MACE collection required')
    reference=read_json(solver_validation);rm=read_json(verify(reference['manifest']))
    if not reference['numerical_checks_pass']:raise InvalidArtifact('GB solver gate failed')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,m['implementation']);result=copy.deepcopy(m)
    for t in result['tasks']:
        t.pop('cache_key');r=c['rows'][t['task_id']]
        t.update(solver='native',platform='CUDA',solvent_dielectric=78.5,
                 source_density=r['density_coefficients'],source_vacuum_energy_eV=r['energy_eV'])
    result.update(schema_version=GB,protocol_id='mace_polar_1m_coupled_core_obc2_v1',
                  agreement=record(agreement),model=GB_MODEL,software=rm['software'],implementation=pins,
                  source_mace_collection=record(collection),solver_validation=record(solver_validation),
                  run_inventory={'new_GB_calls':116,'reused_GB_calls':20,'new_MACE_calls':0,'new_DFT_calls':0})
    seal(out,result);return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest);schema=m['schema_version'];short=schema==SHORT;is_gb=schema==GB
    if schema not in (CORE,GB,SHORT) or m['tolerances']!=TOL:raise InvalidArtifact('unsupported mechanics execution schema')
    p=preparation(verify(m['preparation']));gate,gm=actual_collection(verify(m['short_gate']))
    if not gate['numerical_checks_pass'] or numerical_parent_gate(m)['status']!='pass':raise InvalidArtifact('numerical prerequisite failed')
    software=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),software['python'],software['requirements'],
                *read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    if not short:
        ref,rm=actual_collection(verify(m['reference_mace']));gb,gbm=actual_collection(verify(m['reference_gb']))
        if gbm['source_mace_collection']!=m['reference_mace']:raise InvalidArtifact('reused solvent/model mismatch')
        old={t['task_id']:t for t in rm['tasks']}
        prepared_states={f'{name}_{point}_{metal}':endpoint
                         for name,ref in p['cases'].items()
                         for point,grid in read_json(verify(ref))['grids'].items()
                         for metal,endpoint in grid['endpoints'].items()}
        for key,r in p['reused_GGR'].items():
            t=old[r['source_task_id']]
            state=prepared_states[key]
            if (r['xyz']['sha256']!=t['xyz']['sha256'] or state['xyz']!=r['xyz'] or
                    state['charge']!=t['charge'] or state['spin_multiplicity']!=t['spin_multiplicity']):
                raise InvalidArtifact('cached GGR coordinates/charge/spin differ')
            verify(r['xyz'])
    if is_gb:
        source,source_m=actual_collection(verify(m['source_mace_collection']));validate(verify(source['manifest']))
        validated=read_json(verify(m['solver_validation']));vm=read_json(verify(validated['manifest']))
        if (not validated['numerical_checks_pass'] or m['model']!=GB_MODEL or m['software']!=vm['software'] or
                source_m['preparation']!=m['preparation']):raise InvalidArtifact('GB model/source validation mismatch')
    else:
        parent=gm if short else rm
        expected=copy.deepcopy(parent['model']);expected['preparation_policy']=POLICY
        if m['model']!=expected or m['software']!=parent['software']:raise InvalidArtifact('MACE method/checkpoint differs')
    expected=expected_tasks(p,short)
    if len(m['tasks'])!=len(expected):raise InvalidArtifact('finite task inventory changed')
    for t,base in zip(m['tasks'],expected):
        if any(t.get(k)!=v for k,v in base.items()):raise InvalidArtifact('task differs from physical preparation')
        check_atoms(xyz(verify(t['xyz'])),t['charge'])
        if is_gb:
            r=source['rows'][t['task_id']]
            if (t['source_density']!=r['density_coefficients'] or t['source_vacuum_energy_eV']!=r['energy_eV'] or
                    t['solver']!='native' or t['platform']!='CUDA' or t['solvent_dielectric']!=78.5):
                raise InvalidArtifact('GB source/solver changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('scientific cache key changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def collect_mechanics(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:
                if m['schema_version']==GB:verify(r['serialized_system'])
                valid.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=valid[-1] if valid else {'status':'unavailable','energy_eV':None}
    reuse={}
    if m['schema_version']!=SHORT:
        p=preparation(verify(m['preparation']));source_key='reference_gb' if m['schema_version']==GB else 'reference_mace'
        c,_=actual_collection(verify(m[source_key]))
        reuse={key:{'source_collection':m[source_key],'source_task_id':r['source_task_id'],'result':c['rows'][r['source_task_id']]}
               for key,r in p['reused_GGR'].items()}
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete',
            'manifest':record(mp),'protocol_id':m['protocol_id'],'rows':rows,'reused_rows':reuse,'attempts':attempts,**UNAVAILABLE}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare-core')
    for k in ('prepared','reference-mace','reference-gb','short-gate','agreement','output'):a.add_argument('--'+k,required=True)
    a=sub.add_parser('prepare-short')
    for k in ('prepared','short-gate','agreement','output'):a.add_argument('--'+k,required=True)
    a=sub.add_parser('prepare-gb')
    for k in ('collection','solver-validation','agreement','output'):a.add_argument('--'+k,required=True)
    args=vars(p.parse_args());command=args.pop('command')
    print(json.dumps({'prepare-core':prepare_core,'prepare-short':prepare_short,'prepare-gb':prepare_gb}[command](**args),indent=2))
