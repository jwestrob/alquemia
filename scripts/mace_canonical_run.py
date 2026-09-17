"""Finite canonical PQQ MACE tasks using the existing runner and exact caches."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz,paired
from mace_hybrid import accepted_attempt,check_atoms
from mace_global_benchmark import snapshot,numerical_parent_gate
from mace_curvature import actual_collection
from mace_gb import MODEL as GB_MODEL

MACE='alquemia.mace_canonical.v1'
GB='alquemia.mace_canonical_gb.v1'
POLICY='canonical_frozen_PQQ_core_direct_MACE_obc2_v1'
TOL={'total_charge_e':1e-5,'energy_kcal_mol':.01,'force_max_eV_A':.001}
PENDING={'calibration_status':'not_yet_calibrated','calibrated_class':None,
         'relaxation_correction_kcal_mol':None,'entropy_correction_kcal_mol':None}


def inventory(path):
    data=read_json(path)
    if data['schema_version']!='alquemia.mace_canonical_inventory.v1' or data['status']!='complete':
        raise InvalidArtifact('complete exact canonical input inventory required')
    for ref in [data['implementation'],data['source_verifier'],*data['sources'].values()]:verify(ref)
    cal=read_json(verify(data['sources']['calibration']))
    expected={r['panel_id']:(r['class'],'calibration') for r in cal['scores']}
    expected.update({'1H4I':('Ca','retrospective_structural_transfer'),'4MAE':('La','retrospective_structural_transfer'),
                     '1KB0':('Ca','retrospective_external_class_transfer')})
    if len(data['rows'])!=28 or {r['case_id'] for r in data['rows']}!=set(expected):raise InvalidArtifact('declared28-case inventory changed')
    for row in data['rows']:
        if (row['expected_class'],row['evaluation_role'])!=expected[row['case_id']] or row['prospectively_blind']:
            raise InvalidArtifact('canonical label/calibration/evaluation assignment changed')
        for ref in [row['source_geometry'],row['source_manifest'],*row['preparation_provenance'].values()]:verify(ref)
        if row['explicit_water_inventory']!=[] or row['cofactor']['microstate_id']!='pqq_ox_3minus_v1':
            raise InvalidArtifact('water/cofactor state changed')
        endpoints=row['endpoints'];la,ca=(endpoints[metal] for metal in ('La','Ca'))
        paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        for state in endpoints.values():
            verify(state['baseline_input']);check_atoms(xyz(verify(state['xyz'])),state['charge'])
            if state['spin_multiplicity']!=1:raise InvalidArtifact('spin changed')
    return data


def all_tasks(data):
    return [{'task_id':r['case_id']+'_'+metal,'case_id':r['case_id'],'metal':metal,'kind':'core','variant':'primary',
             'expected_class':r['expected_class'],'evaluation_role':r['evaluation_role'],**e}
            for r in data['rows'] for metal,e in r['endpoints'].items()]


def effective_model(model):
    return {k:v for k,v in model.items() if k!='preparation_policy'}


def reuse_records(data,parent,reference_mace,reference_gb):
    c,m=actual_collection(reference_mace);g,gm=actual_collection(reference_gb)
    if (effective_model(parent['model'])!=effective_model(m['model']) or parent['software']!=m['software'] or
            gm['model']!=GB_MODEL or gm['source_mace_collection']!=record(reference_mace)):
        raise InvalidArtifact('existing crystal cache is not the same MACE/GB method')
    old={t['task_id']:t for t in m['tasks']};rows={r['case_id']:r for r in data['rows']};reuse={}
    for name in ('1H4I','4MAE'):
        for metal in ('La','Ca'):
            key=f'PQQ_{name}_archived_{metal}';t=old[key];e=rows[name]['endpoints'][metal]
            if (t['xyz']['sha256']!=e['xyz']['sha256'] or t['charge']!=e['charge'] or
                    t['spin_multiplicity']!=e['spin_multiplicity']):
                raise InvalidArtifact('exact archived crystal coordinate/charge/spin cache mismatch')
            reuse[name+'_'+metal]={'source_task_id':key,'source_mace':record(reference_mace),'source_gb':record(reference_gb),
                                  'xyz_sha256':e['xyz']['sha256'],'charge':e['charge'],'spin_multiplicity':e['spin_multiplicity']}
    return reuse


def seal(out,m):
    for t in m['tasks']:
        t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
    write_new(out/'manifest.json',m)


def prepare_mace(source_inventory,parent,checkpoint,reference_mace,reference_gb,agreement,output):
    data=inventory(source_inventory);pc,pm=actual_collection(parent)
    if pm['checkpoint_label']!=checkpoint or checkpoint not in ('medium','large'):
        raise InvalidArtifact('declared checkpoint parent required')
    reuse=reuse_records(data,pm,reference_mace,reference_gb) if checkpoint=='medium' else {}
    if (checkpoint=='large' and (reference_mace or reference_gb)) or len(reuse)!=(4 if checkpoint=='medium' else 0):
        raise InvalidArtifact('unexpected cache inventory')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,(*pm['implementation'],'mace_canonical_run.py','mace_canonical.py',
                       'mace_curvature.py','affordable_response.py'))
    model=copy.deepcopy(pm['model']);model['preparation_policy']=POLICY
    m={'schema_version':MACE,'protocol_id':f'mace_polar_1{checkpoint[0]}_canonical_PQQ_direct_v1',
       'inventory':record(source_inventory),'agreement':record(agreement),'parent_collection':record(parent),
       'checkpoint_label':checkpoint,'model':model,'software':pm['software'],'numerical_reference':pm['numerical_reference'],
       'implementation':pins,'tolerances':TOL,'reused':reuse,
       'tasks':[t for t in all_tasks(data) if t['task_id'] not in reuse],**PENDING}
    seal(out,m);return validate(out/'manifest.json')


def prepare_gb(collection,solver_validation,agreement,output):
    c,m=actual_collection(collection);validate(verify(c['manifest']))
    if m['schema_version']!=MACE:raise InvalidArtifact('completed canonical MACE source required')
    v=read_json(solver_validation);vm=read_json(verify(v['manifest']))
    if not v['numerical_checks_pass']:raise InvalidArtifact('validated solvent backend required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,m['implementation']);result=copy.deepcopy(m)
    for t in result['tasks']:
        t.pop('cache_key');r=c['rows'][t['task_id']]
        t.update(solver='native',platform='CUDA',solvent_dielectric=78.5,
                 source_density=r['density_coefficients'],source_vacuum_energy_eV=r['energy_eV'])
    result.update(schema_version=GB,protocol_id=f'mace_polar_1{m["checkpoint_label"][0]}_canonical_PQQ_obc2_v1',
                  agreement=record(agreement),model=GB_MODEL,software=vm['software'],implementation=pins,
                  source_mace_collection=record(collection),solver_validation=record(solver_validation))
    seal(out,result);return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest);is_gb=m['schema_version']==GB
    if m['schema_version'] not in (MACE,GB) or m['tolerances']!=TOL:raise InvalidArtifact('unsupported canonical execution protocol')
    data=inventory(verify(m['inventory']));pc,pm=actual_collection(verify(m['parent_collection']))
    if m['checkpoint_label']!=pm['checkpoint_label'] or numerical_parent_gate(m)['status']!='pass':
        raise InvalidArtifact('checkpoint/numerical parent differs')
    software=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),software['python'],software['requirements'],
                *read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    if m['checkpoint_label']=='medium':
        entry=next(iter(m['reused'].values()))
        if m['reused']!=reuse_records(data,pm,verify(entry['source_mace']),verify(entry['source_gb'])):
            raise InvalidArtifact('reused canonical endpoints differ')
    elif m['reused']:raise InvalidArtifact('large has no declared matching cache')
    if is_gb:
        c,cm=actual_collection(verify(m['source_mace_collection']));validate(verify(c['manifest']))
        v=read_json(verify(m['solver_validation']));vm=read_json(verify(v['manifest']))
        if (not v['numerical_checks_pass'] or m['model']!=GB_MODEL or m['software']!=vm['software'] or
                cm['inventory']!=m['inventory'] or cm['checkpoint_label']!=m['checkpoint_label']):
            raise InvalidArtifact('solvent model/source changed')
    else:
        model=copy.deepcopy(pm['model']);model['preparation_policy']=POLICY
        if m['model']!=model or m['software']!=pm['software']:raise InvalidArtifact('MACE Hamiltonian/checkpoint changed')
    expected=[t for t in all_tasks(data) if t['task_id'] not in m['reused']]
    if len(expected)!=(52 if m['checkpoint_label']=='medium' else 56) or len(m['tasks'])!=len(expected):
        raise InvalidArtifact('finite canonical task count changed')
    for t,e in zip(m['tasks'],expected):
        if any(t.get(k)!=v for k,v in e.items()):raise InvalidArtifact('task differs from exact input inventory')
        if is_gb:
            r=c['rows'][t['task_id']]
            if (t['source_density']!=r['density_coefficients'] or t['source_vacuum_energy_eV']!=r['energy_eV'] or
                    t['solver']!='native' or t['platform']!='CUDA' or t['solvent_dielectric']!=78.5):
                raise InvalidArtifact('solvent charge/radius/cavity recipe changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('scientific cache identity changed')
    return {'status':'pass','tasks':len(m['tasks']),'reused':len(m['reused']),'manifest':record(manifest)}


def collect(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];reused={}
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
    key='source_gb' if m['schema_version']==GB else 'source_mace'
    for task_id,source in m['reused'].items():
        c,_=actual_collection(verify(source[key]))
        reused[task_id]={'source_collection':source[key],'source_task_id':source['source_task_id'],
                         'result':c['rows'][source['source_task_id']]}
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete',
            'manifest':record(mp),'rows':rows,'reused_rows':reused,'attempts':attempts,**PENDING}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare-mace')
    for key in ('source-inventory','parent','checkpoint','agreement','output'):a.add_argument('--'+key,required=True)
    for key in ('reference-mace','reference-gb'):a.add_argument('--'+key)
    a=sub.add_parser('prepare-gb')
    for key in ('collection','solver-validation','agreement','output'):a.add_argument('--'+key,required=True)
    args=vars(p.parse_args());command=args.pop('command')
    print(json.dumps({'prepare-mace':prepare_mace,'prepare-gb':prepare_gb}[command](**args),indent=2))
