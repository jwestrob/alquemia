"""Native quantum response to a frozen physical protein field; opt-in research."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz,paired
from mace_file_checks import cached_file_checks
from mace_omol_vacuum import METHOD,embedded_input,parse_endpoint

PROTOCOL='normalized_r2scan3c_fixed_protein_field_response_v1'
SOURCE_SHA='3e9e6bc11464157d1bb1b2df9151f0db97bcdfec17e8df4b8432267e583fca23'
TOL={'variational_response_kcal_mol':.05,'algebra_kcal_mol':1e-7}
NULLS={'aqueous_affinity_score':None,'reference':None,'calibrated_class':None,'combined_gradient':None,
       'relaxation_correction':None,'response_status':'response_model_not_validated','baseline_changed':False}


@cached_file_checks
def sources(source):
    from mace_explicit_field_short import collect
    if record(source)['sha256']!=SOURCE_SHA:raise InvalidArtifact('declared frozen-density source differs')
    r=read_json(source);mp=verify(r['manifest']);actual=collect(mp)
    for value in (r,actual):verify(value['collection_implementation'])
    if ({k:v for k,v in r.items() if k!='collection_implementation'}!={k:v for k,v in actual.items() if k!='collection_implementation'}
            or r['collection_implementation']['sha256']!=actual['collection_implementation']['sha256']
            or not r['numerical_gate_pass'] or not r['charge_representation_gate_pass']):
        raise InvalidArtifact('actual frozen source does not replay')
    fm=read_json(mp);cm=read_json(verify(fm['charges']));qm=read_json(verify(cm['quantum']))
    return r,fm,cm,qm


def pointcharge_text(task):
    from affordable_common import BOHR_TO_A
    p=np.loadtxt(verify(task['points']),skiprows=1)*BOHR_TO_A
    s=read_json(verify(task['state']));w=read_json(verify(task['weights']))
    # Use full-precision source positions, not rounded utility probe coordinates.
    xyzA=np.array([s['physical_atoms'][i]['xyz_A'] for i in w['physical_indices']])
    if not np.allclose(p,xyzA,atol=1e-11,rtol=0):raise InvalidArtifact('potential/physical positions differ')
    return str(len(xyzA))+'\n'+''.join(' '.join(format(float(v),'.17g') for v in (q,*x))+'\n' for q,x in zip(w['weights_e'],xyzA))


def prepare(source,agreement,accounting,output):
    r,fm,cm,qm=sources(source);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={};paths={k:verify(v) for k,v in fm['implementation'].items()}
    for name in ('mace_responsive_field.py','mace_omol_vacuum.py'):paths[name]=Path(__file__).with_name(name)
    for name,path in paths.items():shutil.copyfile(path,impl/name);pins[name]=record(impl/name)
    tasks=[]
    for old in fm['tasks']:
        d=out/old['task_id'];d.mkdir();xp=d/'core.xyz';ip=d/'endpoint.inp';pc=d/'environment.pc'
        shutil.copyfile(verify(old['xyz']),xp);ip.write_text(embedded_input(old['charge']));pc.write_text(pointcharge_text(old))
        t={'task_id':old['task_id'],'case_id':old['case_id'],'metal':old['metal'],'charge':old['charge'],'multiplicity':1,
            'input':record(ip),'xyz':record(xp),'pointcharges':record(pc),'source_state':old['state'],'source_task':old,
            'output_path':str(d/'endpoint.out'),'engrad_path':str(d/'endpoint.engrad'),'task_type':'analytic_gradient'}
        tasks.append(t)
    tasks.sort(key=lambda t:(-len(xyz(verify(t['xyz']))),t['task_id']))
    m={'protocol_id':PROTOCOL,'source':record(source),'agreement':record(agreement),'accounting':record(accounting),
        'implementation':pins,'tasks':tasks,'orca':qm['orca'],'method':METHOD,'tolerances':TOL,
        'energy_scope':'embedded_permanent_field_endpoint_no_external_self_no_solvent',
        'execution_policy':{k:pins[n] for k,n in [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]},
        'compute_budget':None,'wall_time_limit':None,'new_DFT_calls':8,**NULLS}
    for t in tasks:t['cache_key']=key(t,m)
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def key(t,m):return cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},**{k:m[k] for k in ('protocol_id','method','source','orca','implementation')}})


@cached_file_checks
def validate(manifest):
    from affordable_workflow import dry_run
    m=read_json(manifest);r,fm,cm,qm=sources(verify(m['source']))
    if m['protocol_id']!=PROTOCOL or m['method']!=METHOD or m['tolerances']!=TOL or m['orca']!=qm['orca']:
        raise InvalidArtifact('native method/settings changed')
    for pin in [m['accounting'],*m['implementation'].values()]:verify(pin)
    if len(m['tasks'])!=8 or {t['task_id'] for t in m['tasks']}!={t['task_id'] for t in fm['tasks']}:
        raise InvalidArtifact('embedded endpoint inventory differs')
    for t in m['tasks']:
        old=next(v for v in fm['tasks'] if v['task_id']==t['task_id'])
        for k in ('case_id','metal','charge'):
            if t[k]!=old[k]:raise InvalidArtifact('paired state changed')
        if (t['multiplicity']!=1 or t['source_task']!=old or t['source_state']!=old['state'] or t['xyz']['sha256']!=old['xyz']['sha256']
            or verify(t['input']).read_text()!=embedded_input(old['charge']) or verify(t['pointcharges']).read_text()!=pointcharge_text(old)
            or t['cache_key']!=key(t,m)):raise InvalidArtifact('embedded source/state/field/cache changed')
    for name in sorted({t['case_id'] for t in m['tasks']}):
        a,b=(next(t for t in m['tasks'] if t['case_id']==name and t['metal']==metal) for metal in ('La','Ca'))
        paired(verify(a['xyz']),verify(b['xyz']),0,-1)
        if a['pointcharges']['sha256']!=b['pointcharges']['sha256']:raise InvalidArtifact('paired permanent field differs')
    return dry_run(manifest)


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest);m,rows=executed(manifest);old=read_json(verify(m['source']));checks=[];pairs={}
    for t in m['tasks']:
        r=rows[t['task_id']]
        if r['status']!='complete':continue
        try:
            r.update(parse_endpoint(t,verify(r['output']),t['engrad_path'],permanent_field=True))
            olde=old['cases'][t['case_id']]['endpoints'][t['metal']]
            response=r['energy_hartree']*HA_TO_KCAL-olde['DFT_vacuum_kcal_mol']-olde['direct_exact_kcal_mol']
            r['electronic_response_kcal_mol']=response
            r['variational_check_pass']=response<=TOL['variational_response_kcal_mol']
            r['pointcharges']=t['pointcharges'];r['source_state']=t['source_state']
            checks.append({'name':t['task_id']+'_variational','response_kcal_mol':response,'pass':r['variational_check_pass']})
        except (ValueError,OSError) as exc:r.update(status='invalid',reason=str(exc),energy_hartree=None)
    complete=all(r['status']=='complete' for r in rows.values())
    if complete:
        for name in old['cases']:
            a,b=(rows[name+'_'+metal] for metal in ('Ca','La'))
            pairs[name]={'embedded_R_kcal_mol':(a['energy_hartree']-b['energy_hartree'])*HA_TO_KCAL,
                'electronic_response_R_kcal_mol':a['electronic_response_kcal_mol']-b['electronic_response_kcal_mol']}
    return {'protocol_id':PROTOCOL,'manifest':record(manifest),'rows':rows,'status':'complete' if complete else 'incomplete',
        'checks':checks,'variational_gate_pass':complete and all(c['pass'] for c in checks),'paired':pairs,
        'collection_implementation':record(__file__),'parser_implementation':record(Path(__file__).with_name('mace_omol_vacuum.py')),'new_DFT_calls':8,**NULLS}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('source','agreement','accounting','output'):a.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.source,a.agreement,a.accounting,a.output)
    elif a.op=='dry-run':r=validate(a.manifest)
    elif a.op=='execute':
        from affordable_workflow import execute
        validate(a.manifest);r=execute(a.manifest)
    else:r=collect(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','tasks')},indent=2))
