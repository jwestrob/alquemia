"""Separate component qualification; v1 absolute-energy failures remain immutable."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
import context_factorial_gfn2 as f
from concurrent.futures import ThreadPoolExecutor
import argparse, json, csv, urllib.request

OLD=ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_v1'
BASE=ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2'

def prepare():
    BASE.mkdir(exist_ok=False)
    original=f.read(OLD/'manifest.json')
    source_url='https://raw.githubusercontent.com/tblite/tblite/v0.7.0/src/tblite/disp/d4.f90'
    (BASE/'tblite_d4_v070.f90').write_bytes(urllib.request.urlopen(source_url).read())
    original['protocol_id']='Linux_xtb671_tblite070_component_factorial_v2'
    original['parent_absolute_protocol']=f.pin(OLD/'manifest.json')
    original['driver']=f.pin(__file__)
    original['plan']=f.pin(Path(__file__).with_name('COMPONENT_V2_PLAN.md'))
    original['source_proof']=[{'url':source_url,**f.pin(BASE/'tblite_d4_v070.f90')},
        f.pin(ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_diagnosis_v1/source/src_disp_dftd4.F90')]
    original['settings']['require_independent_parity']=False
    original['settings']['require_independent_component_parity']=True
    original['settings']['component_parity_Eh']=1e-6
    original['settings']['absolute_parity_status']='diagnostic_only_v1_gate_remains_failed'
    original['reuse']['scout_results']={t['id']:f.pin(OLD/'points'/t['id']/'RESULT.json') for t in original['tasks'] if t['scout']}
    original['scope']={'additional_inputs_per_backend':76,'total_numerical_threads':2,'no_retries':True}
    f.write(BASE/'manifest.json',original)
    print(json.dumps({'manifest':str(BASE/'manifest.json'),'new_inputs_per_backend':76}))

def raw_results(manifest):
    results={}
    for task in manifest['tasks']:
        p=BASE/'points'/task['id']/'RESULT.json'
        if task['scout']:
            rec=manifest['reuse']['scout_results'][task['id']];r=f.read(rec);p=Path(rec['path'])
        elif p.exists():r=f.read(p)
        else:r=None
        results[task['id']]=(task,r,f.pin(p) if r else None)
    return results

def collect(name):
    manifest=f.read(BASE/'manifest.json');raw=raw_results(manifest);groups={}
    for tid,(task,result,receipt) in raw.items():
        key='__'.join((task['protein'],task['context'],task['geometry'],task['metal']))
        groups.setdefault(key,{})[task['medium']]=(task,result,receipt)
    rows=[]
    for key,media in groups.items():
        v,vr,vp=media['vacuum'];w,wr,wp=media['water']
        assert v['xyz']==w['xyz'] and v['charge']==w['charge'] and v['uhf']==w['uhf']
        row={'id':key,'protein':v['protein'],'context':v['context'],'geometry':v['geometry'],
             'metal':v['metal'],'charge':v['charge'],'uhf':v['uhf'],'xyz':v['xyz'],
             'inputs':{'vacuum':vp,'water':wp},'status':'unavailable','accepted':False,
             'primary_vacuum_energy_hartree':None,'primary_water_energy_hartree':None,'primary_solvation_hartree':None,
             'independent_solvation_hartree':None,'component_error_hartree':None}
        if vr is None or wr is None:row['failure_reason']='not_executed'
        elif any(b['status']!='available' for r in (vr,wr) for b in r['backends'].values()):
            row['failure_reason']='backend_convergence_or_state_failure'
        else:
            de={b:wr['backends'][b]['energy_hartree']-vr['backends'][b]['energy_hartree'] for b in ('xtb','tblite')}
            err=de['xtb']-de['tblite']
            absolute={m:r['parity']['energy_difference_hartree'] for m,r in [('vacuum',vr),('water',wr)]}
            dq=max(vr['parity']['maximum_population_difference_e'],wr['parity']['maximum_population_difference_e'])
            ref=max(abs(r.get('reference_energy_difference_hartree',0.0)) for r in (vr,wr))
            accepted=abs(err)<=1e-6 and dq<=1e-5 and ref<=1e-6
            row.update(accepted=accepted,status='component_qualified' if accepted else 'component_failed',
                primary_vacuum_energy_hartree=vr['backends']['xtb']['energy_hartree'],
                primary_water_energy_hartree=wr['backends']['xtb']['energy_hartree'],
                primary_solvation_hartree=de['xtb'] if accepted else None,
                independent_solvation_hartree=de['tblite'],component_error_hartree=err,
                absolute_energy_errors_hartree=absolute,max_population_error_e=dq,
                original_absolute_energy_gate_pass=all(abs(x)<=1e-6 for x in absolute.values()),
                maximum_available_archived_primary_error_hartree=ref)
            if not accepted:row['failure_reason']='component_or_population_or_archived_parity_failed'
        rows.append(row)
    result={'protocol_id':manifest['protocol_id'],'manifest':f.pin(BASE/'manifest.json'),
            'declared_inputs':80,'executed_inputs':sum(r is not None for _,r,_ in raw.values()),
            'declared_components':40,'accepted_components':sum(r['accepted'] for r in rows),'rows':rows}
    f.write(BASE/name,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','manifest')}),flush=True)

def execute():
    manifest=f.read(BASE/'manifest.json');runtimes=f.read(manifest['runtime'])
    assert f.pin(__file__)==manifest['driver'] and f.pin(f.__file__)==manifest['implementation']
    for key in ('xtb','xtb_parameter_file','tblite_python'):assert f.pin(runtimes[key]['path'])==runtimes[key]
    tasks=[t for t in manifest['tasks'] if not t['scout']]
    assert len(tasks)==76
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda t:f.execute_point(t,manifest,BASE,runtimes),tasks))
    collect('FINAL_COLLECTION.json')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','execute','collect']);parser.add_argument('--output',default='CURRENT_COLLECTION.json');a=parser.parse_args()
    if a.action=='prepare':prepare()
    elif a.action=='execute':execute()
    else:collect(a.output)
