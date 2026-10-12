"""Four finite causal calls; this deliberately is not a factorial scoring backend."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'scripts'))
import context_factorial_gfn2 as f
import json, re, shutil
from concurrent.futures import ThreadPoolExecutor

BASE = ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_diagnosis_v1'
OLD = ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_v1'

def point(task):
    import numpy as np
    from tblite import library
    from tblite.interface import Calculator
    p=library.new_param();library.export_gfn2_param(p)
    t=library.new_table();library.dump_param(p,t);d=library.table_to_dict(t)
    f.write('ORIGINAL_PARAMS.json',d)
    assert d['dispersion']['d4']['s9']==5.0
    d['dispersion']['d4']['s9']=0.0
    f.write('PARAMS_s9zero.json',d)
    symbols,xyz=f.coordinates('core.xyz')
    c=Calculator(d,np.array([f.Z[s] for s in symbols]),np.array(xyz)/f.BOHR_A,charge=task['charge'],uhf=0)
    c.set('accuracy',.0001);c.set('max-iter',500);c.set('temperature',300*3.166808578545117e-6);c.set('verbosity',2)
    r=c.singlepoint()
    f.write('scientific.json',{'energy_hartree':float(r.get('energy')),'charges_e':r.get('charges').tolist(),
        'orbital_occupations_e':r.get('orbital-occupations').tolist(),'library_returned_without_error':True})

def run(task):
    rt=f.read(OLD/'RUNTIME.json');p=BASE/'points'/task['id'];p.mkdir(parents=True,exist_ok=False)
    f.write(p/'TASK.json',task)
    for backend in ('xtb','tblite'):
        out=p/backend;out.mkdir();shutil.copyfile(task['xyz']['path'],out/'core.xyz')
        assert f.pin(out/'core.xyz')['sha256']==task['xyz']['sha256']
        env=f.clean_env(rt['xtb_parameter_file']['path'])
        if backend=='xtb':
            original=Path(rt['xtb_parameter_file']['path']).read_text()
            changed,n=re.subn(r'(?m)^s9\s+5\.00000\s*$', 's9          0.00000',original);assert n==1
            (out/'param_gfn2-xtb.txt').write_text(changed)
            env['XTBPATH']=str(out)
            cmd=[rt['xtb']['path'],'core.xyz','--gfn','2','--chrg',str(task['charge']),'--uhf','0',
                 '--etemp','300','--acc','.0001','--iterations','500','--norestart','--json','--pop']
        else:
            cmd=[rt['tblite_launcher_path'],str(Path(__file__).resolve()),'point',str(p/'TASK.json')]
        f.call(cmd,out,env)
    parsed={b:f.parse_backend(p/b,b,task) for b in ('xtb','tblite')}
    old=f.read(OLD/'points'/task['id']/'RESULT.json')
    result={'task':task,'backends':parsed,'diagnostic_only':True,'original':f.pin(OLD/'points'/task['id']/'RESULT.json')}
    if all(v['status']=='available' for v in parsed.values()):
        result['energy_difference_Eh']=parsed['xtb']['energy_hartree']-parsed['tblite']['energy_hartree']
        result['max_population_difference_e']=max(abs(x-y) for x,y in zip(parsed['xtb']['charges_e'],parsed['tblite']['charges_e']))
        result['removed_ATM_energy_Eh']={b:old['backends'][b]['energy_hartree']-parsed[b]['energy_hartree'] for b in parsed}
    f.write(p/'RESULT.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ('task','backends')}),flush=True)
    return result

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='point':point(f.read(sys.argv[2]))
    else:
        tasks=[t for t in f.read(OLD/'manifest.json')['tasks'] if t['scout'] and t['medium']=='vacuum']
        assert len(tasks)==2
        f.write(BASE/'PLAN.json',{'tasks':tasks,'plan':f.pin(Path(__file__).with_name('DISPERSION_DIAGNOSTIC_PLAN.md')),
            'driver':f.pin(__file__),'runtime':f.pin(OLD/'RUNTIME.json'),'max_new_calls':4,'max_total_threads':2,
            'source_pins':[f.pin(p) for p in sorted((BASE/'source').glob('*'))]})
        with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(run,tasks))
        f.write(BASE/'RESULT.json',{'rows':rows})
