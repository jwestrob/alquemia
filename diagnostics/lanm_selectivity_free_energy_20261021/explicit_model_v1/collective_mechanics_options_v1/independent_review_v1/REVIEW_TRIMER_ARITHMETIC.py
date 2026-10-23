#!/usr/bin/env python3
"""One saved-only eight-cycle arithmetic review; native admission is separate."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

HA=627.509474
D=0.529177210544*4.803206802
COLLECTOR='afef783b12ec89e0aeb968a5ed9642626b74bc62cbd7f6ed89b3f4d6e1f3fb4a'

def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def read(record):
    assert pin(record['path'])==record, 'Changed input pin'
    return json.loads(Path(record['path']).read_text())

def review(path,sha):
    target={'path':str(Path(path).resolve()),'sha256':sha};actual=read(target)
    assert actual['status']=='complete_saved_trimer_comparison' and len(actual['cycles'])==8
    cfg=actual['inputs'];assert cfg['implementation']['sha256']==COLLECTOR
    assert pin(cfg['implementation']['path'])==cfg['implementation']
    manifest=read(cfg['manifest']);model=read(cfg['model_result'])
    assert len(manifest['reference_cells'])==len(manifest['model_rows'])==len(model['rows'])==16
    assert model['status']=='complete' and model['native_program_invocations']==16
    qmom=actual['quantum_moment_rows'];qnew=actual['reference_native_metadata']
    assert len(qnew)==16 and set(qnew)=={c['cell_id'] for c in manifest['reference_cells']}
    qenergy={k:r['electronic_energy_hartree'] for k,r in qnew.items()}
    qenergy.update({k:r['record']['electronic_energy_hartree'] for k,r in manifest['reused_reference'].items()})
    mrows={k:r['parsed_moments'] for k,r in model['rows'].items()}
    mrows.update({k:r['record']['moments'] for k,r in manifest['reused_model'].items()})
    checks=[];rows={}
    def close(a,b,label,tol=1e-12):
        assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tol,label
        checks.append(label)
    for gid in ('aquo8_central_La','aquo8_central_Dy','aquo9_central_La','aquo9_central_Dy'):
        n=int(gid.split('_')[0][4:]);metal=gid.rsplit('_',1)[1];prefix=f'pbe0__aquo{n}__{metal}'
        for label in ('near','far'):
            cid=gid+'__'+label;saved=actual['cycles'][cid];declared=manifest['cycles'][cid];i,j=declared['selection']['water_indices']
            qt={prefix+'__'+label+'__metal_two_water':1,prefix+'__'+label+'__two_water':-1,
                prefix+f'__water_{i:02d}__metal_water':-1,prefix+f'__water_{j:02d}__metal_water':-1,
                prefix+'__ion_ghost_shell':1,prefix+f'__water_{i:02d}__water_only':1,prefix+f'__water_{j:02d}__water_only':1}
            mt={gid+'__'+label+'__metal_two_water':1,gid+'__'+label+'__two_water':-1,
                gid+f'__dimer_{i:02d}':-1,gid+f'__dimer_{j:02d}':-1,'ion_'+metal:1,
                gid+f'__water_{i:02d}':1,gid+f'__water_{j:02d}':1}
            assert qt==saved['reference_terms']==declared['reference_terms'] and mt==saved['model_terms']==declared['model_terms']
            checks.append(cid+':independent seven-term identities')
            qe=HA*math.fsum(c*qenergy[k] for k,c in qt.items())
            me=math.fsum(c*mrows[k]['energy_kcal_mol'] for k,c in mt.items())
            qmu=[D*math.fsum(c*qmom[k]['total_dipole_au'][a] for k,c in qt.items()) for a in range(3)]
            mmu=[math.fsum(c*math.fsum(site[a] for site in mrows[k]['induced_Debye']) for k,c in mt.items()) for a in range(3)]
            tmu=[math.fsum(c*mrows[k]['total_lab_Debye'][a] for k,c in mt.items()) for a in range(3)]
            components={name:math.fsum(c*mrows[k]['components'].get(name,{}).get('energy_kcal_mol',0.) for k,c in mt.items()) for name in set().union(*(mrows[k]['components'] for k in mt))}
            qb=D*math.fsum(abs(c)*qmom[k]['print_bounds_au']['total_dipole_au'] for k,c in qt.items());mb=14e-6
            for key,val in [('reference_electronic_E3_kcal_mol',qe),('model_E3_kcal_mol',me),('model_minus_reference_E3_kcal_mol',me-qe)]:close(val,saved[key],cid+':'+key,1e-9)
            for key,val in [('reference_mu3_Debye',qmu),('model_mu3_Debye',mmu),('model_total_M_mu3_Debye',tmu),('model_minus_reference_mu3_Debye',[a-b for a,b in zip(mmu,qmu)])]:
                for a in range(3):close(val[a],saved[key][a],cid+':'+key+str(a))
            for name,val in components.items():close(val,saved['model_E3_components_kcal_mol'][name],cid+':component:'+name,1e-9)
            close(qb+mb,saved['combined_print_bound_Debye_per_component'],cid+':print bound')
            assert max(abs(a-b) for a,b in zip(tmu,mmu))<=.0035+mb+2e-9*7
            rows[cid]={'selection':declared['selection'],'reference_E3_kcal_mol':qe,'model_E3_kcal_mol':me,'energy_error_kcal_mol':me-qe,'reference_mu3_Debye':qmu,'model_mu3_Debye':mmu,'dipole_error_Debye':[a-b for a,b in zip(mmu,qmu)],'combined_print_bound_Debye_per_component':qb+mb,'model_energy_components_kcal_mol':components}
    return {'status':'passed_independent_trimer_arithmetic_review','input_result':target,'implementation':pin(__file__),'collector':cfg['implementation'],'manifest':cfg['manifest'],'model_result':cfg['model_result'],'cycles':rows,'checks_passed':len(checks),'checks':checks,'denominator':{'new_reference':16,'new_model':16,'cycles':8},'new_native_calls':0,'raw_native_admission_repeated':False,'physical_acceptance_threshold':None,'limitations':['Raw native admission is a separate review.','Electronic E3 excludes D4; eight selected cycles do not reconstruct full B.','Printing bounds are numerical, not reference uncertainty or a physical gate.']}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--result',required=True);parser.add_argument('--sha256',required=True);parser.add_argument('--out',required=True);args=parser.parse_args()
    out=Path(args.out);assert not out.exists();tick=time.monotonic();cpu=time.process_time()
    try:r=review(args.result,args.sha256)
    except Exception as exc:r={'status':'failed_independent_trimer_arithmetic_review','error':repr(exc),'implementation':pin(__file__),'new_native_calls':0}
    r['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'max_local_threads':1}
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'status':r['status'],'review':pin(out)}))
