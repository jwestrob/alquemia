#!/usr/bin/env python3
"""Original reference run/collection bodies, with one reviewed joint binding gate."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('radial_binding',HERE/'BIND.py');B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B)
require,pin,verify,read=B.require,B.pin,B.verify,B.read
def save(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
    return pin(path)
CELLS=tuple('ecp_pbe0__'+g for g in B.IDS)

def stage_metadata(text):
    inp=read(pin(HERE/'SOURCE_INPUTS.json'));old=B.module(inp['old_reference'],'reference_stage')
    return old.stage_metadata(text)

def gate(config,manifest,release):
    bp=pin(config);b,ctx,r=B.release(bp,pin(release),'reference');m,inp=ctx[:2]
    require(pin(manifest)==m['reference_manifest'],'Unchanged original reference manifest')
    old=B.module(inp['old_reference'],'radial_reference_gate')
    cfg,ref,original,geometry,cells,centers=old.load(verify(inp['old_reference_config']))
    prepared=read(m['reference_manifest']);p=ref.load_panel(verify(prepared['native_panel']))
    require(prepared['native_panel']==m['reference_panel'] and p['cells']==cells and prepared['anchor_reference']==b['anchor_reference'],'Exact original panel and central references')
    evidence=dict(release=pin(release),manifest=m['reference_manifest'],bound=bp,candidate_result=b['candidate_result'],selected_review=b['candidate_review'],selected=b['selected'],parameters=b['parameters'],root_acceptance=b['root_acceptance'],task_output_directories=r['task_output_directories'])
    return cfg,ref,p,geometry,cells,prepared,evidence

def run(config,manifest,release,cell,out):
    cfg,ref,p,geometry,cells,prepared,evidence=gate(config,manifest,release);require(cell in CELLS,'Cell outside four fresh references')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','0'))==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Exact four MPI ranks required')
    out=Path(out).resolve();require(str(out)==evidence['task_output_directories'][cell] and not out.exists(),'Exact released fresh task path required; no retry')
    out.parent.mkdir(parents=True,exist_ok=True);save(Path(str(out)+'.release.json'),evidence|{'cell_id':cell,'job_id':os.environ['SLURM_JOB_ID'],'native_calls_authorized':1})
    return ref.run(prepared['native_panel']['path'],cell,out)

def collect(config,manifest,release,result_inventory,out):
    cfg,ref,p,geometry,cells,prepared,evidence=gate(config,manifest,release);paths=json.loads(Path(result_inventory).read_text());require(set(paths)==set(CELLS),'Exact four-result denominator required')
    rows={};available=0;cost=0.
    for c in cells:
        gid=c['geometry_id'];meta=geometry['rows'][gid];row={'status':'unavailable','geometry':c['geometry'],'input':meta['input'],'bodies':meta['bodies'],'reference':None};rows[gid]=row
        expected=Path(evidence['task_output_directories'][c['cell_id']])/'RESULT.json';require(Path(paths[c['cell_id']]).resolve()==expected,'Result inventory differs from released task path')
        if not expected.exists():row['reason']='missing_result';continue
        authority=pin(expected)
        try:
            saved=read(authority);cost+=saved.get('allocated_cpu_seconds',0.)
            require(saved['status']=='complete' and saved['returncode']==0 and all(saved['checks'].values()) and saved['cell_id']==c['cell_id'] and saved['panel']==prepared['native_panel'] and saved['adapter']==cfg['native_adapter'],'Native reference admission/provenance failed')
            replay=ref.collect_one(p,c,expected.parent,saved['returncode']);require(all(saved[k]==v for k,v in replay.items()),'Raw reference replay differs')
            require(len(saved['forces_kcal_mol_A'])==28 and all(math.isfinite(v) for f in saved['forces_kcal_mol_A'] for v in f),'Finite full28 force inventory required')
            native_start=json.loads((expected.parent/'START.json').read_text());require(native_start['ranks']==4 and native_start['input']==saved['input'] and verify(saved['input']).read_text()==ref.deck(p,c,4,native_start['maxcore_MiB']),'Executed native recipe/resource differs')
            stages=stage_metadata(verify(saved['output']).read_text());center=prepared['anchor_reference'][c['metal']]['energy_stages']
            changed=stages['COSX_grid_settings_signature']!=center['COSX_grid_settings_signature'] or stages['final_exchange_recomputed']!=center['final_exchange_recomputed'] or stages['analytic_gradient_stage']!=center['analytic_gradient_stage']
            row.update(status='complete',reference={'energy':saved['energy_hartree']*ref.HA,'energy_hartree':saved['energy_hartree'],'forces':saved['forces_kcal_mol_A'],'raw_result':authority,'native_checks':saved['checks'],'state':{k:c[k] for k in ('charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')},'energy_stages':stages,'changed_central_energy_stage_signature':changed,'output':saved['output'],'input':saved['input'],'engrad':saved['engrad'],'physical_state_qualified':False,'force_derivative_qualified':False,'method_accuracy_qualified':False});available+=1
        except Exception as exc:row.update(reason=repr(exc),raw_result=authority)
    return save(out,{'status':'complete' if available==4 else 'unavailable','declared':4,'available':available,'rows':rows,'anchor_reference':prepared['anchor_reference'],'binding':evidence,'result_inventory':pin(result_inventory),'allocated_native_CPU_seconds':cost,'native_calls_this_collection':0,'model_scoring_performed':False,'physical_accuracy_qualified':False})

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    for action in ('run','collect'):
        q=s.add_parser(action)
        for name in ('bound','manifest','release','out'):q.add_argument('--'+name,required=True)
        q.add_argument('--cell' if action=='run' else '--result-inventory',required=True)
    a=p.parse_args()
    if a.action=='run':raise SystemExit(run(a.bound,a.manifest,a.release,a.cell,a.out))
    print(json.dumps(collect(a.bound,a.manifest,a.release,a.result_inventory,a.out)))
