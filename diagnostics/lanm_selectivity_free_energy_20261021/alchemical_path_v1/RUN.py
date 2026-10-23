#!/usr/bin/env python3
"""Prepare a frozen path or execute its 20-call numerical check under Slurm."""
import argparse
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import time

PROTOCOL='alchemical_path_v1'
GATES={'energy_absolute_kcal_mol':1e-4,'force_max_component_kcal_mol_A':1e-3,
       'spatial_FD_max_component_kcal_mol_A':1e-3,'lambda_derivative_difference_kcal_mol':1e-3}
FIXTURES=('water_central','acetate_bi_central')
THREAD_KEYS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])==record,'Changed artifact: '+record['path']);return Path(record['path'])


def load(record):return json.loads(verify(record).read_text())


def module(record,name):
    spec=importlib.util.spec_from_file_location(name,verify(record));m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m


def save(path,data):
    tmp=Path(str(path)+'.tmp');tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)


def pair_inputs(controls,fixture):
    a,b=controls[fixture+'_La']['input'],controls[fixture+'_Dy']['input']
    for label in ('key','parameters'):require(a[label]['sha256']==b[label]['sha256'],'Endpoint force fields differ')
    la=verify(a['xyz']).read_text();dy=verify(b['xyz']).read_text()
    al=[v.split()for v in la.splitlines()[1:]if v.strip()];bl=[v.split()for v in dy.splitlines()[1:]if v.strip()]
    require(len(al)==len(bl)==a['atoms']==b['atoms'],'Endpoint dimensions differ')
    require(al[0][0]==bl[0][0]=='1' and al[0][2:5]==bl[0][2:5] and al[0][5:]==['500']and bl[0][5:]==['501'],
            'Endpoint metal position/type/topology differs')
    require(al[1:]==bl[1:] and a['charge']==b['charge'],'Endpoint nonmetals or charge differ')
    return a,la


def inspect(cfg):
    require(cfg['protocol']==PROTOCOL and cfg['gates']==GATES,'Changed protocol/gates')
    require(cfg['fixtures']==list(FIXTURES) and cfg['maximum_native_invocations']==20,'Changed finite inventory')
    require(cfg['implementation']==pin(__file__),'Wrong run adapter')
    verify(cfg['specification']);renderer=module(cfg['renderer'],'lambda_renderer')
    helper=module(cfg['gpu_adapter'],'qualified_gpu_adapter');authority=load(cfg['gpu_authority'])
    require(authority['adapter']==cfg['gpu_adapter'],'GPU adapter authority differs')
    controls=helper.inspect_inputs(authority,bound=True)
    gpu=load(cfg['gpu_parity_result'])
    require(gpu['status']=='gpu_parity_passed'and gpu['config']==cfg['gpu_authority'],'GPU parity not admitted')
    require(cfg['gpu_binary']==authority['executables']['tinker9'],'GPU binary differs');verify(cfg['gpu_binary'])
    cpu_config=load(load(authority['cpu_controls_result'])['config'])
    require(cfg['cpu_binary']==cpu_config['executables']['testgrad'],'CPU binary differs');verify(cfg['cpu_binary'])
    require(set(cfg['native_semantics'])=={'kvdw','kpolar','initprm'},'Native semantics inventory')
    for name,record in cfg['native_semantics'].items():
        require(verify(record)==Path(authority['source']['path'])/'tinker/source'/(name+'.f'),'Native semantics source path')
    require("vdwindex = 'CLASS'" in verify(cfg['native_semantics']['initprm']).read_text(),'Native default CLASS absent')
    require(cfg['lambda_values']==list(renderer.LAMBDAS),'Changed lambda inventory')
    return renderer,helper,controls


def prepare(cfg,out,renderer,controls):
    windows=[];fixtures={}
    for fixture in FIXTURES:
        a,xyz_text=pair_inputs(controls,fixture)
        key_text=verify(a['key']).read_text();prm_text=verify(a['parameters']).read_text()
        new_xyz,atoms=renderer.render_xyz(xyz_text);types=[int(row[5])for row in atoms[1:]]
        fixtures[fixture]={'La':controls[fixture+'_La']['input'],'Dy':controls[fixture+'_Dy']['input']}
        for lam in renderer.LAMBDAS:
            name=f'{fixture}__lambda_{lam:.8f}';directory=out/name;directory.mkdir()
            key,params=renderer.render_key(prm_text,key_text,types,lam)
            (directory/'system.xyz').write_text(new_xyz);(directory/'system.key').write_text(key)
            shutil.copyfile(verify(a['parameters']),directory/'amoeba09.prm')
            require(pin(directory/'amoeba09.prm')['sha256']==a['parameters']['sha256'],'Parameters altered')
            windows.append({'id':name,'fixture':fixture,'lambda':lam,'atoms':a['atoms'],'charge':a['charge'],
                            'spatial_FD':lam==.5,'parameters':params,
                            'input':{k:pin(directory/v)for k,v in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]}})
    manifest={'protocol':PROTOCOL,'implementation':cfg['implementation'],'renderer':cfg['renderer'],
              'specification':cfg['specification'],'gpu_authority':cfg['gpu_authority'],
              'native_semantics':cfg['native_semantics'],'fixtures':fixtures,'windows':windows,
              'native_program_invocations':{'GPU':18,'CPU_midpoints':2,'total':20},
              'physical_accuracy_qualified':False,'free_energy_estimated':False}
    save(out/'MANIFEST.json',manifest)
    return manifest


def validate_manifest(record,cfg,renderer,controls):
    manifest=load(record)
    for field in ('implementation','renderer','specification','gpu_authority','native_semantics'):
        require(manifest[field]==cfg[field],'Manifest provenance differs: '+field)
    require(len(manifest['windows'])==18 and len({x['id']for x in manifest['windows']})==18,'Manifest window count')
    require(manifest['native_program_invocations']=={'GPU':18,'CPU_midpoints':2,'total':20},'Manifest call cap')
    for fixture in FIXTURES:
        a,xyz=pair_inputs(controls,fixture);new_xyz,atoms=renderer.render_xyz(xyz)
        require(manifest['fixtures'][fixture]=={'La':controls[fixture+'_La']['input'],'Dy':controls[fixture+'_Dy']['input']},'Fixture binding differs')
        windows=[x for x in manifest['windows']if x['fixture']==fixture]
        require([x['lambda']for x in windows]==list(renderer.LAMBDAS),'Window ordering differs')
        for row in windows:
            key,params=renderer.render_key(verify(a['parameters']).read_text(),verify(a['key']).read_text(),
                                            [int(v[5])for v in atoms[1:]],row['lambda'])
            require(row['id']==f"{fixture}__lambda_{row['lambda']:.8f}"and row['parameters']==params,'Rendered parameter mismatch')
            require(row['spatial_FD']==(row['lambda']==.5)and row['atoms']==a['atoms']and row['charge']==a['charge'],'Window state mismatch')
            require(verify(row['input']['xyz']).read_text()==new_xyz and verify(row['input']['key']).read_text()==key,'Rendered input mismatch')
            require(verify(row['input']['parameters']).read_bytes()==verify(a['parameters']).read_bytes(),'Nonmetal force field changed')
    return manifest


def parse_analytic(text,atoms,helper):
    """Same native positive evidence as parity parser, with Numer deliberately off."""
    require(not helper.FAILURE.search(text),'Native failure/nonfinite marker')
    energies=re.findall(r'^\s*Total Potential Energy\s*:\s*('+helper.NUMBER+r')\s+Kcal/mole\s*$',text,re.M)
    require(len(energies)==1 and text.count('Total Potential Energy')==1,'Native energy incomplete')
    grad=[]
    pattern=r'\s*Anlyt\s+(\d+)\s+('+helper.NUMBER+r')\s+('+helper.NUMBER+r')\s+('+helper.NUMBER+r')\s+('+helper.NUMBER+r')\s*'
    for line in text.splitlines():
        m=re.fullmatch(pattern,line)
        if m:
            require(int(m[1])==len(grad)+1,'Gradient duplicate/nonsequential');grad.append([helper.number(m[i])for i in (2,3,4)])
    require(len(grad)==atoms and not re.search(r'\bNumer\s+\d+',text),'Incomplete or unexpected gradient')
    summary=re.findall(r'Anlyt\s+Total Gradient Norm Value\s*:?\s*('+helper.NUMBER+')',text)
    require(len(summary)==1 and helper.number(summary[0])>=0,'Native final gradient summary absent')
    return {'total_energy_kcal_mol':helper.number(energies[0]),'gradient_kcal_mol_A':grad,
            'force_kcal_mol_A':[[-v for v in row]for row in grad],'finite_difference':{'status':'not_requested'}}


def compare(native,reference):
    af,bf=native['force_kcal_mol_A'],reference['force_kcal_mol_A']
    require(len(af)==len(bf)and all(len(row)==3 for row in af+bf),'Comparison dimensions')
    errors=[a-b for ar,br in zip(af,bf)for a,b in zip(ar,br)]
    energy=native['total_energy_kcal_mol']-reference['total_energy_kcal_mol']
    require(all(math.isfinite(v)for v in errors+[energy]),'Nonfinite comparison')
    return {'energy_difference_kcal_mol':energy,'force_max_component_difference_kcal_mol_A':max(map(abs,errors)),
            'force_RMS_component_difference_kcal_mol_A':math.sqrt(sum(x*x for x in errors)/len(errors)),
            'force_RMS_denominator':len(errors),'passed':abs(energy)<=GATES['energy_absolute_kcal_mol']and
            max(map(abs,errors))<=GATES['force_max_component_kcal_mol_A']}


def execute(cfg,out,manifest,helper,controls,result):
    require(os.environ.get('SLURM_JOB_ID'),'Native execution requires Slurm')
    require(os.environ.get('SLURM_CPUS_PER_TASK')=='4'and os.environ.get('SLURM_JOB_PARTITION')=='gpu_h200','Requires 4 CPUs/H200 allocation')
    gpu=os.environ.get('CUDA_VISIBLE_DEVICES','').split(',')
    require(len(gpu)==1 and gpu[0]not in ('','-1','NoDevFiles'),'Requires one visible GPU')
    require(all(os.environ.get(k)=='1'for k in THREAD_KEYS),'All local thread caps must be one')
    result['rows']=[]
    for engine,windows in [('GPU',manifest['windows']),('CPU',[r for r in manifest['windows']if r['lambda']==.5])]:
        for window in windows:
            require(len(result['native_invocations'])<20,'Native invocation cap exceeded')
            row={'id':window['id'],'engine':engine,'fixture':window['fixture'],'lambda':window['lambda'],'status':'running','window':window}
            result['rows'].append(row);directory=out/engine/window['id'];directory.mkdir(parents=True)
            for key,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:
                shutil.copyfile(verify(window['input'][key]),directory/name)
            binary=verify(cfg['gpu_binary'if engine=='GPU'else'cpu_binary'])
            argv=[str(binary)]+(['testgrad']if engine=='GPU'else[])+['system.xyz','Y','Y'if window['spatial_FD']else'N']
            if window['spatial_FD']:argv+=['0.00001']
            helper.invoke(argv,directory,directory/'gradient.out',300,result)
            text=(directory/'gradient.out').read_text()
            native=helper.parse_gradient(text,window['atoms'])if window['spatial_FD']else parse_analytic(text,window['atoms'],helper)
            native['gradient_output']=pin(directory/'gradient.out');row['native']=native
            for key,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:
                require(pin(directory/name)['sha256']==window['input'][key]['sha256'],'Native modified inputs')
            if window['spatial_FD']:require(native['finite_difference']['passed'],'Spatial FD failed')
            if engine=='GPU'and window['lambda']in (0.,1.):
                metal='La' if window['lambda']==0 else 'Dy'
                row['endpoint_comparison']=compare(native,controls[window['fixture']+'_'+metal]['native'])
                require(row['endpoint_comparison']['passed'],'Endpoint equivalence failed')
            if engine=='CPU':
                matched=next(x for x in result['rows']if x['engine']=='GPU'and x['id']==window['id'])
                row['CPU_minus_GPU_comparison']=compare(native,matched['native'])
                require(row['CPU_minus_GPU_comparison']['passed'],'Midpoint CPU/GPU parity failed')
            row['status']='complete';save(out/'RESULT.json',result)
    result['lambda_derivatives']=[]
    for fixture in FIXTURES:
        energies={x['lambda']:x['native']['total_energy_kcal_mol']for x in result['rows']if x['fixture']==fixture and x['engine']=='GPU'}
        derivatives=[(energies[.5+h]-energies[.5-h])/(2*h)for h in (1e-4,5e-5)]
        difference=abs(derivatives[0]-derivatives[1])
        result['lambda_derivatives'].append({'fixture':fixture,'steps':[1e-4,5e-5],
                                             'central_FD_dE_dlambda_kcal_mol':derivatives,
                                             'difference_kcal_mol':difference,'passed':difference<=GATES['lambda_derivative_difference_kcal_mol']})
        require(difference<=GATES['lambda_derivative_difference_kcal_mol'],'Lambda derivative consistency failed')
    require(len(result['native_invocations'])==20 and len(result['rows'])==20,'Incomplete numerical inventory')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True)
    args=ap.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    result={'protocol':PROTOCOL,'status':'started','config':pin(args.config),'implementation':pin(__file__),
            'native_invocations':[],'physical_accuracy_qualified':False,'free_energy_estimated':False,
            'production_sampling_released':False,'allocation':{k:os.environ.get(k)for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','SLURM_JOB_PARTITION','CUDA_VISIBLE_DEVICES')}}
    cfg=None
    try:
        cfg=load(result['config']);renderer,helper,controls=inspect(cfg)
        require(cfg['stage']in ('prepare','run','saved_preflight'),'Unknown stage')
        if cfg['stage']=='prepare':
            prepare(cfg,out,renderer,controls);result['manifest']=pin(out/'MANIFEST.json')
            validate_manifest(result['manifest'],cfg,renderer,controls);result['status']='preparation_passed'
        else:
            manifest=validate_manifest(cfg['manifest'],cfg,renderer,controls);result['manifest']=cfg['manifest']
            if cfg['stage']=='saved_preflight':result['status']='saved_preflight_passed'
            else:
                execute(cfg,out,manifest,helper,controls,result)
                validate_manifest(cfg['manifest'],cfg,renderer,controls)
                result['status']='numerical_path_passed'
        verify(result['config']);inspect(cfg)
        return 0
    except Exception as exc:
        result.update(status='failed',error_type=type(exc).__name__,error=str(exc))
        for row in result.get('rows',[]):
            if row['status']=='running':row['status']='failed_execution_or_gate'
        return 1
    finally:
        result.update(stage=cfg.get('stage')if cfg else None,elapsed_seconds=time.monotonic()-start,
                      native_program_invocation_count=len(result['native_invocations']),
                      utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'native_calls':len(result['native_invocations']),'result':str(out/'RESULT.json')}))


if __name__=='__main__':raise SystemExit(main())
