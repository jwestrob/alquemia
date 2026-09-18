"""Actual source-array equivalence and timing for an isolated OpenMP build."""
from pathlib import Path
import argparse
import gc
import json
import os
import resource
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A,InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_ddx_multipole_bridge import preflight as bridge_preflight
from mace_ddx_source_parallel import patch

PROTOCOL='native_ddX_dense_source_openmp_equivalence_scaling_v1'
THREADS=(1,8,64)


def prepare(bridge,collection,software,plan,output):
    p=bridge_preflight(bridge);co=read_json(collection);sw=read_json(software);build=read_json(verify(sw['manifest']))
    if not co['complete'] or not co['checks_pass'] or co['manifest']!=record(bridge) or sw['status']!='built_and_imported':
        raise InvalidArtifact('completed source bridge and isolated build required')
    if build['protocol']!='isolated_pyddx_dense_source_openmp_build_v1':raise InvalidArtifact('source parallel build required')
    env=Path(software).resolve().parent/'venv';module=list((env/'lib/python3.11/site-packages').glob('pyddx*.so'))
    if len(module)!=1:raise InvalidArtifact('expected one isolated module')
    py=env/'bin/python';pypin={**record(py),'path':str(py)}
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in p['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name in ('mace_ddx_source_scaling.py','mace_ddx_source_parallel.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    group=next(g for g in p['groups'] if g['group_id']=='GGR_2FW0_primary')
    m=dict(protocol=PROTOCOL,bridge=record(bridge),collection=record(collection),software=record(software),plan=record(plan),
           python=pypin,module=record(module[0]),group=group,reference=co['groups'][group['group_id']],model=p['model'],grid=p['grid'],
           threads=list(THREADS),requested_potential_calls=12,requested_integral_calls=12,new_solver_calls=0,
           tolerances=dict(potential_au=1e-12,psi=1e-12),implementation={x.name:record(x) for x in impl.glob('*.py')})
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')


def preflight(path):
    m=read_json(path);p=bridge_preflight(verify(m['bridge']));co=read_json(verify(m['collection']));sw=read_json(verify(m['software']))
    b=read_json(verify(sw['manifest']))
    if m['protocol']!=PROTOCOL or m['threads']!=list(THREADS) or m['model']!=p['model'] or m['grid']!=p['grid'] or m['tolerances']!=dict(potential_au=1e-12,psi=1e-12):
        raise InvalidArtifact('source scaling scope changed')
    if m['cache_key']!=cache_key({k:v for k,v in m.items() if k!='cache_key'}) or not co['checks_pass']:
        raise InvalidArtifact('source scaling manifest or prerequisite changed')
    expected=next(g for g in p['groups'] if g['group_id']=='GGR_2FW0_primary')
    if m['group']!=expected or m['reference']!=co['groups'][expected['group_id']]:raise InvalidArtifact('actual source reference changed')
    for pin in (m['plan'],m['python'],m['module'],m['reference']['arrays'],*m['implementation'].values(),*b['pins']):verify(pin)
    import tarfile
    original=b['original_native_source']
    with tarfile.open(verify(original['archive'])) as f:text=f.extractfile(original['member']).read().decode()
    if patch(text)!=verify(b['parallel_native_source']).read_text():raise InvalidArtifact('native arithmetic or loop patch changed')
    return m


def execute(path,output):
    import pyddx
    m=preflight(path);root=Path(path).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64 or record(pyddx.__file__)!=m['module'] or pyddx.__version__!='0.9.0' or np.__version__!='1.26.4':
        raise InvalidArtifact('isolated source scaling runtime differs')
    physical=m['group']['physical'];centers=np.asfortranarray(np.array(physical['coordinates_A']).T/BOHR_TO_A);radii=np.array(physical['radii_A'])/BOHR_TO_A
    start=time.monotonic();cpu=time.process_time()
    with np.load(verify(m['reference']['arrays']),allow_pickle=False) as reference:
        labels=sorted(m['reference']['states'])
        for threads in THREADS:
            d=root/'groups'/str(threads)/'attempt_0001';d.mkdir(parents=True,exist_ok=False);tic=time.monotonic()
            model=pyddx.Model(sphere_centres=centers,sphere_radii=radii,**{**m['model'],'n_proc':threads},**m['grid'])
            cavity=np.array(model.cavity);r=dict(threads=threads,states={},setup_wall_seconds=time.monotonic()-tic,
                cavity_identical=np.array_equal(cavity,reference['cavity_bohr']),manifest=record(path),cache_key=m['cache_key'],
                actual_parameters={k:v for k,v in model.input_parameters.items() if k not in ('sphere_centres','sphere_radii')})
            arrays=dict(cavity_bohr=cavity)
            for label in labels:
                before=time.monotonic();before_cpu=time.process_time();state=dict(status='failed')
                try:
                    if not r['cavity_identical']:raise InvalidArtifact('native cavity changed')
                    coeff=np.asfortranarray(reference[label+'_coefficients']);phi=np.array(model.multipole_electrostatics(coeff,derivative_order=0)['phi'])
                    psi=np.array(model.multipole_psi(coeff));elapsed=time.monotonic()-before;used=time.process_time()-before_cpu
                    pe=float(np.max(abs(phi-reference[label+'_phi'])));se=float(np.max(abs(psi-reference[label+'_psi'])))
                    if not np.isfinite(phi).all() or not np.isfinite(psi).all():raise InvalidArtifact('nonfinite source values')
                    state.update(status='computed',potential_error_au=pe,psi_error=se,checks_pass=pe<=1e-12 and se<=1e-12,
                         potential_bitwise_identical=np.array_equal(phi,reference[label+'_phi']),psi_bitwise_identical=np.array_equal(psi,reference[label+'_psi']),
                         native_property_wall_seconds=elapsed,native_property_CPU_seconds=used)
                    arrays[label+'_phi']=phi;arrays[label+'_psi']=psi
                except Exception as exc:state.update(failure_reason=str(exc))
                r['states'][label]=state;write_new(d/(label+'.json'),state);print(threads,label,state,flush=True)
            np.savez_compressed(d/'arrays.npz',**arrays);r['arrays']=record(d/'arrays.npz');write_new(d/'result.json',r)
            model=None;arrays.clear();gc.collect()
    write_new(root/'execution_receipt.json',dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,slurm_job_id=os.environ['SLURM_JOB_ID'],new_solver_calls=0))
    return collect(path,output)


def collect(path,output):
    m=preflight(path);root=Path(path).resolve().parent;groups={};checks=[];speedup={}
    for threads in THREADS:
        pin=record(root/'groups'/str(threads)/'attempt_0001/result.json');r=read_json(verify(pin));verify(r['arrays'])
        if r['manifest']!=record(path) or r['cache_key']!=m['cache_key']:raise InvalidArtifact('source scaling receipt changed')
        groups[str(threads)]=dict(receipt=pin,**r)
        checks.append(dict(name=str(threads)+'_cavity',pass_=r['cavity_identical']))
        for k,v in {**m['model'],'n_proc':threads,**m['grid']}.items():
            if k not in ('model','enable_force'):checks.append(dict(name=str(threads)+'_'+k,pass_=r['actual_parameters'].get(k)==v))
        for label,state in r['states'].items():checks.append(dict(name=str(threads)+'_'+label,pass_=state.get('checks_pass',False)))
    for threads in THREADS:
        speedup[str(threads)]={}
        for label,state in groups[str(threads)]['states'].items():
            serial=groups['1']['states'][label]
            speedup[str(threads)][label]=serial['native_property_wall_seconds']/state['native_property_wall_seconds'] if serial['status']==state['status']=='computed' else None
    result=dict(protocol=PROTOCOL,manifest=record(path),groups=groups,checks=checks,checks_pass=all(c['pass_'] for c in checks),
                speedup_relative_to_serial=speedup,new_solver_calls=0,new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False)
    write_new(output,result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare')
    for k in ('bridge','collection','software','plan','output'):p.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        p=sub.add_parser(op);p.add_argument('--manifest',required=True)
        if op!='dry-run':p.add_argument('--output',required=True)
    a=parser.parse_args()
    if a.command=='prepare':r=prepare(a.bridge,a.collection,a.software,a.plan,a.output)
    elif a.command=='dry-run':r=preflight(a.manifest)
    elif a.command=='execute':r=execute(a.manifest,a.output)
    else:r=collect(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('groups','implementation')}))
