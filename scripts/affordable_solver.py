"""Real ESP/APBS validation with concurrent charging solves; never substitutes a baseline score on failure."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import fcntl
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from affordable_common import InvalidArtifact, read_json, write_new, record, verify, xyz, BOHR_TO_A, cache_key, snapshot_implementation
from affordable_environment import mbis_charges, prepare, collect


def real_esp_points(atoms, radii):
    """Two fixed Fibonacci shells per real atom; exclude every atomic interior."""
    coords=np.array([a[1:] for a in atoms]); points=[]
    n=32; indices=np.arange(n); z=1-2*(indices+.5)/n
    theta=indices*math.pi*(3-math.sqrt(5)); radial=np.sqrt(1-z*z)
    directions=np.column_stack((radial*np.cos(theta),radial*np.sin(theta),z))
    atom_radii=np.array([radii[a[0]] for a in atoms])
    for origin,r in zip(coords,atom_radii):
        for scale in (1.4,1.8):
            for point in origin+scale*r*directions:
                if np.all(np.linalg.norm(coords-point,axis=1)>=1.4*atom_radii-1e-8): points.append(point)
    if len(points)<32: raise InvalidArtifact('insufficient exterior ESP points')
    return np.array(points)/BOHR_TO_A


def run_command(command,cwd,log,cost):
    start=time.monotonic()
    with Path(log).open('x') as f:
        r=subprocess.run(['/usr/bin/time','-v','-o',str(Path(cost).resolve()),*command],cwd=cwd,stdout=f,stderr=subprocess.STDOUT,
                         env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},check=False)
    result={'command':command,'returncode':r.returncode,'wall_seconds':time.monotonic()-start,
            'log':record(log),'resource_usage':record(cost),'slurm_job_id':os.environ.get('SLURM_JOB_ID')}
    return result


def esp_check(task,orca_dir,radii,workdir=None):
    endpoint=Path(task['output_path']); d=Path(workdir) if workdir else endpoint.parent/'esp_check'
    if (d/'quality.json').exists():
        raise InvalidArtifact('existing ESP attempt requires explicit collection, not overwrite')
    d.mkdir(parents=True,exist_ok=True)
    charges=mbis_charges(endpoint,verify(task['xyz']),task['charge'])
    atoms=xyz(task['xyz']['path']); points=real_esp_points(atoms,radii)
    pp=d/'points_bohr.xyz'; pp.write_text(str(len(points))+'\n'+'\n'.join(' '.join(f'{v:.12f}' for v in row) for row in points)+'\n')
    gbw=endpoint.parent/'endpoint.runtime.gbw'; densities=endpoint.parent/'endpoint.runtime.densities'
    if not gbw.exists() or not densities.exists(): raise InvalidArtifact('wavefunction/density unavailable')
    potential=d/'potential.out'; utility=Path(orca_dir)/'orca_vpot'
    receipt=run_command([str(utility),str(gbw),'endpoint.runtime.scfp',str(pp.resolve()),str(potential.resolve())],
                        endpoint.parent,d/'utility.log',d/'resources.txt')
    write_new(d/'execution.json',receipt)
    if receipt['returncode']!=0: raise InvalidArtifact('orca_vpot failed; see actual utility output')
    values=np.loadtxt(potential,skiprows=1)
    if values.shape!=(len(points),4): raise InvalidArtifact('unexpected potential row count')
    if np.allclose(values[:,:3],points,atol=1e-5,rtol=0): qm=values[:,3]
    elif np.allclose(values[:,1:],points,atol=1e-5,rtol=0): qm=values[:,0]
    else: raise InvalidArtifact('potential coordinate/order mismatch')
    coords=np.array([a[1:] for a in atoms])/BOHR_TO_A
    approximation=np.sum(np.array(charges['charge_e'])[None,:]/np.linalg.norm(points[:,None,:]-coords[None,:,:],axis=2),axis=1)
    rms=float(np.sqrt(np.mean((qm-approximation)**2))); norm=float(np.sqrt(np.mean(qm**2)))
    passed=rms<=.005 or (norm>0 and rms/norm<=.10)
    quality={'status':'passed' if passed else 'failed','RMS_potential_error_au':rms,'relative_RMS':rms/norm if norm else None,
             'point_count':len(points),'actual_quantum_potential':record(potential),'points':record(pp),
             'gbw':record(gbw),'density':record(densities),'utility':record(utility),'charges':charges,
             'execution_receipt':record(d/'execution.json'),'units':'atomic_units_e_per_bohr',
             'rule':'relative_RMS<=0.10 OR absolute_RMS<=0.005_au'}
    write_new(d/'quality.json',quality)
    return quality,record(d/'quality.json')


def settings(physical,spacing=.5,padding=20.,extent=None,center=None):
    coords=np.array([a['xyz_A'] for a in physical])
    if center is None: center=(coords.min(0)+coords.max(0))/2
    if extent is None:
        # Multiples of 64 A keep both 0.5 and 0.4 A multigrid lattices at the
        # exact same extent, independently of endpoint and partition.
        extent=np.ceil((coords.max(0)-coords.min(0)+2*padding)/64)*64
    dims=np.rint(np.asarray(extent)/spacing).astype(int)+1
    return {'solute_dielectric':1.,'solvent_dielectric':78.54,'salt_molar':0.,'temperature_K':298.15,
            'grid_spacing_A':spacing,'grid_dimensions':dims.tolist(),'grid_center_A':np.asarray(center).tolist(),
            'surface':'mol','probe_radius_A':1.4,'charge_discretization':'spl2','boundary_condition':'mdh',
            'radii_policy':'Bondi_CHNOS_common_1p8A_metal_v1'}


def run_state(state,root,label,apbs):
    d=root/label; d.mkdir(parents=True,exist_ok=True)
    sp=d/'state.json'; write_new(sp,state)
    prepare(sp,d/'calculation')
    c=d/'calculation';receipt=run_command([str(apbs),'transfer.in'],c,c/'apbs.out',c/'resources.txt')
    write_new(c/'execution.json',receipt)
    if receipt['returncode']!=0: raise InvalidArtifact('APBS failed; raw output retained')
    result=collect(c/'apbs_manifest.json',c/'apbs.out');result['execution_receipt']=record(c/'execution.json')
    result['execution_cache_key']=cache_key({'preparation':record(c/'apbs_manifest.json'),'solver':record(apbs)})
    write_new(d/'result.json',result)
    return result


def transform(state,rotation=None,translation=None):
    s=copy.deepcopy(state); matrix=np.eye(3) if rotation is None else np.asarray(rotation)
    shift=np.zeros(3) if translation is None else np.asarray(translation)
    center=np.asarray(s['settings']['grid_center_A'])
    for field in ('physical_atoms','core_atoms','environment_atoms'):
        for a in s[field]: a['xyz_A']=(matrix@(np.asarray(a['xyz_A'])-center)+center+shift).tolist()
    # Keep the box fixed under sub-grid translation to measure grid placement
    # error, rather than making the test trivially identical by translating it.
    return s


def prior_allocated_cost(receipts):
    total=0;seen=set();records=[]
    for path in receipts:
        data=read_json(path);job=str(data['job_id'])
        if job in seen: raise InvalidArtifact('duplicate prior solver accounting job')
        seen.add(job)
        lines=data['sacct'].splitlines();header=lines[0].split('|')
        rows=[dict(zip(header,line.split('|'))) for line in lines[1:] if line.startswith(job+'|')]
        if len(rows)!=1 or rows[0]['State'] not in ('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED'):
            raise InvalidArtifact('prior solver terminal accounting unavailable')
        row=rows[0];cost=int(row['AllocCPUS'])*int(row['ElapsedRaw'])
        if cost!=int(row['CPUTimeRAW']): raise InvalidArtifact('inconsistent allocated accounting')
        total+=cost;records.append({'receipt':record(path),'job_id':job,'allocated_core_seconds':cost})
    return total,records


def frozen_schedule(states):
    schedule=[(name+'/primary',s) for name,s in states.items()]
    for name,s in states.items():
        refined=copy.deepcopy(s);refined['settings']=settings(s['physical_atoms'],spacing=.4)
        extended=copy.deepcopy(s);extended['settings']=settings(s['physical_atoms'],padding=30.)
        schedule.extend([(name+'/refined',refined),(name+'/extended',extended)])
        if name.startswith('1h4i_qm33'):
            theta=.37;rot=[[math.cos(theta),-math.sin(theta),0],[math.sin(theta),math.cos(theta),0],[0,0,1]]
            schedule.extend([(name+'/translated',transform(s,translation=[.173,.117,.231])),
                             (name+'/rotated',transform(s,rotation=rot))])
    anchor=states.get('1h4i_qm33_La')
    if anchor:
        identity=copy.deepcopy(anchor);identity['physical_atoms']=[dict(a,charge_e=0.) for a in identity['core_atoms']]
        identity['environment_atoms']=[];identity['expected_environment_charge_e']=0.
        schedule.extend([('identity/1h4i_qm33_La',identity),('repeat/1h4i_qm33_La',copy.deepcopy(anchor))])
    return schedule


def split_charging_input(text):
    """Independent charging solves; each retains the exact read/grid/boundary block."""
    import re
    first=text.index('elec name ');read=text[:first]
    blocks=re.findall(r'^elec name (\w+)\n(.*?)^end\s*$',text[first:],re.M|re.S)
    if len(blocks)!=6: raise InvalidArtifact('expected six charging blocks')
    return {name:read+f'elec name {name}\n'+body+'end\nprint elecEnergy 1 end\nquit\n' for name,body in blocks}


def prepare_completion(pilot_path,previous,output,apbs,agreement):
    from affordable_state import validate_skeleton_pair
    started=time.monotonic();previous=Path(previous).resolve();output=Path(output).resolve()
    m=read_json(pilot_path);old=read_json(previous/'solver_result.json')
    if old['pilot']!=record(pilot_path) or old['apbs']!=record(apbs):
        raise InvalidArtifact('previous result pilot or executable differs')
    states={};quality_records=[]
    for task in m['tasks']:
        qp=previous/'esp'/task['task_id']/'quality.json';q=read_json(qp)
        for key in ('actual_quantum_potential','points','gbw','density','utility','execution_receipt'): verify(q[key])
        if q['charges']['source_output']!=record(task['output_path']) or q['charges']['source_xyz']!=task['xyz']:
            raise InvalidArtifact('ESP source mismatch')
        quality_records.append(record(qp))
        sp=Path(pilot_path).resolve().parent.parent/'environment'/task['task_id']/'skeleton.json'
        if q['status']!='passed' or not sp.exists(): continue
        s=read_json(sp)
        if len(s['core_atoms'])!=len(q['charges']['charge_e']): raise InvalidArtifact('charge atom count mismatch')
        for a,charge in zip(s['core_atoms'],q['charges']['charge_e']): a['charge_e']=charge
        s['charge_quality']={'status':'passed','receipt':record(qp)};s['settings']=settings(s['physical_atoms'])
        states[task['task_id']]=s
    for name in ('1h4i_qm33','1h4i_qm36'): validate_skeleton_pair(states[name+'_La'],states[name+'_Ca'])
    scheduled={e['label']:e for e in read_json(previous/'numerical_schedule.json')['entries']}
    schedule=frozen_schedule(states)
    if set(scheduled)!=set(n for n,s in schedule): raise InvalidArtifact('scheduled state set changed')
    for label,s in schedule:
        if cache_key(s)!=scheduled[label]['state_cache_key']: raise InvalidArtifact('frozen state changed: '+label)
    output.mkdir(parents=True,exist_ok=False);entries=[];tasks=[]
    for label,s in schedule:
        prior=next(r for r in old['numerical_checks'] if r['label']==label)
        if prior['status']=='computed':
            rp=previous/label/'result.json';result=read_json(rp)
            prior_manifest=read_json(verify(result['source_manifest']))
            if cache_key(read_json(verify(prior_manifest['state'])))!=cache_key(s): raise InvalidArtifact('cached state differs')
            receipt=read_json(verify(result['execution_receipt']))
            if receipt['returncode']!=0: raise InvalidArtifact('cached execution failed')
            checked=collect(verify(result['source_manifest']),verify(result['output']))
            if checked['components']!=result['components']: raise InvalidArtifact('cached result differs from raw output')
            entries.append({'label':label,'state_cache_key':cache_key(s),'cached_result':record(rp)})
            continue
        directory=output/label;directory.mkdir(parents=True)
        sp=directory/'state.json';write_new(sp,s);manifest=prepare(sp,directory/'calculation')
        parts=split_charging_input(verify(manifest['input']).read_text())
        task_ids=[]
        for component,inp in parts.items():
            path=directory/'calculation'/f'{component}.in';path.write_text(inp)
            task_id=label+'/'+component;task_ids.append(task_id)
            tasks.append({'task_id':task_id,'label':label,'component':component,'input':record(path),
                          'output':str(path.with_suffix('.out')),'memory_estimate_bytes':math.prod(s['settings']['grid_dimensions'])*256})
        entries.append({'label':label,'state_cache_key':cache_key(s),'source_manifest':record(directory/'calculation/apbs_manifest.json'),'task_ids':task_ids})
    result={'schema_version':'alquemia.apbs_parallel_completion.v1','pilot':record(pilot_path),'previous_result':record(previous/'solver_result.json'),
            'previous_schedule':record(previous/'numerical_schedule.json'),'apbs':record(apbs),'agreement':record(agreement),
            'implementation':record(__file__),'environment_implementation':record(Path(__file__).with_name('affordable_environment.py')),
            'quality_receipts':quality_records,'preparation':old['preparation'],'entries':entries,'tasks':tasks,
            'limits':{'compute_budget':None,'wall_time_limit':None,'scheduling':'all manifested independent charging solves, limited concurrently only by CPUs and memory'},
            'preparation_wall_seconds':time.monotonic()-started}
    write_new(output/'completion_manifest.json',result)
    return result


def collect_charging(output):
    import re
    text=Path(output).read_text()
    values=re.findall(r'Global net ELEC energy\s*=\s*([-+\d.eE]+)\s+kJ/mol',text)
    if len(values)!=1 or 'Thanks for using APBS' not in text: raise InvalidArtifact('charging output incomplete')
    value=float(values[0])
    if not math.isfinite(value): raise InvalidArtifact('nonfinite charging energy')
    return value


def completion_preflight(path):
    m=read_json(path)
    for key in ('pilot','previous_result','previous_schedule','apbs','agreement','implementation','environment_implementation'): verify(m[key])
    for rec in m['quality_receipts']: verify(rec)
    if m['limits']['compute_budget'] is not None or m['limits']['wall_time_limit'] is not None:
        raise InvalidArtifact('completion policy must not impose compute or time budgets')
    for entry in m['entries']:
        if 'cached_result' in entry: verify(entry['cached_result']);continue
        p=read_json(verify(entry['source_manifest']));verify(p['state']);verify(p['input'])
        for rec in p['pqr'].values(): verify(rec)
    for t in m['tasks']: verify(t['input'])
    return m


def execute_completion(path):
    from affordable_environment import transfer_components
    from concurrent.futures import as_completed
    if not os.environ.get('SLURM_JOB_ID'): raise InvalidArtifact('requires SLURM allocation')
    m=completion_preflight(path);root=Path(path).resolve().parent
    campaign=Path(m['pilot']['path']).parent.parent
    lock=(campaign/'solver_campaign.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    cpus=int(os.environ['SLURM_CPUS_ON_NODE']);started=time.monotonic()
    available_kib=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))
    largest=max(t['memory_estimate_bytes'] for t in m['tasks'])
    workers=min(len(m['tasks']),cpus,max(1,int(available_kib*1024*.75)//largest))
    implementation=snapshot_implementation(root/f"implementation_{os.environ['SLURM_JOB_ID']}")
    write_new(root/f"admission_{os.environ['SLURM_JOB_ID']}.json",{'tasks':len(m['tasks']),'workers':workers,'allocated_cpus':cpus,
               'available_memory_bytes':available_kib*1024,'largest_task_memory_estimate_bytes':largest,
               'compute_budget':None,'wall_time_limit':None,'implementation':implementation})
    apbs=verify(m['apbs'])
    def one(t):
        inp=verify(t['input']);op=Path(t['output']);receipt_path=op.with_suffix('.execution.json')
        try:
            if receipt_path.exists():
                r=read_json(receipt_path)
                if r.get('input')!=t['input'] or r.get('solver')!=m['apbs'] or r['returncode']!=0: raise InvalidArtifact('cached charging receipt mismatch')
                verify(r['log']);value=collect_charging(op)
                return {'task_id':t['task_id'],'status':'computed','value_kJ_mol':value,'execution_receipt':record(receipt_path),'cache_reused':True}
            if op.exists(): raise InvalidArtifact('partial charging output retained; no overwrite')
            r=run_command([str(apbs),inp.name],inp.parent,op,op.with_suffix('.resources.txt'))
            r.update(input=t['input'],solver=m['apbs']);write_new(receipt_path,r)
            if r['returncode']!=0: raise InvalidArtifact('APBS charging failed')
            return {'task_id':t['task_id'],'status':'computed','value_kJ_mol':collect_charging(op),'execution_receipt':record(receipt_path),'cache_reused':False}
        except Exception as exc: return {'task_id':t['task_id'],'status':'failed','reason':str(exc)}
    results={};events=root/f"charging_events_{os.environ['SLURM_JOB_ID']}.jsonl"
    with ThreadPoolExecutor(max_workers=workers) as pool,events.open('x') as log:
        for future in as_completed([pool.submit(one,t) for t in m['tasks']]):
            r=future.result();results[r['task_id']]=r;log.write(json.dumps(r)+'\n');log.flush()
    numerical=[]
    for entry in m['entries']:
        if 'cached_result' in entry:
            numerical.append({'label':entry['label'],'status':'computed','cache_reused':True,'result':read_json(verify(entry['cached_result']))});continue
        rows=[results[t] for t in entry['task_ids']]
        if any(r['status']!='computed' for r in rows):
            numerical.append({'label':entry['label'],'status':'failed','reason':'one or more charging solves failed','charging_tasks':rows,'result':None});continue
        prep=read_json(verify(entry['source_manifest']))
        charging={t.rsplit('/',1)[1]:results[t]['value_kJ_mol'] for t in entry['task_ids']}
        result={'protocol_id':prep['protocol_id'],'source_manifest':entry['source_manifest'],'charging_energies_kJ_mol':charging,
                'components':transfer_components(charging,prep['direct_coulomb_kcal_mol']),
                'physical_boundary_hash':prep['physical_boundary_hash'],'charging_tasks':rows,
                'physical_validation_status':'not_yet_validated','decision':'uncalibrated_protocol'}
        rp=root/entry['label']/'result.json';write_new(rp,result)
        numerical.append({'label':entry['label'],'status':'computed','cache_reused':False,'result':result})
    elapsed=time.monotonic()-started
    summary={'schema_version':'alquemia.apbs_physical_pilot.v1','pilot':m['pilot'],'apbs':m['apbs'],'completion_manifest':record(path),
             'preparation':m['preparation'],'numerical_checks':numerical,'slurm_job_id':os.environ['SLURM_JOB_ID'],
             'allocated_core_seconds':elapsed*cpus,'allocated_cpus':cpus,'elapsed_seconds':elapsed,'workers':workers,
             'budget_exceeded':False,'compute_budget':None,'wall_time_limit':None,'prior_costs_retained_in':m['previous_result'],
             'predictive_claim':'none_development_only','actual_implementation':implementation}
    write_new(root/'solver_result.json',summary)
    from affordable_compare import compare,report
    data=compare(Path(__file__).resolve().parents[1],root/'solver_result.json')
    write_new(root/'comparison.json',data)
    with (root/'REPORT.md').open('x') as f:f.write(report(data))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare-completion')
    for name in ('pilot-manifest','previous','output','apbs','agreement'):q.add_argument('--'+name,type=Path,required=True)
    for name in ('dry-run-completion','execute-completion'):
        q=sub.add_parser(name);q.add_argument('--manifest',type=Path,required=True)
    a=p.parse_args()
    if a.op=='prepare-completion':
        r=prepare_completion(a.pilot_manifest,a.previous,a.output,a.apbs.resolve(),a.agreement)
        print(f"{len(r['entries'])} states; {len(r['tasks'])} independent charging solves; no compute/time budget")
    elif a.op=='dry-run-completion':
        r=completion_preflight(a.manifest);print(f"verified {len(r['tasks'])} charging tasks")
    else:
        r=execute_completion(a.manifest);print(json.dumps({'states':len(r['numerical_checks']),'allocated_core_seconds':r['allocated_core_seconds']}))


if __name__=='__main__':main()
