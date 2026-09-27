"""Finite classical additive component checks on pinned real 1H4I.

prepare/dry-run never create a Context. execute requires a Slurm allocation.
No quantum calls, solvent model, optimizer, or combined hybrid qualification.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
import resource
from pathlib import Path
import shutil
import sys
import time
import traceback
import numpy as np
from affordable_common import read_json, record, verify, write_new

PROTOCOL='nikasha_1h4i_classical_component_checks_v1'
GROUPS=('retained_bonded','MM_LJ','MM_Coulomb','QM_MM_LJ')
TOLERANCES={'repeat_energy_kcal_mol':1e-8,'rigid_energy_kcal_mol':1e-6,
            'rigid_force_kcal_mol_A':1e-5,'fd_absolute_kcal_mol_A':.001,'fd_relative':.0001}
CONTEXTS={}


def rotation():
    c,s=np.cos(.37),np.sin(.37)
    return np.array([[c,-s,0],[s,c,0],[0,0,1.]])


def load_ledger(path):
    r=read_json(path)
    if r['protocol_id']!='nikasha_1h4i_additive_mechanics_ledger_v1':
        raise ValueError('unsupported ledger')
    for pin in r['inputs'].values():verify(pin)
    verify(r['native_system']);verify(r['forcefield'])
    a={k:read_json(verify(v)) for k,v in r['artifacts'].items()}
    return r,a


def build_system(ledger, artifacts, metal):
    """Native OpenMM force construction only; no Context or energy call."""
    import openmm as mm
    from openmm import unit
    from openmm.app import element
    parent=mm.XmlSerializer.deserialize(verify(ledger['native_system']).read_text())
    particles=artifacts['particles.json']; n=len(particles)
    system=mm.System()
    for i,a in enumerate(particles):
        system.addParticle(parent.getParticleMass(i) if i<parent.getNumParticles()
                           else element.get_by_symbol(metal if a.get('kind')=='metal' else a['element']).mass)
    selected={}
    for t in artifacts['bonded_terms.json']:
        selected.setdefault(t['force'],set()).add(t['index'])
    counts={}
    for force in parent.getForces():
        name=type(force).__name__
        if name in ('NonbondedForce','CMMotionRemover'):continue
        if name not in ('HarmonicBondForce','HarmonicAngleForce','PeriodicTorsionForce','CMAPTorsionForce'):
            raise ValueError('unexpected native force '+name)
        new=type(force)();keep=selected.get(name,set())
        if name=='CMAPTorsionForce':
            for j in range(force.getNumMaps()):new.addMap(*force.getMapParameters(j))
            for j in sorted(keep):new.addTorsion(*force.getTorsionParameters(j))
        elif name=='HarmonicBondForce':
            for j in sorted(keep):new.addBond(*force.getBondParameters(j))
        elif name=='HarmonicAngleForce':
            for j in sorted(keep):new.addAngle(*force.getAngleParameters(j))
        else:
            for j in sorted(keep):new.addTorsion(*force.getTorsionParameters(j))
        new.setForceGroup(0);new.setName('retained_'+name);system.addForce(new);counts[name]=len(keep)
    lj=mm.NonbondedForce();coul=mm.NonbondedForce()
    for group,f in ((1,lj),(2,coul)):
        f.setForceGroup(group);f.setName(GROUPS[group]);f.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        f.setUseDispersionCorrection(False)
    cross=mm.CustomNonbondedForce('4*epsilon*((sigma/r)^12-(sigma/r)^6);sigma=(sigma1+sigma2)/2;epsilon=sqrt(epsilon1*epsilon2)')
    cross.addPerParticleParameter('sigma');cross.addPerParticleParameter('epsilon')
    cross.setForceGroup(3);cross.setName('QM_MM_LJ_generic');cross.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    qset=set();mset=set()
    for i,a in enumerate(particles):
        pars=a['endpoint_lj'][metal] if a.get('kind')=='metal' else a
        sig,eps=pars['sigma_nm'],pars['epsilon_kJ_mol']
        ismm=a['region']=='MM'
        (mset if ismm else qset).add(i)
        lj.addParticle(0.,sig,eps if ismm else 0.)
        coul.addParticle(a['charge_e'] if ismm else 0.,1.,0.)
        cross.addParticle([sig,eps])
    cross.addInteractionGroup(qset,mset)
    replacements=mm.CustomBondForce('4*epsilon*((sigma/r)^12-(sigma/r)^6)')
    replacements.addPerBondParameter('sigma');replacements.addPerBondParameter('epsilon')
    replacements.setForceGroup(3);replacements.setName('QM_MM_LJ_native_exceptions')
    for e in artifacts['exceptions.json']:
        a,b=e['atoms'];ismm=e['region']=='MM'
        lj.addException(a,b,0.,e['sigma_nm'],e['epsilon_kJ_mol'] if ismm else 0.)
        coul.addException(a,b,e['chargeprod_e2'] if ismm else 0.,1.,0.)
        cross.addExclusion(a,b)
        if e['region']=='cross':replacements.addBond(a,b,[e['sigma_nm'],e['epsilon_kJ_mol']])
    for f in (lj,coul,cross,replacements):system.addForce(f)
    if system.getNumConstraints()!=0:raise ValueError('energy diagnostic must have no constraints')
    return system,dict(particle_count=n,real_QM=len(qset),MM=len(mset),bonded_counts=counts,
                       MM_exception_count=lj.getNumExceptions(),cross_replacements=replacements.getNumBonds(),
                       constraints=0,energy_force_calls=0)


def geometry_pool(ledger, artifacts):
    p=artifacts['particles.json'];ids={a['id']:i for i,a in enumerate(p)}
    x=np.array([a['xyz_A'] for a in p]);y=x.copy()
    for a in read_json(verify(ledger['environments']['B']['atoms'])):
        i=ids[a['id']]
        if p[i]['charge_e']!=a['charge_e']:raise ValueError('A/B charge differs')
        y[i]=a['xyz_A']
    hi,oi,ci=(ids['A/159/ /'+k] for k in ('HG1','OG1','CB'))
    axis=x[oi]-x[ci];axis/=np.linalg.norm(axis)
    tangent=np.cross(axis,x[hi]-x[oi]);radius=np.linalg.norm(tangent);tangent/=radius
    bi,ri=ids['A/177/ /CA'],ids['A/177/ /CB']
    u=x[bi]-x[ri];u/=np.linalg.norm(u)
    bd=np.eye(3)[np.argmin(abs(u))];bd-=u*np.dot(u,bd);bd/=np.linalg.norm(bd)
    mi=ids['metal'];md=x[ids['A/177/ /OE1']]-x[mi];md/=np.linalg.norm(md)
    modes={'hydroxyl':{'physical_index':hi,'direction':tangent.tolist(),'axis':axis.tolist(),'pivot':oi,'radius_A':float(radius)},
           'boundary':{'physical_index':bi,'direction':bd.tolist()},
           'metal':{'physical_index':mi,'direction':md.tolist()}}
    pool={'A':x,'B':y,'repeat':x.copy(),'rigid':x@rotation().T+(.173,.117,.231)}
    for name,m in modes.items():
        for h in (-.001,.001,-.0005,.0005):
            z=x.copy();i=m['physical_index'];direction=np.array(m['direction'])
            if name=='hydroxyl':
                v=x[i]-x[oi];angle=h/radius
                z[i]=x[oi]+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+axis*np.dot(axis,v)*(1-np.cos(angle))
            else:z[i]+=h*direction
            pool[f'{name}_{h:+.4f}']=z
    return pool,modes


def prepare(ledger_path, agreement, output):
    import openmm as mm
    ledger,artifacts=load_ledger(ledger_path)
    pool,modes=geometry_pool(ledger,artifacts)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir()
    for name in (Path(__file__).name,'affordable_common.py'):
        shutil.copyfile(Path(__file__).parent/name,impl/name)
    models={};counts={}
    for metal in ('Ca','La'):
        system,counts[metal]=build_system(ledger,artifacts,metal)
        path=out/(metal+'.xml');path.write_text(mm.XmlSerializer.serialize(system));models[metal]=record(path)
    coords={}
    for name,x in pool.items():
        path=out/(name+'.npy');np.save(path,x);coords[name]=record(path)
    tasks=[dict(task_id=metal+'_'+name,metal=metal,configuration=name,coordinates=coords[name]) for metal in ('Ca','La') for name in pool]
    manifest=dict(protocol_id=PROTOCOL,ledger=record(ledger_path),agreement=record(agreement),platform='Reference',
                  software=dict(openmm=mm.version.full_version,python=sys.version,executable=sys.executable),
                  models=models,modes=modes,groups=list(GROUPS),tolerances=TOLERANCES,tasks=tasks,counts=counts,
                  implementation={p.name:record(p) for p in impl.iterdir()},
                  limitations=['finite dry classical components only','no quantum calls or optimization','full solvent and hybrid qualification unavailable'])
    write_new(out/'manifest.json',manifest)
    return {'manifest':record(out/'manifest.json'),'tasks':len(tasks),'energy_force_calls':0}


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['tolerances']!=TOLERANCES or m['platform']!='Reference':
        raise ValueError('protocol/tolerance/platform mismatch')
    verify(m['ledger']);verify(m['agreement'])
    for p in (*m['models'].values(),*m['implementation'].values()):verify(p)
    if len(m['tasks'])!=32 or len({t['task_id'] for t in m['tasks']})!=32:raise ValueError('requires full 32-cell manifest')
    for t in m['tasks']:verify(t['coordinates'])
    return m


def evaluate_one(args):
    manifest,task=args
    import openmm as mm
    from openmm import unit
    m=read_json(manifest);out=Path(manifest).parent/'results'/task['task_id'];out.mkdir(parents=True,exist_ok=True)
    receipt=out/'receipt.json';pin=record(manifest)
    if receipt.exists():
        old=read_json(receipt)
        if old.get('manifest')!=pin:raise ValueError('cached manifest mismatch')
        if old['status']=='complete':
            verify(old['forces']);return old
        raise ValueError('failed task preserved; use an explicitly new attempt directory')
    start=time.monotonic();cpu_start=time.process_time()
    row=dict(task_id=task['task_id'],manifest=pin,status='failed',platform='Reference',
             job_id=os.environ['SLURM_JOB_ID'],allocated_cpus=os.environ.get('SLURM_CPUS_ON_NODE'),
             openmm=mm.version.full_version)
    try:
        key=m['models'][task['metal']]['sha256']
        reused=key in CONTEXTS
        if not reused:
            system=mm.XmlSerializer.deserialize(verify(m['models'][task['metal']]).read_text())
            integrator=mm.VerletIntegrator(.001)
            context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
            CONTEXTS[key]=(context,integrator)
        context=CONTEXTS[key][0]
        x=np.load(verify(task['coordinates']));context.setPositions(x*.1*unit.nanometer)
        energies={};forces=[];component_wall={}
        for index,name in enumerate(GROUPS):
            component_start=time.monotonic()
            state=context.getState(getEnergy=True,getForces=True,groups=1<<index)
            energies[name]=state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            forces.append(state.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom))
            component_wall[name]=time.monotonic()-component_start
        f=np.asarray(forces)
        if not np.isfinite(f).all() or not np.isfinite(list(energies.values())).all():raise ValueError('nonfinite components')
        fp=out/'forces.npz';np.savez_compressed(fp,forces_kcal_mol_A=f)
        row.update(status='complete',energies_kcal_mol=energies,forces=record(fp),warm_context_reused=reused,component_wall_seconds=component_wall)
    except Exception as exc:row.update(failure=str(exc),traceback=traceback.format_exc())
    row['wall_seconds']=time.monotonic()-start
    row['task_cpu_seconds']=time.process_time()-cpu_start
    row['process_peak_rss_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    row['cpu_seconds_process_lifetime']=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime
    row['task']=task
    row['model']=m['models'][task['metal']]
    write_new(receipt,row)
    return row


def execute(manifest, workers, task_id=None):
    if not os.environ.get('SLURM_JOB_ID'):raise ValueError('molecular evaluations require Slurm allocation')
    m=validate(manifest)
    allocation=int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0])
    if workers<1 or workers>allocation:raise ValueError('workers exceed actual allocated CPUs')
    import openmm as mm
    if mm.version.full_version!=m['software']['openmm']:raise ValueError('OpenMM version differs')
    tasks=[t for t in m['tasks'] if task_id is None or t['task_id']==task_id]
    if not tasks:raise ValueError('unknown task ID')
    args=[(str(Path(manifest).resolve()),t) for t in tasks]
    if workers==1:rows=[evaluate_one(a) for a in args]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:rows=list(pool.map(evaluate_one,args))
    return {'completed':sum(r['status']=='complete' for r in rows),'requested':len(rows)}


def collect(manifest):
    m=validate(manifest);root=Path(manifest).parent;rows={};checks=[]
    for t in m['tasks']:
        rp=root/'results'/t['task_id']/'receipt.json'
        r=read_json(rp) if rp.exists() else {'status':'unavailable'}
        if r['status']=='complete':
            if r['manifest']!=record(manifest):raise ValueError('receipt identity differs')
            verify(r['forces'])
        rows[t['task_id']]=r
    def force(row):return np.load(verify(row['forces']))['forces_kcal_mol_A']
    def check(name,value,limit):checks.append(dict(name=name,residual=float(value),tolerance=float(limit),passed=bool(value<=limit)))
    responses={}
    for metal in ('Ca','La'):
        base=rows[metal+'_A']
        if base['status']!='complete':continue
        f0=force(base)
        for mode in ('repeat','rigid'):
            row=rows[metal+'_'+mode]
            if row['status']!='complete':continue
            ff=force(row);expected=f0@rotation().T if mode=='rigid' else f0
            for gi,g in enumerate(GROUPS):
                check(f'{metal}/{mode}/{g}/energy',abs(row['energies_kcal_mol'][g]-base['energies_kcal_mol'][g]),TOLERANCES[mode+'_energy_kcal_mol'])
                check(f'{metal}/{mode}/{g}/forces',np.max(abs(ff[gi]-expected[gi])),TOLERANCES['rigid_force_kcal_mol_A'])
        for mode,desc in m['modes'].items():
            index=desc['physical_index'];direction=np.array(desc['direction'])
            for gi,g in enumerate(GROUPS):
                analytic=float(-f0[gi,index]@direction);fds={}
                limit=TOLERANCES['fd_absolute_kcal_mol_A']+TOLERANCES['fd_relative']*abs(analytic)
                for h in (.001,.0005):
                    plus=rows[f'{metal}_{mode}_{h:+.4f}'];minus=rows[f'{metal}_{mode}_{-h:+.4f}']
                    if plus['status']!='complete' or minus['status']!='complete':continue
                    fd=(plus['energies_kcal_mol'][g]-minus['energies_kcal_mol'][g])/(2*h);fds[h]=fd
                    check(f'{metal}/{mode}/{g}/FD/{h}',abs(fd-analytic),limit)
                    checks[-1].update(analytic_kcal_mol_A=analytic,finite_difference_kcal_mol_A=fd)
                if len(fds)==2:check(f'{metal}/{mode}/{g}/refinement',abs(fds[.001]-fds[.0005]),limit)
        other=rows[metal+'_B']
        if other['status']=='complete':responses[metal]={g:other['energies_kcal_mol'][g]-base['energies_kcal_mol'][g] for g in GROUPS}
    for name in {t['configuration'] for t in m['tasks']}:
        ca,la=rows['Ca_'+name],rows['La_'+name]
        if ca['status']!='complete' or la['status']!='complete':continue
        fc,fl=force(ca),force(la)
        for gi,g in enumerate(GROUPS[:3]):
            check(f'paired/{name}/{g}/energy',abs(ca['energies_kcal_mol'][g]-la['energies_kcal_mol'][g]),TOLERANCES['repeat_energy_kcal_mol'])
            check(f'paired/{name}/{g}/forces',np.max(abs(fc[gi]-fl[gi])),TOLERANCES['rigid_force_kcal_mol_A'])
    complete=sum(r['status']=='complete' for r in rows.values())
    result=dict(status='passed_classical_component_checks' if complete==32 and len(checks)==200 and all(c['passed'] for c in checks) else 'incomplete_or_failed',
                manifest=record(manifest),rows=rows,checks=checks,responses_kcal_mol=responses,
                completed=complete,required=32,full_hybrid_qualified=False,solvent_present=False)
    if len(responses)==2:result['La_minus_Ca_response_kcal_mol']={g:responses['La'][g]-responses['Ca'][g] for g in GROUPS}
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare')
    for k in ('ledger','agreement','output'):q.add_argument('--'+k,required=True,type=Path)
    q=sub.add_parser('dry-run');q.add_argument('--manifest',required=True,type=Path)
    q=sub.add_parser('execute');q.add_argument('--manifest',required=True,type=Path);q.add_argument('--workers',type=int,required=True);q.add_argument('--task-id')
    q=sub.add_parser('collect');q.add_argument('--manifest',required=True,type=Path);q.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.ledger,a.agreement,a.output)
    elif a.command=='dry-run':
        m=validate(a.manifest);r={'status':'ready_no_energies','task_count':len(m['tasks']),'platform':m['platform']}
    elif a.command=='execute':r=execute(a.manifest,a.workers,a.task_id)
    else:r=collect(a.manifest);write_new(a.output,r);r={k:v for k,v in r.items() if k not in ('rows','checks')}
    print(json.dumps(r,indent=2))

if __name__=='__main__':main()
