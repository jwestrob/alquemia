"""Opt-in frozen-density, AMOEBA-GK proxy, and MACE short-context pilot."""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,energy,read_json,record,verify,write_new
from mace_hybrid import EV_TO_KCAL,rotation,accepted_attempt
from mace_density_multipoles import coupling,parse_moments,validate as validate_density
from mace_density_gk_boundary import CASES,validate as validate_boundary
from mace_native_field_input import derive as derive_input
from mace_tinker_framework_solver import ENERGY_NAMES

PROTOCOL='vacuum_density_direct_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1'
TOL=dict(refinement_kcal=.01,rigid_kcal=.01,identity_kcal=1e-8,algebra_kcal=1e-7,
         radius_relative_kcal=1.,partition_kcal=2.,ordering_kcal=.02,large_dipole_eA=1.)
VARIANTS=dict(primary=1.,rigid=1.,radius_minus=.95,radius_plus=1.05)
NULLS=dict(reference=None,calibrated_class=None,aqueous_affinity_score=None,combined_gradient=None,
           relaxation_correction=None,response_status='response_model_not_validated',baseline_changed=False)


def derive(text):
    original,supplied=derive_input(text)
    out=supplied.replace('subroutine alquemia_induce0c (source_fields)','subroutine alquemia_induce0c (source_fields,eps_final)',1)
    out=out.replace('      real*8 source_fields(3,n,4)','      real*8 source_fields(3,n,4),eps_final',1)
    # The final diagnostic copies the native value without changing any equation.
    i=out.rfind('      return\n')
    if i<0:raise InvalidArtifact('native final return unavailable')
    out=out[:i]+'      eps_final = eps\n'+out[i:]
    return original,out


def build(parent_software,native_source,frontend,output):
    p=read_json(parent_software);lib=verify(p['library']);root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    if p['returncode']:raise InvalidArtifact('parent solver unavailable')
    shutil.copyfile(native_source,root/'original_induce.f');original,derived=derive(Path(native_source).read_text())
    (root/'original_induce0c.f').write_text(original);(root/'supplied_induce0c.f').write_text(derived)
    shutil.copyfile(frontend,root/Path(frontend).name);shutil.copyfile(__file__,root/Path(__file__).name)
    cmd=[p['command'][0],'-O2','-g','-fno-fast-math','-fopenmp','-I'+str(lib.parent/'mod'),str(root/Path(frontend).name),
        str(root/'supplied_induce0c.f'),str(lib),'-lfftw3_threads','-lfftw3','-o',str(root/'density_gk')]
    before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
    with (root/'build.log').open('x') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result=dict(command=cmd,returncode=r.returncode,library=p['library'],parent_software=record(parent_software),
        original_source=record(root/'original_induce.f'),original_routine=record(root/'original_induce0c.f'),
        derived_routine=record(root/'supplied_induce0c.f'),frontend=record(root/Path(frontend).name),
        implementation=record(root/Path(__file__).name),compiler=p['compiler'],log=record(root/'build.log'),
        wall_seconds=time.monotonic()-start,child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        peak_child_RSS_KiB=after.ru_maxrss,new_scientific_calls=0)
    if not r.returncode:result['executable']=record(root/'density_gk')
    write_new(root/'receipt.json',result);return result


def load_terms(path,bm):
    r=read_json(path);m=read_json(verify(r['manifest']))
    if record(path)['sha256']!='3e9e6bc11464157d1bb1b2df9151f0db97bcdfec17e8df4b8432267e583fca23' or r['status']!='complete':
        raise InvalidArtifact('declared vacuum/short-component archive required')
    reuse=read_json(verify(m['reuse_audit']));out={}
    for case in CASES:
        gid='GGR_1GLG' if case.startswith('GGR_') else case
        for metal in ('Ca','La'):
            tid=case+'_'+metal;t=next(t for t in m['tasks'] if t['task_id']==tid)
            bt=next(t for t in bm['tasks'] if t['case_id']==case and t['metal']==metal and t['mode']=='source')
            if t['xyz']!=bt['source_QM_xyz'] or t['source_receipt']!=bt['source_quantum_receipt']:
                raise InvalidArtifact('reused energy geometry/quantum state differs')
            short=r['rows'][tid]['short'];mp=verify(r['manifest'])
            if not any(accepted_attempt(a,t,mp)==short for a in (mp.parent/'execution'/tid).glob('attempt_*')):
                raise InvalidArtifact('missing actual core MACE receipt')
            whole=reuse['whole_reuses'][gid][metal];c=read_json(verify(whole['collection']));cm=read_json(verify(c['manifest']))
            wt=next(x for x in cm['tasks'] if x['task_id']==whole['task_id']);wr=c['rows'][whole['task_id']]
            if not any(accepted_attempt(a,wt,verify(c['manifest']))==wr for a in (Path(c['manifest']['path']).parent/'execution'/whole['task_id']).glob('attempt_*')):
                raise InvalidArtifact('missing actual whole MACE receipt')
            dft=energy(verify(t['source_output']));core=short['energy_eV'];full=whole['short_energy_eV']
            old=r['cases'][case]['endpoints'][metal]
            if dft*HA_TO_KCAL!=old['DFT_vacuum_kcal_mol'] or core!=old['short_core_eV'] or full!=old['short_full_eV']:
                raise InvalidArtifact('component values disagree with actual archive')
            out[tid]=dict(DFT_vacuum_hartree=dft,short_core_eV=core,short_full_eV=full,
                short_context_kcal=(full-core)*EV_TO_KCAL,source_output=t['source_output'],source_receipt=t['source_receipt'],
                core_receipt=short,whole_reuse=whole,evidence=r['cases'][case]['evidence'])
    return out


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},protocol=PROTOCOL,sources=m['sources'],
        software=m['software'],implementation=m['implementation'],plan=m['plan'],tolerances=TOL))


def prepare(boundary,density,short_report,software,plan,output):
    start=time.monotonic();cpu=time.process_time();br=read_json(boundary);bm=validate_boundary(verify(br['manifest']))
    dr=read_json(density);dm=validate_density(verify(dr['manifest']));sw=read_json(software)
    if not br['complete'] or not br['gates_pass'] or not dr['complete'] or not dr['numerical_pass'] or sw['returncode']:
        raise InvalidArtifact('boundary/density/software prerequisite failed')
    if sw['library']!=read_json(verify(bm['software']))['library']:raise InvalidArtifact('native library changed')
    terms=load_terms(short_report,bm);root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    impl=root/'implementation';impl.mkdir()
    for name,pin in bm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name in ('mace_native_field_input.py','mace_density_gk_hybrid.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    write_new(root/'reused_terms.json',terms)
    m=dict(protocol=PROTOCOL,sources=dict(boundary=record(boundary),boundary_manifest=br['manifest'],density=record(density),
        density_manifest=dr['manifest'],short_report=record(short_report)),software=record(software),plan=record(plan),
        implementation={p.name:record(p) for p in sorted(impl.glob('*.py'))},tolerances=TOL,threads=64,terms=record(root/'reused_terms.json'),
        static_tasks=[],response_tasks=[],new_energy_calls=51,new_field_queries=48,new_response_solves=60,new_DFT_calls=0,new_MACE_calls=0,**NULLS)
    def task(case,state,variant,mode,eps):
        metal='Ca' if state=='environment' else state
        old=next(t for t in bm['tasks'] if t['case_id']==case and t['metal']==metal and t['mode']==('zero_source' if state=='environment' else 'source'))
        r=read_json(verify(br['tasks'][old['task_id']]['result']));meta=read_json(verify(old['boundary']))
        rid=mode+'_'+variant+'_'+case+'_'+state+('_standard' if eps==1e-7 else '')
        d=root/'tasks'/rid;d.mkdir(parents=True)
        rotate=variant=='rigid';mat=rotation() if rotate else np.eye(3);shift=np.array([17.,-23.,9.]) if rotate else np.zeros(3)
        coords=np.array(meta['positions_A'])@mat.T+shift
        lines=verify(old['xyz']).read_text().splitlines();rows=[lines[0]]
        for line,p in zip(lines[1:],coords):
            f=line.split();rows.append(' '.join(f[:2]+[format(v,'.17g') for v in p]+f[5:]))
        (d/'framework.xyz').write_text('\n'.join(rows)+'\n')
        radius=VARIANTS.get(variant,1.);diameter=format(3.6497*radius,'.17g');opts=[]
        for line in verify(old['key']).read_text().splitlines():
            f=line.split()
            if f[0]=='polar-eps':line='polar-eps '+str(eps)
            if f[0]=='solute' and f[1] in ('4998','4999'):
                f[4]=diameter;line=' '.join(f)
            opts.append(line)
        (d/'framework.key').write_text('\n'.join(opts)+'\n')
        for k in ('mask','overrides'):shutil.copyfile(verify(old[k]),d/Path(old[k]['path']).name)
        dt=next(t for t in dm['tasks'] if t['task_id']==case+'_'+metal)
        t=dict(task_id=rid,case_id=case,state=state,metal=metal,variant=variant,mode=mode,poleps=eps,
            rotation=mat.tolist(),translation_A=shift.tolist(),metal_GK_diameter=diameter,
            xyz=record(d/'framework.xyz'),key=record(d/'framework.key'),mask=record(d/'framework.freeze'),overrides=record(d/'framework.charges'),
            prepared_result=br['tasks'][old['task_id']]['result'],boundary=old['boundary'],frozen_indices=old['frozen_indices'],
            density_task=dt,density_arrays=dr['rows'][dt['task_id']]['arrays'])
        t['cache_key']=key(t,m);return t
    for var in VARIANTS:
        for case in CASES:
            for state in ('environment','Ca','La'):
                m['static_tasks'].append(task(case,state,var,'static',1e-9))
                m['response_tasks'].append(task(case,state,var,'solve',1e-9))
                if var=='primary':m['response_tasks'].append(task(case,state,var,'solve',1e-7))
    for state in ('environment','Ca','La'):m['static_tasks'].append(task('GGR_extended',state,'identity','identity',1e-9))
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['tolerances']!=TOL or len(m['static_tasks'])!=51 or len(m['response_tasks'])!=60:
        raise InvalidArtifact('changed hybrid pilot inventory/settings')
    for pin in [m['software'],m['plan'],m['terms'],*m['sources'].values(),*m['implementation'].values()]:verify(pin)
    sw=read_json(verify(m['software']))
    for k in ('executable','library','original_source','original_routine','derived_routine','frontend','implementation','log'):verify(sw[k])
    a,b=derive(verify(sw['original_source']).read_text())
    if a!=verify(sw['original_routine']).read_text() or b!=verify(sw['derived_routine']).read_text():raise InvalidArtifact('native solver equations changed')
    for t in m['static_tasks']+m['response_tasks']:
        if t['cache_key']!=key(t,m):raise InvalidArtifact('changed scientific task key')
        for k in ('xyz','key','mask','overrides','prepared_result','boundary','density_arrays'):verify(t[k])
        for k in ('probes','multipole_inputs','field_arrays'):verify(t['density_task'][k])
    return m


def parse(text,t):
    mode=t['mode']
    if 'ALQUEMIA_DENSITY_GK_COMPLETE '+mode not in text or 'not Converged' in text:raise InvalidArtifact('incomplete/nonconverged hybrid native call')
    marker='ALQUEMIA_MOMENTS_COMPLETE';prefix=text[:text.index(marker)]+marker
    moments=parse_moments(prefix,None)
    if not moments['pass_']:raise InvalidArtifact('native moment rotation failed')
    p=moments['parameters'];old=read_json(verify(t['prepared_result']));expect=copy.deepcopy(old['parameters'])
    expect['inventory'][-1]=64;expect['global'][5]=t['poleps']
    coords=[[float(v) for v in line.split()[2:5]] for line in verify(t['xyz']).read_text().splitlines()[1:]]
    for atom,pos in zip(expect['atoms'],coords):atom['xyz_A']=pos
    expect['atoms'][-1]['radius_A']=.5*float(t['metal_GK_diameter'])
    if mode=='identity':
        expect['inventory'][2]=0
        for i,atom in enumerate(expect['atoms'],1):
            atom.update(polarizability_A3=0.,damping_A_half=0.,response_allowed=False)
            if i not in t['frozen_indices']:atom['charge_e']=0.
    if p!=expect:raise InvalidArtifact('runtime native parameters differ from declared geometry/state')
    mat=np.array(t['rotation']);moment_error=0.
    for i,(actual,base) in enumerate(zip(moments['atoms'],old['atoms']),1):
        v=np.array(base['global_']);expected=np.r_[v[0],mat@v[1:4],(mat@v[4:].reshape(3,3)@mat.T).ravel()]
        if mode=='identity' and i not in t['frozen_indices']:expected[:]=0.
        moment_error=max(moment_error,float(np.max(abs(expected-np.array(actual['global_'])))))
    if moment_error>1e-10:raise InvalidArtifact('permanent source/environment moments changed')
    fields={};mu={};energies=None;eps=None;timing=None;constants=None;active=None
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0] in ('FIELD','RESPONSE','VACUUM_RESPONSE'):
            dest=fields if f[0]=='FIELD' else mu;i=int(f[1]);width=14 if f[0]=='VACUUM_RESPONSE' else 15
            if i in dest or len(f)!=width:raise InvalidArtifact('native field/response columns or duplicate')
            dest[i]=list(map(float,f[2:]))
        elif f[0]=='ENERGIES':energies=dict(zip(ENERGY_NAMES,map(float,f[1:])))
        elif f[0]=='EXACT_FINAL_RMS_DEBYE':eps=float(f[1])
        elif f[0]=='KERNEL_TIMING':timing=list(map(float,f[1:]))
        elif f[0]=='CONSTANTS':constants=list(map(float,f[1:]))
        elif f[0]=='ACTIVE_TERMS':active=f[1:]
    n=p['inventory'][0];ids=set(range(1,n+1))
    if set(mu)!=ids or constants is None or timing is None or active is None:raise InvalidArtifact('missing native calculation evidence')
    if mode=='solve':
        if eps is None or eps>t['poleps'] or eps<0 or energies is not None or fields:raise InvalidArtifact('invalid/nonconverged native response')
        if any(v!=0 for i in t['frozen_indices'] for v in mu[i][1:]):raise InvalidArtifact('source induced response not frozen')
    else:
        if energies is None or set(energies)!=set(ENERGY_NAMES) or energies['polarization']!=0 or eps is not None:
            raise InvalidArtifact('invalid static energy')
        if any(energies[k]!=0 for k in ENERGY_NAMES[4:]) or abs(energies['total']-sum(energies[k] for k in ENERGY_NAMES[1:4]))>TOL['algebra_kcal']:
            raise InvalidArtifact('native static component accounting failed')
        if mode=='static' and set(fields)!=ids:raise InvalidArtifact('missing static driving fields')
        if any(v!=0 for a in mu.values() for v in (a if mode=='identity' else a[1:])):raise InvalidArtifact('static induced state not zero')
    if active!=[mode,*(['F','F','F'] if mode=='identity' else ['T','T','T'])]:raise InvalidArtifact('native enabled term mismatch')
    if mode=='identity' and (fields or energies['solvation_with_nonpolar']!=0):raise InvalidArtifact('identity retains solvent')
    result=dict(status='computed_native_hybrid_component',parameters=p,moments=moments['atoms'],max_moment_error=moment_error,
        fields=fields,response=mu,energies_kcal_mol=energies,exact_RMS_Debye=eps,electric=constants[0],dielec=constants[1],
        kernel_wall_seconds=timing[0],kernel_CPU_seconds=timing[1],active_terms=active)
    cache_key(result);return result


def accepted(t,directory):
    for p in sorted(Path(directory).glob('attempt_*/result.json')):
        r=read_json(p)
        if r['cache_key']!=t['cache_key']:raise InvalidArtifact('incompatible cached native component')
        if r.get('receipt'):verify(r['receipt']['log'])
        if r['status']=='computed_native_hybrid_component':
            if not r.get('receipt'):raise InvalidArtifact('successful component has no execution receipt')
            if r.get('supplied_fields'):verify(r['supplied_fields'])
            for pin in r.get('dependencies',[]):verify(pin)
            return r,p
    return None


def supplied(t,m,root):
    static=next(x for x in m['static_tasks'] if x['variant']==t['variant'] and x['case_id']==t['case_id'] and x['state']==t['state'])
    env=next(x for x in m['static_tasks'] if x['variant']==t['variant'] and x['case_id']==t['case_id'] and x['state']=='environment')
    a=accepted(static,root/'tasks'/static['task_id']);b=accepted(env,root/'tasks'/env['task_id'])
    if not a or not b:raise InvalidArtifact('required static source/environment output unavailable')
    ar,ap=a;br,bp=b;n=ar['parameters']['inventory'][0]
    Q=np.array([ar['fields'][str(i)] for i in range(1,n+1)]);E=np.array([br['fields'][str(i)] for i in range(1,n+1)])
    if np.max(abs(Q[:,0]-E[:,0]))>1e-12:raise InvalidArtifact('paired source/environment Born radii differ')
    result=E.copy();rf=np.zeros((n,2,3))
    if t['state']!='environment':
        meta=read_json(verify(t['boundary']));ids={pid:i for i,pid in enumerate(meta['physical_ids'])}
        data=read_json(verify(t['density_task']['probes']));indices=np.array([ids[i] for i in data['physical_ids']])
        expected=set(range(n))-{i-1 for i in t['frozen_indices']}
        if set(indices)!=expected:raise InvalidArtifact('density field does not cover every active site')
        with np.load(verify(t['density_task']['field_arrays'])) as arr:exact=arr['exact']@np.array(t['rotation']).T/BOHR_TO_A**2
        for start in (1,4,7,10):result[indices,start:start+3]+=exact
        rf[:,0,:]=(Q[:,7:10]-Q[:,1:4])-(E[:,7:10]-E[:,1:4])
        rf[:,1,:]=(Q[:,10:13]-Q[:,4:7])-(E[:,10:13]-E[:,4:7])
        result[:,7:10]+=rf[:,0,:];result[:,10:13]+=rf[:,1,:]
    d=root/'supplied_fields'/static['task_id'];d.mkdir(parents=True,exist_ok=True);target=d/'fields.dat'
    text=str(n)+'\n'+''.join(str(i)+' '+' '.join(format(v,'.17g') for v in row)+'\n' for i,row in enumerate(result,1))
    if target.exists():
        if target.read_text()!=text:raise InvalidArtifact('derived supplied field changed on recovery')
    else:target.write_text(text)
    return record(target),[record(ap),record(bp)],dict(max_d_p_RF_difference_e_A2=float(np.max(abs(rf[:,0]-rf[:,1]))))


def execute(path,retry=False):
    m=validate(path);root=Path(path).resolve().parent;sw=read_json(verify(m['software']))
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=64:raise InvalidArtifact('64CPU native allocation required')
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for t in m['static_tasks']+m['response_tasks']:
            directory=root/'tasks'/t['task_id']
            if accepted(t,directory):continue
            attempts=sorted(directory.glob('attempt_*'))
            if attempts and not retry:
                print(json.dumps(dict(task_id=t['task_id'],status='failed_attempt_requires_explicit_retry')),flush=True);continue
            out=directory/f'attempt_{len(attempts)+1:04d}';out.mkdir()
            result=dict(task_id=t['task_id'],cache_key=t['cache_key'],status='failed',job_id=os.environ['SLURM_JOB_ID'],
                energy_attempted=False,field_attempted=False,response_attempted=False,receipt=None)
            try:
                cmd=[str(verify(sw['executable'])),str(verify(t['xyz'])),str(verify(t['mask'])),str(verify(t['overrides'])),t['mode']]
                if t['mode']=='solve':
                    pin,deps,diagnostic=supplied(t,m,root);cmd.append(str(verify(pin)))
                    result.update(supplied_fields=pin,dependencies=deps,field_diagnostic=diagnostic)
                env={**os.environ,'OMP_NUM_THREADS':'64','OMP_MAX_ACTIVE_LEVELS':'1','OMP_PROC_BIND':'spread',
                     'OMP_PLACES':'cores','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
                start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
                with (out/'native.log').open('x') as f:r=subprocess.run(cmd,cwd=directory,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,env=env)
                after=resource.getrusage(resource.RUSAGE_CHILDREN)
                result['receipt']=dict(command=cmd,returncode=r.returncode,wall_seconds=time.monotonic()-start,
                    child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,peak_child_RSS_KiB=after.ru_maxrss,log=record(out/'native.log'))
                text=verify(result['receipt']['log']).read_text()
                result.update(energy_attempted='ALQUEMIA_ENERGY_BEGIN' in text,field_attempted='ALQUEMIA_FIELDS_BEGIN' in text,response_attempted='ALQUEMIA_RESPONSE_BEGIN' in text)
                if r.returncode:raise InvalidArtifact('native execution failed')
                result.update(parse(text,t))
            except Exception as error:result['failure_reason']=f'{type(error).__name__}: {error}'
            write_new(out/'result.json',result)
            print(json.dumps({k:result.get(k) for k in ('task_id','status','failure_reason')}),flush=True)


def collect(path,output):
    m=validate(path);root=Path(path).resolve().parent;terms=read_json(verify(m['terms']));rows={};parsed={};checks=[]
    for t in m['static_tasks']+m['response_tasks']:
        item=accepted(t,root/'tasks'/t['task_id'])
        if item:parsed[t['task_id']]=item[0];rows[t['task_id']]=dict(status=item[0]['status'],result=record(item[1]))
        else:rows[t['task_id']]=dict(status='unavailable')
    def fetch(case,state,variant,mode,eps=1e-9):
        ts=m['response_tasks'] if mode=='solve' else m['static_tasks']
        t=next(t for t in ts if t['case_id']==case and t['state']==state and t['variant']==variant and t['mode']==mode and t['poleps']==eps)
        return t,parsed.get(t['task_id'])
    def induction(t,r):
        vals=np.loadtxt(verify(r['supplied_fields']),skiprows=1)[:,1:];n=len(vals)
        mu=np.array([r['response'][str(i)] for i in range(1,n+1)])
        I=-.5*r['electric']/r['dielec']*math.fsum((mu[:,7:10]*vals[:,10:13]).ravel())
        return I,float(np.max(np.linalg.norm(mu[:,7:10],axis=1)))
    variants={}
    for label,var,eps in [('standard','primary',1e-7)]+[(v,v,1e-9) for v in VARIANTS]:
        cases={}
        for case in CASES:
            te,se=fetch(case,'environment',var,'static');tr,er=fetch(case,'environment',var,'solve',eps)
            if not se or not er:cases[case]=dict(status='unavailable');continue
            I0,mu0=induction(tr,er);endpoints={}
            for metal in ('Ca','La'):
                ts,s=fetch(case,metal,var,'static');tt,r=fetch(case,metal,var,'solve',eps)
                if not s or not r:continue
                I,mu=induction(tt,r);meta=read_json(verify(tt['boundary']));ids={p:i for i,p in enumerate(meta['physical_ids'])}
                data=read_json(verify(tt['density_task']['probes']));ix=np.array([ids[i] for i in data['physical_ids']])
                moments=np.array([a['global_'] for a in r['moments']])[ix];q=moments[:,0];dip=moments[:,1:4]/BOHR_TO_A;Q=moments[:,4:].reshape(-1,3,3)/BOHR_TO_A**2
                with np.load(verify(tt['density_arrays'])) as arr:phi=arr['phi'];H=arr['H']
                with np.load(verify(tt['density_task']['field_arrays'])) as arr:E=arr['exact']
                mat=np.array(tt['rotation']);E=E@mat.T;H=np.einsum('ij,njk,lk->nil',mat,H,mat)
                direct=coupling(q,dip,Q,phi,E,H);V=sum(direct[k] for k in ('charge','dipole','quadrupole'))
                G=s['energies_kcal_mol']['solvation_with_nonpolar']-se['energies_kcal_mol']['solvation_with_nonpolar']
                term=terms[case+'_'+metal];parts=dict(DFT_vacuum_kcal=term['DFT_vacuum_hartree']*HA_TO_KCAL,direct_density_kcal=V,
                    GK_permanent_transfer_kcal=G,environment_induction_transfer_kcal=I-I0,short_context_kcal=term['short_context_kcal'])
                endpoints[metal]=dict(components=parts,environment_correction_kcal=V+G+I-I0,total_kcal=math.fsum(parts.values()),
                    direct_components_kcal=direct,induction_total_kcal=I,induction_environment_only_kcal=I0,
                    max_induced_dipole_eA=mu,large_response_flag=mu>TOL['large_dipole_eA'],native_static_components=s['energies_kcal_mol'],
                    native_environment_static_components=se['energies_kcal_mol'],exact_RMS_Debye=r['exact_RMS_Debye'])
            if len(endpoints)!=2:cases[case]=dict(status='unavailable');continue
            components={k:endpoints['Ca']['components'][k]-endpoints['La']['components'][k] for k in endpoints['Ca']['components']}
            R=math.fsum(components.values());error=R-(endpoints['Ca']['total_kcal']-endpoints['La']['total_kcal'])
            checks.append(dict(name=label+'_'+case+'_algebra',error_kcal=error,pass_=abs(error)<=TOL['algebra_kcal']))
            cases[case]=dict(status='computed_development_contrast',endpoints=endpoints,components_R_kcal=components,R_kcal=R,evidence=terms[case+'_Ca']['evidence'])
        if all(c['status']!='unavailable' for c in cases.values()):
            partition=cases['GGR_connected']['R_kcal']-cases['GGR_extended']['R_kcal']
            contrasts=[dict(alpha=a,GGR=g,difference_kcal=cases[a]['R_kcal']-cases[g]['R_kcal'],pass_=cases[a]['R_kcal']-cases[g]['R_kcal']>TOL['ordering_kcal'])
                for a in ('ALPHA_1F6S','ALPHA_6IP9') for g in ('GGR_extended','GGR_connected')]
            variants[label]=dict(status='complete',cases=cases,partition_kcal=partition,partition_pass=abs(partition)<=TOL['partition_kcal'],contrasts=contrasts,ordering_pass=all(c['pass_'] for c in contrasts))
        else:variants[label]=dict(status='unavailable',cases=cases)
    primary=variants['primary']
    if primary['status']=='complete':
        for var in ('standard','rigid'):
            if variants[var]['status']!='complete':continue
            for case in CASES:
                a=primary['cases'][case];b=variants[var]['cases'][case]
                errors={'R':b['R_kcal']-a['R_kcal']}
                for metal in ('Ca','La'):
                    errors[metal+'_environment']=b['endpoints'][metal]['environment_correction_kcal']-a['endpoints'][metal]['environment_correction_kcal']
                    if var=='rigid':
                        errors.update({metal+'_'+k:b['endpoints'][metal]['components'][k]-a['endpoints'][metal]['components'][k] for k in a['endpoints'][metal]['components']})
                checks.append(dict(name=case+'_'+var,errors_kcal=errors,pass_=max(abs(v) for v in errors.values())<=TOL['refinement_kcal' if var=='standard' else 'rigid_kcal']))
        for var in ('radius_minus','radius_plus'):
            if variants[var]['status']!='complete':continue
            diffs={'partition':variants[var]['partition_kcal']-primary['partition_kcal']}
            diffs.update({x['alpha']+'_'+x['GGR']:x['difference_kcal']-y['difference_kcal'] for x,y in zip(variants[var]['contrasts'],primary['contrasts'])})
            checks.append(dict(name=var+'_relative_sensitivity',errors_kcal=diffs,pass_=max(abs(v) for v in diffs.values())<=TOL['radius_relative_kcal']))
    identity={}
    for state in ('environment','Ca','La'):
        t,r=fetch('GGR_extended',state,'identity','identity')
        if r:identity[state]=dict(native_static_components=r['energies_kcal_mol'],active_terms=r['active_terms'],
            environment_moments_zero=all(all(v==0 for v in a['global_']) for i,a in enumerate(r['moments'],1) if i not in t['frozen_indices']),
            response_all_zero=all(all(v==0 for v in row) for row in r['response'].values()))
    if len(identity)==3:
        for metal in ('Ca','La'):
            value=identity[metal]['native_static_components']['solvation_with_nonpolar']-identity['environment']['native_static_components']['solvation_with_nonpolar']
            identity[metal]['environment_correction_kcal']=value
            checks.append(dict(name=metal+'_vacuum_identity',error_kcal=value,pass_=abs(value)<=TOL['identity_kcal'] and all(x['environment_moments_zero'] and x['response_all_zero'] for x in identity.values())))
    attempts=[read_json(p) for p in root.glob('tasks/*/attempt_*/result.json')];complete=len(parsed)==111
    if complete and (len(checks)!=32 or any(v['status']!='complete' for v in variants.values())):
        raise InvalidArtifact('complete pilot lacks required comparisons')
    result=dict(protocol=PROTOCOL,manifest=record(path),complete=complete,variants=variants,identity=identity,checks=checks,tasks=rows,
        numerical_pass=complete and all(c['pass_'] for c in checks if 'relative_sensitivity' not in c['name']),
        radius_sensitivity_pass=complete and all(c['pass_'] for c in checks if 'relative_sensitivity' in c['name']),
        partition_pass=primary.get('partition_pass',False),ordering_pass=primary.get('ordering_pass',False),
        actual_energy_calls=sum(a['energy_attempted'] for a in attempts),actual_field_queries=sum(a['field_attempted'] for a in attempts),
        actual_response_solves=sum(a['response_attempted'] for a in attempts),
        native_process_wall_seconds=sum(a['receipt']['wall_seconds'] for a in attempts if a['receipt']),
        native_process_CPU_seconds=sum(a['receipt']['child_CPU_seconds'] for a in attempts if a['receipt']),
        native_kernel_wall_seconds=sum(a.get('kernel_wall_seconds',0) for a in attempts),**NULLS)
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('build')
    for name in ('parent-software','native-source','frontend','output'):s.add_argument('--'+name,required=True)
    s=sub.add_parser('prepare')
    for name in ('boundary','density','short-report','software','plan','output'):s.add_argument('--'+name,required=True)
    for name in ('dry-run','execute','collect'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name=='collect':s.add_argument('--output',required=True)
        if name=='execute':s.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.command=='build':
        r=build(a.parent_software,a.native_source,a.frontend,a.output);print(json.dumps({'returncode':r['returncode']}));raise SystemExit(r['returncode'])
    elif a.command=='prepare':
        m=prepare(a.boundary,a.density,a.short_report,a.software,a.plan,a.output);print(json.dumps({'static_tasks':len(m['static_tasks']),'response_tasks':len(m['response_tasks'])}))
    elif a.command=='dry-run':
        m=validate(a.manifest);print(json.dumps({'static_tasks':len(m['static_tasks']),'response_tasks':len(m['response_tasks']),'new_native_calls':0}))
    elif a.command=='execute':execute(a.manifest,a.retry_failed)
    else:
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','numerical_pass','radius_sensitivity_pass','partition_pass','ordering_pass')}))


if __name__=='__main__':main()
