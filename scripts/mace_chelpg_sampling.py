"""Native fixed-orbital CHELPG replay; separate from converged SCF scoring."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
import re
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz

PROTOCOL='native_NoIter_CHELPG_same_density_sampling_v1'
METHOD='r2SCAN-3c NoAutostart DefGrid3 TightSCF SP MORead CalcGuessEnergy NoIter CHELPG'


def scientific_input(charge,grid,extent):
    return (f'! {METHOD}\n%moinp "seed.gbw"\n%method\n DoEQ false\nend\n'
        '%pointcharges "environment.pc"\n'+f'%chelpg\n GRID {grid}\n RMAX {extent}\n VDWRADII COSMO\n DIPOLE false\nend\n'
        +f'* xyzfile {charge} 1 core.xyz\n')


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},protocol=m['protocol_id'],config=m['config'],orca=m['orca'],implementation=m['implementation']))


def prepare(config,output):
    c=read_json(config);verify(c['plan'])
    if c['protocol']!=PROTOCOL or c['phase']!='qualification' or c['case_ids']!=['GGR_2FVY'] or c['settings']!={'original':{'grid_A':.3,'extent_A':2.8}}:
        raise InvalidArtifact('only declared first-stage identity qualification implemented')
    obs=read_json(verify(c['observations_manifest']));dr=read_json(verify(c['density_report']))
    qr=read_json(verify(obs['quantum']));qm=read_json(verify(qr['manifest']))
    if dr['manifest']!=c['observations_manifest'] or dr['status']!='complete' or not dr['projection_gate_pass'] or qr['status']!='complete':
        raise InvalidArtifact('actual density/quantum observations unavailable')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in qm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    shutil.copyfile(__file__,impl/Path(__file__).name);pins={p.name:record(p) for p in impl.glob('*.py')}
    m=dict(protocol_id=PROTOCOL,config=record(config),phase=c['phase'],agreement=c['plan'],implementation=pins,method=METHOD,orca=qm['orca'],
        execution_policy={k:pins[n] for k,n in [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]},
        compute_budget=None,wall_time_limit=None,tasks=[],planned_SCF_optimizations=0,requested_property_starts=2,
        generic_runner_SCF_status='retain raw flag; NoIter qualification is a separate explicit contract',new_score=None,baseline_changed=False)
    for case in c['case_ids']:
        for metal in ('Ca','La'):
            tid=case+'_'+metal;src=next(t for t in obs['tasks'] if t['task_id']==tid)
            qt=next(t for t in qm['tasks'] if t['task_id']==tid);d=root/(tid+'_original');d.mkdir()
            for pin,name in [(qt['xyz'],'core.xyz'),(qt['pointcharges'],'environment.pc'),(src['files']['gbw'],'seed.gbw')]:shutil.copyfile(verify(pin),d/name)
            (d/'endpoint.inp').write_text(scientific_input(qt['charge'],.3,2.8))
            t=dict(task_id=d.name,case_id=case,metal=metal,charge=qt['charge'],multiplicity=1,task_type='fixed_orbital_population',
                sampling={'grid_A':.3,'extent_A':2.8},input=record(d/'endpoint.inp'),xyz=record(d/'core.xyz'),pointcharges=record(d/'environment.pc'),seed_gbw=record(d/'seed.gbw'),
                source_gbw=src['files']['gbw'],source_quantum_receipt=qr['rows'][tid]['receipt'],source_quantum_output=qr['rows'][tid]['output'],
                source_energy_hartree=qr['rows'][tid]['energy_hartree'],source_density_row=dr['rows'][tid],source_observation_task=src,
                output_path=str(d/'endpoint.out'))
            t['cache_key']=key(t,m);m['tasks'].append(t)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    from affordable_workflow import dry_run
    m=read_json(path);c=read_json(verify(m['config']))
    if m['protocol_id']!=PROTOCOL or m['method']!=METHOD or m['phase']!='qualification' or len(m['tasks'])!=2 or m['planned_SCF_optimizations']!=0:
        raise InvalidArtifact('fixed-orbital qualification scope differs')
    for pin in [m['agreement'],m['orca'],*m['implementation'].values(),c['observations_manifest'],c['density_report']]:verify(pin)
    obs=read_json(verify(c['observations_manifest']));qm=read_json(verify(read_json(verify(obs['quantum']))['manifest']))
    if {t['task_id'] for t in m['tasks']}!={'GGR_2FVY_Ca_original','GGR_2FVY_La_original'}:raise InvalidArtifact('paired qualification inventory differs')
    for t in m['tasks']:
        if t['cache_key']!=key(t,m) or t['sampling']!=c['settings']['original'] or verify(t['input']).read_text()!=scientific_input(t['charge'],.3,2.8):
            raise InvalidArtifact('changed method or sampling input')
        qt=next(x for x in qm['tasks'] if x['task_id']==t['case_id']+'_'+t['metal'])
        if t['charge']!=qt['charge'] or t['multiplicity']!=qt['multiplicity'] or t['xyz']['sha256']!=qt['xyz']['sha256'] or t['pointcharges']['sha256']!=qt['pointcharges']['sha256']:
            raise InvalidArtifact('geometry, generating field or electronic state differs')
        for name in ('seed_gbw','source_gbw','source_quantum_receipt','source_quantum_output','pointcharges'):verify(t[name])
        if t['seed_gbw']['sha256']!=t['source_gbw']['sha256']:raise InvalidArtifact('imported orbital file differs')
    return dry_run(path)


def native_replay(t,manifest,charge_parser=None):
    """NoIter has no convergence claim: verify its separate native contract."""
    from density_embedding import parse_chelpg
    from mace_omol_vacuum import NUMBERS
    p=Path(t['output_path']);receipt=read_json(str(p)+'.execution.json');text=p.read_text()
    if receipt['manifest']!=record(manifest) or receipt['returncode']!=0 or not receipt['normal_termination'] or receipt['orca_version']!='6.1.1':
        raise InvalidArtifact('native replay did not terminate normally under the pinned manifest')
    m=read_json(manifest)
    for pin in receipt['artifacts'].values():verify(pin)
    if (receipt['artifacts']['output']!=record(p) or receipt['artifacts']['template_input']!=t['input']
        or receipt['artifacts']['xyz']!=t['xyz'] or receipt['orca_executable']!=m['orca']
        or receipt['task_runner']!=m['execution_policy']['task_runner']):
        raise InvalidArtifact('native execution provenance differs')
    required=('Input Geometry matches current geometry (good)','Input basis set matches current basis set (good)',
              'Occupation numbers from the input are taken because of NoIter','MaxIter==0 has been chosen - skipping the SCF',
              'CalcGuessEnergy requested: will build one Fock matrix to compute the energy')
    if any(s not in text for s in required) or re.search(r'SCF NOT CONVERGED|SCF CONVERGED AFTER|^\s*SCF ITERATIONS\s*$|ORCA SCF GRADIENT CALCULATION',text,re.M):
        raise InvalidArtifact('not a supported fixed-orbital population replay')
    atoms=xyz(verify(t['xyz']));electrons=sum(NUMBERS[a[0]] for a in atoms)-t['charge']-(46 if t['metal']=='La' else 0)
    for pattern,wanted in [(r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)',t['charge']),
                           (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)',1),
                           (r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)',electrons)]:
        if re.findall(pattern,text)!=[str(wanted)]:raise InvalidArtifact('native replay electronic state differs')
    ecp=re.findall(r'Type\s+(\w+)\s+ECP\s+(\S+)\s+\(replacing\s+(\d+)\s+core electrons',text)
    if ecp!=([('La','Def2-ECP','46')] if t['metal']=='La' else []):raise InvalidArtifact('native replay ECP differs')
    counts=re.findall(r'Reading point charge file\s+\.{2,}\s+ok\s+\((\d+) point charges\)',text)
    count=int(verify(t['pointcharges']).read_text().splitlines()[0])
    if not counts or any(int(n)!=count for n in counts):raise InvalidArtifact('generating point-charge inventory differs')
    es=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+\d.EeDd]+)',text)
    if len(es)!=1 or 'DFTD4' not in text or not re.search(r'gCP correction\s+[-+0-9.]',text):raise InvalidArtifact('missing native replay energy/components')
    energy=float(es[0].replace('D','E'));delta=energy-t['source_energy_hartree']
    if text.count('CHELPG CHARGES GENERATION')!=1:raise InvalidArtifact('CHELPG population section missing or duplicated')
    section=text.split('CHELPG CHARGES GENERATION',1)[1].split('CHELPG charges calculated...',1)[0]+'CHELPG charges calculated...'
    if 'CHELPG charges calculated...' not in text:raise InvalidArtifact('charge fit incomplete')
    parser=parse_chelpg if charge_parser is None else charge_parser
    q=parser(section,[a[0] for a in atoms],t['charge']);old=np.array(t['source_density_row']['charge_e'])
    qerr=float(np.max(abs(q-old)))
    files={}
    for kind,suffix in (('gbw','gbw'),('density','densities'),('density_info','densitiesinfo')):
        f=p.parent/('endpoint.runtime.'+suffix)
        if not f.exists():raise InvalidArtifact('replay saved density unavailable')
        files[kind]=record(f)
    return dict(status='native_NoIter_population_completed',execution_receipt=record(str(p)+'.execution.json'),output=record(p),
        original_generic_SCF_converged=receipt['scf_converged'],SCF_optimization_performed=False,property_energy_hartree=energy,
        energy_identity_error_hartree=delta,default_charge_max_error_e=qerr,charge_e=q.tolist(),files=files,
        preliminary_identity_pass=abs(delta)<=1e-7 and qerr<=2e-6,density_potential_identity_pass=None,
        qualification_pass=None,qualification_status='requires_actual_saved_density_potential_comparison')


def collect_native(path,output):
    validate(path);m=read_json(path);rows={}
    for t in m['tasks']:
        try:rows[t['task_id']]=native_replay(t,path)
        except (ValueError,OSError) as exc:rows[t['task_id']]=dict(status='invalid_or_unavailable',reason=str(exc),preliminary_identity_pass=False,qualification_pass=None)
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,root/Path(__file__).name)
    result=dict(protocol_id=PROTOCOL,manifest=record(path),implementation=record(root/Path(__file__).name),rows=rows,
        native_population_complete=all(r['status']=='native_NoIter_population_completed' for r in rows.values()),
        preliminary_identity_pass=all(r['preliminary_identity_pass'] for r in rows.values()),qualification_pass=None,
        new_score=None,baseline_changed=False,generic_runner_note='SCF convergence flag retained; intentionally absent for documented NoIter property replay')
    write_new(root/'result.json',result);return result


def prepare_identity(population,output):
    from density_embedding import parse_potential
    pop=read_json(population);mp=verify(pop['manifest']);validate(mp);parent=read_json(mp)
    if not pop['preliminary_identity_pass']:raise InvalidArtifact('native default replay identity failed')
    cfg=read_json(verify(parent['config']));obs=read_json(verify(cfg['observations_manifest']))
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for n,pin in obs['implementation'].items():shutil.copyfile(verify(pin),impl/n)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    m=dict(protocol=PROTOCOL+'_density_identity',population=record(population),plan=parent['agreement'],
        utility=obs['utilities']['orca_vpot'],implementation={p.name:record(p) for p in impl.glob('*.py')},tasks=[],
        potential_tolerance_au=1e-8,workers=2,requested_queries=2,requested_DFT_calls=0)
    for t in parent['tasks']:
        row=native_replay(t,mp)
        if row!=pop['rows'][t['task_id']]:raise InvalidArtifact('native population collection changed')
        source=t['source_observation_task'];data=read_json(verify(source['probes']));coords=np.loadtxt(verify(source['points']),skiprows=1)
        values=parse_potential(verify(t['source_density_row']['potential']),coords);ix=slice(data['exterior_start'],None)
        expected=values[ix];points=coords[ix]
        if len(points)!=data['exterior_count']:raise InvalidArtifact('fixed exterior validation inventory differs')
        d=root/'tasks'/t['task_id'];d.mkdir(parents=True);pp=d/'points_bohr.xyz'
        pp.write_text(str(len(points))+'\n'+''.join(' '.join(format(float(x),'.14f') for x in v)+'\n' for v in points))
        write_new(d/'expected.json',dict(potential_au=expected.tolist(),source_potential=t['source_density_row']['potential'],source_points=source['points'],exterior_start=data['exterior_start']))
        task=dict(task_id=t['task_id'],files=row['files'],points=record(pp),expected=record(d/'expected.json'),source_replay_receipt=row['execution_receipt'])
        task['cache_key']=cache_key(dict(task=task,population=m['population'],utility=m['utility'],implementation=m['implementation'],tolerance=m['potential_tolerance_au']))
        m['tasks'].append(task)
    write_new(root/'manifest.json',m);return validate_identity(root/'manifest.json')


def validate_identity(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL+'_density_identity' or m['potential_tolerance_au']!=1e-8 or len(m['tasks'])!=2 or m['requested_DFT_calls']:
        raise InvalidArtifact('density identity scope changed')
    for pin in [m['population'],m['plan'],m['utility'],*m['implementation'].values()]:verify(pin)
    for t in m['tasks']:
        original={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key(dict(task=original,population=m['population'],utility=m['utility'],implementation=m['implementation'],tolerance=m['potential_tolerance_au'])):raise InvalidArtifact('density identity cache changed')
        for pin in [t['points'],t['expected'],t['source_replay_receipt'],*t['files'].values()]:verify(pin)
        expected=read_json(verify(t['expected']));verify(expected['source_potential']);verify(expected['source_points'])
    return m


def execute_identity(path,output):
    from affordable_solver import run_command
    from density_embedding import parse_potential
    m=validate_identity(path);root=Path(path).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<2:raise InvalidArtifact('declared utility allocation required')
    def work(t):
        d=root/'tasks'/t['task_id'];attempt=d/'attempt_0001'
        if attempt.exists():raise InvalidArtifact('retain existing attempt; collect it rather than rerunning')
        attempt.mkdir();r=dict(task_id=t['task_id'],cache_key=t['cache_key'],status='failed',receipt=None)
        try:
            for name,pin in t['files'].items():shutil.copyfile(verify(pin),attempt/Path(pin['path']).name)
            gbw=attempt/Path(t['files']['gbw']['path']).name;potential=attempt/'potential.out'
            r['receipt']=run_command([str(verify(m['utility'])),str(gbw),'endpoint.runtime.scfp',str(verify(t['points'])),str(potential)],attempt,attempt/'vpot.log',attempt/'resources.txt')
            if r['receipt']['returncode']:raise InvalidArtifact('actual density query failed')
            values=parse_potential(potential,np.loadtxt(verify(t['points']),skiprows=1));expected=np.array(read_json(verify(t['expected']))['potential_au'])
            error=float(np.max(abs(values-expected)))
            r.update(status='computed',potential=record(potential),maximum_identity_error_au=error,pass_=error<=m['potential_tolerance_au'])
        except Exception as exc:r.update(failure_reason=str(exc),pass_=False)
        write_new(attempt/'result.json',r);return r
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(work,m['tasks']))
    result=dict(protocol=m['protocol'],manifest=record(path),rows={r['task_id']:r for r in rows},complete=all(r['status']=='computed' for r in rows),qualification_pass=all(r['pass_'] for r in rows),new_score=None)
    write_new(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('prepare-identity');q.add_argument('--population',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute','collect-native'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op!='dry-run':q.add_argument('--output',required=True)
    for op in ('dry-run-identity','execute-identity'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op=='execute-identity':q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output)
    elif a.command=='prepare-identity':r=prepare_identity(a.population,a.output);r={'tasks':len(r['tasks']),'status':'prepared'}
    elif a.command=='dry-run':r=validate(a.manifest)
    elif a.command=='dry-run-identity':r=validate_identity(a.manifest);r={'tasks':len(r['tasks']),'status':'pass'}
    elif a.command=='execute-identity':r=execute_identity(a.manifest,a.output)
    elif a.command=='collect-native':r=collect_native(a.manifest,a.output)
    else:
        from affordable_workflow import execute
        validate(a.manifest);r=execute(a.manifest);write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k!='rows'}))
