"""Frozen-density field + full GB + qualified local MACE readout experiment."""
from __future__ import annotations
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
from pathlib import Path
import shutil
import socket
import time
import numpy as np
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from affordable_solver import run_command
from density_embedding import parse_potential
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms

SCHEMA='alquemia.mace_explicit_field_short.v1'
PROTOCOL='normalized_vacuum_density_ff19SB_full_OBC2_POLAR_short_hybrid_v1'
GB_SHA='69f88af3655ea6d9784d682e528f7cf3a780410b93dda999671f7446e425e446'
QM_SHA='d63801217ed7b6092cf67f4224b34b28574199aefc1ca291abde34c07fee6275'
CHECKPOINT='fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a'
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')
TOL={'coupling_error_kcal_mol':1.,'algebra_kcal_scale':1e-7,'partition_kcal_scale':2.,'ordering_kcal_scale':.02}
NULLS={'aqueous_affinity_score':None,'reference':None,'calibrated_class':None,'combined_gradient':None,
       'response_status':'response_model_not_validated','relaxation_correction':None,'baseline_changed':False}


def accepted_source(collection,task_id):
    c=read_json(verify(collection));mp=verify(c['manifest']);m=read_json(mp)
    t=next(t for t in m['tasks'] if t['task_id']==task_id);r=c['rows'][task_id]
    if not any(accepted_attempt(a,t,mp)==r for a in (mp.parent/'execution'/task_id).glob('attempt_*')):
        raise InvalidArtifact('source lacks accepted receipt: '+task_id)
    return t,r,m


def source_data(gb,charge_manifest,reuse,short_reference):
    g=read_json(gb);gm=read_json(verify(g['manifest']));cm=read_json(charge_manifest);a=read_json(reuse)
    if (g['manifest']['sha256']!=GB_SHA or g['status']!='complete' or not g['numerical_gate_pass']
            or cm['quantum']['sha256']!=QM_SHA or a['status']!='complete' or a['checkpoint']['sha256']!=CHECKPOINT):
        raise InvalidArtifact('declared completed sources required')
    h=read_json(verify(gm['vacuum_hybrid']));q=read_json(verify(cm['quantum']))
    sr=read_json(short_reference);sm=read_json(verify(sr['manifest']))
    if sr['status']!='complete' or not sr['numerical_checks_pass'] or sm['model']['checkpoint']!=a['checkpoint']:
        raise InvalidArtifact('qualified short adapter required')
    if sm['software']!=a['software']:raise InvalidArtifact('short/full software differs')
    for name in CASES:
        state=read_json(verify(gm['states'][name]));gid=read_json(verify(state['normalized_preparation']))['case_id']
        for metal in ('Ca','La'):
            t,r,fm=accepted_source(a['full_source'],gid+'_'+metal+'_primary')
            reuse_row=a['whole_reuses'][gid][metal]
            if (reuse_row['result']!=r or reuse_row['short_energy_eV']!=r['energy_components_eV']['interaction_energy']
                    or fm['model']['checkpoint']!=a['checkpoint'] or fm['software']!=a['software']
                    or xyz(verify(t['xyz']))!=xyz(verify(state['endpoints'][metal]['xyz']))
                    or t['charge']!=state['endpoints'][metal]['charge'] or t['spin_multiplicity']!=1):
                raise InvalidArtifact('whole short reuse changed')
            _,br,_=accepted_source(record(gb),name+'_'+metal+'_full_primary')
            if br['GB_reaction_kcal_mol']!=g['cases'][name]['endpoints'][metal]['GB_total_kcal_mol']:
                raise InvalidArtifact('GB scalar/source differs')
            ct=next(t for t in cm['tasks'] if t['task_id']==name+'_'+metal)
            qt=next(t for t in q['tasks'] if t['task_id']==ct['task_id'])
            receipt=read_json(verify(ct['source_receipt']))
            if (receipt['manifest']!=cm['quantum'] or receipt['returncode'] or not receipt['normal_termination']
                    or receipt['artifacts']['output']!=ct['source_output'] or ct['xyz']!=qt['xyz'] or ct['charge']!=qt['charge']):
                raise InvalidArtifact('quantum receipt/state differs')
            verify(ct['source_output'])
            for key,pin in ct['files'].items():
                verify(pin);verify(ct['source_wavefunctions'][key])
                if pin['sha256']!=ct['source_wavefunctions'][key]['sha256']:raise InvalidArtifact('saved density changed')
    return g,gm,cm,a,sm,h


def probes(state,atoms):
    env=np.array(state['environment_charges_e']);ids=np.flatnonzero(env)
    positions=np.array([a['xyz_A'] for a in state['physical_atoms']])[ids]
    nearest=float(np.min(np.linalg.norm(positions[:,None,:]-np.array([a[1:] for a in atoms])[None,:,:],axis=2)))
    if nearest<1e-6:raise InvalidArtifact('environment probe coincides with QM nucleus/cap')
    return positions/BOHR_TO_A,{'physical_indices':ids.tolist(),'weights_e':env[ids].tolist(),
        'physical_ids':[state['physical_atoms'][i]['id'] for i in ids],'nearest_QM_nucleus_or_cap_A':nearest}


def prepare(gb,charges,reuse,short_reference,agreement,output):
    g,gm,cm,a,sm,h=source_data(gb,charges,reuse,short_reference)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    sources={k:verify(v) for parent in (cm,sm) for k,v in parent['implementation'].items()}
    for name in ('mace_hybrid.py','mace_explicit_field_short.py','mace_omol_coupling_audit.py'):sources[name]=Path(__file__).with_name(name)
    pins={}
    for name,p in sources.items():shutil.copyfile(p,impl/name);pins[name]=record(impl/name)
    model=copy.deepcopy(sm['model']);model['preparation_policy']='exact_matched_normalized_H_QM_core_v1'
    tasks=[]
    for old in cm['tasks']:
        name,metal=old['case_id'],old['metal'];s=read_json(verify(gm['states'][name]));atoms=xyz(verify(old['xyz']))
        pp,weights=probes(s,atoms);d=out/'inputs'/old['task_id'];d.mkdir(parents=True)
        p=d/'points_bohr.xyz';p.write_text(str(len(pp))+'\n'+''.join(' '.join(format(x,'.12f') for x in row)+'\n' for row in pp))
        wp=d/'weights.json';write_new(wp,weights)
        t={'task_id':old['task_id'],'case_id':name,'metal':metal,'kind':'core','xyz':old['xyz'],
           'charge':old['charge'],'spin_multiplicity':1,'energy_component':'interaction_energy',
           'state':gm['states'][name],'points':record(p),'weights':record(wp),'files':old['files'],
           'source_receipt':old['source_receipt'],'source_output':old['source_output']}
        tasks.append(t)
    m={'schema_version':SCHEMA,'protocol_id':PROTOCOL,'agreement':record(agreement),'GB':record(gb),
       'charges':record(charges),'reuse_audit':record(reuse),'short_reference':record(short_reference),
       'software':sm['software'],'model':model,'implementation':pins,'tasks':tasks,'tolerances':TOL,
       'utilities':{'orca_vpot':cm['utilities']['orca_vpot']},'utility_python':cm['python'],
       'numerical_reference':sm['numerical_reference'],'compute_budget':None,'wall_time_limit':None,
       'run_inventory':{'new_DFT':0,'new_GB':0,'new_charge_fit':0,'new_short_core':8,'new_vpot':8,'whole_short_reused':6},**NULLS}
    for t in tasks:t['cache_key']=key(t,m)
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def key(task,m):
    return cache_key({'task':{k:v for k,v in task.items() if k!='cache_key'},
        **{k:m[k] for k in ('model','software','implementation','utilities','GB','charges','reuse_audit','short_reference','tolerances')}})


def validate(manifest):
    m=read_json(manifest)
    if m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['tolerances']!=TOL:raise InvalidArtifact('protocol changed')
    for pin in [m['agreement'],m['software'],m['utility_python'],m['numerical_reference'],*m['implementation'].values(),*m['utilities'].values()]:verify(pin)
    software=read_json(verify(m['software']))
    for pin in [software['python'],software['requirements'],*read_json(verify(software['backend_source_inventory']))['files']]:verify(pin)
    g,gm,cm,a,sm,h=source_data(*(verify(m[k]) for k in ('GB','charges','reuse_audit','short_reference')))
    expected=copy.deepcopy(sm['model']);expected['preparation_policy']='exact_matched_normalized_H_QM_core_v1'
    if m['model']!=expected or m['software']!=sm['software'] or m['utilities']['orca_vpot']!=cm['utilities']['orca_vpot']:
        raise InvalidArtifact('model/software/utility changed')
    verify(m['model']['checkpoint'])
    if len(m['tasks'])!=8 or {t['task_id'] for t in m['tasks']}!={t['task_id'] for t in cm['tasks']}:
        raise InvalidArtifact('eight exact endpoints required')
    paired={}
    for t in m['tasks']:
        old=next(x for x in cm['tasks'] if x['task_id']==t['task_id'])
        for k in ('case_id','metal','xyz','charge','files','source_receipt','source_output'):
            if t[k]!=old[k]:raise InvalidArtifact('endpoint source changed: '+k)
        if (t['state']!=gm['states'][t['case_id']] or t['kind']!='core' or t['spin_multiplicity']!=1 or t['energy_component']!='interaction_energy'):
            raise InvalidArtifact('short state/component changed')
        atoms=xyz(verify(t['xyz']));check_atoms(atoms,t['charge']);s=read_json(verify(t['state']))
        pp,w=probes(s,atoms);saved=np.loadtxt(verify(t['points']),skiprows=1)
        if saved.shape!=pp.shape or not np.allclose(saved,pp,atol=1e-11,rtol=0) or read_json(verify(t['weights']))!=w:
            raise InvalidArtifact('environment probes/weights changed')
        if key(t,m)!=t['cache_key']:raise InvalidArtifact('scientific cache key differs')
        paired.setdefault(t['case_id'],[]).append((saved,w,atoms))
    for (a,aw,ax),(b,bw,bx) in paired.values():
        if not np.array_equal(a,b) or aw!=bw or [x[1:] for x in ax]!=[x[1:] for x in bx] or [x[0] for x in ax[1:]]!=[x[0] for x in bx[1:]]:
            raise InvalidArtifact('paired coordinates/weights differ')
    return {'status':'pass','tasks':8,'potential_calls':8,'manifest':record(manifest)}


def accepted_potential(attempt,t,manifest):
    try:
        r=read_json(attempt/'execution.json')
        if r['task']!=t or r['manifest']!=record(manifest) or r['status']!='complete' or not r['slurm_job_id']:return None
        u=r['vpot']
        if u['returncode'] or u['slurm_job_id']!=r['slurm_job_id']:return None
        for pin in (u['log'],u['resource_usage'],r['potential']):verify(pin)
        return r
    except (OSError,ValueError,KeyError):return None


def execute_potential(manifest,retry_failed=False):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);root=mp.parent/'potential_execution';root.mkdir(exist_ok=True)
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('potential utility requires allocation')
    def one(t):
        td=root/t['task_id'];td.mkdir(exist_ok=True);attempts=sorted(td.glob('attempt_*'))
        for a in reversed(attempts):
            if accepted_potential(a,t,mp):return {'task_id':t['task_id'],'status':'reused','receipt':record(a/'execution.json')}
        if attempts and not retry_failed:return {'task_id':t['task_id'],'status':'failed_attempt_requires_explicit_retry'}
        d=td/f'attempt_{len(attempts)+1:04d}';d.mkdir();start=time.monotonic()
        r={'task':t,'manifest':record(mp),'host':socket.gethostname(),'slurm_job_id':os.environ['SLURM_JOB_ID'],'status':'failed'}
        try:
            for pin in t['files'].values():shutil.copyfile(verify(pin),d/Path(pin['path']).name)
            r['vpot']=run_command([str(verify(m['utilities']['orca_vpot'])),str(d/'endpoint.runtime.gbw'),
                'endpoint.runtime.scfp',str(verify(t['points'])),str(d/'potential.out')],d,d/'vpot.log',d/'vpot.resources.txt')
            if r['vpot']['returncode']:raise InvalidArtifact('native potential utility failed')
            points=np.loadtxt(verify(t['points']),skiprows=1);parse_potential(d/'potential.out',points)
            for pin in t['files'].values():
                if record(d/Path(pin['path']).name)['sha256']!=pin['sha256']:raise InvalidArtifact('utility mutated density')
            r.update(status='complete',potential=record(d/'potential.out'))
        except Exception as exc:r['reason']=str(exc)
        r['wall_seconds']=time.monotonic()-start;write_new(d/'execution.json',r)
        return {'task_id':t['task_id'],'status':r['status'],'receipt':record(d/'execution.json')}
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ThreadPoolExecutor(max_workers=min(8,int(os.environ['SLURM_CPUS_ON_NODE']))) as pool:rows=list(pool.map(one,m['tasks']))
    return {'status':'complete' if all(r['status'] in ('complete','reused') for r in rows) else 'incomplete','rows':rows,'manifest':record(mp)}


def collect(manifest):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        good=[];pots=[]
        for folder,accept,acc in [('execution',accepted_attempt,good),('potential_execution',accepted_potential,pots)]:
            for a in sorted((mp.parent/folder/t['task_id']).glob('attempt_*')):
                r=accept(a,t,mp);rp=a/('receipt.json' if folder=='execution' else 'execution.json')
                attempts.append({'task_id':t['task_id'],'kind':folder,'directory':str(a),'accepted':r is not None,'receipt':record(rp) if rp.exists() else None})
                if r:acc.append(r)
        row={'status':'unavailable','short':good[-1] if good else None,'potential_receipt':pots[-1] if pots else None}
        if good and pots:
            pp=np.loadtxt(verify(t['points']),skiprows=1);v=parse_potential(verify(pots[-1]['potential']),pp);weights=read_json(verify(t['weights']))
            from mace_omol_coupling_audit import cross_terms
            s=read_json(verify(t['state']));terms=cross_terms(s['physical_atoms'],np.array(s['endpoints'][t['metal']]['QM_charges_e']),np.array(s['environment_charges_e']))
            exact=v*np.array(weights['weights_e'])*HA_TO_KCAL
            row.update(status='computed',direct_exact_kcal_mol=float(exact.sum()),direct_projected_kcal_mol=float(terms.sum()),
                exact_per_probe_kcal_mol=exact.tolist(),potential_au=v.tolist(),weights=t['weights'])
        rows[t['task_id']]=row
    complete=all(r['status']=='computed' for r in rows.values());cases={};checks=[];quality=[];contrasts=[];partition=None
    def check(name,error):checks.append({'name':name,'error_kcal_scale':error,'pass':abs(error)<=TOL['algebra_kcal_scale']})
    if complete:
        g=read_json(verify(m['GB']));gm=read_json(verify(g['manifest']));h=read_json(verify(gm['vacuum_hybrid']));reuse=read_json(verify(m['reuse_audit']))
        for name in CASES:
            ends={};gid=h['cases'][name]['global_id']
            for metal in ('Ca','La'):
                row=rows[name+'_'+metal];core=row['short']['energy_eV'];full=reuse['whole_reuses'][gid][metal]['short_energy_eV']
                e={'DFT_vacuum_kcal_mol':h['cases'][name]['endpoints'][metal]['DFT_vacuum_hartree']*HA_TO_KCAL,
                   'direct_exact_kcal_mol':row['direct_exact_kcal_mol'],'GB_kcal_mol':g['cases'][name]['endpoints'][metal]['GB_total_kcal_mol'],
                   'short_context_model_kcal':(full-core)*EV_TO_KCAL}
                e['hybrid_kcal_scale']=sum(e.values());e.update(short_full_eV=full,short_core_eV=core,direct_projected_kcal_mol=row['direct_projected_kcal_mol']);ends[metal]=e
            components={k:ends['Ca'][k]-ends['La'][k] for k in ('DFT_vacuum_kcal_mol','direct_exact_kcal_mol','GB_kcal_mol','short_context_model_kcal')}
            r=sum(components.values());check(name+'_direct_algebra',r-(ends['Ca']['hybrid_kcal_scale']-ends['La']['hybrid_kcal_scale']))
            error=(ends['Ca']['direct_projected_kcal_mol']-ends['La']['direct_projected_kcal_mol'])-components['direct_exact_kcal_mol']
            quality.append({'name':name,'paired_projected_minus_exact_kcal_mol':error,'pass':abs(error)<=TOL['coupling_error_kcal_mol']})
            cases[name]={'endpoints':ends,'components_R':components,'hybrid_R_kcal_scale':r,'coupling_error_kcal_mol':error,'evidence':h['cases'][name]['evidence']}
        error=cases['GGR_connected']['coupling_error_kcal_mol']-cases['GGR_extended']['coupling_error_kcal_mol']
        quality.append({'name':'GGR_partition','paired_projected_minus_exact_kcal_mol':error,'pass':abs(error)<=TOL['coupling_error_kcal_mol']})
        partition=cases['GGR_connected']['hybrid_R_kcal_scale']-cases['GGR_extended']['hybrid_R_kcal_scale']
        for a in ('ALPHA_1F6S','ALPHA_6IP9'):
            for b in ('GGR_extended','GGR_connected'):
                delta=cases[a]['hybrid_R_kcal_scale']-cases[b]['hybrid_R_kcal_scale']
                contrasts.append({'alpha':a,'GGR':b,'difference_kcal_scale':delta,'pass':delta>TOL['ordering_kcal_scale']})
    return {'protocol_id':PROTOCOL,'manifest':record(mp),'status':'complete' if complete else 'incomplete','rows':rows,'attempts':attempts,
        'cases':cases,'checks':checks,'numerical_gate_pass':complete and all(c['pass'] for c in checks),'coupling_quality':quality,
        'charge_representation_gate_pass':complete and all(c['pass'] for c in quality),'partition_shift_kcal_scale':partition,
        'partition_gate_pass':partition is not None and abs(partition)<=TOL['partition_kcal_scale'],
        'contrasts':contrasts,'ordering_gate_pass':complete and all(c['pass'] for c in contrasts),'tolerances':TOL,
        'collection_implementation':record(__file__),**NULLS}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('gb','charges','reuse','short-reference','agreement','output'):a.add_argument('--'+k,required=True)
    for op in ('dry-run','execute-potential','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
        if op=='execute-potential':a.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.gb,a.charges,a.reuse,a.short_reference,a.agreement,a.output)
    elif a.op=='dry-run':r=validate(a.manifest)
    elif a.op=='execute-potential':r=execute_potential(a.manifest,a.retry_failed)
    else:r=collect(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','attempts','cases')},indent=2))
