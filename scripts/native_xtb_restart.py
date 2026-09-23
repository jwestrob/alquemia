#!/usr/bin/env python3
"""Eight pinned native-xTB self/cross restarts; existing ORCA executor unchanged."""
from __future__ import annotations
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import time

from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, record, verify, write_new
from run_orca_task_manifest import load_manifest_tasks, run_manifest
from compact_solvation import completed, diagnostics
from structure_informed_starts import scf_details

PROTOCOL='native_GFN2_Q88_La_xtbw_self_cross_restart_v1'
GEOMETRIES=('old_adaptive','new_template')
MEDIA=('vacuum','alpb')
SEEDS=('self','opposite')
SETTINGS={'charge':-2,'multiplicity':1,'electrons':438,'temperature_K':300,
          'MaxIter':500,'mpi_ranks':8,'workers':8,'maxcore_MiB':2000,
          'endpoint_tolerance_kcal_mol':.1,'differential_tolerance_kcal_mol':.2,
          'restart_marker':'INITIAL GUESS: XTBRESTART'}


def recipe(medium, restart=True):
    return ('! Native-GFN2-xTB'+('' if restart else ' NoAutostart')+
            (' ALPB(Water)' if medium=='alpb' else '')+
            '\n%maxcore 2000\n%method\n WriteXTBParam true\n ReadXTBParam false\nend\n'
            '%scf\n MaxIter 500\n SmearTemp 300\n UseXTBMixer true\nend\n'
            '* xyzfile -2 1 core.xyz\n')


def source_rows(sources):
    s=read_json(sources);verify(s['collection'])
    rows={(r['geometry'],r['medium']):r for r in s['rows']}
    if len(s['rows'])!=4 or set(rows)!={(g,m) for g in GEOMETRIES for m in MEDIA}:
        raise InvalidArtifact('exact four seed sources required')
    for (g,m),r in rows.items():
        if (r['charge'],r['multiplicity'],r['electron_count'],r['xtbw_bytes'])!=(-2,1,438,11864):
            raise InvalidArtifact('source electronic state or seed size differs')
        for k in ('xyz','xtbw','output','parameters','execution_receipt'):verify(r[k])
        receipt=read_json(verify(r['execution_receipt']))
        if (receipt['returncode']!=0 or not receipt['normal_termination'] or
                not receipt['scf_converged'] or receipt['orca_version']!='6.1.1'):
            raise InvalidArtifact('source receipt incomplete')
        if verify(receipt['artifacts']['template_input']).read_text()!=recipe(m,False):
            raise InvalidArtifact('source primary recipe differs')
        if 'INFO: Using special xTB SCF mixer' not in verify(r['output']).read_text():
            raise InvalidArtifact('source native mixer absent')
    if len({r['parameters']['sha256'] for r in rows.values()})!=1:
        raise InvalidArtifact('source parameter sets differ')
    for g in GEOMETRIES:
        if rows[g,'vacuum']['xyz']['sha256']!=rows[g,'alpb']['xyz']['sha256']:
            raise InvalidArtifact('same-geometry solvent coordinates differ')
    return rows


def prepare(sources,agreement,authorization,output,activation='seed_only'):
    rows=source_rows(sources);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    tasks=[]
    for g in GEOMETRIES:
        for medium in MEDIA:
            for kind in SEEDS:
                dest=rows[g,medium];sg=g if kind=='self' else next(x for x in GEOMETRIES if x!=g)
                seed=rows[sg,medium];tid=f'q88jh5__La__{g}__{medium}__{kind}'
                d=out/'tasks'/tid;d.mkdir(parents=True)
                shutil.copyfile(verify(dest['xyz']),d/'core.xyz')
                shutil.copyfile(verify(seed['xtbw']),d/'seed.immutable.xtbw')
                if activation=='matched_gbw_auto':shutil.copyfile(verify(seed['gbw']),d/'seed.immutable.gbw')
                (d/'endpoint.inp').write_text(recipe(medium))
                tasks.append(dict(task_id=tid,case_id='q88jh5-pqq-la_model',metal='La',
                    geometry=g,medium=medium,seed_kind=kind,seed_geometry=sg,
                    charge=-2,multiplicity=1,source=dest,seed_source=seed,
                    xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),
                    immutable_seed=record(d/'seed.immutable.xtbw'),
                    active_seed_path=str(d/'endpoint.runtime.xtbw'),output_path=str(d/'endpoint.out')))
                if activation=='matched_gbw_auto':tasks[-1].update(immutable_gbw=record(d/'seed.immutable.gbw'),active_gbw_path=str(d/'endpoint.runtime.gbw'))
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        dst=impl/p.name;shutil.copyfile(p,dst);pins[p.name]=record(dst)
    receipt=read_json(verify(next(iter(rows.values()))['execution_receipt']))
    m=dict(protocol_id=PROTOCOL+('_matched_gbw_activation_v2' if activation=='matched_gbw_auto' else ''),activation=activation,settings=SETTINGS,sources=record(sources),agreement=record(agreement),authorization=record(authorization),
           tasks=tasks,orca=receipt['orca_executable'],implementation=pins,
           execution_resources={'mpi_ranks':8,'concurrent_tasks':8},
           execution_policy={'task_runner':pins['run_orca_task_manifest.py'],
                             'runtime_renderer':pins['render_orca_runtime_input.py']},
           new_molecular_calls=8,new_DFT_MACE_calls=0,baseline_changed=False)
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,fresh=True);write_new(out/'PREFLIGHT.json',v)
    return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp);rows=source_rows(verify(m['sources']))
    mode=m.get('activation','seed_only')
    if mode not in ('seed_only','matched_gbw_auto'):raise InvalidArtifact('unknown activation route')
    if m['protocol_id']!=PROTOCOL+('_matched_gbw_activation_v2' if mode=='matched_gbw_auto' else '') or m['settings']!=SETTINGS or m['new_molecular_calls']!=8:
        raise InvalidArtifact('fixed restart scope differs')
    verify(m['agreement']);verify(m['authorization']);verify(m['orca'])
    for pin in m['implementation'].values():verify(pin)
    if len(m['tasks'])!=8 or {(t['geometry'],t['medium'],t['seed_kind']) for t in m['tasks']}!={
        (g,s,k) for g in GEOMETRIES for s in MEDIA for k in SEEDS}:
        raise InvalidArtifact('exact eight restart tasks required')
    _,normalized=load_manifest_tasks(mp)
    for t in m['tasks']:
        g,s,k=t['geometry'],t['medium'],t['seed_kind'];sg=g if k=='self' else next(x for x in GEOMETRIES if x!=g)
        if t['source']!=rows[g,s] or t['seed_source']!=rows[sg,s] or t['seed_geometry']!=sg:
            raise InvalidArtifact('seed/destination mapping changed')
        if (t['charge'],t['multiplicity'],t['metal'])!=(-2,1,'La'):
            raise InvalidArtifact('state changed')
        if verify(t['input']).read_text()!=recipe(s) or verify(t['xyz']).read_bytes()!=verify(rows[g,s]['xyz']).read_bytes():
            raise InvalidArtifact('source coordinates or exact recipe changed')
        if verify(t['immutable_seed']).read_bytes()!=verify(rows[sg,s]['xtbw']).read_bytes():
            raise InvalidArtifact('immutable seed changed')
        d=Path(t['output_path']).parent
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw':raise InvalidArtifact('runtime seed basename differs')
        expected={'core.xyz','endpoint.inp','seed.immutable.xtbw'}
        if mode=='matched_gbw_auto':
            if verify(t['immutable_gbw']).read_bytes()!=verify(rows[sg,s]['gbw']).read_bytes():
                raise InvalidArtifact('matched GBW seed changed')
            if Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('runtime GBW basename differs')
            expected.add('seed.immutable.gbw')
        if fresh and set(p.name for p in d.iterdir())!=expected:
            raise InvalidArtifact('nonfresh task directory: '+t['task_id'])
    return dict(status='validated',manifest=record(mp),tasks=len(normalized),new_molecular_calls_in_validation=0,
                native_restart_confirmation='requires_actual_output')


def execute(manifest):
    mp=Path(manifest).resolve();m=read_json(mp)
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('allocation required')
    if int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64:raise InvalidArtifact('64 allocated CPUs required')
    with (mp.parent/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);validate(mp,fresh=True)
        before=[]
        for t in m['tasks']:
            shutil.copyfile(verify(t['immutable_seed']),t['active_seed_path'])
            if m.get('activation')=='matched_gbw_auto':shutil.copyfile(verify(t['immutable_gbw']),t['active_gbw_path'])
            before.append({'task_id':t['task_id'],'active_seed':record(t['active_seed_path']),
                           'immutable_seed':t['immutable_seed'],'gbw_present':m.get('activation')=='matched_gbw_auto',
                           'active_gbw':record(t['active_gbw_path']) if m.get('activation')=='matched_gbw_auto' else None})
        write_new(mp.parent/'SEEDS_BEFORE.json',{'manifest':record(mp),'rows':before})
        start=time.time();error=None;results=None
        try:
            results=run_manifest(mp,orca_path=verify(m['orca']),workers=8,nprocs=8)
        except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        finally:
            after=[]
            for t in m['tasks']:
                p=Path(t['active_seed_path']);d=Path(t['output_path']).parent
                saved=d/'seed.after.xtbw'
                if p.exists():shutil.copyfile(p,saved)
                after.append({'task_id':t['task_id'],'active_seed_after':record(p) if p.exists() else None,
                              'preserved_after':record(saved) if saved.exists() else None})
                if m.get('activation')=='matched_gbw_auto':
                    gbw=Path(t['active_gbw_path']);gbw_after=d/'seed.after.gbw'
                    if gbw.exists():shutil.copyfile(gbw,gbw_after)
                    ges=d/'endpoint.runtime.ges'
                    after[-1].update(active_gbw_after=record(gbw) if gbw.exists() else None,
                                     preserved_gbw_after=record(gbw_after) if gbw_after.exists() else None,
                                     autostart_ges=record(ges) if ges.exists() else None)
            write_new(mp.parent/'SEEDS_AFTER.json',{'manifest':record(mp),'rows':after})
            write_new(mp.parent/'EXECUTION.json',dict(manifest=record(mp),results=results,error=error,
                wall_seconds=time.time()-start,slurm_job_id=os.environ['SLURM_JOB_ID'],
                allocated_cpus=int(os.environ['SLURM_CPUS_ON_NODE']),completed_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
        if error:raise InvalidArtifact(error)
    return {'status':'executed','tasks':8}


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows=[]
    before_path=mp.parent/'SEEDS_BEFORE.json';after_path=mp.parent/'SEEDS_AFTER.json'
    before={r['task_id']:r for r in read_json(before_path)['rows']} if before_path.exists() else {}
    after={r['task_id']:r for r in read_json(after_path)['rows']} if after_path.exists() else {}
    for t in m['tasks']:
        r={k:t[k] for k in ('task_id','geometry','medium','seed_kind','seed_geometry')}
        r.update(status='unavailable',energy_hartree=None,observed_energy_hartree=None,source=t['source'],seed_source=t['seed_source'])
        pin=completed(mp,t['task_id'])
        if pin:
            try:
                text=verify(pin['output']).read_text();detail=scf_details(text);audit=diagnostics(pin,t)
                r.update(actual=pin,observed_energy_hartree=pin['energy_hartree'],details=detail,audit=audit,
                         restart_marker_observed=SETTINGS['restart_marker'] in text,
                         seed_before=before.get(t['task_id']),seed_after=after.get(t['task_id']))
                if not r['restart_marker_observed']:
                    r.update(status='restart_not_confirmed',reason='no positive XTBRESTART output marker')
                else:
                    b=before[t['task_id']];a=after[t['task_id']]
                    wants_gbw=m.get('activation')=='matched_gbw_auto'
                    if b['active_seed']['sha256']!=t['immutable_seed']['sha256'] or b['gbw_present']!=wants_gbw:
                        raise InvalidArtifact('seed preexecution audit differs')
                    if wants_gbw and b['active_gbw']['sha256']!=t['seed_source']['gbw']['sha256']:
                        raise InvalidArtifact('matched GBW preexecution audit differs')
                    verify(a['preserved_after'])
                    if (not detail['native_mixer_observed'] or (detail['charge'],detail['multiplicity'],detail['electrons'])!=(-2,1,438)
                            or audit['parameter_export']['sha256']!=t['source']['parameters']['sha256']
                            or audit['charge_sanity_status']!='pass'):
                        raise InvalidArtifact('native method/state/charge not confirmed')
                    r.update(status='confirmed_restart',energy_hartree=pin['energy_hartree'],
                             change_from_unseeded_kcal_mol=(pin['energy_hartree']-t['source']['energy_hartree'])*HA_TO_KCAL)
            except Exception as exc:r.update(status='audit_failed',reason=str(exc))
        else:
            rp=Path(t['output_path']+'.execution.json')
            if rp.exists():r.update(status='execution_failed',receipt=record(rp))
        rows.append(r)
    transfers={};continuity={}
    for g in GEOMETRIES:
        transfers[g]={}
        for k in SEEDS:
            rr={r['medium']:r['energy_hartree'] for r in rows if r['geometry']==g and r['seed_kind']==k}
            transfers[g][k]=(rr['alpb']-rr['vacuum'])*HA_TO_KCAL if all(x is not None for x in rr.values()) else None
    for s in MEDIA:
        values=[r['energy_hartree'] for r in rows if r['medium']==s]
        continuity[s]=(max(values)-min(values))*HA_TO_KCAL if all(v is not None for v in values) else None
    result=dict(protocol_id=m['protocol_id'],manifest=record(mp),rows=rows,denominator=8,
                confirmed=sum(r['status']=='confirmed_restart' for r in rows),
                matched_La_transfer_kcal_mol=transfers,all_four_initialized_energy_range_kcal_mol=continuity,
                full_Ca_La_score=None,new_molecular_calls_in_collection=0,baseline_changed=False,
                interpretation='technical initialization diagnostic; no ground-state policy or new reference')
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='op',required=True)
    s=sp.add_parser('prepare');s.add_argument('--sources',type=Path,required=True);s.add_argument('--agreement',type=Path,required=True);s.add_argument('--authorization',type=Path,required=True);s.add_argument('--output',type=Path,required=True)
    s.add_argument('--activation',choices=['seed_only','matched_gbw_auto'],default='seed_only')
    for name in ('dry-run','execute','collect'):
        s=sp.add_parser(name);s.add_argument('--manifest',type=Path,required=True)
        if name=='collect':s.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.sources,a.agreement,a.authorization,a.output,a.activation)
    elif a.op=='dry-run':r=validate(a.manifest,fresh=True)
    elif a.op=='execute':r=execute(a.manifest)
    else:r=collect(a.manifest,a.output)
    print(json.dumps(r,indent=2,sort_keys=True))


if __name__=='__main__':main()
