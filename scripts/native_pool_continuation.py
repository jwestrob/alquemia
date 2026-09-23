"""Two fixed native continuations of the archived four PQQ geometry pools."""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import time

from affordable_common import InvalidArtifact, HA_TO_KCAL, BOHR_TO_A, read_json, record, verify, write_new, xyz
from compact_solvation import input_text, completed, diagnostics
from run_orca_task_manifest import load_manifest_tasks, run_manifest
from structure_informed_starts import scf_details

PROTOCOL='native_PQQ_fixed_pool_two_continuations_v1'
CASES=('1H4I','4MAE','q88jh5-pqq-la_model','a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1')
CANDIDATES=('origin','proposal_Ca','proposal_La','adaptive_Ca','adaptive_La')
SETTINGS={'MaxIter':500,'temperature_K':300,'mpi_ranks':8,'workers':8,'maxcore_MiB':2000,
          'stages':2,'reported_stage':2,'cell_tolerance_kcal_mol':.1,
          'contrast_tolerance_kcal_mol':.2,'pool_tolerance_kcal_mol':.2,
          'origin_selection_tolerance_kcal_mol':.1,'maximum_logical_calls':160}


def key(t): return (t['case_id'],t['candidate'],t['metal'],t['medium'])
def tid(t): return '__'.join(key(t))


def recipe(charge,multiplicity,medium,gradient=False):
    body=input_text(charge,multiplicity,medium,'native').replace(' NoAutostart','')
    body=body.replace('%scf\n','%scf\n MaxIter 500\n')
    if gradient: body=body.replace('! Native-GFN2-xTB','! Native-GFN2-xTB EnGrad')
    return body


def want_gradient(stage,row):
    return stage==2 and row['candidate']=='origin' and row['case_id'] in ('4MAE','q88jh5-pqq-la_model')


def source_rows(inventory):
    inv=read_json(inventory)
    expected={(c,q,z,s) for c in CASES for q in CANDIDATES for z in ('Ca','La') for s in ('vacuum','alpb')}
    rows=inv['rows']
    if len(rows)!=80 or {key(t) for t in rows}!=expected or inv['available']!=80:
        raise InvalidArtifact('exact verified 80-cell archive required')
    for r in rows:
        if r['status']!='compatible_seed_pair_available':raise InvalidArtifact('source seed unavailable')
        for name in ('gbw','xtbw','xyz','input','output','receipt','parameter_export'):verify(r[name])
        if Path(r['gbw']['path']).parent!=Path(r['xtbw']['path']).parent:raise InvalidArtifact('unmatched seed pair')
        body=re.sub(r'^\s*MaxIter 500\s*\n','',verify(r['input']).read_text(),flags=re.M)
        if body!=input_text(r['charge'],r['multiplicity'],r['medium'],'native'):
            raise InvalidArtifact('archive primary recipe differs')
    return inv,rows


def prepare(inventory,agreement,output,previous=None):
    inv,rows=source_rows(inventory);prev=read_json(previous) if previous else None
    stage=2 if previous else 1
    if prev:
        pm=read_json(verify(prev['manifest']))
        if prev['stage']!=1 or pm['inventory']!=record(inventory) or pm['agreement']!=record(agreement):
            raise InvalidArtifact('stage 2 requires exact stage 1 collection')
        previous_rows={key(r):r for r in prev['rows']}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];missing=[]
    for r in rows:
        seed={'gbw':r['gbw'],'xtbw':r['xtbw']}
        if prev:
            pr=previous_rows[key(r)]
            if pr['status']!='confirmed_restart':
                missing.append({'task_id':tid(r),**{k:r[k] for k in ('case_id','candidate','metal','medium')},
                                'status':'unavailable','reason':'stage1_restart_unavailable','source_row':pr})
                continue
            seed={'gbw':pr['seed_after']['preserved_gbw_after'],'xtbw':pr['seed_after']['preserved_after']}
        d=out/'tasks'/tid(r);d.mkdir(parents=True)
        shutil.copyfile(verify(r['xyz']),d/'core.xyz')
        for name in ('gbw','xtbw'):shutil.copyfile(verify(seed[name]),d/('seed.immutable.'+name))
        gradient=want_gradient(stage,r);(d/'endpoint.inp').write_text(recipe(r['charge'],r['multiplicity'],r['medium'],gradient))
        t={k:r[k] for k in ('case_id','candidate','metal','medium','charge','multiplicity')}
        t.update(task_id=tid(r),source=r,seed_source=seed,gradient_requested=gradient,
                 xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),output_path=str(d/'endpoint.out'),
                 immutable_seed=record(d/'seed.immutable.xtbw'),immutable_gbw=record(d/'seed.immutable.gbw'),
                 active_seed_path=str(d/'endpoint.runtime.xtbw'),active_gbw_path=str(d/'endpoint.runtime.gbw'))
        if gradient:t['engrad_path']=str(d/'endpoint.engrad')
        tasks.append(t)
    impl=out/'implementation';impl.mkdir();pins={}
    source_impl=verify(pm['implementation']['native_pool_continuation.py']).parent if prev else Path(__file__).parent
    for p in source_impl.glob('*.py'):
        dst=impl/p.name;shutil.copyfile(p,dst);pins[p.name]=record(dst)
    orca=read_json(verify(rows[0]['receipt']))['orca_executable']
    inputs=read_json(verify(inv['inputs']))
    m={'protocol_id':PROTOCOL,'stage':stage,'settings':SETTINGS,'tasks':tasks,'missing':missing,
       'cell_denominator':80,'inventory':record(inventory),'agreement':record(agreement),
       'previous':record(previous) if prev else None,'orca':orca,'implementation':pins,
       'reference':inputs['adaptive_reference'],'inputs':inv['inputs'],
       'execution_resources':{'mpi_ranks':8,'concurrent_tasks':8},
       'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
       'new_MACE_DFT_calls':0,'baseline_changed':False}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp);inv,rows=source_rows(verify(m['inventory']))
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['stage'] not in (1,2):raise InvalidArtifact('fixed protocol differs')
    verify(m['agreement']);verify(m['orca']);verify(m['reference']);verify(m['inputs'])
    for p in m['implementation'].values():verify(p)
    original={key(r):r for r in rows};actual=m['tasks']+m['missing']
    if len(actual)!=80 or {key(t) for t in actual}!=set(original):raise InvalidArtifact('fixed denominator differs')
    previous={}
    if m['stage']==1:
        if m['previous'] or m['missing']:raise InvalidArtifact('first stage must include all sources')
    else:
        p=read_json(verify(m['previous']));pm=read_json(verify(p['manifest']))
        if p['stage']!=1 or pm['inventory']!=m['inventory'] or pm['agreement']!=m['agreement']:raise InvalidArtifact('stage linkage differs')
        previous={key(r):r for r in p['rows']}
        if {key(t) for t in m['missing']}!={k for k,r in previous.items() if r['status']!='confirmed_restart'}:
            raise InvalidArtifact('missing previous cell removed or successful cell dropped')
    _,normalized=load_manifest_tasks(mp)
    for t in m['tasks']:
        r=original[key(t)]
        seed={'gbw':r['gbw'],'xtbw':r['xtbw']} if m['stage']==1 else {
            'gbw':previous[key(t)]['seed_after']['preserved_gbw_after'],
            'xtbw':previous[key(t)]['seed_after']['preserved_after']}
        if t['source']!=r or t['seed_source']!=seed or t['task_id']!=tid(r):raise InvalidArtifact('source/seed mapping changed')
        if any(t[k]!=r[k] for k in ('charge','multiplicity')):raise InvalidArtifact('state changed')
        grad=want_gradient(m['stage'],r)
        if t['gradient_requested']!=grad or verify(t['input']).read_text()!=recipe(t['charge'],t['multiplicity'],t['medium'],grad):
            raise InvalidArtifact('recipe/gradient request differs')
        if verify(t['xyz']).read_bytes()!=verify(r['xyz']).read_bytes():raise InvalidArtifact('coordinates changed')
        for field,name in (('immutable_seed','xtbw'),('immutable_gbw','gbw')):
            if verify(t[field]).read_bytes()!=verify(seed[name]).read_bytes():raise InvalidArtifact('seed changed')
        d=Path(t['output_path']).parent
        if (Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw'):
            raise InvalidArtifact('runtime basename differs')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp','seed.immutable.xtbw','seed.immutable.gbw'}:
            raise InvalidArtifact('nonfresh task directory')
    return {'status':'validated','manifest':record(mp),'stage':m['stage'],'tasks':len(normalized),'missing':len(m['missing']),
            'gradient_requests':sum(t['gradient_requested'] for t in m['tasks']),'new_calls_in_validation':0}


def execute_tasks(manifest,validator,seeded=True):
    """Common seed activation/receipt wrapper; caller supplies its finite-scope validator."""
    mp=Path(manifest).resolve();m=read_json(mp)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64:
        raise InvalidArtifact('64-CPU allocation required')
    with (mp.parent/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);validator(mp,fresh=True)
        before=[]
        for t in m['tasks']:
            if seeded:
                shutil.copyfile(verify(t['immutable_seed']),t['active_seed_path'])
                shutil.copyfile(verify(t['immutable_gbw']),t['active_gbw_path'])
            else:
                if 'NoAutostart' not in verify(t['input']).read_text():raise InvalidArtifact('seedless caller requires explicit NoAutostart')
                if any(Path(t[k]).exists() for k in ('active_seed_path','active_gbw_path')):raise InvalidArtifact('seedless caller has existing restart files')
            before.append({'task_id':t['task_id'],'active_seed':record(t['active_seed_path']) if seeded else None,
                           'active_gbw':record(t['active_gbw_path']) if seeded else None,'seeded':seeded})
        write_new(mp.parent/'SEEDS_BEFORE.json',{'manifest':record(mp),'rows':before})
        start=time.monotonic();error=None;results=None
        try:
            results=run_manifest(mp,orca_path=verify(m['orca']),workers=8,nprocs=8) if m['tasks'] else []
        except Exception as exc:error=repr(exc)
        finally:
            after=[]
            for t in m['tasks']:
                d=Path(t['output_path']).parent;row={'task_id':t['task_id']}
                for active,name,field in (('active_seed_path','xtbw','preserved_after'),('active_gbw_path','gbw','preserved_gbw_after')):
                    p=Path(t[active]);dest=d/('seed.after.'+name)
                    if p.exists():shutil.copyfile(p,dest)
                    row[field]=record(dest) if dest.exists() else None
                ges=d/'endpoint.runtime.ges';row['autostart_ges']=record(ges) if ges.exists() else None;after.append(row)
            write_new(mp.parent/'SEEDS_AFTER.json',{'manifest':record(mp),'rows':after})
            elapsed=time.monotonic()-start
            write_new(mp.parent/'EXECUTION.json',{'manifest':record(mp),'results':results,'error':error,'wall_seconds':elapsed,
                'allocated_core_seconds':int(os.environ['SLURM_CPUS_ON_NODE'])*elapsed,
                'slurm_job_id':os.environ['SLURM_JOB_ID'],'allocated_cpus':int(os.environ['SLURM_CPUS_ON_NODE'])})
        if error:raise InvalidArtifact(error)
    return {'status':'executed','tasks':len(m['tasks'])}


def gradient_result(task,pin):
    import gemmi
    import numpy as np
    from affordable_response import read_engrad
    receipt=read_json(verify(pin['receipt']));path=verify(receipt['artifacts']['engrad'])
    text=verify(pin['output']).read_text()
    if re.search('numerical gradient|numerical differentiation',text,re.I):raise InvalidArtifact('numerical gradient forbidden')
    if 'XTB SCF gradient' not in text or 'XTB CN gradient' not in text:raise InvalidArtifact('native analytic driver absent')
    if ('ALPB gradient' in text)!=(task['medium']=='alpb'):raise InvalidArtifact('solvent gradient differs')
    g=read_engrad(path);atoms=xyz(verify(task['xyz']))
    if abs(g['energy_Ha']-pin['energy_hartree'])>1e-8 or g['atom_count']!=len(atoms):raise InvalidArtifact('gradient energy/count mismatch')
    if not np.array_equal(g['atomic_numbers'],[gemmi.Element(a[0]).atomic_number for a in atoms]):raise InvalidArtifact('gradient order differs')
    if not np.allclose(g['coordinates_bohr']*BOHR_TO_A,[a[1:] for a in atoms],atol=1e-6,rtol=0):raise InvalidArtifact('gradient coordinates differ')
    return {'engrad':record(path),'gradient_kcal_mol_A':(g['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A).tolist(),
            'quantity':'gradient_not_force','numerical_gradient_used':False}


def collect_task(manifest,t,before,after):
    row={k:t[k] for k in ('task_id','case_id','candidate','metal','medium') if k in t}
    row.update(status='unavailable',energy_hartree=None,observed_energy_hartree=None,gradient=None,
               source=t['source'],seed_source=t['seed_source'],reason=None)
    try:
        pin=completed(manifest,t['task_id'])
        if not pin:
            op=Path(t['output_path']);rp=Path(str(op)+'.execution.json')
            row.update(reason='execution_failed_or_not_complete',artifacts=[record(p) for p in (op,rp) if p.exists()]);return row
        text=verify(pin['output']).read_text();detail=scf_details(text);audit=diagnostics(pin,t)
        row.update(actual=pin,observed_energy_hartree=pin['energy_hartree'],details=detail,audit=audit,seed_before=before,seed_after=after)
        if 'INITIAL GUESS: XTBRESTART' not in text:
            row.update(status='restart_not_confirmed',reason='actual XTBRESTART absent');return row
        if (before['active_seed']['sha256']!=t['immutable_seed']['sha256'] or before['active_gbw']['sha256']!=t['immutable_gbw']['sha256']):
            raise InvalidArtifact('preexecution seed audit differs')
        verify(after['preserved_after']);verify(after['preserved_gbw_after'])
        source=t['source']
        if (not detail['native_mixer_observed'] or (detail['charge'],detail['multiplicity'],detail['electrons'])!=(t['charge'],t['multiplicity'],source['electron_count'])
                or audit['parameter_export']['sha256']!=source['parameter_export']['sha256'] or audit['charge_sanity_status']!='pass'):
            raise InvalidArtifact('native method/state/parameter/charge audit failed')
        if t['gradient_requested']:row['gradient']=gradient_result(t,pin)
        row.update(status='confirmed_restart',energy_hartree=pin['energy_hartree'],
                   change_from_archive_kcal_mol=(pin['energy_hartree']-source['source_energy_hartree'])*HA_TO_KCAL)
    except Exception as exc:row.update(status='audit_failed',reason=str(exc))
    return row


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def index(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']} if p.exists() else {}
    before=index('SEEDS_BEFORE.json');after=index('SEEDS_AFTER.json')
    rows=[collect_task(mp,t,before.get(t['task_id']),after.get(t['task_id'])) for t in m['tasks']]
    rows.extend({**t,'energy_hartree':None,'observed_energy_hartree':None,'gradient':None} for t in m['missing'])
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'stage':m['stage'],'rows':rows,'denominator':80,
            'confirmed':sum(r['status']=='confirmed_restart' for r in rows),
            'gradients_available':sum(r.get('gradient') is not None for r in rows),'new_calls_in_collection':0}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


def matrices(inventory,collection):
    from nikasha_pool import choose_rows
    inv=read_json(inventory);new={key(r):r for r in read_json(collection)['rows']};cases={}
    for group in inv['groups']:
        case=next(c for c in read_json(verify(group['collection']))['cases'] if c['case_id']==group['case_id'])
        matrix=copy.deepcopy(case['matrix'])
        for z in ('Ca','La'):
            for q in CANDIDATES:
                cell=matrix[z][q];rr={s:new[group['case_id'],q,z,s] for s in ('vacuum','alpb')}
                cell['status']='complete' if all(r['status']=='confirmed_restart' for r in rr.values()) else 'unavailable'
                cell['archived_low']=cell.pop('low')
                cell['low']={medium:({'status':'complete',**r['actual']} if r['status']=='confirmed_restart'
                                    else {'status':'unavailable','energy_hartree':None,'reason':r.get('reason')})
                             for medium,r in rr.items()}
                cell['MACE_reused']=True;cell['reused']=False
                cell['continuations']=rr
                if cell['status']=='complete':cell['components']={
                    'MACE_eV':case['matrix'][z][q]['components']['MACE_eV'],
                    'GFN2_vacuum_hartree':rr['vacuum']['energy_hartree'],
                    'GFN2_ALPB_hartree':rr['alpb']['energy_hartree']}
                else:cell['components']=None
        cases[group['case_id']]={'old':case,'matrix':matrix,'pool':choose_rows(matrix,list(CANDIDATES))}
    return cases


def compare(stage1,stage2,output):
    from nikasha_pool import choose_rows
    from accommodation_nonlinear import contrast_components, relative_components
    from accommodation_folds_compare import decision
    from accommodation_fold_proposals import outcome
    a=read_json(stage1);b=read_json(stage2);am=read_json(verify(a['manifest']));bm=read_json(verify(b['manifest']))
    if a['stage']!=1 or b['stage']!=2 or bm['previous']!=record(stage1):raise InvalidArtifact('stage comparison differs')
    for field in ('inventory','agreement','settings','reference'):
        if am[field]!=bm[field]:raise InvalidArtifact('stage method/source differs')
    first=matrices(verify(am['inventory']),stage1);second=matrices(verify(am['inventory']),stage2)
    ar={key(r):r for r in a['rows']};br={key(r):r for r in b['rows']}
    ref=read_json(verify(am['reference']));inputs=read_json(verify(am['inputs']))
    labels={r['case_id']:r['known_class'] for r in inputs['cases']}
    bands={'released':ref['old_frozen_bands'],'adaptive':ref['variants']['operational']['bands']}
    rows=[]
    for cid in CASES:
        f=first[cid];s=second[cid];old=f['old'];cell_changes=[];geometry=[]
        for q in CANDIDATES:
            for z in ('Ca','La'):
                for medium in ('vacuum','alpb'):
                    k=(cid,q,z,medium);v1=ar[k]['energy_hartree'];v2=br[k]['energy_hartree']
                    delta=None if v1 is None or v2 is None else (v2-v1)*HA_TO_KCAL
                    cell_changes.append({'candidate':q,'metal':z,'medium':medium,'delta_kcal_mol':delta,
                                         'pass':delta is not None and abs(delta)<=SETTINGS['cell_tolerance_kcal_mol']})
            values={}
            for name,mat in (('old',old['matrix']),('stage1',f['matrix']),('stage2',s['matrix'])):
                values[name]=contrast_components(mat['Ca'][q]['components'],mat['La'][q]['components']) if all(mat[z][q]['status']=='complete' for z in ('Ca','La')) else None
            delta=(values['stage2']['composite_R_model_kcal_mol']-values['stage1']['composite_R_model_kcal_mol']) if values['stage1'] and values['stage2'] else None
            geometry.append({'candidate':q,'contrasts':values,'stage2_minus_stage1_R':delta,
                             'pass':delta is not None and abs(delta)<=SETTINGS['contrast_tolerance_kcal_mol']})
        pool_delta={mode:(s['pool'][mode]['composite_R_model_kcal_mol']-f['pool'][mode]['composite_R_model_kcal_mol']
                         if s['pool']['status']==f['pool']['status']=='available' else None) for mode in ('mathematical','operational')}
        qualified=all(r['pass'] for r in cell_changes+geometry) and all(v is not None and abs(v)<=SETTINGS['pool_tolerance_kcal_mol'] for v in pool_delta.values())
        pools={'old':old['pool'],'stage1':f['pool'],'stage2':s['pool']};calls={}
        for band,limits in bands.items():
            calls[band]={}
            for name,pool in pools.items():
                value=pool['operational']['composite_R_model_kcal_mol'] if pool['status']=='available' else None
                call=decision(value,limits);calls[band][name]={'decision':call,'outcome':outcome(call,labels[cid])}
        delta_old=None;work=None
        if s['pool']['status']=='available':
            delta_old=s['pool']['operational']['composite_R_model_kcal_mol']-old['pool']['operational']['composite_R_model_kcal_mol']
            work={z:relative_components(s['matrix'][z][s['pool']['rows'][z]['operational_candidate']]['components'],
                  old['matrix'][z][old['pool']['rows'][z]['operational_candidate']]['components']) for z in ('Ca','La')}
            if abs(work['Ca']['composite_kcal_mol']-work['La']['composite_kcal_mol']-delta_old)>1e-6:raise InvalidArtifact('work sign/unit mismatch')
        rows.append({'case_id':cid,'expected_class_for_reporting':labels[cid],
             'status':'qualified_within_declared_test' if qualified else 'numerically_unqualified',
             'qualified_R':s['pool']['operational']['composite_R_model_kcal_mol'] if qualified else None,
             'diagnostic_pools':pools,'cell_settling':cell_changes,'same_geometry_contrasts':geometry,
             'pool_stage2_minus_stage1':pool_delta,'stage2_minus_old_R':delta_old,
             'endpoint_work_from_old_pool':work,'old_band_transfer_only':calls,
             'matrix_stage2':s['matrix']})
    result={'protocol_id':PROTOCOL,'stage1':record(stage1),'stage2':record(stage2),'settings':SETTINGS,
            'reference':am['reference'],'own_calibration':None,'case_denominator':4,'cell_denominator_per_stage':80,
            'stage1_confirmed':a['confirmed'],'stage2_confirmed':b['confirmed'],
            'qualified_cases':sum(r['qualified_R'] is not None for r in rows),'cases':rows,
            'no_new_reference_or_affinity_claim':True,'production_changed':False}
    write_new(output,result);return {k:v for k,v in result.items() if k!='cases'}


def main():
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='op',required=True)
    a=sp.add_parser('prepare')
    for name in ('inventory','agreement','output'):a.add_argument('--'+name,required=True,type=Path)
    a.add_argument('--previous',type=Path)
    for name in ('dry-run','execute','collect'):
        a=sp.add_parser(name);a.add_argument('--manifest',required=True,type=Path)
        if name=='collect':a.add_argument('--output',required=True,type=Path)
    a=sp.add_parser('compare')
    for name in ('stage1','stage2','output'):a.add_argument('--'+name,required=True,type=Path)
    args=vars(p.parse_args());op=args.pop('op')
    if op=='dry-run':r=validate(args['manifest'],True)
    elif op=='execute':r=execute_tasks(args['manifest'],validate)
    else:r=globals()[op](**args)
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
