"""Promoted, versioned water preparation around the existing ORCA workflow.

Entry point: affordable_workflow.py baseline. Original inputs/results stay intact.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from affordable_common import (InvalidArtifact, HA_TO_KCAL, cache_key, classify_raw,
    digest, energy, read_json, record, verify, write_new, xyz)
import contextual_water_prepare as water
from hydration_core_transfer import transfer_water_hydrogens, PROTOCOL as PREPARED_PROTOCOL
from hydration_network import write_xyz

WORKFLOW = 'baseline_contextual_water_v1'
CANONICAL = 'pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'
DEFAULT_RELEASE = Path(__file__).resolve().parents[1]/'params/baseline_water_v1.json'


def freeze_code(output):
    directory=Path(output)/'implementation';directory.mkdir()
    result={}
    for name in ('baseline_water.py','affordable_workflow.py','hydration_core_transfer.py',
                 'hydration_network.py','contextual_water_prepare.py','contextual_water_context.py'):
        target=directory/name;shutil.copyfile(Path(__file__).with_name(name),target);result[name]=record(target)
    return result


def release(path):
    r=read_json(path)
    if r['workflow_id']!=WORKFLOW or r['default_water_policy']!='contextual_if_supported':
        raise InvalidArtifact('unsupported baseline water release')
    for pin in r['artifacts'].values(): verify(pin)
    return r


def request(preparation, output, release_path=DEFAULT_RELEASE, context_group=None,
            reuse_request=None, water_policy='contextual_if_supported'):
    r=release(release_path); refs=r['artifacts']
    water.request(preparation,verify(refs['policy']),verify(refs['proposals']),
                  verify(refs['scorer']),verify(refs['agreement']),output,context_group)
    p=Path(output); req=read_json(p)
    if reuse_request:
        old=read_json(reuse_request); old_cases={c['case_id']:c for c in old['cases']}
        for c in req['cases']:
            prior=old_cases.get(c['case_id'])
            if prior is None or prior['preparation']!=c['preparation']:
                raise InvalidArtifact('reuse request must contain the identical prepared case')
            if prior.get('context_group',prior['case_id'])!=c.get('context_group',c['case_id']):
                raise InvalidArtifact('reuse request changes context grouping')
            for key in ('reuse_proposals','reuse_scores','disconnected_reference'):
                if key in prior:c[key]=prior[key]
        req['reuse_request']=record(reuse_request)
    req.update(workflow_id=WORKFLOW,release=record(release_path),water_policy=water_policy)
    # This request is not exposed to an executor until the complete file is written.
    p.write_text(json.dumps(req,indent=2,sort_keys=True)+'\n')
    return {'status':'request_prepared','request':record(p),'water_policy':water_policy}


def prepare(request, output):
    req=read_json(request); r=release(verify(req['release']))
    if req['workflow_id']!=WORKFLOW or req['water_policy'] not in ('contextual_if_supported','original'):
        raise InvalidArtifact('unsupported requested workflow')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    cases=[]
    for c in req['cases']:
        p=read_json(verify(c['preparation']))
        cases.append({'case_id':c['case_id'],'preparation':c['preparation'],
            'source_core':p['source_preparation'],'water_count':len(p['explicit_waters']),
            'operation':'original_reference' if req['water_policy']=='original' else
                'contextual_water_H' if p['explicit_waters'] else 'exact_dry_identity'})
    pin=None
    if req['water_policy']=='contextual_if_supported':
        water.prepare(request,out/'water')
        pin=record(out/'water/manifest.json')
    else:
        from mace_omol_prepared import audit_preparation
        for c in cases:
            (audit_preparation if c['water_count'] else water.identity_audit)(verify(c['preparation']))
    plan={'workflow_id':WORKFLOW,'request':record(request),'release':req['release'],
          'water_policy':req['water_policy'],'cases':cases,'water_manifest':pin,
          'original_records_preserved':True,'implementation':freeze_code(out)}
    write_new(out/'plan.json',plan)
    return {'status':'prepared','plan':record(out/'plan.json'),'cases':len(cases),
            'water_policy':req['water_policy'],'next_stage':'water-execute' if pin else 'prepare-dft'}


def checked_plan(path):
    p=read_json(path);req=read_json(verify(p['request']));release(verify(p['release']))
    for pin in p['implementation'].values():verify(pin)
    if p['workflow_id']!=WORKFLOW or p['release']!=req['release'] or p['water_policy']!=req['water_policy']:
        raise InvalidArtifact('baseline plan differs from its request')
    expected=[]
    for c in req['cases']:
        w=read_json(verify(c['preparation']))
        expected.append({'case_id':c['case_id'],'preparation':c['preparation'],
            'source_core':w['source_preparation'],'water_count':len(w['explicit_waters']),
            'operation':'original_reference' if req['water_policy']=='original' else
                'contextual_water_H' if w['explicit_waters'] else 'exact_dry_identity'})
    if p['cases']!=expected:
        raise InvalidArtifact('baseline prepared cases changed')
    for c in p['cases']:verify(c['source_core']);verify(c['preparation'])
    if p['water_manifest']:water.validate(verify(p['water_manifest']))
    return p


def water_execute(plan):
    p=checked_plan(plan)
    if not p['water_manifest']:return {'status':'original_reference_no_water_tasks'}
    manifest=verify(p['water_manifest']);m=read_json(manifest)
    if not m['tasks']:return {'status':'water_preparation_reused_or_dry','new_model_calls':0}
    subprocess.run([sys.executable,str(Path(water.__file__).resolve()),'execute','--manifest',str(manifest)],check=True)
    return {'status':'water_executor_finished_collect_required'}


def recipe(path):
    body='\n'.join(line.split('#',1)[0] for line in Path(path).read_text().splitlines())
    body=re.sub(r'(?mi)^\s*%maxcore\s+\d+\s*$','',body)
    rest=re.sub(r'(?mi)^\s*\*\s+xyzfile\s+[^\n]+$','',body).strip()
    if rest.lower().split()!=['!','r2scan-3c','noautostart','cpcm(water)','defgrid3']:
        raise InvalidArtifact('baseline requires unchanged native r2SCAN-3c/CPCM SP recipe')
    return rest


def completed(manifest, task_id):
    from run_orca_task_manifest import load_manifest_tasks,_completed_attempt_is_valid
    mp=Path(manifest);m,ts=load_manifest_tasks(mp)
    t=next(t for t in ts if t['task_id']==task_id)
    op=t['output'];rp=Path(str(op)+'.execution.json')
    if not op.exists() or not rp.exists():return None
    if not _completed_attempt_is_valid(rp,op,manifest_sha256=digest(mp),task=t,
            runner_identity=m['execution_policy']['task_runner'],
            runtime_renderer_identity=m['execution_policy']['runtime_renderer']):return None
    if read_json(rp)['orca_version']!='6.1.1':return None
    return {'energy_hartree':energy(op),'output':record(op),'receipt':record(rp),
            'manifest':record(mp),'task_id':task_id}


def reusable(task, manifests):
    from run_orca_task_manifest import load_manifest_tasks
    for path in manifests:
        m,old_tasks=load_manifest_tasks(Path(path))
        for old in old_tasks:
            xp=verify({'path':str(old['xyz']),'sha256':old['xyz_sha256']})
            ip=verify({'path':str(old['input']),'sha256':old['input_sha256']})
            coords=xyz(xp)
            if coords[0][0]!=task['metal'] or coords!=xyz(verify(task['xyz'])):continue
            state=re.search(r'(?mi)^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+',ip.read_text())
            if state is None or (int(state[1]),int(state[2]))!=(task['charge'],task['multiplicity']):continue
            if recipe(ip).lower().split()!=recipe(verify(task['input'])).lower().split():continue
            result=completed(path,old['task_id'])
            if result:return result
    return None


def prepare_dft(plan, output, reuse_manifest=None):
    p=checked_plan(plan);r=release(verify(p['release']));w=None
    if p['water_manifest']:
        w=water.collect(verify(p['water_manifest']))
        if w['status']!='complete':raise InvalidArtifact('water preparation incomplete; no baseline substitution')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    all_tasks=[];reused={};case_records=[]
    for c in p['cases']:
        whole=read_json(verify(c['preparation']));core=read_json(verify(c['source_core']))
        eps=whole['source_audit_row']['endpoints']
        wet=c['operation']=='contextual_water_H'
        wc=next(x for x in w['cases'] if x['case_id']==c['case_id']) if wet else None
        contexts=read_json(verify(wc['context'])) if wet else None
        row={**c,'original_protocol_id':core['protocol_id'],
             'selected_arm':'water_prepared' if wet else 'original',
             'prepared_protocol_id':PREPARED_PROTOCOL if wet else core['protocol_id']}
        for arm in (('original','water_prepared') if wet else ('original',)):
            for z in ('Ca','La'):
                e=eps[z];body=recipe(verify(e['input']));atoms=xyz(verify(e['xyz']));mobile=[]
                if e.get('multiplicity',1)!=1:raise InvalidArtifact('unsupported baseline multiplicity')
                if wet and arm=='water_prepared':
                    proposal=wc['endpoints'][z]['proposal']['proposal_xyz']
                    atoms,mobile=transfer_water_hydrogens(core,atoms,contexts,xyz(verify(proposal)))
                d=out/'tasks'/(c['case_id']+'__'+arm+'__'+z);d.mkdir(parents=True)
                xp=d/'core.xyz';write_xyz(xp,atoms,PREPARED_PROTOCOL if arm=='water_prepared' else core['protocol_id'])
                ip=d/'endpoint.inp';ip.write_text(body+f'\n* xyzfile {e["charge"]} 1 core.xyz\n')
                task={'task_id':d.name,'case':c['case_id'],'arm':arm,'metal':z,'charge':e['charge'],
                    'multiplicity':1,'input':record(ip),'xyz':record(xp),'output_path':str(d/'endpoint.out'),
                    'source_core':c['source_core'],'original_xyz':e['xyz'],'original_input':e['input'],
                    'water_H_indices':mobile,'protocol_id':PREPARED_PROTOCOL if arm=='water_prepared' else core['protocol_id']}
                task['cache_key']=cache_key({'task':task,'release':p['release']})
                all_tasks.append(task)
                prior=reusable(task,reuse_manifest or [])
                if prior:reused[task['task_id']]=prior
        case_records.append(row)
    manifest={'workflow_id':WORKFLOW,'protocol_id':WORKFLOW,'plan':record(plan),'release':p['release'],
        'agreement':r['artifacts']['agreement'],'orca':r['artifacts']['orca'],
        'execution_policy':{'task_runner':r['artifacts']['task_runner'],'runtime_renderer':r['artifacts']['runtime_renderer']},
        'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},
        'cases':case_records,'all_tasks':all_tasks,'reused':reused,
        'tasks':[t for t in all_tasks if t['task_id'] not in reused],
        'water_collection':w,'implementation':freeze_code(out),'original_records_preserved':True}
    write_new(out/'manifest.json',manifest)
    return dry_run(out/'manifest.json')


def dry_run(manifest):
    m=read_json(manifest);p=checked_plan(verify(m['plan']));r=release(verify(m['release']))
    for pin in m['implementation'].values():verify(pin)
    if m['workflow_id']!=WORKFLOW:raise InvalidArtifact('unrecognized baseline manifest')
    if m['release']!=p['release'] or m['orca']!=r['artifacts']['orca']:
        raise InvalidArtifact('baseline release/executable changed')
    if m['water_collection']!=(water.collect(verify(p['water_manifest'])) if p['water_manifest'] else None):
        raise InvalidArtifact('baseline water preparation collection changed')
    for c,source in zip(m['cases'],p['cases']):
        core=read_json(verify(source['source_core']))
        wet=source['operation']=='contextual_water_H'
        expected_case={**source,'original_protocol_id':core['protocol_id'],
            'selected_arm':'water_prepared' if wet else 'original',
            'prepared_protocol_id':PREPARED_PROTOCOL if wet else core['protocol_id']}
        if c!=expected_case:raise InvalidArtifact('baseline case routing changed')
    if len(m['cases'])!=len(p['cases']):raise InvalidArtifact('baseline case denominator changed')
    expected={c['case_id']+'__'+a+'__'+z for c in m['cases']
        for a in (('original','water_prepared') if c['selected_arm']=='water_prepared' else ('original',)) for z in ('Ca','La')}
    if len(m['all_tasks'])!=len(expected) or {t['task_id'] for t in m['all_tasks']}!=expected:
        raise InvalidArtifact('baseline original/prepared endpoint coverage changed')
    if m['tasks']!=[t for t in m['all_tasks'] if t['task_id'] not in m['reused']]:
        raise InvalidArtifact('baseline cache/task partition changed')
    for t in m['all_tasks']:
        atoms=xyz(verify(t['xyz']));old=xyz(verify(t['original_xyz']));recipe(verify(t['input']))
        case=next(c for c in m['cases'] if c['case_id']==t['case'])
        if t['task_id']!=t['case']+'__'+t['arm']+'__'+t['metal'] or t['arm'] not in ('original','water_prepared'):
            raise InvalidArtifact('baseline endpoint identity changed')
        protocol=case['original_protocol_id'] if t['arm']=='original' else case['prepared_protocol_id']
        if t['protocol_id']!=protocol or t['source_core']!=case['source_core']:
            raise InvalidArtifact('baseline protocol or core changed')
        whole=read_json(verify(case['preparation']));e=whole['source_audit_row']['endpoints'][t['metal']]
        if t['original_xyz']!=e['xyz'] or t['original_input']!=e['input'] or t['charge']!=e['charge'] or t['multiplicity']!=1:
            raise InvalidArtifact('baseline source or electronic state changed')
        state=re.findall(r'(?mi)^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)\s*$',verify(t['input']).read_text())
        if state!=[(str(t['charge']),str(t['multiplicity']),'core.xyz')]:
            raise InvalidArtifact('baseline ORCA input state changed')
        if recipe(verify(t['input'])).lower().split()!=recipe(verify(e['input'])).lower().split():
            raise InvalidArtifact('baseline Hamiltonian changed')
        expected_atoms=old;expected_mobile=[]
        if t['arm']=='water_prepared':
            wc=next(c for c in m['water_collection']['cases'] if c['case_id']==t['case'])
            expected_atoms,expected_mobile=transfer_water_hydrogens(read_json(verify(case['source_core'])),old,
                read_json(verify(wc['context'])),xyz(verify(wc['endpoints'][t['metal']]['proposal']['proposal_xyz'])))
        expected_atoms=[(a[0],*(float(f'{x:.10f}') for x in a[1:])) for a in expected_atoms]
        if atoms!=expected_atoms or t['water_H_indices']!=expected_mobile:
            raise InvalidArtifact('baseline coordinates do not match the verified water transfer')
        if len(atoms)!=len(old) or [a[0] for a in atoms]!=[a[0] for a in old]:raise InvalidArtifact('core composition changed')
        mobile=set(t['water_H_indices'])
        if any(atoms[i]!=old[i] for i in range(len(old)) if i not in mobile):raise InvalidArtifact('nonwater-H source coordinates changed')
        if any(atoms[i][0]!='H' for i in mobile):raise InvalidArtifact('nonhydrogen mobile atom')
        if t['cache_key']!=cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},'release':m['release']}):
            raise InvalidArtifact('baseline scientific cache changed')
        pin=m['reused'].get(t['task_id'])
        if pin and (completed(verify(pin['manifest']),pin['task_id'])!=pin or
                    reusable(t,[verify(pin['manifest'])])!=pin):raise InvalidArtifact('invalid archived quantum reuse')
    if m['tasks']:
        from affordable_workflow import dry_run as native_dry_run
        native_dry_run(manifest)
    return {'status':'dry_run_pass','manifest':record(manifest),'new_DFT_endpoints':len(m['tasks']),
            'reused_DFT_endpoints':len(m['reused']),'cases':len(m['cases'])}


def execute(manifest):
    result=dry_run(manifest)
    if not result['new_DFT_endpoints']:return {'status':'all_DFT_reused','new_DFT_endpoints':0}
    from affordable_workflow import execute as native_execute
    return native_execute(manifest)


def prepare_retry(manifest, output):
    """Fresh workspace for missing/partial attempts; reuse verified successes."""
    dry_run(manifest);m=read_json(manifest)
    sources=sorted({str(verify(pin['manifest'])) for pin in m['reused'].values()})
    if m['tasks']:sources.append(str(Path(manifest).resolve()))
    return prepare_dft(verify(m['plan']),output,sources)


def collect(manifest, output):
    dry_run(manifest);m=read_json(manifest);r=release(verify(m['release']));rows={};cases=[]
    for t in m['all_tasks']:
        result=m['reused'].get(t['task_id']) or completed(manifest,t['task_id'])
        rows[t['task_id']]={'status':'complete',**result} if result else {'status':'unavailable','energy_hartree':None}
    for c in m['cases']:
        item={'case_id':c['case_id'],'operation':c['operation'],'selected_arm':c['selected_arm'],
              'original':None,'water_prepared':None,'selected_score':None,'preparation_delta_R_kcal_mol':None}
        for arm in ('original','water_prepared'):
            ids={z:c['case_id']+'__'+arm+'__'+z for z in ('Ca','La')}
            if not all(i in rows for i in ids.values()):continue
            eps={z:rows[i] for z,i in ids.items()}
            proto=c['original_protocol_id'] if arm=='original' else c['prepared_protocol_id']
            score={'status':'unavailable','protocol_id':proto,'endpoints':eps,'R_kcal_mol':None,
                   'S_kcal_mol':None,'decision':'unavailable','reference_status':'unavailable_for_this_protocol'}
            if all(e['status']=='complete' for e in eps.values()):
                contrast=(eps['Ca']['energy_hartree']-eps['La']['energy_hartree'])*HA_TO_KCAL
                score.update(status='complete',R_kcal_mol=contrast,decision='uncalibrated_protocol')
                if proto==CANONICAL:
                    ref=read_json(verify(r['artifacts']['canonical_release']));g=ref['aquo_reporting_gauge']
                    verify({'path':g['path'],'sha256':g['sha256']})
                    score.update(S_kcal_mol=contrast-g['A_kcal_mol'],decision=classify_raw(contrast,ref,CANONICAL),
                                 reference_status='unchanged_released_canonical_protocol')
            item[arm]=score
        selected=item[c['selected_arm']]
        item['status']=selected['status'];item['selected_score']=selected
        if item['original'] and item['water_prepared'] and all(item[a]['status']=='complete' for a in ('original','water_prepared')):
            item['preparation_delta_R_kcal_mol']=item['water_prepared']['R_kcal_mol']-item['original']['R_kcal_mol']
        cases.append(item)
    result={'workflow_id':WORKFLOW,'manifest':record(manifest),'status':'complete' if all(x['status']=='complete' for x in rows.values()) else 'incomplete',
            'cases':cases,'original_records_preserved':True,'occupancy_probabilities':None,'entropy_correction':None}
    write_new(output,result);return result


def report(collection, output):
    c=read_json(collection);verify(c['manifest'])
    if c['workflow_id']!=WORKFLOW:raise InvalidArtifact('not a promoted baseline collection')
    lines=['# Baseline with contextual water preparation','',
        'Original and prepared protocols remain separate. Changed wet sites have no inherited absolute bands.','',
        '| Case | Selected arm | Original R (kcal/mol) | Prepared R (kcal/mol) | Status |',
        '|---|---|---:|---:|---|']
    for row in c['cases']:
        original=(row['original'] or {}).get('R_kcal_mol')
        prepared=(row['water_prepared'] or {}).get('R_kcal_mol')
        lines.append(f'| {row["case_id"]} | {row["selected_arm"]} | {original} | {prepared} | {row["status"]} |')
    destination=Path(output);destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('x') as f:f.write('\n'.join(lines)+'\n')
    return {'status':c['status'],'report':record(destination),'collection':record(collection)}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('request');q.add_argument('--preparation',action='append',required=True)
    q.add_argument('--context-group',action='append');q.add_argument('--reuse-request')
    q.add_argument('--release-path',default=str(DEFAULT_RELEASE));q.add_argument('--output',required=True)
    q.add_argument('--water-policy',choices=('contextual_if_supported','original'),default='contextual_if_supported')
    for op,fields in [('prepare',('request','output')),('water-execute',('plan',)),
                     ('prepare-dft',('plan','output')),('prepare-retry',('manifest','output')),('dry-run',('manifest',)),
                     ('execute',('manifest',)),('collect',('manifest','output')),('report',('collection','output'))]:
        q=s.add_parser(op)
        for f in fields:q.add_argument('--'+f,required=True)
        if op=='prepare-dft':q.add_argument('--reuse-manifest',action='append')
    a=vars(p.parse_args(argv));op=a.pop('op').replace('-','_');result=globals()[op](**a)
    print(json.dumps({k:v for k,v in result.items() if k!='cases' or not isinstance(v,list)},indent=2))


if __name__=='__main__':main()
