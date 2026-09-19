"""Reusable contextual rigid-water-H preparation; production defaults stay unchanged."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import accepted_attempt,check_atoms,write_xyz
from mace_file_checks import cached_file_checks
from hydration_proposal_opt import SETTINGS as OPTIMIZER
from contextual_water_context import policy,discover,transfer,verify_proposal_geometry

STAGE='contextual_water_orientation'
PROTOCOL='source_context_rigid_water_H_native_omol_preparation_v1'
REQUEST_SCHEMA='alquemia.contextual_water_request.v1'


def request(preparation,policy_source,proposal_source_manifest,scorer_source_manifest,agreement,output,context_group=None):
    """Build hash-pinned request JSON from paths; no source-code editing required."""
    groups={}
    for item in context_group or []:
        if '=' not in item:raise InvalidArtifact('context group must be CASE_ID=GROUP_ID')
        case,group=item.split('=',1)
        if not case or not group or case in groups:raise InvalidArtifact('invalid/duplicate context group')
        groups[case]=group
    cases=[]
    for path in preparation:
        p=read_json(path);case=p['case_id'];cases.append({'case_id':case,'preparation':record(path),
            'context_group':groups.get(case,case)})
    if set(groups)-{c['case_id'] for c in cases}:raise InvalidArtifact('context group names an absent case')
    r={'schema_version':REQUEST_SCHEMA,'cases':cases,'policy_source':record(policy_source),
       'proposal_source_manifest':record(proposal_source_manifest),'scorer_source_manifest':record(scorer_source_manifest),
       'agreement':record(agreement)}
    write_new(output,r);return {'status':'request_prepared','cases':len(cases),'request':record(output)}


def identity_audit(path):
    """An exact identity needs no repeated force-field/protonation construction."""
    p=read_json(path)
    if (p.get('status')!='prepared' or p['explicit_waters']
            or any(a['kind']=='retained_site_water' for a in p['physical_atoms'])):
        raise InvalidArtifact('zero-water identity preparation is inconsistent')
    for key in ('source','source_preparation'):verify(p[key])
    for metal in ('Ca','La'):
        expected=[(metal if a['element']=='M' else a['element'],*a['xyz_A']) for a in p['physical_atoms']]
        if xyz(verify(p['endpoints'][metal]['xyz']))!=expected:
            raise InvalidArtifact('identity endpoints disagree with prepared physical atoms')
    return {'status':'pass','preparation':record(path),'operation':'exact_zero_water_identity',
            'forcefield_reconstruction_performed':False,'new_model_forwards':0}


def completed_result(pin):
    r=read_json(verify(pin));mp=verify(r['manifest']);m=read_json(mp)
    matches=[t for t in m['tasks'] if t['task_id']==r['task_id']]
    if len(matches)!=1 or accepted_attempt(Path(pin['path']).parent,matches[0],mp)!=r:
        raise InvalidArtifact('archived result lacks a compatible successful execution receipt')
    return r,matches[0],m


def reusable_proposal(pin,task,model,software,implementation):
    r,t,m=completed_result(pin)
    if m['model']!=model or m['software']!=software or m.get('optimization')!=OPTIMIZER:
        raise InvalidArtifact('proposal reuse changes model/software/optimizer')
    for name in ('hydration_proposal_opt.py','hydration_mace.py'):
        if m['implementation'][name]['sha256']!=implementation[name]['sha256']:
            raise InvalidArtifact('proposal optimizer implementation changed')
    for key in ('charge','spin_multiplicity','water_groups','mobile_indices','state'):
        if t[key]!=task[key]:raise InvalidArtifact('proposal reuse changes '+key)
    if set(t['starts'])!=set(task['starts']):raise InvalidArtifact('proposal starts differ')
    for seed in task['starts']:
        if xyz(verify(t['starts'][seed]))!=xyz(verify(task['starts'][seed])):
            raise InvalidArtifact('proposal source coordinates differ')
    verify_proposal_geometry(task,r)
    return r


def isolated_dispatch(out,m):
    names=('contextual_water_context.py','contextual_water_prepare.py','contextual_water_score.py',
        'hydration_proposal_opt.py','hydration_mace.py','hydration_network.py','hydration_square.py',
        'hydration_core_transfer.py','hydration_scanner.py','affordable_peptide.py','carve_generic.py','coordination_policy.py')
    for name in names:
        path=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),path);m['implementation'][name]=record(path)
    path=out/'implementation/mace_omol.py';text=path.read_text()
    for name in ('validate','collect'):
        header=f'def {name}(manifest):\n'
        body=f"    if read_json(manifest).get('stage')=='{STAGE}':\n        from contextual_water_prepare import {name} as operation\n        return operation(manifest)\n"
        body+=f"    if read_json(manifest).get('stage')=='contextual_water_scoring':\n        from contextual_water_score import {name} as operation\n        return operation(manifest)\n"
        if text.count(header)!=1:raise InvalidArtifact('unexpected runner dispatch')
        text=text.replace(header,header+body)
    header='def worker(manifest, task_id, output, memory_mode):\n'
    body=f"    if read_json(manifest).get('stage')=='{STAGE}':\n        from hydration_proposal_opt import worker as operation\n        return operation(manifest,task_id,output,memory_mode)\n"
    text=text.replace(header,header+body);path.write_text(text);m['implementation']['mace_omol.py']=record(path)


@cached_file_checks
def prepare(request,output):
    import mace_omol as omol
    from mace_omol_prepared import audit_preparation
    from hydration_network import write_xyz as write_source_xyz
    req=read_json(request)
    if req.get('schema_version')!=REQUEST_SCHEMA or not req.get('cases'):raise InvalidArtifact('explicit prepared-case request required')
    ids=[c['case_id'] for c in req['cases']]
    if len(ids)!=len(set(ids)) or any(not re.fullmatch(r'[A-Za-z0-9_.-]+',i) for i in ids):raise InvalidArtifact('invalid/duplicate case identifiers')
    settings=policy(verify(req['policy_source']));parent=read_json(verify(req['proposal_source_manifest']))
    if parent['optimization']!=OPTIMIZER or parent['model']!=omol.model(verify(parent['software'])):
        raise InvalidArtifact('unsupported native orientation qualification source')
    audits={}
    for c in req['cases']:
        path=verify(c['preparation'])
        audits[c['case_id']]=(audit_preparation(path) if read_json(path)['explicit_waters'] else identity_audit(path))
    contexts=discover(req['cases'],settings)
    _,out,m=omol.common(verify(parent['inventory']),verify(parent['software']),verify(req['agreement']),output,STAGE)
    isolated_dispatch(out,m)
    # Preserve the actual demonstrated kernel; generic request handling lives
    # outside it. Later source-file edits need not force new scientific calls.
    for name in ('hydration_proposal_opt.py','hydration_mace.py'):
        path=out/'implementation'/name;shutil.copyfile(verify(parent['implementation'][name]),path);m['implementation'][name]=record(path)
    all_tasks=[];reused={};case_records=[]
    for case in req['cases']:
        name=case['case_id'];prep=read_json(verify(case['preparation']));entry=copy.deepcopy(case)
        entry['preparation_audit']=audits[name]
        if not prep['explicit_waters']:
            entry.update(operation='exact_zero_water_identity',context=None);case_records.append(entry);continue
        context=contexts[name];directory=out/'contexts'/name;directory.mkdir(parents=True)
        context_path=directory/'context.json';write_new(context_path,{k:v for k,v in context.items() if k!='starts'})
        entry.update(operation='native_context_rigid_water_H',context=record(context_path));case_records.append(entry)
        for metal in ('Ca','La'):
            starts={}
            for seed,rows in context['starts'].items():
                path=directory/(metal+'_'+seed+'.xyz');write_source_xyz(path,[(metal,*rows[0][1:])]+rows[1:],PROTOCOL)
                starts[seed]=record(path)
            task={'task_id':name+'__'+metal,'case_id':name,'metal':metal,'metal_index':0,'kind':'core',
                'variant':'primary','energy_component':omol.COMPONENT,'energy_only':False,
                'charge':context['charges'][metal],'spin_multiplicity':1,'xyz':starts['source'],
                'state':check_atoms(xyz(verify(starts['source'])),context['charges'][metal]),
                'water_groups':context['water_groups'],'mobile_indices':sorted(i for w in context['water_groups'] if w['role']=='variable' for i in w['hydrogen_indices']),
                'starts':starts,'water_orientation_optimization':True,'context':record(context_path),
                'target_preparation':case['preparation']}
            task['cache_key']=cache_key({'task':task,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
            if metal in case.get('reuse_proposals',{}):
                pin=case['reuse_proposals'][metal];reusable_proposal(pin,task,m['model'],m['software'],m['implementation']);reused[task['task_id']]=pin
            all_tasks.append(task)
    m.update(protocol_id=PROTOCOL,request=record(request),settings=settings,optimization=OPTIMIZER,
        cases=case_records,proposal_tasks=all_tasks,reused=reused,tasks=[t for t in all_tasks if t['task_id'] not in reused],
        baseline_changed=False,new_quantum_calls=0,occupancy_probabilities=None)
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


@cached_file_checks
def validate(manifest):
    import mace_omol as omol
    m=read_json(manifest);req=read_json(verify(m['request']))
    if m['stage']!=STAGE or m['protocol_id']!=PROTOCOL or m['optimization']!=OPTIMIZER or m['model']!=omol.model(verify(m['software'])):
        raise InvalidArtifact('context preparation method changed')
    if m['settings']!=policy(verify(req['policy_source'])):raise InvalidArtifact('context policy changed')
    for pin in [m['agreement'],req['proposal_source_manifest'],*m['implementation'].values()]:verify(pin)
    expected=set();lookup={c['case_id']:c for c in m['cases']};definitions={c['case_id']:c for c in req['cases']}
    if len(lookup)!=len(m['cases']) or set(lookup)!={c['case_id'] for c in req['cases']}:raise InvalidArtifact('case coverage changed')
    for c in m['cases']:
        if any(c.get(k)!=v for k,v in definitions[c['case_id']].items()):raise InvalidArtifact('explicit case request changed')
        p=read_json(verify(c['preparation']))
        for z in ('Ca','La'):verify(p['endpoints'][z]['xyz'])
        if c['operation']=='exact_zero_water_identity':
            if p['explicit_waters'] or c['context'] is not None:raise InvalidArtifact('zero-water identity is false')
        else:
            context=read_json(verify(c['context']))
            if context['parent']!=p['source_preparation']:raise InvalidArtifact('different source context')
            expected.update(c['case_id']+'__'+z for z in ('Ca','La'))
    if len(m['proposal_tasks'])!=len(expected) or {t['task_id'] for t in m['proposal_tasks']}!=expected:raise InvalidArtifact('paired proposal coverage changed')
    if len(m['tasks'])!=len(expected-set(m['reused'])) or {t['task_id'] for t in m['tasks']}!=expected-set(m['reused']):raise InvalidArtifact('partial/reused task partition changed')
    for t in m['proposal_tasks']:
        old=xyz(verify(t['xyz']));fixed=set(range(len(old)))-set(t['mobile_indices'])
        for seed,pin in t['starts'].items():
            rows=xyz(verify(pin))
            if check_atoms(rows,t['charge'])!=t['state'] or any(rows[i]!=old[i] for i in fixed):raise InvalidArtifact('seed changed fixed coordinates/state')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('preparation cache changed')
        if t['task_id'] in m['reused']:reusable_proposal(m['reused'][t['task_id']],t,m['model'],m['software'],m['implementation'])
    return {'status':'pass','tasks':len(m['tasks']),'reused_proposals':len(m['reused']),
        'zero_water_identities':sum(c['operation']=='exact_zero_water_identity' for c in m['cases']),
        'new_quantum_calls':0,'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest);m=read_json(manifest);mp=Path(manifest);results={};cases=[]
    for t in m['proposal_tasks']:
        pin=m['reused'].get(t['task_id']);r=None
        if pin:r=reusable_proposal(pin,t,m['model'],m['software'],m['implementation'])
        else:
            for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
                possible=accepted_attempt(a,t,mp)
                if possible is not None:r=possible;pin=record(a/'result.json')
        if r is None:results[t['task_id']]={'status':'unavailable','reason':'no successful proposal receipt'};continue
        try:verify_proposal_geometry(t,r)
        except InvalidArtifact as exc:results[t['task_id']]={'status':'unsupported','reason':str(exc),'result':pin};continue
        results[t['task_id']]={'status':'prepared','result':pin,'proposal_xyz':r['proposed_xyz'],'reused':t['task_id'] in m['reused']}
    for c in m['cases']:
        prep=read_json(verify(c['preparation']));row={'case_id':c['case_id'],'source_preparation':c['preparation'],
            'operation':c['operation'],'explicit_water_count':len(prep['explicit_waters']),'endpoints':{}}
        if c['operation']=='exact_zero_water_identity':
            row.update(status='prepared',endpoints=prep['endpoints'],common_paired_geometry=True)
        else:
            ts={t['metal']:t for t in m['proposal_tasks'] if t['case_id']==c['case_id']}
            if any(results[t['task_id']]['status']!='prepared' for t in ts.values()):
                row.update(status='unavailable',failure='both converged endpoint water preparations required')
            else:
                context=read_json(verify(c['context']));paired=[]
                for metal,t in ts.items():
                    r=results[t['task_id']];coords,mapping=transfer(prep,context,xyz(verify(r['proposal_xyz'])),metal)
                    destination=mp.parent/'prepared'/c['case_id']/(metal+'.xyz');destination.parent.mkdir(parents=True,exist_ok=True)
                    if destination.exists():
                        if xyz(destination)!=coords:raise InvalidArtifact('existing transferred coordinates differ')
                    else:write_xyz(destination,coords)
                    row['endpoints'][metal]={**prep['endpoints'][metal],'xyz':record(destination),'water_transfer':mapping,'proposal':r}
                    paired.append([a[1:] for a in coords])
                row.update(status='prepared',context=c['context'],common_paired_geometry=paired[0]==paired[1],
                    required_scoring_states=['Ca_bound','Ca_detached','La_bound','La_detached'])
        cases.append(row)
    return {'status':'complete' if all(c['status']=='prepared' for c in cases) else 'incomplete','manifest':record(manifest),
        'cases':cases,'proposal_results':results,'baseline_changed':False,'new_quantum_calls':0,
        'occupancy_probabilities':None,'calibrated_decision':None}


def report(collection,output):
    c=read_json(collection)
    if c!=collect(verify(c['manifest'])):raise InvalidArtifact('preparation collection changed')
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    write_new(out/'result.json',c)
    lines=['# Contextual water preparation','',f"Status: {c['status']}. No baseline, water inventory or occupancy change.",'',
        '| Case | Waters | Operation | Status | Common Ca/La geometry |','|---|---:|---|---|---|']
    for r in c['cases']:lines.append(f"| {r['case_id']} | {r['explicit_water_count']} | {r['operation']} | {r['status']} | {r.get('common_paired_geometry')} |")
    (out/'TABLES.md').write_text('\n'.join(lines)+'\n');return {'status':c['status'],'cases':len(c['cases'])}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('request');q.add_argument('--preparation',action='append',required=True);q.add_argument('--context-group',action='append')
    for field in ('policy-source','proposal-source-manifest','scorer-source-manifest','agreement','output'):q.add_argument('--'+field,required=True)
    for op,fields in [('prepare',('request','output')),('dry-run',('manifest',)),('execute',('manifest',)),('collect',('manifest','output')),('report',('collection','output'))]:
        q=sub.add_parser(op)
        for name in fields:q.add_argument('--'+name,required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='dry-run':r=validate(**a)
    elif op=='execute':
        validate(a['manifest']);m=read_json(a['manifest'])
        if m['tasks']:
            run=subprocess.run([sys.executable,str(verify(m['implementation']['mace_hybrid.py'])),'execute','--manifest',str(Path(a['manifest']).resolve()),'--memory-mode','native'],check=True)
            r={'status':'executor_returned_collect_required','returncode':run.returncode}
        else:r={'status':'all_proposals_reused_or_zero_water_identity','new_model_calls':0}
    elif op=='collect':
        dest=a.pop('output');r=collect(**a);write_new(dest,r);r={'status':r['status'],'cases':len(r['cases'])}
    else:r=globals()[op](**a)
    print(json.dumps(r,indent=2))
