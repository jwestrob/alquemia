"""Finite original-state endpoints for the approved three-La-membership pilot."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import time
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from adaptive_origin_recovery import snapshot
import consistent_context as context
import compact_solvation as solvent
from compact_solvation_compare import native_endpoint
from affordable_workflow import dry_run,execute as run_native
from slsqp_precision_transfer import qualification
from union_triple_preparation import STRESS

PROTOCOL='Nikasha_three_La_membership_canonical28_stress_origins_v1'


def prepare(canonical,canonical_reuse,primary,primary_reuse,pool_reuse,template,rank_qualification,agreement,output):
    cp=read_json(canonical);cr=read_json(canonical_reuse);pp=read_json(primary);pr=read_json(primary_reuse);reuse=read_json(pool_reuse)
    if (cp['denominator'],cp['supported'],reuse['complete_pool_reuses'],reuse['new_pools'])!=(28,28,19,9):
        raise InvalidArtifact('approved canonical preparation and19 full-pool reuses required')
    if cr['preparation']!=record(canonical) or pr['preparation']!=record(primary) or reuse['canonical_preparation']!=record(canonical):
        raise InvalidArtifact('origin/full-pool audit source differs')
    rank=qualification(rank_qualification);base=read_json(template)
    selected=cp['cases']+[next(c for c in pp['cases'] if c['case_id']==STRESS and c['is_separate_stress_probe'])]
    audit={r['case_id']:r for r in cr['rows']};audit[STRESS]=next(r for r in pr['rows'] if r['case_id']==STRESS and r['is_separate_stress_probe'])
    primary_selection=read_json(verify(cp['selection']));cfg=primary_selection['config']
    if cfg['model']!=base['model'] or cfg['software']!=base['software']:raise InvalidArtifact('actual electronic runtime differs')
    new=[c for c in selected if audit[c['case_id']]['missing_native_endpoints']]
    if len(new)!=6 or any(audit[c['case_id']]['missing_native_endpoints']!=['Ca','La'] or len(audit[c['case_id']]['missing_solvent_cells'])!=4 for c in new):
        raise InvalidArtifact('expected exactly six wholly missing origin pairs')
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces/')
    out.mkdir(parents=True,exist_ok=False);impl=snapshot(Path(__file__).parent,out/'implementation')
    inputs={'protocol_id':PROTOCOL,'canonical':record(canonical),'canonical_reuse':record(canonical_reuse),'primary':record(primary),
        'primary_reuse':record(primary_reuse),'pool_reuse':record(pool_reuse),'template':record(template),'rank_qualification':rank,
        'agreement':record(agreement),'cases':selected,'reuse_rows':audit,'denominator':29,'new_origin_sources':[c['case_id'] for c in new],
        'config':cfg,'new_origin_MACE_calls':12,'new_origin_GFN2_calls':24,'implementation':impl,'production_changed':False}
    ip=out/'INPUTS.json';write_new(ip,inputs)
    mm={'protocol_id':'mace_omol_0_100m_vacuum_descriptor_v1','stage':PROTOCOL,'preparation':record(ip),
        'agreement':record(agreement),'software':cfg['software'],'model':cfg['model'],'implementation':impl,
        'tasks':[],'reused':{},'verification_references':{}}
    tasks=[]
    for c in new:
        for z,ep in c['representations']['context']['endpoints'].items():
            t={'task_id':c['case_id']+'__context__'+z,'case_id':c['case_id'],'representation':'context','metal':z,'metal_index':0,
                'kind':'core','xyz':ep['xyz'],'charge':ep['charge'],'spin_multiplicity':ep['multiplicity'],
                'energy_component':'MACE_OMOL_total_vacuum_energy'}
            t['cache_key']=cache_key({'task':t,'model':mm['model'],'software':mm['software'],'implementation':impl});mm['tasks'].append(t)
            for medium in ('vacuum','alpb'):
                tid=t['task_id']+'__'+medium;td=out/'solvent'/'tasks'/tid;td.mkdir(parents=True)
                xp=td/'core.xyz';shutil.copyfile(verify(ep['xyz']),xp);ipath=td/'endpoint.inp'
                body=solvent.input_text(ep['charge'],ep['multiplicity'],medium,'native').replace('%scf\n','%scf\n MaxIter 500\n');ipath.write_text(body)
                task={'task_id':tid,'case_id':c['case_id'],'case':c['case_id'],'metal':z,'medium':medium,'solver':'native',
                    'representation':'context','charge':ep['charge'],'multiplicity':ep['multiplicity'],'xyz':record(xp),
                    'input':record(ipath),'output_path':str(td/'endpoint.out')}
                task['scientific_key']=cache_key({'xyz_sha256':task['xyz']['sha256'],'charge':ep['charge'],'multiplicity':ep['multiplicity'],
                    'medium':medium,'solver':'native','method':solvent.METHOD,'input_body':body,'orca':base['orca']});tasks.append(task)
    sm={'protocol_id':solvent.PROTOCOL,'method_id':solvent.METHOD,'pilot_protocol_id':PROTOCOL,'inputs':record(ip),
        'agreement':record(agreement),'orca':base['orca'],'implementation':impl,'rank_qualification':rank,
        'tasks':tasks,'all_tasks':tasks,'reused':{},'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},
        'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
        'GFN2_maxiter':500,'allocated_cpus':32,'allocated_memory_MiB':65536}
    write_new(out/'mace'/'manifest.json',mm);write_new(out/'solvent'/'manifest.json',sm)
    check=validate(out/'solvent'/'manifest.json');write_new(out/'PREFLIGHT.json',check)
    ready={'inputs':record(ip),'MACE_manifest':record(out/'mace'/'manifest.json'),'GFN2_manifest':record(out/'solvent'/'manifest.json'),
        'new_origin_MACE_calls':12,'new_origin_GFN2_calls':24,'cpu_python':base['cpu_python'],'gpu_python':base['gpu_python']}
    write_new(out/'READY.json',ready);return ready


def validate(manifest):
    m=read_json(manifest);i=read_json(verify(m['inputs']));qualification(verify(m['rank_qualification']))
    if (m['pilot_protocol_id']!=PROTOCOL or len(m['tasks'])!=24 or m['tasks']!=m['all_tasks'] or
            m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32} or i['new_origin_MACE_calls']!=12):
        raise InvalidArtifact('finite origin task scope differs')
    for p in m['implementation'].values():verify(p)
    for t in m['tasks']:
        c=next(c for c in i['cases'] if c['case_id']==t['case_id']);ep=c['representations']['context']['endpoints'][t['metal']]
        if t['case_id'] not in i['new_origin_sources'] or not context.reusable_state(ep,t):raise InvalidArtifact('origin source/state differs')
        if verify(t['input']).read_text()!=solvent.input_text(t['charge'],t['multiplicity'],t['medium'],'native').replace('%scf\n','%scf\n MaxIter 500\n'):
            raise InvalidArtifact('native origin recipe differs')
    return dry_run(manifest)


def execute(manifest):
    validate(manifest)
    if int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=65536:
        raise InvalidArtifact('declared32CPU/64GiB allocation required')
    return run_native(manifest)


def collect(stage,output):
    start=time.monotonic();stage=Path(stage);ready=read_json(stage/'READY.json');i=read_json(verify(ready['inputs']))
    mm=read_json(verify(ready['MACE_manifest']));sm=read_json(verify(ready['GFN2_manifest']));validate(verify(ready['GFN2_manifest']))
    rows=[]
    for c in i['cases']:
        cid=c['case_id'];prior=i['reuse_rows'][cid];r={'case_id':cid,'source':c,'native_endpoints':{},'solvent_endpoints':{}}
        for z in ('Ca','La'):
            ep=c['representations']['context']['endpoints'][z]
            try:
                pin=prior['native_endpoints'][z]['receipt'] if z in prior['native_endpoints'] else record(verify(ready['MACE_manifest']).parent/'results'/(cid+'__context__'+z)/'result.json')
                native=native_endpoint(read_json(verify(pin)),mm['model'])
                if not context.reusable_state(ep,native):raise InvalidArtifact('actual origin native state differs')
                r['native_endpoints'][z]={**native,'status':'complete','reused':z in prior['native_endpoints']}
            except (OSError,InvalidArtifact,KeyError) as exc:r['native_endpoints'][z]={'status':'unavailable','reason':str(exc)}
            r['solvent_endpoints'][z]={}
            for medium in ('vacuum','alpb'):
                try:
                    pin=prior['solvent_endpoints'].get(z,{}).get(medium)
                    reused=bool(pin)
                    if not pin:pin=solvent.completed(verify(ready['GFN2_manifest']),cid+'__context__'+z+'__'+medium)
                    if pin is None:raise InvalidArtifact('origin native solvent unavailable')
                    lm=read_json(verify(pin['manifest']));lt=next(t for t in lm.get('all_tasks',lm['tasks']) if t['task_id']==pin['task_id'])
                    if not context.reusable_state(ep,lt):raise InvalidArtifact('actual low-level origin geometry/state differs')
                    audit=solvent.diagnostics(pin,lt)
                    if audit['charge_sanity_status']!='pass':raise InvalidArtifact('native state audit failed')
                    r['solvent_endpoints'][z][medium]={**pin,'status':'complete','component_audit':audit,'reused':reused}
                except (OSError,InvalidArtifact,KeyError) as exc:r['solvent_endpoints'][z][medium]={'status':'unavailable','reason':str(exc)}
        r['status']='complete' if all(r['native_endpoints'][z]['status']=='complete' and all(r['solvent_endpoints'][z][s]['status']=='complete' for s in ('vacuum','alpb')) for z in ('Ca','La')) else 'unavailable'
        r['source_preparation']=c['representations']['context']['preparation'];rows.append(r)
    result={'protocol_id':PROTOCOL,'ready':record(stage/'READY.json'),'inputs':ready['inputs'],'rows':rows,'denominator':29,
        'complete':sum(r['status']=='complete' for r in rows),'new_calls_in_collection':0,'wall_seconds':time.monotonic()-start}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    commands={'prepare':('canonical','canonical_reuse','primary','primary_reuse','pool_reuse','template','rank_qualification','agreement','output'),
        'validate':('manifest',),'execute':('manifest',),'collect':('stage','output')}
    for op,fields in commands.items():
        q=s.add_parser(op)
        for f in fields:q.add_argument('--'+f.replace('_','-'),required=True)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
