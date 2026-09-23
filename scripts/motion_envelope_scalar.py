"""Finite strict-native scalar cells for fixed envelope origin/candidate geometry."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from affordable_workflow import dry_run,execute as run_native
from adaptive_origin_recovery import snapshot
from strict_native_pool import recipe
import compact_solvation as solvent
import consistent_context as context
from structure_informed_starts import scf_details

PROFILE='Nikasha_native_GFN2_ALPB_vacuum_TolE1e10_fresh_MaxIter500_rank1_v1'


def qualified(path):
    q=read_json(path)
    if (q['case_denominator'],q['cell_denominator'],q['agreeing_cells'])!=(32,384,384):raise InvalidArtifact('strict32 native scalar qualification unavailable')
    if q['settings']['TolE_hartree']!=1e-10 or q['counts']['all32']['qualified']!=32:raise InvalidArtifact('strict settings differ')
    return record(path)


def prepare_origins(origins,template,qualification,agreement,output):
    d=read_json(origins);base=read_json(template)
    if d['denominator']!=34:raise InvalidArtifact('fixed34 origin sources required')
    entries=[]
    for r in d['rows']:
        for z,ep in r['native_endpoints'].items():
            if ep['status']!='complete':continue
            entries.append({'case_id':r['case_id'],'source_case_id':r['source_case_id'],'selection_id':r['selection_id'],
                'candidate':'origin','metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity']})
    return prepare(entries,base,record(origins),qualification,agreement,output,'origins',136)


def prepare(entries,base,source,qualification,agreement,output,stage,maximum):
    qualification=qualified(qualification);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=snapshot(Path(__file__).parent,out/'implementation');tasks=[]
    for ep in entries:
        if ep['multiplicity']!=1:raise InvalidArtifact('only actual closed-shell state admitted')
        for medium in ('vacuum','alpb'):
            tid=ep['case_id']+'__'+ep['candidate']+'__'+ep['metal']+'__'+medium;td=out/'tasks'/tid;td.mkdir(parents=True)
            xp=td/'core.xyz';shutil.copyfile(verify(ep['xyz']),xp);ip=td/'endpoint.inp';ip.write_text(recipe(ep['charge'],medium,'fresh'))
            task={**ep,'task_id':tid,'case':ep['case_id'],'representation':'context','medium':medium,'solver':'native',
                'source_xyz':ep['xyz'],'xyz':record(xp),'input':record(ip),'output_path':str(td/'endpoint.out')}
            task['scientific_key']=cache_key({'xyz_sha256':task['xyz']['sha256'],'charge':ep['charge'],'multiplicity':1,
                'medium':medium,'solver':'native','method':solvent.METHOD,'input_body':ip.read_text(),'orca':base['orca']});tasks.append(task)
    if len(tasks)>maximum or len({t['task_id']for t in tasks})!=len(tasks):raise InvalidArtifact('finite cell scope exceeded')
    m={'protocol_id':solvent.PROTOCOL,'method_id':solvent.METHOD,'numerical_policy_id':PROFILE,'envelope_stage':stage,
        'source':source,'qualification':qualification,'agreement':record(agreement),'orca':base['orca'],'implementation':impl,
        'tasks':tasks,'all_tasks':tasks,'reused':{},'maximum_tasks':maximum,'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},
        'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
        'allocated_cpus':32,'allocated_memory_MiB':65536,'GFN2_maxiter':500,'TolE_hartree':1e-10}
    write_new(out/'manifest.json',m);v=validate(out/'manifest.json');write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest):
    m=read_json(manifest);qualified(verify(m['qualification']));verify(m['source']);verify(m['agreement'])
    if (m['numerical_policy_id']!=PROFILE or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}
        or m['tasks']!=m['all_tasks'] or len(m['tasks'])>m['maximum_tasks'] or m['reused']):raise InvalidArtifact('strict finite scalar profile differs')
    source=read_json(verify(m['source']))
    if m['envelope_stage']=='origins':
        expected={(r['case_id'],z):ep for r in source['rows']for z,ep in r['native_endpoints'].items()if ep['status']=='complete'}
        if len(m['tasks'])!=2*len(expected) or m['maximum_tasks']!=136:raise InvalidArtifact('actual origin denominator differs')
    elif m['envelope_stage']=='candidates':
        expected={(t['case_id'],t['metal'],t['candidate']):t for t in source['tasks']}
        if len(m['tasks'])!=2*len(expected) or m['maximum_tasks']!=272:raise InvalidArtifact('actual candidate denominator differs')
    else:raise InvalidArtifact('unknown strict scalar stage')
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        ep=expected[(t['case_id'],t['metal'])if m['envelope_stage']=='origins'else(t['case_id'],t['metal'],t['candidate'])]
        if not context.reusable_state(t,ep) or verify(t['input']).read_text()!=recipe(t['charge'],t['medium'],'fresh'):
            raise InvalidArtifact('strict scalar source/state/recipe differs')
    return dry_run(manifest)


def execute(manifest):
    validate(manifest)
    if int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=65536:raise InvalidArtifact('declared32CPU/64GiB allocation required')
    m=read_json(manifest)
    for t in m['tasks']:
        td=Path(t['output_path']).parent
        if not Path(t['output_path']).exists() and any((td/('endpoint.runtime.'+ext)).exists()for ext in ('gbw','xtbw')):raise InvalidArtifact('fresh task has hidden restart')
    return run_native(manifest)


def endpoint(manifest,task):
    pin=solvent.completed(manifest,task['task_id'])
    if pin is None:return {'status':'unavailable','energy_hartree':None,'reason':'no valid completed actual scalar receipt'}
    try:
        audit=solvent.diagnostics(pin,task);text=verify(pin['output']).read_text();details=scf_details(text)
        tolerance=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
        receipt=read_json(verify(pin['receipt']))
        if tolerance!=['1.0000e-10'] and not(len(tolerance)==1 and float(tolerance[0])==1e-10):raise InvalidArtifact('actual native TolE differs')
        if receipt['parallelism']['nprocs']!=1 or not details['native_mixer_observed'] or audit['charge_sanity_status']!='pass':raise InvalidArtifact('actual strict state/rank differs')
        return {**pin,'status':'complete','component_audit':audit,'SCF':details,'observed_TolE_hartree':1e-10}
    except (OSError,InvalidArtifact,KeyError,ValueError) as exc:return {'status':'unavailable','energy_hartree':None,'reason':str(exc),'actual':pin}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);rows=[{k:t[k]for k in ('task_id','case_id','candidate','metal','medium')}|endpoint(manifest,t)for t in m['tasks']]
    d={'protocol_id':PROFILE,'manifest':record(manifest),'rows':rows,'denominator':len(rows),'complete':sum(r['status']=='complete'for r in rows),'new_calls_in_collection':0}
    write_new(output,d);return {'collection':record(output),'complete':d['complete'],'denominator':len(rows)}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    for op,keys in {'prepare_origins':('origins','template','qualification','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        q=s.add_parser(op)
        for k in keys:q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());cmd=a.pop('command');print(json.dumps(globals()[cmd](**a),indent=2))
