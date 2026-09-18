"""Explicit conductor-like model on the frozen real protein source inventory."""
from pathlib import Path
import argparse
import json
import os
import shutil
import time
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_ddx_recovery import preflight as parent_preflight,execute_group

PROTOCOL='fixed_source_full_protein_ddCPCM_component_v1'
REFINED_PROTOCOL='fixed_source_full_protein_ddCPCM_component_refinement_v2'
REFINED_GRIDS=dict(coarse=dict(lmax=12,n_lebedev=590),primary=dict(lmax=18,n_lebedev=974),refined=dict(lmax=24,n_lebedev=2030))


def reuse_inventory(g,collection,pin):
    if g['variant']!='primary' or g['grid']!='coarse':return {}
    old=collection['groups'][g['case_id']+'_primary_refined']
    if old['source']!=g['source'] or any(old['rows'][s]['status']!='computed' for s in ('Ca','La')):raise InvalidArtifact('completed matching coarse-reference states required')
    return {s:dict(kind='original',receipt=pin,group=old['group_id'],state=s,arrays=old['arrays'],maxiter=1200) for s in ('Ca','La')}


def group_key(g,m):
    return cache_key(dict(group={k:v for k,v in g.items() if k!='cache_key'},config=m['config'],model=m['model'],
                         energy_prefactor=m['energy_prefactor'],implementation=m['implementation']))


def prepare(config,output):
    c=read_json(config);protocol=c.get('protocol',PROTOCOL)
    if protocol==PROTOCOL:parent=parent_preflight(verify(c['parent']));prior=None
    elif protocol==REFINED_PROTOCOL:
        parent=preflight(verify(c['parent']));prior=read_json(verify(c['prior_collection']))
        if parent['protocol']!=PROTOCOL or prior['manifest']!=c['parent'] or not prior['complete']:raise InvalidArtifact('matching first conductor collection required')
    else:raise InvalidArtifact('unknown conductor protocol')
    verify(c['plan'])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in parent['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name in ('mace_ddx_cpcm.py','mace_ddx_recovery.py','mace_ddx_source_report.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    model={**parent['model'],'model':'cosmo'};eps=model['solvent_epsilon']
    m=dict(protocol=protocol,config=record(config),parent=c['parent'],plan=c['plan'],python=parent['python'],module=parent['module'],
           model=model,energy_prefactor=(eps-1)/eps,grids=REFINED_GRIDS if prior else parent['grids'],tolerances=parent['tolerances'],
           implementation={p.name:record(p) for p in impl.glob('*.py')},groups=[],requested_state_roles=20,requested_new_forward_solves=16 if prior else 20,
           reused_state_roles=4 if prior else 0,new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False)
    if prior:m['prior_collection']=c['prior_collection']
    for g in parent['groups']:
        new={k:v for k,v in g.items() if k not in ('cache_key','reuse','new_states')};new.update(reuse={},new_states=list(g['states']))
        if prior:
            new['reuse']=reuse_inventory(new,prior,m['prior_collection']);new['new_states']=[s for s in new['states'] if s not in new['reuse']]
        new['cache_key']=group_key(new,m);m['groups'].append(new)
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')


def preflight(path):
    m=read_json(path);prior=None
    if m['protocol']==PROTOCOL:parent=parent_preflight(verify(m['parent']));grids=parent['grids']
    elif m['protocol']==REFINED_PROTOCOL:
        parent=preflight(verify(m['parent']));prior=read_json(verify(m['prior_collection']));grids=REFINED_GRIDS
        if parent['protocol']!=PROTOCOL or prior['manifest']!=m['parent'] or not prior['complete']:raise InvalidArtifact('first conductor collection changed')
    else:raise InvalidArtifact('unknown conductor protocol')
    eps=parent['model']['solvent_epsilon']
    if m['cache_key']!=cache_key({k:v for k,v in m.items() if k!='cache_key'}) or m['model']!={**parent['model'],'model':'cosmo'} or m['energy_prefactor']!=(eps-1)/eps or m['grids']!=grids or m['tolerances']!=parent['tolerances'] or m['python']!=parent['python'] or m['module']!=parent['module'] or len(m['groups'])!=8:
        raise InvalidArtifact('conductor model, source software or frozen numerical settings differ')
    for pin in (m['config'],m['plan'],m['python'],m['module'],*m['implementation'].values()):verify(pin)
    for g,p in zip(m['groups'],parent['groups']):
        reuse=reuse_inventory(g,prior,m['prior_collection']) if prior else {}
        if any(g[k]!=v for k,v in p.items() if k not in ('cache_key','reuse','new_states')) or g['reuse']!=reuse or g['new_states']!=[s for s in g['states'] if s not in reuse] or g['cache_key']!=group_key(g,m):raise InvalidArtifact('source identity, execution inventory or cache changed')
        for item in reuse.values():verify(item['arrays'])
    if m['requested_new_forward_solves']!=(16 if prior else 20) or m['reused_state_roles']!=(4 if prior else 0) or m['new_DFT_calls'] or m['new_MACE_calls']:raise InvalidArtifact('execution count differs')
    return m


def collect(path,output):
    m=preflight(path);root=Path(path).resolve().parent;groups={};receipts={}
    for g in m['groups']:
        p=root/'groups'/g['group_id']/'attempt_0001/result.json';r=read_json(p)
        if r['manifest']!=record(path) or r['cache_key']!=g['cache_key'] or set(r['rows'])!=set(g['states']):raise InvalidArtifact('group receipt mismatch')
        verify(r['arrays']);groups[g['group_id']]=r;receipts[g['group_id']]=record(p)
    r=dict(protocol=m['protocol'],manifest=record(path),groups=groups,receipts=receipts,complete=all(g['status']=='complete' for g in groups.values()),
           requested_new_forward_solves=m['requested_new_forward_solves'],reused_state_roles=m['reused_state_roles'],actual_new_forward_solve_starts=sum(g['actual_new_forward_solve_starts'] for g in groups.values()),
           energy_prefactor=m['energy_prefactor'],new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False)
    write_new(output,r);return r


def execute(path,output):
    m=preflight(path);root=Path(path).resolve().parent;start=time.monotonic()
    for g in m['groups']:execute_group(m,path,g['group_id'],root/('group_'+g['group_id']+'.json'))
    r=collect(path,output)
    write_new(root/'execution_receipt.json',dict(manifest=record(path),collection=record(output),wall_seconds=time.monotonic()-start,slurm_job_id=os.environ.get('SLURM_JOB_ID')))
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute','collect'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op in ('execute','collect'):q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output)
    elif a.command=='dry-run':r=preflight(a.manifest)
    elif a.command=='execute':r=execute(a.manifest,a.output)
    else:r=collect(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('groups','implementation')}))
