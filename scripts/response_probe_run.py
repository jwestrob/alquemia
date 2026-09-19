#!/usr/bin/env python3
"""Narrow manifest adapter around the preserved native MACE runner and worker."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz

SCHEMA='alquemia.response_probe_core.v1'


def validate(manifest):
    from response_probe_backend import check_atoms
    m=read_json(manifest);p=read_json(verify(m['preparation']))
    parent=read_json(verify(m['parent_manifest']))
    effective=lambda model:{k:v for k,v in model.items() if k!='preparation_policy'}
    if m['schema_version']!=SCHEMA or effective(m['model'])!=effective(parent['model']) or m['software']!=parent['software']:
        raise InvalidArtifact('changed native model/software')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    wanted={(c,metal) for c in p['cases'] for metal in ('Ca','La')}
    if {(t['case_id'],t['metal']) for t in m['tasks']}!=wanted or len(m['tasks'])!=len(wanted):
        raise InvalidArtifact('changed finite endpoint inventory')
    for t in m['tasks']:
        state=p['cases'][t['case_id']]['states'][t['metal']]
        if any(t[k]!=state[k] for k in ('xyz','charge','spin_multiplicity')):
            raise InvalidArtifact('task differs from exact preparation')
        check_atoms(xyz(verify(t['xyz'])),t['charge'])
        base={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':base,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('scientific cache identity changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def collect(manifest):
    from response_probe_backend import accepted_attempt
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
            if r is not None:valid.append(r)
        rows[t['task_id']]=valid[-1] if valid else {'status':'unavailable','energy_eV':None}
    return {'manifest':record(mp),'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete',
            'rows':rows,'attempts':attempts,'baseline_changed':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    for name in ('dry-run','execute','collect','worker'):
        a=s.add_parser(name);a.add_argument('--manifest',required=True)
        if name in ('execute','worker'):a.add_argument('--memory-mode',default='native')
        if name=='execute':a.add_argument('--task-id',action='append')
        if name in ('worker','collect'):a.add_argument('--output',required=True)
        if name=='worker':a.add_argument('--task-id',required=True)
    args=p.parse_args()
    if args.operation=='dry-run':result=validate(args.manifest)
    elif args.operation=='collect':result=collect(args.manifest);write_new(args.output,result)
    else:
        import response_probe_backend as backend
        # Existing execution/caching/locking/native worker are unchanged. Only
        # this explicitly named six-endpoint manifest validator is substituted.
        backend.dry_run=validate
        if args.operation=='execute':result=backend.execute(args.manifest,args.memory_mode,args.task_id)
        else:result=backend.worker(args.manifest,args.task_id,args.output,args.memory_mode)
    print(json.dumps(result,indent=2))
    if result.get('status') in ('failed','partial_failure'):
        raise SystemExit(1)
