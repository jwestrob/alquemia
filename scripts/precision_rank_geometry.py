"""Complete only the declared four missing native geometry/rank cells."""
from __future__ import annotations
import argparse,json,re,shutil
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,HA_TO_KCAL
from affordable_workflow import dry_run
from compact_solvation import completed,diagnostics
from union_adaptive import snapshot

PROTOCOL='native_GFN2_two_source_geometry_rank_2x2_v1'
GATE=.1


def prepare(inventory,agreement,output):
    inv=read_json(inventory)
    if inv['new_molecular_calls']!=0 or inv['maximum_new_scalar_calls_if_authorized']!=4 or len(inv['targets'])!=2 or len(inv['missing_cells'])!=4:
        raise InvalidArtifact('four-cell read-only source inventory required')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=snapshot(root/'implementation');manifests=[]
    backend=inv['targets'][0]['versions']['old']['orca']
    for ranks in (1,8):
        dest=root/('rank_'+str(ranks));tasks=[]
        for cell in inv['missing_cells']:
            if cell['requested_ranks']!=ranks:continue
            source=cell['source_task'];sm=read_json(verify(cell['source_manifest']))
            if sm['orca']!=backend or (source['charge'],source['multiplicity'],source['medium'],source['metal'])!=(-1,1,'vacuum','La'):
                raise InvalidArtifact('target chemistry/backend differs')
            text=verify(source['input']).read_text()
            if not all(s in text for s in ('Native-GFN2-xTB NoAutostart','MaxIter 500','SmearTemp 300','UseXTBMixer true')):
                raise InvalidArtifact('source numerical recipe differs')
            tid=cell['case_id']+'__'+cell['geometry']+'_geometry__rank_'+str(ranks);td=dest/'tasks'/tid;td.mkdir(parents=True)
            for key,name in [('input','endpoint.inp'),('xyz','core.xyz')]:shutil.copyfile(verify(source[key]),td/name)
            tasks.append({'task_id':tid,'case_id':cell['case_id'],'candidate':cell['candidate'],'geometry_version':cell['geometry'],
                'metal':'La','medium':'vacuum','charge':-1,'multiplicity':1,'input':record(td/'endpoint.inp'),'xyz':record(td/'core.xyz'),
                'output_path':str(td/'endpoint.out'),'actual_source_task':source,'actual_source_manifest':cell['source_manifest']})
        low={'protocol_id':PROTOCOL,'inventory':record(inventory),'agreement':record(agreement),'orca':backend,'tasks':tasks,'all_tasks':tasks,
             'rank_count':ranks,'component_gate_kcal_mol':GATE,'execution_resources':{'mpi_ranks':ranks,'concurrent_tasks':2},
             'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
             'implementation':impl,'allocated_cpus_shared':18,'allocated_memory_MiB_shared':65536,'no_seed_files':True}
        mp=dest/'manifest.json';write_new(mp,low);v=validate(mp);write_new(dest/'PREFLIGHT.json',v);manifests.append(record(mp))
    top={'protocol_id':PROTOCOL,'inventory':record(inventory),'agreement':record(agreement),'manifests':manifests,'new_calls':4,
         'reused_cells':4,'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0,'implementation':impl,'production_changed':False}
    write_new(root/'READY.json',top);return {k:v for k,v in top.items() if k!='implementation'}


def validate(manifest):
    m=read_json(manifest);inv=read_json(verify(m['inventory']));r=m['rank_count']
    if m['protocol_id']!=PROTOCOL or r not in (1,8) or len(m['tasks'])!=2 or m['execution_resources']!={'mpi_ranks':r,'concurrent_tasks':2} or m['component_gate_kcal_mol']!=GATE:
        raise InvalidArtifact('fixed parallelism/population differs')
    verify(m['agreement'])
    for t in m['tasks']:
        expected=next(c for c in inv['missing_cells'] if (c['case_id'],c['geometry'],c['requested_ranks'])==(t['case_id'],t['geometry_version'],r))
        if t['actual_source_task']!=expected['source_task'] or t['actual_source_manifest']!=expected['source_manifest']:raise InvalidArtifact('source task changed')
        for k in ('input','xyz'):
            if verify(t[k]).read_bytes()!=verify(expected['source_task'][k]).read_bytes():raise InvalidArtifact('exact source input changed')
        if (t['charge'],t['multiplicity'],t['metal'],t['medium'])!=(-1,1,'La','vacuum'):raise InvalidArtifact('state changed')
    return {'status':'validated','manifest':record(manifest),'dry_run':dry_run(manifest),'new_calls':2,'new_physical_geometries':0}


def collect(ready,output):
    top=read_json(ready);inv=read_json(verify(top['inventory']));fresh=[]
    for pin in top['manifests']:
        mp=verify(pin);validate(mp);m=read_json(mp)
        for t in m['tasks']:
            row={'case_id':t['case_id'],'geometry':t['geometry_version'],'ranks':m['rank_count'],'status':'unavailable',
                 'reason':None,'manifest':pin,'task_id':t['task_id'],'energy_hartree':None,'receipt':None}
            try:
                actual=completed(mp,t['task_id'])
                if actual is None:raise InvalidArtifact('scientific output incomplete or nonconverged')
                audit=diagnostics(actual,t);receipt=read_json(verify(actual['receipt']));text=verify(actual['output']).read_text()
                original=next(x for x in inv['targets'] if x['case_id']==t['case_id'])['versions'][t['geometry_version']]
                if audit['parameter_export']['sha256']!=original['parameter_export']['sha256'] or receipt['parallelism']['nprocs']!=m['rank_count']:
                    raise InvalidArtifact('actual parameters/ranks differ')
                cycles=re.findall(r'SCF CONVERGED AFTER\s+(\d+)\s+CYCLES',text)
                if len(cycles)!=1:raise InvalidArtifact('ambiguous SCF termination')
                row.update(status='complete',energy_hartree=actual['energy_hartree'],receipt=actual['receipt'],output=actual['output'],
                    audit=audit,SCF_cycles=int(cycles[0]),runtime_line=next((l for l in text.splitlines() if l.startswith('TOTAL RUN TIME:')),None))
            except (InvalidArtifact,KeyError,OSError) as exc:
                row['reason']=str(exc);rp=Path(t['output_path']+'.execution.json');row['receipt']=record(rp) if rp.exists() else None
            fresh.append(row)
    targets=[]
    for target in inv['targets']:
        matrix={g:{} for g in ('old','new')}
        for g,r in target['versions'].items():matrix[g][str(r['actual_ranks'])]={'status':'complete','energy_hartree':r['energy_hartree'],'receipt':r['receipt'],'output':r['output'],'SCF_cycles':r['SCF_cycles'],'reused':True}
        for row in fresh:
            if row['case_id']==target['case_id']:matrix[row['geometry']][str(row['ranks'])]={**row,'reused':False}
        ok=all(matrix[g].get(str(r),{}).get('status')=='complete' for g in ('old','new') for r in (1,8));effects=None
        if ok:
            E=lambda g,r:matrix[g][str(r)]['energy_hartree']
            rank={g:(E(g,1)-E(g,8))*HA_TO_KCAL for g in ('old','new')}
            geom={str(r):(E('new',r)-E('old',r))*HA_TO_KCAL for r in (1,8)}
            effects={'rank1_minus_rank8_kcal_mol':rank,'new_minus_old_geometry_kcal_mol':geom,
                'rank_gate_pass':{g:abs(x)<=GATE for g,x in rank.items()},'interaction_kcal_mol':rank['new']-rank['old']}
        targets.append({'case_id':target['case_id'],'candidate':target['candidate'],'metal':'La','medium':'vacuum',
                        'status':'complete' if ok else 'unavailable','matrix':matrix,'effects':effects})
    result={'protocol_id':PROTOCOL,'ready':record(ready),'inventory':top['inventory'],'fresh_cells':fresh,'fresh_denominator':4,
        'fresh_complete':sum(r['status']=='complete' for r in fresh),'targets':targets,'component_gate_kcal_mol':GATE,
        'production_changed':False,'new_pool_or_calibration':False,'new_molecular_calls_in_collection':0,'implementation':record(__file__)}
    write_new(output,result);return {'collection':record(output),'fresh_complete':result['fresh_complete'],'effects':[t['effects'] for t in targets]}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('inventory','agreement','output'),'validate':('manifest',),'collect':('ready','output')}.items():
        q=s.add_parser(op)
        for field in fields:q.add_argument('--'+field.replace('_','-'),required=True,type=Path)
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
