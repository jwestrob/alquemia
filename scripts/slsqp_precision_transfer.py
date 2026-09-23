"""Four disjoint precision-policy transfer shards; same scientific engine."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import numpy as np
import scipy
import slsqp_precision as pilot
import union_adaptive as original
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz, paired
from adaptive_force_diagnostic import project, preview
from mace_site_kinematics import Kinematics
from compact_solvation import input_text

REFERENCE_HASH = '1f8470bbdf7a056098ca261aa10b72bdc940cf71ab860fbbc54090bda1e73250'
POPULATION = 'uniform_precision_primary225'
LOW_RESOURCES = {'mpi_ranks': 1, 'concurrent_tasks': 32}


def qualification(path):
    r = read_json(path)
    if (r['protocol_id'] != 'Nikasha_union_reference28_native_GFN2_rank1_scalar_qualification_v1'
        or (r['complete_cells'], r['cell_denominator'], r['complete'], r['case_denominator']) != (336,336,28,28)
        or r['all_numerical_gates_pass'] is not True or not all(c['all_numerical_gates_pass'] for c in r['cases'])):
        raise InvalidArtifact('complete scalar rank28 qualification required')
    return record(path)


def population(inventory, rank_qualification):
    inv=read_json(inventory); old=read_json(verify(inv['original_selection'])); ref=read_json(verify(inv['new_reference']))
    qualification(rank_qualification)
    if inv['new_reference']['sha256']!=REFERENCE_HASH or ref['optimizer_settings']!=pilot.SETTINGS or ref['settings']!=original.POOL_SETTINGS:
        raise InvalidArtifact('frozen precision reference differs')
    if ref['calibration_denominator']!=25 or ref['crystals_or_noncanonical_used_for_fit']:
        raise InvalidArtifact('reference population differs')
    rows=inv['rows']; ready=[r for r in rows if r['union_origin_status']=='complete']
    reuse=[r for r in ready if r['precision34_reuse']]; fresh=[r for r in ready if not r['precision34_reuse']]
    if (len(rows),len(ready),len(reuse),len(fresh))!=(225,208,6,202):raise InvalidArtifact('225/208/6/202 denominator differs')
    if len({r['case_id'] for r in rows})!=225 or [{k:r[k] for k in o} for r,o in zip(rows,old['all_cases'])]!=old['all_cases']:
        raise InvalidArtifact('original source identity/order changed')
    comparison=read_json(verify(inv['precision34_comparison'])); by={r['case_id']:r for r in comparison['rows']}
    for row in reuse:
        c=by[row['case_id']]
        if row['precision34_reuse']!=c['collection'] or row['reused_score']!=c['scores'] or c['status']!='available':
            raise InvalidArtifact('actual six precision reuses differ')
    return inv, rows, fresh


def prepare(inventory, rank_qualification, sources, agreement, output):
    inv, rows, fresh=population(inventory,rank_qualification)
    if len(sources)!=4:raise InvalidArtifact('four original union transfer manifests required')
    parents=[read_json(p) for p in sources]; lookup={}
    for path,p in zip(sources,parents):
        if p['protocol_id']!=original.PROTOCOL or p['settings']!=pilot.original_completion.SETTINGS:raise InvalidArtifact('old optimizer protocol differs')
        for c in p['cases']:
            if c['case_id'] in lookup:raise InvalidArtifact('duplicate original source')
            lookup[c['case_id']]=(path,p,c)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    top={'population':POPULATION,'inventory':record(inventory),'rank_qualification':qualification(rank_qualification),
         'reference':inv['new_reference'],'agreement':record(agreement),'all_cases':rows,
         'reuse_case_ids':[r['case_id'] for r in rows if r['precision34_reuse']],
         'total_denominator':225,'prepared_sources':208,'preparation_exclusions':17,'new_sources':202,'new_searches':404,
         'maximum_cross_MACE_calls':404,'maximum_GFN2_calls':1616,'new_q0_calls':0,'new_DFT_calls':0,
         'shard_count':4,'selection_UTC':datetime.now(timezone.utc).isoformat(),'production_changed':False}
    tp=out/'SELECTION.json';write_new(tp,top);pins=[]
    for shard in range(4):
        chosen=fresh[shard::4];cases=[];tasks=[];records=[]
        for row in chosen:
            cid=row['case_id'];path,parent,c=lookup[cid];cases.append(c)
            tasks.extend(next(t for t in parent['tasks'] if (t['case_id'],t['metal'])==(cid,z)) for z in ('Ca','La'))
            old_pool=Path(path).parent.parent/'pool/collection_final.json'
            records.append({**row,'parent':record(path),'old_pool':record(old_pool),
                'old_receipts':{z:record(Path(path).parent/'proposals'/(cid+'__'+z)/'result.json') for z in ('Ca','La')}})
        p=parents[0];dest=out/('shard_'+str(shard))/'proposals';dest.mkdir(parents=True)
        m={k:p[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources','optimizer_software')}
        m.update(protocol_id=pilot.PROTOCOL,settings=pilot.SETTINGS,gates=pilot.GATES,agreement=record(agreement),
            reference=top['reference'],selection=record(tp),rank_qualification=top['rank_qualification'],shard=shard,
            cases=cases,tasks=tasks,sources=records,declared_case_ids=[r['case_id'] for r in chosen],
            population=POPULATION+'_shard_'+str(shard),maximum_optimizer_starts=2*len(cases),
            maximum_cross_MACE_calls=2*len(cases),maximum_GFN2_calls=8*len(cases),new_origin_calls=0,new_DFT_calls=0,
            production_changed=False,implementation=original.snapshot(dest/'implementation'))
        mp=dest/'manifest.json';write_new(mp,m);v=validate(mp);write_new(dest/'PREFLIGHT.json',v);pins.append(record(mp))
        print(json.dumps(v),flush=True)
    r={'selection':record(tp),'shards':pins,'sources_per_shard':[len(fresh[i::4]) for i in range(4)],'new_searches':404,'GFN2_maximum':1616,'new_molecular_calls':0}
    write_new(out/'READY.json',r);return r


def validate(manifest):
    m=read_json(manifest);top=read_json(verify(m['selection']));inv,rows,fresh=population(verify(top['inventory']),verify(top['rank_qualification']))
    ids=[r['case_id'] for r in fresh[m['shard']::4]];n=len(ids)
    if (top['population']!=POPULATION or m['shard'] not in range(4) or n not in (50,51) or
        m['protocol_id']!=pilot.PROTOCOL or m['settings']!=pilot.SETTINGS or m['gates']!=pilot.GATES or
        m['declared_case_ids']!=ids or len(m['cases'])!=n or len(m['tasks'])!=2*n or
        m['maximum_optimizer_starts']!=2*n or m['maximum_cross_MACE_calls']!=2*n or m['maximum_GFN2_calls']!=8*n or
        m['reference']!=inv['new_reference'] or top['all_cases']!=rows or m['rank_qualification']!=top['rank_qualification']):
        raise InvalidArtifact('frozen transfer shard differs')
    for key in ('agreement','reference','software','orca','cpu_executable','gpu_executable'):verify(m[key])
    for pin in m['implementation'].values():verify(pin)
    for key in ('wrapper','kernel'):verify(m['optimizer_software'][key])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer software changed')
    if {k for k in pilot.SETTINGS if pilot.SETTINGS[k]!=pilot.original_completion.SETTINGS[k]}!={'optimizer_ftol'}:raise InvalidArtifact('additional optimizer change')
    parents={};expected=[]
    for cid,source,c in zip(ids,m['sources'],m['cases']):
        if source['case_id']!=cid or c['case_id']!=cid:raise InvalidArtifact('source order differs')
        path=verify(source['parent'])
        if str(path) not in parents:parents[str(path)]=read_json(path)
        parent=parents[str(path)]
        if parent['settings']!=pilot.original_completion.SETTINGS or next(x for x in parent['cases'] if x['case_id']==cid)!=c:
            raise InvalidArtifact('source context changed')
        for key in ('model','software','orca','cpu_python','gpu_python','resources'):
            if m[key]!=parent[key]:raise InvalidArtifact('electronic/runtime method differs')
        pair=[];data=[]
        for z in ('Ca','La'):
            t=next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(cid,z));pair.append(t);expected.append(t['task_id'])
            if t!=next(q for q in parent['tasks'] if q['task_id']==t['task_id']):raise InvalidArtifact('original task changed')
            native,forces=original.origin(c['origin_row']['native_endpoints'][z],m['model'])
            kin=Kinematics(read_json(verify(t['mapping']))['context']);coords,v,raw,normed,_=project(kin.data,np.zeros(len(kin.modes)),forces)
            data.append((kin,v,normed))
            if not np.array_equal(raw[t['active_indices']],t['origin_reuse']['point']['gradient_kcal_mol_rad']) or native['energy_eV']!=t['q0']['components']['MACE_eV']:
                raise InvalidArtifact('q0 force/energy changed')
            if not np.allclose(coords,[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('q0 map changed')
            pilot.angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
        paired(verify(pair[1]['xyz']),verify(pair[0]['xyz']),pair[1]['charge'],pair[0]['charge'])
        if read_json(verify(pair[0]['mapping']))!=read_json(verify(pair[1]['mapping'])):raise InvalidArtifact('paired map differs')
        choice=preview(data[0][0].modes,data[0][1],data[0][2],data[1][2])
        if choice!=c['selection'] or any(t['selector']!=choice for t in pair):raise InvalidArtifact('paired selector differs')
    if [t['task_id'] for t in m['tasks']]!=expected:raise InvalidArtifact('task order changed')
    return {'status':'validated','manifest':record(manifest),'shard':m['shard'],'sources':n,'total_denominator':225,
            'new_searches':2*n,'cross_MACE_maximum':2*n,'GFN2_maximum':8*n,'new_q0_calls':0,'new_DFT_calls':0}


def low_prepare(manifest):
    """Same native recipe and runner; qualified scalar parallelism only."""
    from affordable_workflow import dry_run
    m=read_json(manifest);sm=read_json(verify(m['source_manifest']));qualification(verify(sm['rank_qualification']))
    if m['shard_count']!=1:raise InvalidArtifact('one contained solvent shard per source shard required')
    root=Path(manifest).parent/'solvent';sd=root/'shard_0';tasks=[]
    for cell in m['tasks']:
        for medium in ('vacuum','alpb'):
            tid=cell['task_id']+'__'+medium;td=sd/'tasks'/tid;td.mkdir(parents=True)
            xp=td/'core.xyz';shutil.copyfile(verify(cell['xyz']),xp)
            ip=td/'endpoint.inp';ip.write_text(input_text(cell['charge'],cell['multiplicity'],medium,'native').replace('%scf\n','%scf\n MaxIter 500\n'))
            tasks.append({**cell,'task_id':tid,'cell_id':cell['task_id'],'case':cell['case_id'],'medium':medium,
                          'xyz':record(xp),'input':record(ip),'output_path':str(td/'endpoint.out')})
    low={'protocol_id':m['protocol_id'],'pool_manifest':record(manifest),'agreement':m['agreement'],'orca':m['orca'],
         'tasks':tasks,'all_tasks':tasks,'shard':0,'rank_qualification':sm['rank_qualification'],
         'execution_policy':{'task_runner':m['implementation']['run_orca_task_manifest.py'],'runtime_renderer':m['implementation']['render_orca_runtime_input.py']},
         'execution_resources':LOW_RESOURCES,'allocated_cpus':32,'allocated_memory_MiB':65536,'cached_origin_ranks':8}
    lp=sd/'manifest.json';write_new(lp,low);check=dry_run(lp) if tasks else {'status':'no_new_tasks'}
    write_new(sd/'PREFLIGHT.json',check);write_new(root/'INDEX.json',{'pool_manifest':record(manifest),'shards':[record(lp)]})


def engine():
    p=pilot.private('slsqp_precision');p.validate=validate;u,pool=p.engine();pool.low_prepare=low_prepare;return u,pool
def execute(manifest):return engine()[0].execute(manifest)
def collect(manifest,output):return engine()[0].collect(manifest,output)
def prepare_pool(proposals,agreement,output):return engine()[0].prepare_pool(proposals,agreement,output)
def validate_pool(manifest):return engine()[0].validate_pool(manifest)
def execute_mace(manifest):return engine()[1].execute_mace(manifest)
def collect_pool(manifest,output):
    validate_pool(manifest);return engine()[1].collect(manifest,output)


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    specs={'prepare':('inventory','rank_qualification','sources','agreement','output'),'validate':('manifest',),
           'execute':('manifest',),'collect':('manifest','output'),'prepare_pool':('proposals','agreement','output'),
           'validate_pool':('manifest',),'execute_mace':('manifest',),'collect_pool':('manifest','output')}
    for op,fields in specs.items():
        q=s.add_parser(op)
        for field in fields:q.add_argument('--'+field.replace('_','-'),type=Path,required=True,**({'nargs':'+'} if field=='sources' else {}))
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
if __name__=='__main__':main()
