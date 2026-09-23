"""Thin immutable native-search/strict-origin join and existing common-pool execution."""
import argparse
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from mace_site_kinematics import Kinematics
import adaptive_angular_proposals as angular
import consistent_context as context
import motion_envelope_run as proposal
import motion_envelope_scalar as scalar
import slsqp_precision as precision
from union_adaptive import snapshot,POOL_SETTINGS
pool=precision.private('nikasha_pool')
POOL_PROTOCOL='Nikasha_motion_envelope4p3_adaptive_minimal_native_OMOL_strictGFN2_v1'


def joined_origins(collection,source):
    d=read_json(collection);mp=verify(d['manifest']);m=read_json(mp)
    if d['protocol_id']!=scalar.PROFILE or m['source']!=source['origin_collection'] or m['envelope_stage']!='origins':
        raise InvalidArtifact('strict scalar origins do not belong to these native sources')
    scalar.validate(mp);rows={r['task_id']:r for r in d['rows']}
    if set(rows)!={t['task_id']for t in m['tasks']}:raise InvalidArtifact('strict origin collection population changed')
    lows={}
    for t in m['tasks']:
        row=scalar.endpoint(mp,t);saved=rows[t['task_id']]
        if {k:v for k,v in saved.items()if k not in ('task_id','case_id','candidate','metal','medium')}!={k:v for k,v in row.items()if k!='task_id'}:
            raise InvalidArtifact('strict origin collection differs from actual output')
        source_task=next(x for x in source['tasks']if(x['case_id'],x['metal'])==(t['case_id'],t['metal']))
        if not context.reusable_state(t,source_task):raise InvalidArtifact('strict origin physical state differs')
        lows.setdefault((t['case_id'],t['metal']),{})[t['medium']]=row
    joined={}
    for t in source['tasks']:
        low=lows.get((t['case_id'],t['metal']),{})
        complete=len(low)==2 and all(v['status']=='complete'for v in low.values())
        cell={'status':'complete'if complete else'unavailable','components':None,'xyz':t['xyz'],
            'MACE':t['q0']['native_MACE_receipt'],'low':low,'reused':True,'strict_collection':record(collection)}
        if complete:
            params=[read_json(verify(low[k]['component_audit']['parameter_export']))for k in ('vacuum','alpb')]
            if params[0]!=params[1]:raise InvalidArtifact('origin medium parameters differ')
            cell['components']={**t['q0']['components'],'GFN2_vacuum_hartree':low['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':low['alpb']['energy_hartree']}
        joined.setdefault(t['case_id'],{})[t['metal']]=cell
    return joined


def prepare_pool(proposals,scalar_origins,agreement,output):
    """Only assemble cells; unchanged shared executors compute/collect them."""
    result=read_json(proposals);mp=verify(result['manifest']);proposal.validate(mp);m=read_json(mp)
    if result['protocol_id']!=proposal.PROTOCOL or len(result['endpoints'])!=len(m['tasks']):raise InvalidArtifact('proposal population differs')
    joins=joined_origins(scalar_origins,m)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);cases=[];tasks=[]
    for source in m['cases']:
        cid=source['case_id'];ts={z:next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(cid,z)) for z in ('Ca','La')}
        c={'case_id':cid,'status':'unavailable','reason':None,'source':source,'aliases':{},
           'candidates':[{'id':'origin','xyz':ts['Ca']['xyz']}],'matrix':{z:{'origin':joins[cid][z]} for z in ('Ca','La')}}
        pending=[]
        try:
            if any(joins[cid][z]['status']!='complete' for z in ('Ca','La')):raise InvalidArtifact('actual strict scalar origin unavailable')
            points={};atoms_by_name={'origin':xyz(verify(ts['Ca']['xyz']))}
            for z in ('Ca','La'):
                r=next(r for r in result['endpoints'] if (r['case_id'],r['metal'])==(cid,z))
                if r['status']!='candidate_available':raise InvalidArtifact('required adaptive '+z+' unavailable: '+str(r.get('reason')))
                actual=read_json(verify(r['proposal_receipt']));point=r['candidate'];q=np.asarray(point['full_q']);k=Kinematics(read_json(verify(ts[z]['mapping']))['context'])
                if actual['manifest']!=record(mp) or actual['status']!='proposal_available' or actual['proposal']!=point:raise InvalidArtifact('candidate receipt differs')
                angular.final_geometry(k,ts[z],q[ts[z]['active_indices']],[a[0] for a in xyz(verify(ts[z]['xyz']))])
                atoms=xyz(verify(point['coordinate']))
                if not np.allclose(k.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('candidate coordinates differ')
                name='adaptive_'+z;match=next((n for n,a in atoms_by_name.items() if pool.same_geometry(a,atoms)),None)
                if match is None:
                    match=name;atoms_by_name[name]=atoms;c['candidates'].append({'id':name,'xyz':point['coordinate']})
                c['aliases'][name]={'representative':match,'source_xyz':point['coordinate'],'proposal_receipt':r['proposal_receipt']};points[z]=point
            for z in ('Ca','La'):
                t=ts[z]
                for candidate in c['candidates'][1:]:
                    name=candidate['id'];atoms=atoms_by_name[name];cross=[(z,*atoms[0][1:]),*atoms[1:]]
                    tid=cid+'__'+z+'__at_'+name;td=out/'cells'/tid;td.mkdir(parents=True);xp=td/'context.xyz'
                    pool.write_xyz(xp,cross)
                    if xyz(xp)!=cross:raise InvalidArtifact('cross coordinate serialization changed')
                    task={'task_id':tid,'case_id':cid,'metal':z,'candidate':name,'xyz':record(xp),
                        'charge':t['charge'],'multiplicity':t['multiplicity'],'source_mapping':t['mapping'],'source_preparation':t['source_preparation']}
                    if name==c['aliases']['adaptive_'+z]['representative'] and xyz(verify(points[z]['coordinate']))==cross:
                        task['native_reuse']=points[z]['MACE'];pool.native_reuse(task,m)
                    pending.append(task);c['matrix'][z][name]={'status':'pending','task_id':tid,'xyz':record(xp),'reused':False}
            c['status']='prepared';tasks.extend(pending)
        except (InvalidArtifact,KeyError,OSError,StopIteration) as exc:c['reason']=str(exc)
        cases.append(c)
    pm={k:m[k] for k in ('model','software','orca','cpu_python','gpu_python','resources')}
    pm.update(protocol_id=POOL_PROTOCOL,settings=POOL_SETTINGS,scalar_origins=record(scalar_origins),source=record(proposals),source_manifest=record(mp),agreement=record(agreement),
        declared_case_ids=m['declared_case_ids'],cases=cases,tasks=tasks,implementation=snapshot(out/'implementation'),
        population='canonical25_crystal3_probe6_envelope_adaptive_minimal',shard_count=1,new_MACE_cells=sum(not t.get('native_reuse') for t in tasks),
        maximum_new_GFN2_calls=2*len(tasks),new_DFT_calls=0,new_optimizations=0,GFN2_maxiter=500,
        numerical_policy_id=scalar.PROFILE,reference=m.get('reference'),production_changed=False)
    path=out/'manifest.json';write_new(path,pm);check=validate_pool(path);low_prepare(path);write_new(out/'PREFLIGHT.json',check);return check


def validate_pool(manifest):
    m=read_json(manifest);source=read_json(verify(m['source']));sm=read_json(verify(m['source_manifest']));proposal.validate(verify(m['source_manifest']))
    joins=joined_origins(verify(m['scalar_origins']),sm)
    count=len(sm['cases'])
    if (m['protocol_id']!=POOL_PROTOCOL or m['settings']!=POOL_SETTINGS or m['GFN2_maxiter']!=500 or m['numerical_policy_id']!=scalar.PROFILE or count!=34 or
        source['manifest']!=m['source_manifest'] or m['declared_case_ids']!=sm['declared_case_ids'] or m.get('reference')!=sm.get('reference') or len(m['cases'])!=count or
        m['maximum_new_GFN2_calls']!=2*len(m['tasks']) or len(m['tasks'])>4*count or m['new_MACE_cells']>2*count):
        raise InvalidArtifact('fixed minimal union pool differs')
    for key in ('model','software','orca','cpu_python','gpu_python','resources'):
        if m[key]!=sm[key]:raise InvalidArtifact('union pool scientific runtime differs')
    verify(m['agreement'])
    for pin in m['implementation'].values():verify(pin)
    if len({t['task_id'] for t in m['tasks']})!=len(m['tasks']):raise InvalidArtifact('duplicate cell')
    required={cell['task_id'] for c in m['cases'] if c['status']=='prepared' for cells in c['matrix'].values() for cell in cells.values() if not cell['reused']}
    if required!={t['task_id']for t in m['tasks']} or m['new_MACE_cells']!=sum(not t.get('native_reuse')for t in m['tasks']):raise InvalidArtifact('candidate task count or membership differs')
    for c,original in zip(m['cases'],sm['cases']):
        if c['case_id']!=original['case_id'] or c['source']!=original:raise InvalidArtifact('union pool source changed')
        cid=c['case_id'];names=[q['id'] for q in c['candidates']]
        if names[0]!='origin' or len(set(names))!=len(names) or not set(names)<=set(POOL_SETTINGS['candidate_order']):raise InvalidArtifact('candidate membership differs')
        for z in ('Ca','La'):
            t=next(t for t in sm['tasks'] if (t['case_id'],t['metal'])==(cid,z))
            origin=c['matrix'][z]['origin']
            if origin!=joins[cid][z]:raise InvalidArtifact('original source cell differs')
        if c['status']!='prepared':
            if any(t['case_id']==cid for t in m['tasks']):raise InvalidArtifact('unavailable pool gained cells')
            continue
        if set(c['aliases'])!={'adaptive_Ca','adaptive_La'}:raise InvalidArtifact('required endpoint proposal omitted')
        for z in ('Ca','La'):
            e=next(e for e in source['endpoints'] if (e['case_id'],e['metal'])==(cid,z))
            if e['status']!='candidate_available' or c['aliases']['adaptive_'+z]['proposal_receipt']!=e['proposal_receipt']:raise InvalidArtifact('successful proposal identity differs')
            receipt=read_json(verify(e['proposal_receipt']));point=e['candidate'];t=next(t for t in sm['tasks'] if (t['case_id'],t['metal'])==(cid,z))
            if receipt['proposal']!=point or receipt['manifest']!=m['source_manifest']:raise InvalidArtifact('physical proposal changed')
            q=np.asarray(point['full_q']);kin=Kinematics(read_json(verify(t['mapping']))['context'])
            if any(q[i]!=0 for i in set(range(len(q)))-set(t['active_indices'])):raise InvalidArtifact('unselected mode moved')
            angular.final_geometry(kin,t,q[t['active_indices']],[a[0] for a in xyz(verify(t['xyz']))])
            if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in xyz(verify(point['coordinate']))],atol=1e-12,rtol=0):raise InvalidArtifact('candidate map drift')
            alias=c['aliases']['adaptive_'+z]
            representative=next(v for v in c['candidates'] if v['id']==alias['representative'])
            if alias['source_xyz']!=point['coordinate'] or not pool.same_geometry(xyz(verify(representative['xyz'])),xyz(verify(point['coordinate']))):raise InvalidArtifact('candidate alias not numerical copy')
            if set(c['matrix'][z])!=set(names):raise InvalidArtifact('metals do not share same candidates')
            for name in names[1:]:
                cell=c['matrix'][z][name];task=next(t for t in m['tasks'] if t['task_id']==cell['task_id'])
                coords=xyz(verify(next(v['xyz'] for v in c['candidates'] if v['id']==name)))
                want=[(z,*coords[0][1:]),*coords[1:]]
                if (task['case_id'],task['metal'],task['candidate'],task['source_mapping'],task['source_preparation'])!=(cid,z,name,t['mapping'],t['source_preparation']):raise InvalidArtifact('cross source identity differs')
                if xyz(verify(task['xyz']))!=want or (task['charge'],task['multiplicity'])!=(t['charge'],t['multiplicity']):raise InvalidArtifact('cross geometry/state differs')
                if task.get('native_reuse'):pool.native_reuse(task,m)
    return {'status':'validated','manifest':record(manifest),'denominator':count,'prepared':sum(c['status']=='prepared' for c in m['cases']),
            'new_MACE_cells':m['new_MACE_cells'],'new_GFN2_calls':m['maximum_new_GFN2_calls'],'new_DFT_calls':0}


def low_prepare(manifest):
    m=read_json(manifest);origin=read_json(verify(m['scalar_origins']));om=read_json(verify(origin['manifest']));root=Path(manifest).parent/'solvent'
    entries=[{**t,'cell_id':t['task_id']}for t in m['tasks']]
    if entries:
        scalar.prepare(entries,m,record(manifest),verify(om['qualification']),verify(m['agreement']),root/'shard_0','candidates',272)
        shards=[record(root/'shard_0/manifest.json')]
    else:root.mkdir();shards=[]
    write_new(root/'INDEX.json',{'pool_manifest':record(manifest),'shards':shards})


def strict_completed(manifest,task_id):
    m=read_json(manifest);task=next(t for t in m['tasks']if t['task_id']==task_id)
    result=scalar.endpoint(manifest,task)
    if result['status']!='complete':
        if result.get('actual'):raise InvalidArtifact('strict scalar qualification failed: '+result['reason'])
        return None
    return result


def execute_mace(manifest):
    pool.validate=validate_pool
    return pool.execute_mace(manifest)


def collect(manifest,output):
    validate_pool(manifest)
    index=read_json(Path(manifest).parent/'solvent/INDEX.json')
    for pin in index['shards']:scalar.validate(verify(pin))
    pool.completed=strict_completed
    return pool.collect(manifest,output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    for op,keys in {'prepare_pool':('proposals','scalar_origins','agreement','output'),'validate_pool':('manifest',),'execute_mace':('manifest',),'collect':('manifest','output')}.items():
        q=s.add_parser(op)
        for k in keys:q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());cmd=a.pop('command');print(json.dumps(globals()[cmd](**a),indent=2))
