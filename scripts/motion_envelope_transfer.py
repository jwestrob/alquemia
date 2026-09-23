"""Finite full100 envelope adapter around existing native, search and pool runners."""
import argparse
import json
import os
from pathlib import Path
import scipy
import scipy.optimize._slsqp_py as slsqp_wrapper
import scipy.optimize._slsqplib as slsqp_kernel
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,cache_key,energy
from affordable_workflow import dry_run
from compact_solvation_compare import native_endpoint
import consistent_context as context
import union_adaptive as union
import slsqp_precision as precision
import motion_envelope_scalar as scalar
from union_triple_transfer_run import build_tasks,selector_equivalence
from motion_envelope_pool import strict_completed
from motion_envelope_transfer_inventory import PROTOCOL as SCOPE_PROTOCOL

SEARCH_PROTOCOL='Nikasha_motion_envelope4p3_all100_ftol1e8_native_search_v1'
POOL_PROTOCOL='Nikasha_motion_envelope4p3_all100_minimal_strict_native_v1'


def scope(inventory,shard):
    a=read_json(inventory)
    if a['protocol_id']!=SCOPE_PROTOCOL or a['counts']['new_pools']!=101 or shard not in (0,1):raise InvalidArtifact('finite transfer scope differs')
    s=read_json(verify(a['shards'][shard]));ref=read_json(verify(a['reference']))
    if len(s['cases'])!=(51 if shard==0 else 50) or ref['scalar_profile']!=scalar.PROFILE or ref['optimizer_settings']!=precision.SETTINGS:
        raise InvalidArtifact('frozen reference/profile/shard differs')
    verify(a['agreement']);return a,s,ref


def compatible_completed(manifest,task_id):
    """Preserve actual receipt identity; permit only exact runner/renderer byte copies."""
    from run_orca_task_manifest import load_manifest_tasks,_completed_attempt_is_valid
    m,tasks=load_manifest_tasks(Path(manifest));t=next(t for t in tasks if t['task_id']==task_id)
    op=t['output'];rp=Path(str(op)+'.execution.json')
    if not op.exists()or not rp.exists():return None
    r=read_json(rp)
    for key in ('task_runner','runtime_renderer'):
        actual=r[key];declared=m['execution_policy'][key]
        verify(actual);verify(declared)
        if actual['sha256']!=declared['sha256']:raise InvalidArtifact('actual runner/renderer bytes differ')
    if not _completed_attempt_is_valid(rp,op,manifest_sha256=record(manifest)['sha256'],task=t,
        runner_identity=r['task_runner'],runtime_renderer_identity=r['runtime_renderer']):return None
    return {'energy_hartree':energy(op),'output':record(op),'receipt':record(rp),'manifest':record(manifest),'task_id':task_id}


def scalar_engine():
    s=precision.private('motion_envelope_scalar');s.validate=validate_scalar
    solvent=precision.private('compact_solvation');solvent.completed=compatible_completed;s.solvent=solvent
    return s


def strict_completed(manifest,task_id):
    m=read_json(manifest);task=next(t for t in m['tasks']if t['task_id']==task_id)
    r=scalar_engine().endpoint(manifest,task)
    if r['status']!='complete':
        if r.get('actual'):raise InvalidArtifact('actual strict endpoint failed: '+r['reason'])
        return None
    return r


def prepare_origins(inventory,shard,template,qualification,output):
    shard=int(shard);a,s,ref=scope(inventory,shard);base=read_json(template)
    if base['model']!=ref['model']:raise InvalidArtifact('native model changed')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=union.snapshot(out/'implementation');tasks=[];entries=[]
    for r in s['cases']:
        c=r['prepared'];cid=c['case_id']+'__envelope_'+c['selection_id'].rsplit('_',1)[-1]
        for z,ep in c['representations']['context']['endpoints'].items():
            entry={'case_id':cid,'source_case_id':c['case_id'],'selection_id':c['selection_id'],'candidate':'origin',
                'metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity']};entries.append(entry)
            t={**entry,'task_id':cid+'__context__'+z,'representation':'context','metal_index':0,'kind':'core',
                'spin_multiplicity':ep['multiplicity'],'energy_component':'MACE_OMOL_total_vacuum_energy'}
            t['cache_key']=cache_key({'task':t,'model':base['model'],'software':base['software'],'implementation':impl});tasks.append(t)
    inputs={'protocol_id':SCOPE_PROTOCOL,'inventory':record(inventory),'shard':shard,'scope':a['shards'][shard],
        'cases':s['cases'],'entries':entries,'model':base['model'],'software':base['software'],'template':record(template),
        'qualification':scalar.qualified(qualification),'reference':a['reference'],'agreement':a['agreement'],'implementation':impl}
    write_new(out/'INPUTS.json',inputs)
    mm={'protocol_id':'mace_omol_0_100m_vacuum_descriptor_v1','stage':SEARCH_PROTOCOL,'preparation':record(out/'INPUTS.json'),
        'agreement':a['agreement'],'software':base['software'],'model':base['model'],'implementation':impl,
        'tasks':tasks,'reused':{},'verification_references':{}}
    write_new(out/'mace/manifest.json',mm)
    from mace_omol import check_atoms
    for t in tasks:check_atoms(xyz(verify(t['xyz'])),t['charge'])
    scalar_engine().prepare(entries,base,record(out/'INPUTS.json'),qualification,verify(a['agreement']),out/'scalar','transfer_origins',4*len(s['cases']))
    ready={'inputs':record(out/'INPUTS.json'),'MACE_manifest':record(out/'mace/manifest.json'),'scalar_manifest':record(out/'scalar/manifest.json'),
        'origin_MACE':len(tasks),'origin_GFN2':2*len(entries),'new_molecular_calls':0}
    write_new(out/'READY.json',ready);return ready


def validate_scalar(manifest):
    m=read_json(manifest);scalar.qualified(verify(m['qualification']));source=read_json(verify(m['source']));verify(m['agreement'])
    if m['envelope_stage']=='transfer_origins':
        a,s,_=scope(verify(source['inventory']),source['shard']);entries=source['entries'];maximum=4*len(s['cases'])
        expected=[]
        for r in s['cases']:
            c=r['prepared'];cid=c['case_id']+'__envelope_'+c['selection_id'].rsplit('_',1)[-1]
            for z,ep in c['representations']['context']['endpoints'].items():expected.append({'case_id':cid,'source_case_id':c['case_id'],'selection_id':c['selection_id'],
                'candidate':'origin','metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity']})
        if entries!=expected:raise InvalidArtifact('origin physical source entries differ')
    elif m['envelope_stage']=='transfer_candidates':
        sm=read_json(verify(source['source_manifest']));a,s,_=scope(verify(sm['inventory']),sm['shard'])
        entries=source['tasks'];maximum=8*len(s['cases'])
        if source['protocol_id']!=POOL_PROTOCOL:raise InvalidArtifact('candidate source protocol differs')
    else:raise InvalidArtifact('unrecognized finite transfer stage')
    if(m['numerical_policy_id']!=scalar.PROFILE or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}
       or m['tasks']!=m['all_tasks'] or len(m['tasks'])!=2*len(entries) or m['maximum_tasks']!=maximum or len(m['tasks'])>maximum or m['reused']):
        raise InvalidArtifact('finite strict scalar scope/profile differs')
    for p in m['implementation'].values():verify(p)
    by={(e['case_id'],e['metal'],e['candidate']):e for e in entries}
    for t in m['tasks']:
        ep=by[t['case_id'],t['metal'],t['candidate']]
        if not context.reusable_state(t,ep) or verify(t['input']).read_text()!=scalar.recipe(t['charge'],t['medium'],'fresh'):
            raise InvalidArtifact('strict physical input or recipe changed')
    return dry_run(manifest)


def collect_origins(stage,output):
    stage=Path(stage);ready=read_json(stage/'READY.json');inp=read_json(verify(ready['inputs']));mm=read_json(verify(ready['MACE_manifest']))
    lp=verify(ready['scalar_manifest']);validate_scalar(lp);lm=read_json(lp);rows=[];strict=scalar_engine()
    for item in inp['cases']:
        c=item['prepared'];cid=c['case_id']+'__envelope_'+c['selection_id'].rsplit('_',1)[-1]
        r={'case_id':cid,'source_case_id':c['case_id'],'selection_id':c['selection_id'],'source':c,'mapping':c['maps'],
            'source_preparation':c['representations']['context']['preparation'],'native_endpoints':{},'solvent_endpoints':{}}
        for z,ep in c['representations']['context']['endpoints'].items():
            try:
                receipt=verify(ready['MACE_manifest']).parent/'results'/(cid+'__context__'+z)/'result.json'
                native=native_endpoint(read_json(receipt),mm['model']);union.origin(native,mm['model'])
                if not context.reusable_state(ep,native):raise InvalidArtifact('actual native source differs')
                r['native_endpoints'][z]={**native,'status':'complete'}
            except (OSError,KeyError,InvalidArtifact)as exc:r['native_endpoints'][z]={'status':'unavailable','reason':str(exc)}
            r['solvent_endpoints'][z]={}
            for medium in ('vacuum','alpb'):
                t=next(t for t in lm['tasks']if(t['case_id'],t['metal'],t['medium'])==(cid,z,medium))
                r['solvent_endpoints'][z][medium]=strict.endpoint(lp,t)
        r['status']='complete'if all(r['native_endpoints'][z]['status']=='complete'and all(q['status']=='complete'for q in r['solvent_endpoints'][z].values())for z in ('Ca','La'))else'unavailable'
        rows.append(r)
    result={'protocol_id':SCOPE_PROTOCOL,'inputs':ready['inputs'],'inventory':inp['inventory'],'shard':inp['shard'],'model':mm['model'],
        'rows':rows,'denominator':len(rows),'complete':sum(r['status']=='complete'for r in rows),'new_molecular_calls':0}
    write_new(output,result);return {'collection':record(output),'complete':result['complete'],'denominator':len(rows)}


def prepare_searches(origins,output):
    d=read_json(origins);inp=read_json(verify(d['inputs']));base=read_json(verify(inp['template']));a,s,ref=scope(verify(d['inventory']),d['shard'])
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];cases=[];missing=[]
    for r in d['rows']:
        if r['status']!='complete':missing.append({'case_id':r['case_id'],'reason':'actual origin incomplete'});continue
        try:ts,choice=build_tasks(r,base['model'])
        except InvalidArtifact as exc:missing.append({'case_id':r['case_id'],'reason':str(exc)});continue
        tasks.extend(ts);cases.append({'case_id':r['case_id'],'status':'prepared','origin_row':r,'union':r['source'],'selection':choice})
    m={k:base[k]for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')}
    m.update(protocol_id=SEARCH_PROTOCOL,settings=precision.SETTINGS,agreement=a['agreement'],reference=a['reference'],inventory=d['inventory'],shard=d['shard'],
        origin_collection=record(origins),inputs=d['inputs'],cases=cases,tasks=tasks,unavailable=missing,
        declared_case_ids=[c['case_id']for c in cases],declared_all_sources=[r['case_id']for r in d['rows']],
        maximum_optimizer_starts=len(tasks),declared_maximum_optimizer_starts=2*len(s['cases']),maximum_cross_MACE_calls=2*len(cases),maximum_GFN2_calls=8*len(cases),
        new_origin_calls=0,new_DFT_calls=0,production_changed=False,population='all100_envelope_shard_'+str(d['shard']),
        optimizer_software={'version':scipy.__version__,'wrapper':record(slsqp_wrapper.__file__),'kernel':record(slsqp_kernel.__file__)},
        execution_adapter=record(__file__),implementation=union.snapshot(out/'implementation'))
    write_new(out/'manifest.json',m);v=validate_searches(out/'manifest.json');write_new(out/'PREFLIGHT.json',v);return v


def validate_searches(manifest):
    m=read_json(manifest);a,s,ref=scope(verify(m['inventory']),m['shard']);d=read_json(verify(m['origin_collection']));inp=read_json(verify(m['inputs']));base=read_json(verify(inp['template']))
    ids=[r['prepared']['case_id']+'__envelope_'+r['selection_id'].rsplit('_',1)[-1]for r in s['cases']]
    if(m['protocol_id']!=SEARCH_PROTOCOL or m['settings']!=ref['optimizer_settings']or m['reference']!=a['reference']
       or m['declared_all_sources']!=ids or [r['case_id']for r in d['rows']]!=ids or len(m['tasks'])!=2*len(m['cases'])
       or set(m['declared_case_ids'])|{r['case_id']for r in m['unavailable']}!=set(ids) or m['maximum_optimizer_starts']!=len(m['tasks'])):
        raise InvalidArtifact('finite search scope/profile differs')
    for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources'):
        if m[k]!=base[k]:raise InvalidArtifact('scientific runtime differs')
    verify(m['execution_adapter'])
    for p in m['implementation'].values():verify(p)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer runtime changed')
    actual=[];by={r['case_id']:r for r in d['rows']}
    for c in m['cases']:
        r=by[c['case_id']]
        if c['origin_row']!=r or c['union']!=r['source']or r['status']!='complete':raise InvalidArtifact('actual source row differs')
        ts,choice=build_tasks(r,m['model']);selector_equivalence(c['selection'],choice);actual.extend(ts)
    for old,new in zip(m['tasks'],actual):selector_equivalence(old['selector'],new['selector']);new['selector']=old['selector']
    if actual!=m['tasks']:raise InvalidArtifact('actual native forces/state/task differs')
    return {'status':'validated','manifest':record(manifest),'declared_sources':len(ids),'searches':len(actual),'missing':len(m['unavailable'])}


def low_prepare(manifest):
    m=read_json(manifest);sm=read_json(verify(m['source_manifest']));inp=read_json(verify(sm['inputs']));a,s,_=scope(verify(sm['inventory']),sm['shard'])
    root=Path(manifest).parent/'solvent';dest=root/'shard_0'
    if m['tasks']:
        scalar_engine().prepare(m['tasks'],m,record(manifest),verify(inp['qualification']),verify(a['agreement']),dest,'transfer_candidates',8*len(s['cases']))
        shards=[record(dest/'manifest.json')]
    else:
        root.mkdir();shards=[]
    write_new(root/'INDEX.json',{'pool_manifest':record(manifest),'shards':shards})


def engine():
    u,p=precision.engine();u.PROTOCOL=SEARCH_PROTOCOL;u.POOL_PROTOCOL=POOL_PROTOCOL;u.validate=validate_searches
    original_write=u.write_new
    def strict_write(path,value):
        if isinstance(value,dict)and value.get('protocol_id')==POOL_PROTOCOL:value={**value,'numerical_policy_id':scalar.PROFILE}
        return original_write(path,value)
    u.write_new=strict_write;original_validate=u.validate_pool
    def strict_validate(manifest):
        if read_json(manifest)['numerical_policy_id']!=scalar.PROFILE:raise InvalidArtifact('pool strict numerical profile differs')
        return original_validate(manifest)
    u.validate_pool=strict_validate;p.validate=strict_validate;p.low_prepare=low_prepare;p.completed=strict_completed
    return u,p


def execute_searches(manifest):return engine()[0].execute(manifest)
def collect_searches(manifest,output):return engine()[0].collect(manifest,output)
def prepare_pool(proposals,agreement,output):return engine()[0].prepare_pool(proposals,agreement,output)
def validate_pool(manifest):return engine()[0].validate_pool(manifest)
def execute_mace(manifest):return engine()[1].execute_mace(manifest)
def collect_pool(manifest,output):
    validate_pool(manifest);return engine()[1].collect(manifest,output)
def execute_scalar(manifest):return scalar_engine().execute(manifest)
def collect_scalar(manifest,output):return scalar_engine().collect(manifest,output)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    specs={'prepare_origins':('inventory','shard','template','qualification','output'),'validate_scalar':('manifest',),
        'collect_origins':('stage','output'),'prepare_searches':('origins','output'),'validate_searches':('manifest',),
        'execute_searches':('manifest',),'collect_searches':('manifest','output'),'prepare_pool':('proposals','agreement','output'),
        'validate_pool':('manifest',),'execute_mace':('manifest',),'collect_pool':('manifest','output'),
        'execute_scalar':('manifest',),'collect_scalar':('manifest','output')}
    for op,keys in specs.items():
        q=s.add_parser(op)
        for k in keys:q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('command')](**a),indent=2))
