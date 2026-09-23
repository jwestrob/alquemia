"""Strict fresh native scalar replay of the complete225 tenfold-precision panel."""
from __future__ import annotations
import argparse,copy,importlib.util,json,re,shutil,time
from collections import Counter
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new,energy,xyz
from precision_pool_continuation import data,normalized,CANDIDATES,FOLDS
from strict_native_comparator import CASES as COMPARATOR_CASES
from strict_native_pool import recipe,collect_fresh
from compact_solvation import diagnostics,input_text
from structure_informed_starts import scf_details
from run_orca_task_manifest import _completed_attempt_is_valid,load_manifest_tasks
from native_gfn2_rank_panel import native_reuse
from union_adaptive import snapshot
from nikasha_pool import choose_rows,score
from accommodation_nonlinear import relative_components
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome,strict_summary,summarize_rows

PROTOCOL='Nikasha_tenfold_precision225_strict_native_fresh_transfer_v1'
REFERENCE_HASH='171bd31d466ff97ef6073ce23286af8519ef23690f7d6ad23bb444cab30ecd5c'
REUSED=tuple(FOLDS)+tuple(COMPARATOR_CASES)
SETTINGS={'TolE_hartree':1e-10,'MaxIter':500,'electronic_temperature_K':300,'native_mixer':True,
          'mpi_ranks':1,'workers':32,'memory_MiB':65536,'origin_selection_tolerance_kcal_mol':.1,
          'cell_comparison_kcal_mol':.1,'R_comparison_kcal_mol':.2,'candidate_ids':list(CANDIDATES),
          'declared_sources':225,'prepared_sources':208,'reused_pools':8,'fresh_pools':200,'new_scalar_calls':2400}
VARIANTS=('operational','mathematical')
METHODS=('union_precision_strict','union_precision_strict_mathematical')


def tid(r):return '__'.join(r[k] for k in ('case_id','candidate','metal','medium'))


def checked_reference(reference,qualification):
    if record(reference)['sha256']!=REFERENCE_HASH:raise InvalidArtifact('frozen strict32 reference differs')
    ref=read_json(reference);q=read_json(qualification)
    if q['reference']!=record(reference) or q['agreeing_cells']!=384 or q['counts']['all32']['qualified']!=32:
        raise InvalidArtifact('completed strict32 qualification required')
    fresh=ref['branches']['fresh']
    if fresh['calibration_denominator']!=25 or fresh['crystals_or_noncanonical_used_for_fit'] or any(v['available_calibration']!=25 or not v['bands'] for v in fresh['variants'].values()):
        raise InvalidArtifact('canonical-only fresh reference incomplete')
    return fresh


def source_cell(cid,q,z,medium,cell,model):
    low=cell['low'][medium];m,ts=normalized(low['manifest']['path'],low['manifest']['sha256'])
    t=next(t for t in m.get('all_tasks',m['tasks']) if t['task_id']==low['task_id']);receipt=data(low['receipt']);op=verify(low['output'])
    if not _completed_attempt_is_valid(verify(low['receipt']),op,manifest_sha256=low['manifest']['sha256'],task=ts[t['task_id']],runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
        raise InvalidArtifact('actual old receipt incompatible: '+cid)
    body=re.sub(r'^\s*MaxIter 500\s*\n','',verify(t['input']).read_text(),flags=re.M)
    if body!=input_text(t['charge'],1,medium,'native') or energy(op)!=low['energy_hartree'] or (t['metal'],t['medium'],t['multiplicity'])!=(z,medium,1):
        raise InvalidArtifact('old recipe/state differs')
    a=xyz(verify(cell['xyz']));b=xyz(verify(t['xyz']))
    delta=float(np.max(np.abs(np.asarray([x[1:] for x in a])-np.asarray([x[1:] for x in b])))) if len(a)==len(b) else float('inf')
    if len(a)!=len(b) or [x[0] for x in a]!=[x[0] for x in b] or delta>1e-12:raise InvalidArtifact('pool/cell coordinates differ')
    audit=low.get('audit') or low.get('component_audit') or diagnostics(low,t);verify(audit['parameter_export'])
    detail=scf_details(op.read_text());native=native_reuse(cell,t,model)
    if not detail['native_mixer_observed'] or (detail['charge'],detail['multiplicity'])!=(t['charge'],1):raise InvalidArtifact('source native method/state differs')
    return {'case_id':cid,'candidate':q,'metal':z,'medium':medium,'xyz':t['xyz'],'input':t['input'],
        'output':low['output'],'receipt':low['receipt'],'source_manifest':low['manifest'],'source_task_id':low['task_id'],
        'source_energy_hartree':low['energy_hartree'],'charge':t['charge'],'multiplicity':1,
        'electron_count':detail['electrons'],'parameter_export':audit['parameter_export'],'orca':receipt['orca_executable'],
        'MACE_reuse':native,'pool_coordinate_max_delta_A':delta}


def compare_reuse_source(src,task):
    old=task['source']
    if any(task[k]!=src[k] for k in ('case_id','candidate','metal','medium','charge','multiplicity')):
        raise InvalidArtifact('reuse identity/state differs')
    if verify(task['input']).read_text()!=recipe(src['charge'],src['medium'],'fresh'):
        raise InvalidArtifact('reuse not strict fresh recipe')
    a=xyz(verify(src['xyz']));b=xyz(verify(task['xyz']))
    if len(a)!=len(b) or [x[0] for x in a]!=[x[0] for x in b] or max(abs(x[i]-y[i]) for x,y in zip(a,b) for i in (1,2,3))>1e-12:
        raise InvalidArtifact('reuse coordinates differ')
    if old['parameter_export']['sha256']!=src['parameter_export']['sha256'] or old['electron_count']!=src['electron_count']:
        raise InvalidArtifact('reuse native Hamiltonian differs')


def prepare(precision225,reference,qualification,comparator,agreement,output):
    began=time.monotonic();ref=checked_reference(reference,qualification);prior=read_json(precision225)
    if len(prior['rows'])!=225 or len({r['case_id'] for r in prior['rows']})!=225 or Counter(r['precision_status'] for r in prior['rows'])!={'available':208,'preparation_unavailable':17}:
        raise InvalidArtifact('original225/208/17 scope differs')
    comp=read_json(comparator)
    if comp['reference']!=record(reference) or comp['precision225']!=record(precision225) or [c['case_id'] for c in comp['cases']]!=list(COMPARATOR_CASES):
        raise InvalidArtifact('four separately owned comparator pools differ')
    fresh=data(ref['collections']['fresh']);fm=data(fresh['manifest'])
    if fresh['complete']!=384 or fresh['branch']!='fresh':raise InvalidArtifact('strict32 fresh reuse incomplete')
    reuse_tasks={tid(t):t for t in fm['tasks'] if t['case_id'] in FOLDS}
    reuse_tasks.update({tid(t):t for t in comp['tasks']})
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('output must be under workspaces')
    out.mkdir(parents=True,exist_ok=False);rows=[];cases=[];orca=None;model=None
    for old in prior['rows']:
        cid=old['case_id'];case={'case_id':cid,'expected_class':old['expected_class'],'root_case_id':old['root_case_id'],
            'source_conditioning_metal':old['source_conditioning_metal'],'status':old['precision_status'],
            'reason':old.get('precision_reason'),'collection':old.get('precision_collection'),'reuse':cid in REUSED}
        if old['precision_status']=='available':
            coll=data(old['precision_collection']);poolcase=next(c for c in coll['cases'] if c['case_id']==cid);pm=data(coll['manifest'])
            if poolcase['pool']['status']!='available' or choose_rows(poolcase['matrix'],list(CANDIDATES))!=poolcase['pool']:
                raise InvalidArtifact('actual completed common3 pool differs')
            if model is None:model=pm['model']
            if pm['model']!=model or pm['settings']['candidate_order']!=list(CANDIDATES) or pm['settings']['origin_selection_tolerance_kcal_mol']!=.1:
                raise InvalidArtifact('old model/pool policy differs')
            for q in CANDIDATES:
                for z in ('Ca','La'):
                    for medium in ('vacuum','alpb'):
                        src=source_cell(cid,q,z,medium,poolcase['matrix'][z][q],model)
                        if orca is None:orca=src['orca']
                        if src['orca']!=orca:raise InvalidArtifact('mixed executable')
                        if cid in REUSED:compare_reuse_source(src,reuse_tasks[tid(src)])
                        rows.append(src)
        cases.append(case);print(json.dumps({'source':cid,'status':case['status'],'reuse':case['reuse']}),flush=True)
    if len(rows)!=2496 or {r['case_id'] for r in rows if r['case_id'] in REUSED}!=set(REUSED):raise InvalidArtifact('208x12 cells/eight reuses required')
    impl=snapshot(out/'implementation')
    inventory={'protocol_id':PROTOCOL,'settings':SETTINGS,'precision225':record(precision225),'reference':record(reference),
        'qualification':record(qualification),'agreement':record(agreement),'cases':cases,'rows':rows,'model':model,'orca':orca,
        'reuse':{'strict32':{'collection':ref['collections']['fresh'],'manifest':fresh['manifest'],'case_ids':list(FOLDS)},
                 'comparator':{'manifest':record(comparator),'collection_path':str(Path(comparator).parent/'COLLECTION.json'),'case_ids':list(COMPARATOR_CASES)}},
        'implementation':impl,'preparation_wall_seconds':time.monotonic()-began,'new_molecular_calls':0,'production_changed':False}
    ip=out/'INVENTORY.json';write_new(ip,inventory);freshids=[c['case_id'] for c in cases if c['status']=='available' and not c['reuse']]
    shards=[]
    for i in range(4):
        ids=freshids[i::4];tasks=[];base=out/('shard_'+str(i))
        for src in rows:
            if src['case_id'] not in ids:continue
            d=base/'tasks'/tid(src);d.mkdir(parents=True);shutil.copyfile(verify(src['xyz']),d/'core.xyz')
            (d/'endpoint.inp').write_text(recipe(src['charge'],src['medium'],'fresh'))
            t={k:src[k] for k in ('case_id','candidate','metal','medium','charge','multiplicity')}
            t.update(task_id=tid(src),source=src,seed_source=None,gradient_requested=False,xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),
                output_path=str(d/'endpoint.out'),active_seed_path=str(d/'endpoint.runtime.xtbw'),active_gbw_path=str(d/'endpoint.runtime.gbw'));tasks.append(t)
        m={'protocol_id':PROTOCOL,'branch':'fresh','settings':SETTINGS,'inventory':record(ip),'shard':i,'case_ids':ids,'tasks':tasks,
            'orca':orca,'agreement':record(agreement),'implementation':impl,'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},
            'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
            'new_MACE_DFT_optimization_calls':0,'production_changed':False}
        mp=base/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(base/'PREFLIGHT.json',v);shards.append(v)
    result={'inventory':record(ip),'shards':shards,'fresh_pools':200,'fresh_calls':2400,'whole_pool_reuses':8,'old_unavailable':17,'preparation_wall_seconds':time.monotonic()-began}
    write_new(out/'PREFLIGHT.json',result);return result


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp);inv=data(m['inventory']);allnew=[c['case_id'] for c in inv['cases'] if c['status']=='available' and not c['reuse']]
    if m['protocol_id']!=PROTOCOL or inv['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or inv['settings']!=SETTINGS or m['branch']!='fresh':raise InvalidArtifact('fixed policy differs')
    if m['shard'] not in range(4) or m['case_ids']!=allnew[m['shard']::4] or len(m['case_ids'])!=50 or len(m['tasks'])!=600:
        raise InvalidArtifact('finite shard membership differs')
    expected={(c,q,z,s) for c in m['case_ids'] for q in CANDIDATES for z in ('Ca','La') for s in ('vacuum','alpb')}
    if {(t['case_id'],t['candidate'],t['metal'],t['medium']) for t in m['tasks']}!=expected or len({t['task_id'] for t in m['tasks']})!=600:raise InvalidArtifact('paired full-matrix tasks differ')
    for k in ('reference','qualification','precision225','agreement'):verify(inv[k])
    if m['agreement']!=inv['agreement'] or m['orca']!=inv['orca']:raise InvalidArtifact('agreement/executable changed')
    verify(m['orca'])
    for pin in m['implementation'].values():verify(pin)
    if m['implementation']!=inv['implementation'] or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}:raise InvalidArtifact('snapshot/resource profile changed')
    src={tid(r):r for r in inv['rows']};load_manifest_tasks(mp)
    for t in m['tasks']:
        if t['source']!=src[t['task_id']] or any(t[k]!=t['source'][k] for k in ('case_id','candidate','metal','medium','charge','multiplicity')) or t['multiplicity']!=1 or t['seed_source'] is not None or t['gradient_requested']:
            raise InvalidArtifact('state/source differs')
        if verify(t['xyz']).read_bytes()!=verify(t['source']['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(t['charge'],t['medium'],'fresh'):raise InvalidArtifact('coordinates/recipe differ')
        d=mp.parent/'tasks'/t['task_id']
        if Path(t['output_path'])!=d/'endpoint.out' or Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('task outside disjoint shard directory')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp'}:raise InvalidArtifact('nonfresh directory; collect existing attempt, no automatic retry')
    return {'status':'validated','manifest':record(mp),'shard':m['shard'],'cases':50,'fresh_calls':600,'new_calls_in_validation':0}


def execute(manifest):
    path=Path(__file__).with_name('strict_native_pool.py');spec=importlib.util.spec_from_file_location('_strict225_executor',path);engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
    engine.validate=validate;return engine.execute(manifest)


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def idx(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']} if p.exists() else {}
    before,after=idx('SEEDS_BEFORE.json'),idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=collect_fresh(mp,t,before.get(t['task_id']),after.get(t['task_id']));r['observed_TolE_hartree']=None
        if r['status']=='complete':
            text=verify(r['actual']['output']).read_text();tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text);r['observed_TolE_hartree']=float(tol[0]) if len(tol)==1 else None
            if r['observed_TolE_hartree']!=1e-10 or data(r['actual']['receipt'])['parallelism']['nprocs']!=1:r.update(status='audit_failed',reason='effective tolerance/rank differs',energy_hartree=None)
        rows.append(r)
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'inventory':m['inventory'],'shard':m['shard'],'rows':rows,
        'cell_denominator':600,'complete_cells':sum(r['status']=='complete' for r in rows),'new_molecular_calls_in_collection':0}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


def matrix_from_rows(case,low):
    matrix={z:{} for z in ('Ca','La')};cid=case['case_id']
    for z in matrix:
        for q in CANDIDATES:
            rr={s:low['__'.join((cid,q,z,s))] for s in ('vacuum','alpb')};ok=all(r['status']=='complete' for r in rr.values());old=case['matrix'][z][q]
            matrix[z][q]={'status':'complete' if ok else 'unavailable','xyz':old['xyz'],'MACE':old['MACE'],
                'components':{'MACE_eV':old['components']['MACE_eV'],'GFN2_vacuum_hartree':rr['vacuum'].get('energy_hartree'),'GFN2_ALPB_hartree':rr['alpb'].get('energy_hartree')} if ok else None,
                'low':{s:r.get('actual') for s,r in rr.items()},'missing_reasons':{s:r.get('reason') for s,r in rr.items() if r['status']!='complete'}}
    return matrix


def compare(inventory,collections,output):
    inv=read_json(inventory);ref=checked_reference(verify(inv['reference']),verify(inv['qualification']));prior=data(inv['precision225'])
    low={};pins=[];shards=set();source={tid(r):r for r in inv['rows']}
    def add(coll,pin,allowed):
        m=data(coll['manifest']);ts={tid(t):t for t in m['tasks']}
        for r in coll['rows']:
            if r['case_id'] not in allowed:continue
            k=tid(r)
            if k in low:raise InvalidArtifact('duplicate actual scalar cell')
            compare_reuse_source(source[k],ts[k])
            if r['status']=='complete':
                if r.get('observed_TolE_hartree')!=1e-10:raise InvalidArtifact('strict tolerance not confirmed')
                for field in ('output','receipt'):verify(r['actual'][field])
            low[k]=r
        pins.append(pin)
    for path in collections:
        c=read_json(path);m=data(c['manifest']);validate(verify(c['manifest']))
        if c['inventory']!=record(inventory) or c['shard']!=m['shard'] or c['shard'] in shards:raise InvalidArtifact('duplicate/changed shard collection')
        shards.add(c['shard']);add(c,record(path),set(m['case_ids']))
    if shards!={0,1,2,3}:raise InvalidArtifact('all four terminal/partial collections required')
    rr=inv['reuse']['strict32'];add(data(rr['collection']),rr['collection'],set(rr['case_ids']))
    rr=inv['reuse']['comparator'];cp=Path(rr['collection_path']);cc=read_json(cp)
    if cc['manifest']!=rr['manifest']:raise InvalidArtifact('separately owned comparator changed')
    add(cc,record(cp),set(rr['case_ids']))
    if set(low)!=set(source):raise InvalidArtifact('all2496 cells, including failures, required')
    bands={**prior['bands'],**{name:ref['variants'][v]['bands'] for name,v in zip(METHODS,VARIANTS)}};rows=[];full=[]
    for old in prior['rows']:
        cid=old['case_id'];supported=old['precision_status']=='available';new={};work={};numbers=[];selected={};pool=None
        if supported:
            case=next(c for c in data(old['precision_collection'])['cases'] if c['case_id']==cid);matrix=matrix_from_rows(case,low);pool=choose_rows(matrix,list(CANDIDATES));ok=pool['status']=='available'
            for q in CANDIDATES:
                for z in ('Ca','La'):
                    a=case['matrix'][z][q]['components'];b=matrix[z][q]['components']
                    numbers.append({'candidate':q,'metal':z,'vacuum_delta_kcal_mol':(b['GFN2_vacuum_hartree']-a['GFN2_vacuum_hartree'])*HA_TO_KCAL if b else None,
                        'ALPB_delta_kcal_mol':(b['GFN2_ALPB_hartree']-a['GFN2_ALPB_hartree'])*HA_TO_KCAL if b else None,
                        'MACE_delta_kcal_mol':0.0 if b else None})
            full.append({'case_id':cid,'source_collection':old['precision_collection'],'matrix':matrix,'pool':pool,'reused':cid in REUSED})
        else:ok=False
        for name,v in zip(METHODS,VARIANTS):
            val=pool[v]['composite_R_model_kcal_mol'] if ok else None;call=decision(val,bands[name]);new[name]={'R':val,'decision':call,'outcome':outcome(call,old['expected_class'])}
            if ok:
                selected[v]={z:pool['rows'][z][v+'_candidate'] for z in ('Ca','La')}
                work[v]={z:relative_components(matrix[z][selected[v][z]]['components'],matrix[z]['origin']['components']) for z in ('Ca','La')}
                origin=score(matrix['Ca']['origin']['components'],matrix['La']['origin']['components'])['composite_R_model_kcal_mol']
                if abs(val-origin-(work[v]['Ca']['composite_kcal_mol']-work[v]['La']['composite_kcal_mol']))>1e-6:raise InvalidArtifact('score work/sign algebra differs')
        rows.append({**old,'methods':{**old['methods'],**new},'strict_status':pool['status'] if pool else 'prior_preparation_unavailable',
            'strict_reused':cid in REUSED,'strict_components_delta':numbers,'strict_selected':selected or None,'strict_work':work or None,
            'strict_delta_R':{v:new[name]['R']-old['methods'][prev]['R'] if ok else None for name,v,prev in zip(METHODS,VARIANTS,('union_precision','union_precision_mathematical'))},
            'strict_old_reference_transfer':{v:decision(new[name]['R'],prior['bands'][prev]) for name,v,prev in zip(METHODS,VARIANTS,('union_precision','union_precision_mathematical'))}})
    idx={r['case_id']:r for r in rows};groups=[]
    for old in prior['pools']:
        new={}
        for name in METHODS:
            if old['pool'] in ('La4','Ca5'):value=strict_summary(old['members'],idx,name,bands[name],old['expected_class'])
            else:
                arms=[next(p for p in groups if p['root_case_id']==old['root_case_id'] and p['pool']==a)['methods'][name] for a in ('La4','Ca5')]
                val=None if any(a['R'] is None for a in arms) else sum(a['R'] for a in arms)/2;call=decision(val,bands[name])
                value={'R':val,'decision':call,'outcome':outcome(call,old['expected_class']),'required_members':len(old['members']),
                    'missing_members':sorted({x for a in arms for x in a['missing_members']}),'within_protein_range':None,'equal_mean_of_complete_arm_medians':True}
            new[name]=value
        groups.append({**old,'methods':{**old['methods'],**new}})
    triples=[{**t,'methods':{**t['methods'],**{n:strict_summary(t['members'],idx,n,bands[n],t['expected_class']) for n in METHODS}}} for t in prior['triples']]
    subsets={'all225':rows,'La100':[r for r in rows if r['source_conditioning_metal']=='La'],'Ca125':[r for r in rows if r['source_conditioning_metal']=='Ca'],
        **{k:[r for r in groups if r['pool']==k] for k in ('La4','Ca5','balanced')},'La100_triples':triples}
    names=list(rows[0]['methods']);counts={k:summarize_rows(rs,names) for k,rs in subsets.items()};matched={}
    for k,rs in subsets.items():
        matched[k]={}
        for other in names:
            if other in METHODS:continue
            common=[r for r in rs if r['methods'][METHODS[0]]['R'] is not None and r['methods'][other]['R'] is not None]
            matched[k][other]={'declared':len(rs),'common':len(common),'counts':summarize_rows(common,(METHODS[0],other)),
                'transitions':dict(Counter(r['methods'][other]['outcome']+'->'+r['methods'][METHODS[0]]['outcome'] for r in common))}
    result={'protocol_id':PROTOCOL,'inventory':record(inventory),'reference':inv['reference'],'reference_branch':'fresh','qualification':inv['qualification'],
        'prior':inv['precision225'],'collections':pins,'rows':rows,'pools':groups,'triples':triples,'actual_pools':full,'bands':bands,
        'counts':counts,'matched':matched,'complete_cells':sum(r['status']=='complete' for r in low.values()),'cell_denominator':2496,
        'new_thresholds_fitted':False,'production_changed':False,'new_MACE_DFT_optimization_calls':0,'implementation':record(__file__)}
    write_new(output,result);return counts['all225']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('precision225','reference','qualification','comparator','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'compare':('inventory','output')}.items():
        s=sub.add_parser(op)
        for k in fields:s.add_argument('--'+k.replace('_','-'),required=True,type=Path)
        if op=='compare':s.add_argument('--collections',required=True,type=Path,nargs='+')
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
