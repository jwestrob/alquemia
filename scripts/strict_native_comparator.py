"""Four consumed tenfold-precision pools under the qualified fresh strict scalar recipe."""
from __future__ import annotations
import argparse,importlib.util,json,re,shutil
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,energy,xyz,HA_TO_KCAL
from precision_pool_continuation import data,normalized,CANDIDATES
from compact_solvation import diagnostics,input_text
from structure_informed_starts import scf_details
from strict_native_pool import recipe,collect_fresh
from run_orca_task_manifest import _completed_attempt_is_valid,load_manifest_tasks
from union_adaptive import snapshot
from nikasha_pool import choose_rows
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome

PROTOCOL='native_GFN2_strict_fresh_four_tenfold_comparators_v1'
CASES=('a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-3',
       'a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4',
       'a8r3s4-pqq-la_model__conditioned_La__seed-1_sample-1',
       'a8r3s4-pqq-la_model__conditioned_La__seed-1_sample-3')


def tid(r):return '__'.join(r[k] for k in ('case_id','candidate','metal','medium'))


def prepare(precision225,reference,qualification,agreement,output):
    prior=read_json(precision225);ref=read_json(reference);qual=read_json(qualification)
    if qual['reference']!=record(reference) or qual['agreeing_cells']!=384 or qual['counts']['all32']['qualified']!=32:raise InvalidArtifact('strict32 qualification required')
    if any(not r['bands'] for r in ref['branches']['fresh']['variants'].values()):raise InvalidArtifact('frozen fresh reference unavailable')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];groups=[]
    for cid in CASES:
        row=next(r for r in prior['rows'] if r['case_id']==cid);collection=row['precision_collection'];case=next(c for c in data(collection)['cases'] if c['case_id']==cid)
        if case['pool']['status']!='available':raise InvalidArtifact('declared source pool unavailable')
        groups.append({'case_id':cid,'expected_class':row['expected_class'],'collection':collection,'role':'consumed_noncanonical_development'})
        for q in CANDIDATES:
            for z in ('Ca','La'):
                cell=case['matrix'][z][q]
                for medium in ('vacuum','alpb'):
                    low=cell['low'][medium];m,ts=normalized(low['manifest']['path'],low['manifest']['sha256']);t=next(t for t in m['tasks'] if t['task_id']==low['task_id']);receipt=data(low['receipt'])
                    op=verify(low['output'])
                    if not _completed_attempt_is_valid(verify(low['receipt']),op,manifest_sha256=low['manifest']['sha256'],task=ts[t['task_id']],runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']):raise InvalidArtifact('actual old receipt incompatible')
                    body=re.sub(r'^\s*MaxIter 500\s*\n','',verify(t['input']).read_text(),flags=re.M)
                    if body!=input_text(t['charge'],1,medium,'native') or energy(op)!=low['energy_hartree'] or (t['metal'],t['medium'],t['multiplicity'])!=(z,medium,1):raise InvalidArtifact('old recipe/state differs')
                    a=xyz(verify(cell['xyz']));b=xyz(verify(t['xyz']))
                    if [x[0] for x in a]!=[x[0] for x in b] or np.max(np.abs(np.asarray([x[1:] for x in a])-np.asarray([x[1:] for x in b])))>1e-12:raise InvalidArtifact('source/paired coordinates differ')
                    audit=low.get('audit') or diagnostics(low,t);details=scf_details(op.read_text())
                    src={'xyz':t['xyz'],'input':t['input'],'output':low['output'],'receipt':low['receipt'],'source_cell':low,'source_energy_hartree':low['energy_hartree'],
                        'charge':t['charge'],'multiplicity':1,'electron_count':details['electrons'],'parameter_export':audit['parameter_export'],'orca':receipt['orca_executable']}
                    r={'case_id':cid,'candidate':q,'metal':z,'medium':medium,'charge':t['charge'],'multiplicity':1,'source':src,'seed_source':None,'gradient_requested':False}
                    r['task_id']=tid(r);d=out/'tasks'/r['task_id'];d.mkdir(parents=True);shutil.copyfile(verify(t['xyz']),d/'core.xyz');(d/'endpoint.inp').write_text(recipe(t['charge'],medium,'fresh'))
                    r.update(xyz=record(d/'core.xyz'),input=record(d/'endpoint.inp'),output_path=str(d/'endpoint.out'),active_seed_path=str(d/'endpoint.runtime.xtbw'),active_gbw_path=str(d/'endpoint.runtime.gbw'));tasks.append(r)
    impl=snapshot(out/'implementation');m={'protocol_id':PROTOCOL,'branch':'fresh','tasks':tasks,'cases':groups,'reference':record(reference),'qualification':record(qualification),'precision225':record(precision225),
        'agreement':record(agreement),'orca':tasks[0]['source']['orca'],'implementation':impl,'case_denominator':4,'cell_denominator':48,
        'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
        'new_MACE_DFT_optimization_calls':0,'production_changed':False}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp)
    if m['protocol_id']!=PROTOCOL or m['branch']!='fresh' or [r['case_id'] for r in m['cases']]!=list(CASES) or len(m['tasks'])!=48:raise InvalidArtifact('fixed scope differs')
    if {(r['case_id'],r['candidate'],r['metal'],r['medium']) for r in m['tasks']}!={(cid,q,z,s) for cid in CASES for q in CANDIDATES for z in ('Ca','La') for s in ('vacuum','alpb')}:raise InvalidArtifact('matrix denominator differs')
    for name in ('reference','qualification','precision225','agreement','orca'):verify(m[name])
    for pin in m['implementation'].values():verify(pin)
    if m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}:raise InvalidArtifact('resource profile differs')
    load_manifest_tasks(mp)
    for t in m['tasks']:
        if t['seed_source'] is not None or t['gradient_requested'] or t['charge']!=t['source']['charge'] or t['multiplicity']!=1 or verify(t['xyz']).read_bytes()!=verify(t['source']['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(t['charge'],t['medium'],'fresh'):raise InvalidArtifact('state/recipe/source differs')
        d=Path(t['output_path']).parent
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('basename differs')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp'}:raise InvalidArtifact('nonfresh task directory')
    return {'status':'validated','manifest':record(mp),'fresh_calls':48,'cases':4,'new_molecular_calls_in_validation':0}


def execute(manifest):
    # Private instance only: retain tested32worker seedless execution and receipts.
    path=Path(__file__).with_name('strict_native_pool.py');spec=importlib.util.spec_from_file_location('_strict_four_executor',path);engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
    engine.validate=validate
    return engine.execute(manifest)


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def idx(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']} if p.exists() else {}
    before=idx('SEEDS_BEFORE.json');after=idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=collect_fresh(mp,t,before.get(t['task_id']),after.get(t['task_id']));r['observed_TolE_hartree']=None
        if r['status']=='complete':
            text=verify(r['actual']['output']).read_text();tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text);r['observed_TolE_hartree']=float(tol[0]) if len(tol)==1 else None
            if r['observed_TolE_hartree']!=1e-10 or data(r['actual']['receipt'])['parallelism']['nprocs']!=1:r.update(status='audit_failed',reason='effective tolerance/rank differs',energy_hartree=None)
        rows.append(r)
    low={tid(r):r for r in rows};cases=[];ref=data(m['reference'])['branches']['fresh']
    for group in m['cases']:
        cid=group['case_id'];old=next(c for c in data(group['collection'])['cases'] if c['case_id']==cid);matrix={z:{} for z in ('Ca','La')}
        for z in matrix:
            for q in CANDIDATES:
                rr={s:low['__'.join((cid,q,z,s))] for s in ('vacuum','alpb')};ok=all(r['status']=='complete' for r in rr.values());original=old['matrix'][z][q]
                matrix[z][q]={'status':'complete' if ok else 'unavailable','xyz':original['xyz'],
                    'components':{'MACE_eV':original['components']['MACE_eV'],'GFN2_vacuum_hartree':rr['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':rr['alpb']['energy_hartree']} if ok else None,
                    'low':{s:r.get('actual') for s,r in rr.items()}}
        pool=choose_rows(matrix,list(CANDIDATES));oldpool=choose_rows(old['matrix'],list(CANDIDATES));calls={}
        for v in ('mathematical','operational'):
            val=pool[v]['composite_R_model_kcal_mol'] if pool['status']=='available' else None;call=decision(val,ref['variants'][v]['bands']);calls[v]={'decision':call,'outcome':outcome(call,group['expected_class'])}
        cases.append({**group,'status':pool['status'],'matrix':matrix,'pool':pool,'old_pool':oldpool,'calls':calls})
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'reference':m['reference'],'reference_branch':'fresh','rows':rows,'cases':cases,'cell_denominator':48,'complete_cells':sum(r['status']=='complete' for r in rows),
        'case_denominator':4,'complete_cases':sum(c['status']=='available' for c in cases),'new_molecular_calls_in_collection':0,'new_calibration':False,'production_changed':False}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','cases')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('precision225','reference','qualification','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        s=sub.add_parser(op)
        for k in fields:s.add_argument('--'+k,required=True,type=Path)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
