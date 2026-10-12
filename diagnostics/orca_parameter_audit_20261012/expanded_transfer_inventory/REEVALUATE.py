"""Finite archived-pool corrected reevaluation. Preparation never calls a model."""
from pathlib import Path
import argparse, csv, fcntl, hashlib, json, math, os, statistics, sys, time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
import context_factorial_gfn2 as f

ORDER=['origin','adaptive_Ca','adaptive_La']
EV=23.06054783061903
HA=627.509474

@lru_cache(None)
def cached(path,digest):
    assert f.pin(path)['sha256']==digest,('artifact changed',path)
    return f.read(path)

def checked(p):return cached(p['path'],p['sha256'])

def key(t):return (t['xyz']['sha256'],t['charge'],t['uhf'],t['medium'])

def choose(matrix):
    """Published choose_rows algebra/order; origin diagnostic is independent."""
    complete=all(matrix[z][q] is not None for z in ('Ca','La') for q in ORDER)
    origin_ok=all(matrix[z]['origin'] is not None for z in ('Ca','La'))
    chosen={'operational':None,'mathematical':None,'origin':None}
    for policy in chosen:
        if not (origin_ok if policy=='origin' else complete):continue
        mode={}
        for z in ('Ca','La'):
            zero=matrix[z]['origin']
            works={q:(v['MACE_eV']-zero['MACE_eV'])*EV+
                ((v['GFN2_ALPB_hartree']-zero['GFN2_ALPB_hartree'])-(v['GFN2_vacuum_hartree']-zero['GFN2_vacuum_hartree']))*HA
                for q,v in matrix[z].items() if v is not None}
            best=min(ORDER,key=lambda q:works[q]) if complete else 'origin'
            q='origin' if policy=='origin' else (best if policy=='mathematical' or works[best]<-.1 else 'origin')
            mode[z]={'candidate':q,'work_kcal_mol':works[q],'components':matrix[z][q]}
        ca,la=mode['Ca']['components'],mode['La']['components']
        mode['R_model_kcal_mol']=(ca['MACE_eV']-la['MACE_eV'])*EV+((ca['GFN2_ALPB_hartree']-ca['GFN2_vacuum_hartree'])-(la['GFN2_ALPB_hartree']-la['GFN2_vacuum_hartree']))*HA
        chosen[policy]=mode
    return {'status':'available' if complete else 'unavailable',**chosen}

SCIENTIFIC_SETTINGS=('accuracy','maxiter','temperature_K','temperature_Eh','population_parity_e',
    'charge_sum_tolerance_e','bohr_A','numerical_threads_per_worker','solvation','primary','independent','no_restart','retry_policy')

def verify_raw(r,t,m):
    """Bind source protocol, parsed science and immutable process output evidence."""
    assert key(t)==key(r['task'])
    for k in ('multiplicity','valence_electrons','atom_count','metal'):assert t[k]==r['task'][k]
    source=checked(r['manifest'])
    assert source['runtime']==m['runtime'] and r['implementation']==m['implementation']
    assert source['implementation']==m['implementation']
    assert all(source['settings'][k]==m['settings'][k] for k in SCIENTIFIC_SETTINGS)
    if not t.get('reuse'):assert r['manifest']==f.pin(Path(m['_base'])/'manifest.json')
    evidence={}
    for backend,stored in r['backends'].items():
        ep=stored['execution'];receipt=checked(ep);directory=Path(ep['path']).parent
        runtime=checked(m['runtime'])
        start=f.read(directory/'START.json');assert start['command']==receipt['command']
        assert all(v=='1' for v in start['thread_environment'].values())
        if backend=='xtb':
            expected=[runtime['xtb']['path'],'core.xyz','--gfn','2','--chrg',str(t['charge']),'--uhf','0','--etemp','300','--acc','0.0001','--iterations','500','--norestart','--json','--pop']
            if t['medium']=='water':expected+=['--alpb','water','gsolv']
        else:expected=[runtime['tblite_launcher_path'],m['implementation']['path'],'tblite-point','--task',str(directory.parent/'TASK.json')]
        assert receipt['command']==expected,('command mismatch',ep['path'])
        pins=[receipt['stdout'],receipt['stderr'],f.pin(directory/'START.json')]
        for pin in pins:assert f.pin(pin['path'])==pin
        assert f.pin(directory/'core.xyz')['sha256']==t['xyz']['sha256']
        assert f.read(directory.parent/'TASK.json')==r['task']
        parsed=f.parse_backend(directory,backend,r['task'])
        # Python's sum algorithm changed in 3.12; energies/populations remain exact.
        for field in ('charge_sum_e','charge_sum_residual_e'):
            if field in parsed and field in stored:
                assert abs(parsed[field]-stored[field])<=1e-10;parsed[field]=stored[field]
        if backend=='tblite' and parsed['status']=='available':
            a=parsed['solver_evidence'];b=stored['solver_evidence']
            assert abs(a['orbital_occupation_sum']-b['orbital_occupation_sum'])<=1e-10
            a['orbital_occupation_sum']=b['orbital_occupation_sum']
        assert parsed==stored,('raw parser mismatch',ep['path'])
        extra=['core.xyz']+(['charges'] if backend=='xtb' else ['scientific.json'])
        evidence[backend]={'execution':ep,'outputs':pins+[f.pin(directory/n) for n in extra if (directory/n).exists()]}
    return evidence

def prepare(config,output):
    cfg=f.read(config);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    reference=f.read(cfg['reference']);labels={r['case_id']:r for r in reference['rows']}
    inventory=f.read(cfg['transfer_inventory'])
    runtime_source=f.read(cfg['runtime_manifest']);runtime=runtime_source['runtime'];checked(runtime)
    reuse_collection=f.read(cfg['corrected_reuse']);reuse={}
    for row in reuse_collection['rows']:
        if not row['accepted']:continue
        for p in row['inputs'].values():
            r=checked(p);t=r['task']
            assert checked(r['manifest'])['runtime']==runtime
            assert all(checked(r['manifest'])['settings'][k]==runtime_source['settings'][k] for k in SCIENTIFIC_SETTINGS)
            reuse[key(t)]=p
    cases=[];tasks=[];model=None;algebra_errors=[]
    for path in cfg['collections']:
        collection=f.read(path);manifest=checked(collection['manifest'])
        if model is None:model=manifest['model']
        assert manifest['model']==model
        assert manifest['settings']['candidate_order']==ORDER
        assert manifest['settings']['origin_selection_tolerance_kcal_mol']==.1
        for old in collection['cases']:
            cid=old['case_id'];origin=old['source']['origin_row'];prep=origin['source'];source=prep['source']
            assert [q['id'] for q in old['candidates']]==ORDER
            lab=labels.get(cid);role=lab['role'] if lab else 'structural_robustness_development'
            expected=lab['expected_class'] if lab else source['expected_class']
            assert expected in ('Ca','La')
            if role=='canonical_calibration':assert source['canonical_coordinate_match']
            for p in [origin['source_preparation'],prep['state_signature'],*prep['maps'].values()]:checked(p)
            case={'case_id':cid,'source_case_id':origin['source_case_id'],'selection_id':origin['selection_id'],
                  'protein':lab['root_case_id'] if lab else source['root_case_id'],
                  'biological_group':source['biological_group'],'expected_class':expected,'role':role,
                  'source_conditioning_metal':source.get('source_conditioning_metal'),
                  'collection':f.pin(path),'state_signature':prep['state_signature'],'maps':prep['maps'],
                  'source_preparation':origin['source_preparation'],'matrix':{'Ca':{},'La':{}}}
            old_matrix={'Ca':{},'La':{}}
            paired={}
            for z in ('Ca','La'):
                for q in ORDER:
                    cell=old['matrix'][z][q];symbols,positions=f.coordinates(cell['xyz'])
                    assert symbols[0]==z and all(s in f.Z for s in symbols)
                    mr=checked(cell['MACE']);state=mr.get('state',mr.get('input_state_check'))
                    assert mr['status'] in ('computed','complete') and state['head']=='omol'
                    charge=int(state['charge']);assert charge==state['charge'] and state['spin_multiplicity']==1
                    assert state['atoms']==len(symbols) and math.isfinite(mr['energy_eV'])
                    if mr.get('request'):
                        request=checked(mr['request']);assert mr['model']==model
                    else:
                        mm=checked(mr['manifest']);assert mm['model']==model
                        request=next(t for t in mm['tasks'] if t['task_id']==mr['task_id'])
                    f.coordinates(request['xyz'])
                    assert request['xyz']['sha256']==cell['xyz']['sha256'] and request['charge']==charge
                    assert request.get('multiplicity',request.get('spin_multiplicity'))==1
                    electrons=sum(f.VALENCE[s] for s in symbols)-charge;assert electrons%2==0
                    paired[z,q]=(symbols,positions,charge)
                    low={}
                    for medium in ('vacuum','water'):
                        tid='__'.join((cid,z,q,medium))
                        t={'id':tid,'protein':case['protein'],'context':cid,'geometry':q,'metal':z,'medium':medium,
                           'xyz':cell['xyz'],'charge':charge,'multiplicity':1,'uhf':0,'valence_electrons':electrons,
                           'atom_count':len(symbols),'reference':None,'scout':False,
                           'stage':'canonical' if role=='canonical_calibration' else 'transfer'}
                        t['reuse']=reuse.get(key(t));tasks.append(t);low[medium]=tid
                    case['matrix'][z][q]={'xyz':cell['xyz'],'charge':charge,'multiplicity':1,
                        'MACE':cell['MACE'],'MACE_eV':mr['energy_eV'],'inputs':low}
                    old_matrix[z][q]=cell.get('components')
            for q in ORDER:
                ca,la=paired['Ca',q],paired['La',q]
                assert ca[0][1:]==la[0][1:] and ca[1]==la[1] and la[2]-ca[2]==1
            old_replay=choose(old_matrix)
            if old_replay['status']=='available':
                for mode in ('mathematical','operational'):
                    residual=old_replay[mode]['R_model_kcal_mol']-old['pool'][mode]['composite_R_model_kcal_mol']
                    assert abs(residual)<1e-6;algebra_errors.append(residual)
            cases.append(case)
    assert len(cases)==135 and len({c['case_id'] for c in cases})==135
    assert len(tasks)==1620 and sum(t['reuse'] is not None for t in tasks)==24
    assert Counter(c['role'] for c in cases)=={'canonical_calibration':25,'crystal_transfer':3,'separate_probe':6,'structural_robustness_development':101}
    assert len({key(t) for t in tasks})==1620,'additional exact duplicates require explicit alias plan'
    known_pairs={(c['source_case_id'],c['selection_id']):c['case_id'] for c in cases}
    triples=[]
    for t in inventory['triples']:
        triples.append({**t,'pool_case_ids':[known_pairs.get((s,t['selection_id'])) for s in t['members']]})
    m={'protocol_id':'Nikasha_expanded_ownthree_corrected_component_transfer_v1','runtime':runtime,
       'implementation':f.pin(f.__file__),'driver':f.pin(__file__),'config':f.pin(config),
       'plan':f.pin(Path(__file__).with_name('PLAN.md')),'historical_reference':f.pin(cfg['reference']),
       'historical_inventory':f.pin(cfg['transfer_inventory']),'corrected_reuse':f.pin(cfg['corrected_reuse']),
       'model':model,'tasks':tasks,'cases':cases,'triples':triples,
       'settings':runtime_source['settings']|{'component_parity_Eh':1e-6,'absolute_parity_is_diagnostic':True,
           'candidate_order':ORDER,'origin_selection_tolerance_kcal_mol':.1,'minimum_calibration_gap_kcal_mol':.02,
           'EV_TO_KCAL':EV,'HA_TO_KCAL':HA,'workers':'required_explicit_at_launch','maximum_workers':16},
       'historical_algebra_check':{'available_pools':len(algebra_errors)//2,'maximum_absolute_R_error_kcal_mol':max(map(abs,algebra_errors))},
       'scientific_calls_during_preparation':0,'production_changed':False}
    f.write(out/'manifest.json',m)
    dry_run(out)

def dry_run(base):
    m=f.read(base/'manifest.json');assert f.pin(__file__)==m['driver'] and f.pin(f.__file__)==m['implementation']
    runtime=checked(m['runtime'])
    for name in ('xtb','xtb_parameter_file','tblite_python'):
        assert f.pin(runtime[name]['path'])==runtime[name],('runtime changed',name)
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:
        assert f.pin(p['path'])==p,('runtime dependency changed',p['path'])
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']:assert f.pin(p['path'])==p,('dynamic dependency changed',p['path'])
    for t in m['tasks']:
        f.coordinates(t['xyz'])
        if t['reuse']:
            r=checked(t['reuse']);assert key(t)==key(r['task'])
            assert all(x['status']=='available' for x in r['backends'].values())
            verify_raw(r,t,m|{'_base':str(base)})
    summary={'manifest':f.pin(base/'manifest.json'),'status':'prepared_not_executed',
        'pools':len(m['cases']),'MACE_cells_reused':sum(len(c['matrix'][z]) for c in m['cases'] for z in ('Ca','La')),
        'inputs_per_backend':len(m['tasks']),'reused_inputs_per_backend':sum(t['reuse'] is not None for t in m['tasks']),
        'fresh_by_stage_per_backend':{s:sum(t['stage']==s and t['reuse'] is None for t in m['tasks']) for s in ('canonical','transfer')},
        'roles':dict(Counter(c['role'] for c in m['cases'])),'historical_algebra_check':m['historical_algebra_check'],
        'triples':len(m['triples']),'prior_unavailable_triples':sum(t['status']!='prepared' for t in m['triples'])}
    if not (base/'PREPARATION.json').exists():f.write(base/'PREPARATION.json',summary)
    print(json.dumps(summary,indent=2));return summary

@lru_cache(None)
def manifest_data(base):return f.read(Path(base)/'manifest.json')

def raw(base,t):
    if t['reuse']:r=checked(t['reuse']);p=t['reuse']
    else:
        path=base/'points'/t['id']/'RESULT.json'
        if not path.exists():return None,None
        r=f.read(path);p=f.pin(path)
    manifest=manifest_data(str(base))
    verify_raw(r,t,manifest|{'_base':str(base)})
    return r,p

def component(base,v,w):
    assert v['xyz']==w['xyz'] and v['charge']==w['charge'] and v['uhf']==w['uhf']
    vr,vp=raw(base,v);wr,wp=raw(base,w)
    row={'accepted':False,'inputs':{'vacuum':vp,'water':wp},'primary_vacuum_energy_hartree':None,
         'primary_water_energy_hartree':None,'primary_solvation_hartree':None,'component_error_hartree':None}
    if vr is None or wr is None:row['reason']='not_executed';return row
    if any(x['status']!='available' for r in (vr,wr) for x in r['backends'].values()):
        row['reason']='backend_or_state_failure';return row
    de={b:wr['backends'][b]['energy_hartree']-vr['backends'][b]['energy_hartree'] for b in ('xtb','tblite')}
    error=de['xtb']-de['tblite'];dq=max(r['parity']['maximum_population_difference_e'] for r in (vr,wr))
    accepted=abs(error)<=1e-6 and dq<=1e-5
    row.update(accepted=accepted,component_error_hartree=error,maximum_population_error_e=dq,
        verified_output_evidence={q:verify_raw(r,t,manifest_data(str(base))|{'_base':str(base)}) for q,r,t in [('vacuum',vr,v),('water',wr,w)]},
        primary_vacuum_energy_hartree=vr['backends']['xtb']['energy_hartree'],
        primary_water_energy_hartree=wr['backends']['xtb']['energy_hartree'],
        primary_solvation_hartree=de['xtb'] if accepted else None,
        absolute_energy_errors_hartree={q:r['parity']['energy_difference_hartree'] for q,r in [('vacuum',vr),('water',wr)]})
    if not accepted:row['reason']='component_or_population_parity_failed'
    return row

def collect(base,name):
    m=f.read(base/'manifest.json');tasks={t['id']:t for t in m['tasks']};rows=[]
    for c in m['cases']:
        matrix={'Ca':{},'La':{}};components={}
        for z in ('Ca','La'):
            for q,cell in c['matrix'][z].items():
                low=component(base,*[tasks[cell['inputs'][med]] for med in ('vacuum','water')]);components[z+'__'+q]=low
                matrix[z][q]={'MACE_eV':cell['MACE_eV'],'GFN2_vacuum_hartree':low['primary_vacuum_energy_hartree'],
                             'GFN2_ALPB_hartree':low['primary_water_energy_hartree']} if low['accepted'] else None
        rows.append({k:c[k] for k in ['case_id','source_case_id','selection_id','protein','biological_group','role','expected_class','source_conditioning_metal']}|{'components':components,'pool':choose(matrix)})
    canonical=[r for r in rows if r['role']=='canonical_calibration'];assert len(canonical)==25
    calibration={}
    for policy in ('operational','mathematical','origin'):
        available=[r for r in canonical if r['pool'][policy] is not None];bands=None;gap=None
        status='unavailable_calibration_member'
        if len(available)==25:
            extrema={z:(max if z=='Ca' else min)(r['pool'][policy]['R_model_kcal_mol'] for r in canonical if r['expected_class']==z) for z in ('Ca','La')}
            gap=extrema['La']-extrema['Ca'];status='available' if gap>.02 else 'unsupported_class_separation'
            if status=='available':bands={'Ca_max':extrema['Ca'],'La_min':extrema['La']}
        calibration[policy]={'status':status,'available':len(available),'gap_kcal_mol':gap,'bands':bands}
    def call(score,bands,expected):
        if score is None:return {'decision':None,'outcome':'unavailable'}
        if bands is None:return {'decision':None,'outcome':'calibration_unavailable'}
        decision='Ca' if score<=bands['Ca_max'] else ('La' if score>=bands['La_min'] else 'inconclusive')
        return {'decision':decision,'outcome':'inconclusive' if decision=='inconclusive' else ('correct' if decision==expected else 'wrong')}
    for r in rows:
        r['decisions']={p:call(r['pool'][p]['R_model_kcal_mol'] if r['pool'][p] is not None else None,calibration[p]['bands'],r['expected_class']) for p in calibration}
    lookup={r['case_id']:r for r in rows};triples=[]
    for t in m['triples']:
        members=[lookup.get(k) for k in t['pool_case_ids']]
        available=t['status']=='prepared' and all(r is not None and r['pool']['status']=='available' for r in members)
        available_by_policy={p:t['status']=='prepared' and all(r is not None and r['pool'][p] is not None for r in members) for p in calibration}
        medians={p:statistics.median(r['pool'][p]['R_model_kcal_mol'] for r in members) if available_by_policy[p] else None for p in calibration}
        expected=next((r['expected_class'] for r in members if r is not None),None)
        triples.append({'triple_id':t['triple_id'],'protein':t['protein_id'],'members':t['members'],'pool_case_ids':t['pool_case_ids'],
            'status':'available' if available else 'unavailable','available_by_policy':available_by_policy,'prior_preparation_status':t['status'],
            'median_R_kcal_mol':medians,'decisions':{p:call(medians[p],calibration[p]['bands'],expected) for p in calibration}})
    result={'protocol_id':m['protocol_id'],'manifest':f.pin(base/'manifest.json'),'calibration':calibration,'rows':rows,'triples':triples,
        'completed_inputs':sum(raw(base,t)[0] is not None for t in m['tasks']),'declared_inputs':len(m['tasks']),
        'available_pools':sum(r['pool']['status']=='available' for r in rows),'declared_pools':135,
        'production_changed':False,'interpretation':'consumed PQQ functional-class development; no affinity or independent-replica claim'}
    f.write(base/name,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','triples','manifest')},indent=2));return result

def pending_tasks(base,tasks):
    pending=[];completed=[];interrupted=[]
    for t in tasks:
        cell=base/'points'/t['id']
        if (cell/'RESULT.json').exists():raw(base,t);completed.append(t['id']);continue
        if not cell.exists():pending.append(t);continue
        # No automatic repeat of an input that may already have consumed compute.
        for start in cell.glob('*/START.json'):
            if (start.parent/'EXECUTION.json').exists():continue
            pid=f.read(start)['process_group']
            try:os.killpg(pid,0)
            except ProcessLookupError:pass
            else:raise RuntimeError(f'owned unfinished process group {pid} remains; inspect before resumption')
        interrupted.append(t['id'])
    return pending,completed,interrupted

def execute(base,stage,workers):
    assert 1<=workers<=16
    with (base/'EXECUTION.lock').open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        dry_run(base);m=f.read(base/'manifest.json');runtime=checked(m['runtime'])
        if stage=='transfer':
            prior=f.read(base/'CANONICAL_COLLECTION.json');assert prior['manifest']==f.pin(base/'manifest.json')
            assert prior['calibration']['operational']['status']=='available','canonical gate failed: transfer is not executable'
        tasks=[t for t in m['tasks'] if t['stage']==stage and t['reuse'] is None]
        pending,completed,interrupted=pending_tasks(base,tasks)
        attempt=1+len(list(base.glob(stage.upper()+'_EXECUTION_ATTEMPT_*.json')))
        name=f'{stage.upper()}_EXECUTION_ATTEMPT_{attempt:03d}'
        f.write(base/(name+'.json'),{'stage':stage,'workers':workers,'numerical_threads_per_worker':1,
            'pending_inputs':len(pending),'completed_inputs_reused':completed,'interrupted_unavailable_inputs':interrupted,
            'started_epoch':time.time(),'manifest':f.pin(base/'manifest.json')})
        failures=[]
        def run(t):
            try:f.execute_point(t,m,base,runtime)
            except Exception as exc:return {'id':t['id'],'reason':str(exc)}
            return None
        with ThreadPoolExecutor(max_workers=workers) as pool:
            failures=[r for r in pool.map(run,pending) if r is not None]
        f.write(base/(name+'_FINISH.json'),{'finished_epoch':time.time(),'exceptions':failures})
        output='CANONICAL_COLLECTION.json' if stage=='canonical' else 'FINAL_COLLECTION.json'
        if (base/output).exists():output=name+'_COLLECTION.json'
        result=collect(base,output)
        if failures:raise RuntimeError(f'{len(failures)} task exceptions; preserved in {name}_FINISH.json')
        return result

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('dry-run');q.add_argument('--directory',required=True)
    q=sub.add_parser('execute');q.add_argument('--directory',required=True);q.add_argument('--stage',choices=['canonical','transfer'],required=True);q.add_argument('--workers',type=int,required=True)
    q=sub.add_parser('collect');q.add_argument('--directory',required=True);q.add_argument('--output-name',required=True)
    a=p.parse_args()
    if a.op=='prepare':prepare(a.config,a.output)
    elif a.op=='dry-run':dry_run(Path(a.directory).resolve())
    elif a.op=='execute':execute(Path(a.directory).resolve(),a.stage,a.workers)
    else:collect(Path(a.directory).resolve(),a.output_name)
