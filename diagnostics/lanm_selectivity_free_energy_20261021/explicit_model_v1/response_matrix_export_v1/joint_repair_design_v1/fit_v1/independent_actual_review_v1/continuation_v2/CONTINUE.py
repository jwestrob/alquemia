#!/usr/bin/env python3
"""Complete only the unfinished saved-cache/raw phase of the pinned failed audit."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
CHECKS=0
VERIFIED={}

def require(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def verify(p):
    k=(p['path'],p['sha256'])
    if k not in VERIFIED:
        require(pin(p['path'])==p,'Changed authority: '+p['path']);VERIFIED[k]=True
    return Path(p['path'])

def read(p):return json.loads(verify(p).read_text())
def token(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def module(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def save(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

def derived(a,b,path,dim,arithmetic,differences):
    """Only explicitly selected derived metadata enters this arithmetic check."""
    if isinstance(a,dict):
        require(isinstance(b,dict) and set(a)==set(b),'Derived dictionary shape '+path)
        for k in a:derived(a[k],b[k],path+'/'+k,dim,arithmetic,differences)
    elif isinstance(a,list):
        require(isinstance(b,list) and len(a)==len(b),'Derived list shape '+path)
        for i,(x,y) in enumerate(zip(a,b)):derived(x,y,path+'/'+str(i),dim,arithmetic,differences)
    elif isinstance(a,(float,int)) and not isinstance(a,bool):
        require(isinstance(b,(float,int)) and not isinstance(b,bool) and math.isfinite(a) and math.isfinite(b),'Finite derived value '+path)
        bound=arithmetic(dim,a,b);delta=abs(a-b)
        require(delta<=bound,'Existing arithmetic allowance '+path)
        if delta:differences.append(dict(path=path,stored=a,replayed=b,absolute_difference=delta,existing_arithmetic_bound=bound))
    else:require(a==b,'Exact nonnumerical derived metadata '+path)

def compare_qualification(old,new,row,scout):
    require(set(old)==set(new),'Export report denominator')
    require(old['passed'] is True and new['passed'] is True,'Both original export qualifications pass')
    dim=3*row['atoms'];diff=[]
    for k in set(old)-{'checks','exports','polarization_energy'}:
        require(old[k]==new[k],'Exact auxiliary export field '+k)
    require(set(old['checks'])==set(new['checks']),'Exact original numerical gate inventory')
    for name,a in old['checks'].items():
        b=new['checks'][name]
        require(set(a)==set(b)=={'error','bound','passed'},'Original gate schema')
        require(a['passed'] is True and b['passed'] is True and 0<=a['error']<=a['bound'] and 0<=b['error']<=b['bound'],'Unchanged original gate '+name)
        derived(a,b,'checks/'+name,dim,scout.arithmetic,diff)
    require(len(old['exports'])==len(new['exports'])==2,'Exactly two raw exports')
    for i,(a,b) in enumerate(zip(old['exports'],new['exports'])):
        require(set(a)==set(b),'Export dictionary shape')
        # Native arrays, source hashes, site map, restoration and printed checks stay exact.
        require({k:v for k,v in a.items() if k not in ('spectrum','numerical_bounds')}=={k:v for k,v in b.items() if k not in ('spectrum','numerical_bounds')},'Exact native exported tokens and maps')
        verify(a['source'])
        for name in ('spectrum','numerical_bounds'):derived(a[name],b[name],f'exports/{i}/{name}',dim,scout.arithmetic,diff)
    a,b=old['polarization_energy'],new['polarization_energy']
    require(set(a)==set(b)=={'native_contraction','dense_contraction','printed','print_half_unit'},'Polarization diagnostic schema')
    require(a['printed']==b['printed'] and a['print_half_unit']==b['print_half_unit'],'Exact printed polarization authority')
    for k in ('native_contraction','dense_contraction'):derived(a[k],b[k],'polarization_energy/'+k,dim,scout.arithmetic,diff)
    return dict(passed=True,arithmetic_helper=scout.pin(scout.__file__),differences=diff,
                maximum_difference=max((x['absolute_difference'] for x in diff),default=0.),
                original_numerical_gates_unchanged=True,raw_native_arrays_and_tokens_exact=True)

def run(c,out):
    failed=read(c['failed_review']);cfg=read(c['original_configuration']);diag=read(c['diagnosis'])
    strict_failure=read(c['strict_continuation_failure']);all_diag=read(c['all_selected_diagnosis'])
    require(strict_failure['status']=='unavailable_terminal_joint_fit_review_continuation' and strict_failure['error']=="ValueError('Exact selected energy/moment/u-v/chi/force/raw replay aquo9_central_Dy__dimer_01')",'Preserved strict derived-chi equality failure')
    require(all_diag['failed_review']==c['strict_continuation_failure'] and len(all_diag['rows'])==46,'Complete selected-only diagnosis')
    require(failed['status']=='unavailable_terminal_joint_fit_review' and failed['error']=="ValueError('Selected native raw replay')" and failed['checks']==327825,'Exact preserved failure boundary')
    require(failed['configuration']==c['original_configuration'] and failed['implementation']==cfg['implementation'],'Unchanged original review/configuration')
    require(failed['fit_result']==cfg['fit_result']==diag['fit_result'] and diag['failed_review']==c['failed_review'],'Same terminal result and diagnosis')
    verify(cfg['implementation']);r=read(cfg['fit_result'])
    require(r['status']=='bounded_joint_response_fit_complete' and r['stage']=='run' and r['error'] is None and r['candidate_frozen'],'Terminal candidate only')
    dc=read(cfg['driver_config']);ec=read(dc['evaluator_config']);release=read(dc['root_release']);m=read(ec['manifest']);domain=read(ec['domain'])
    require(r['configuration']==cfg['driver_config'] and ec['implementation']==diag['evaluator']==dc['evaluator'],'Same frozen scientific sources')
    e=module(ec['implementation'],'continuation_evaluator');h={k:module(ec[k],'continuation_'+k) for k in e.HELPERS}
    ctx=SimpleNamespace(h=h,m=m,config=ec)
    selected=read(r['selected']);require(r['selected']==diag['selected'],'Same selected artifact')
    boundary=diag['first_difference']['cache_key'];cache=Path(release['cache_directory']);ledger=read(r['evaluator']['ledger'])
    require(boundary==c['boundary_cache_key'] and boundary==diag['checked_selected_rows'][-1]['cache_key'],'Exact interrupted cache key')
    require(r['limits']=={'outer_steps':24,'new_vectors':672,'native_programs':24000,'workers':12},'Original finite limits retained')
    # The completed driver/score/selection phases are inherited, never executed here.
    selected_rows={}
    for cid,s in sorted(selected['cells'].items()):
        row=m['rows'][cid]
        if s['inherited_native']:
            selected_rows[cid]=dict(status='inherited_from_completed_initial_baseline_audit',cache_key=s['cache_key']);continue
        cp=pin(cache/'cells'/s['cache_key']/'RESULT.json');cell=read(cp)
        require({k:v for k,v in s.items() if k!='request_provenance'}==cell,'Exact selected cache authority '+cid)
        parsed=e.Evaluator._parse(ctx,cell,row)
        derived_fields={'export_qualification','chi_a0_cubed','chi_solve_residual'}
        require({k:cell[k] for k in parsed if k not in derived_fields}=={k:v for k,v in parsed.items() if k not in derived_fields},'Exact selected energy/moment/u-v/force/raw replay '+cid)
        q=compare_qualification(cell['export_qualification'],parsed['export_qualification'],row,h['scout_helper']) if 'export_qualification' in parsed else None
        chi_differences=[]
        if q is not None:
            # Exact native arrays; only host-side derived scalar replay may round.
            for report in (cell['export_qualification'],parsed['export_qualification']):
                for xp in report['exports']:
                    lam=xp['spectrum']['A'];scaled=xp['spectrum']['normalized_A']
                    require(0<lam[0] and scaled[0]>1e-6*scaled[-1]>0,'Both stored/replayed original SPD criteria')
            dim=3*row['atoms'];A=np.asarray(parsed['export_qualification']['exports'][1]['A']);z=np.tile([0.,0.,1.],row['atoms']);x=np.linalg.solve(A,z)
            residual_bound=h['scout_helper'].arithmetic(dim,A,x)*max(1.,float(np.max(abs(A))),float(np.max(abs(x))))
            require(0<=cell['chi_solve_residual']<=residual_bound and 0<=parsed['chi_solve_residual']<=residual_bound,'Both original susceptibility solve residual gates')
            for k in ('chi_a0_cubed','chi_solve_residual'):derived(cell[k],parsed[k],k,dim,h['scout_helper'].arithmetic,chi_differences)
        selected_rows[cid]=dict(status='passed_selected_raw_replay',cache_key=s['cache_key'],result=cp,export_arithmetic=q,chi_replay_arithmetic=chi_differences,
            stored_objective_chi_unchanged=True,exact_fields=[k for k in parsed if k not in derived_fields],regression_energy_moment_subset_exact=True)
    require(len(selected_rows)==46,'All selected46 retained')
    save(out/'SELECTED_RAW_REPLAY.json',dict(status='passed_selected_raw_replay',selected=r['selected'],rows=selected_rows,new_native_calls=0))
    print('Selected raw46 complete; beginning unfinished unique-cache receipt phase.',flush=True)

    # Ledger keys were already joined to every vector/cell before the failure.
    required_native={};required_inherited=set();requests={x['vector_key']:x for x in ledger['requests']}
    for req in ledger['requests']:
        for cid,x in req['cells'].items():
            if x['inherited_native']:required_inherited.add(x['cache_key'])
            elif x['cache_key'] not in required_native:required_native[x['cache_key']]=dict(cell_id=cid,vector_key=req['vector_key'])
    files={p.parent.name:p for p in (cache/'cells').glob('*/RESULT.json')}
    require(set(files)==set(required_native) and boundary in files,'Complete unique native cache denominator; no omitted/extra cache')
    suffix={k for k in files if k>boundary};prefix=set(files)-suffix
    print(f'Cache inventory: {len(prefix)} previously audited receipts; {len(suffix)} remaining unique cells.',flush=True)
    # Load vector artifacts only to obtain one exact authority per unique cache.
    # No objective, model, response assembly, selection, Jacobian or oracle rerun.
    # The failed audit did not persist its prefix manifest, so bind that prefix
    # back to its immutable vector records without repeating its receipt checks.
    by_vector={}
    for ck in files:
        rec=required_native[ck];by_vector.setdefault(rec['vector_key'],[]).append(ck)
    signatures={}
    for j,(vk,cks) in enumerate(by_vector.items(),1):
        v=read(r['vectors'][vk]['artifact'])
        for ck in cks:
            s=v['cells'][required_native[ck]['cell_id']]
            require(s['cache_key']==ck,'Unique cache vector authority')
            signatures[ck]=token({k:x for k,x in s.items() if k!='request_provenance'})
        if j%100==0:print(f'Unique cache vector authorities {j}/{len(by_vector)}.',flush=True)
    rendered={}
    def render(values,cid,exports=None):
        row=m['rows'][cid];enabled=row['needs_export'] if exports is None else exports
        key=(row['source_input']['key']['sha256'],row['source_input']['parameters']['sha256'],row['metal'],tuple(e.effective(values,row['metal'],domain)),enabled)
        if key not in rendered:rendered[key]=e.render(values,row,domain,h['resolver'],enabled)
        return rendered[key]
    native=[];receipts=[];failed_cells=[]
    for j,ck in enumerate(sorted(files),1):
        path=files[ck];cp=pin(path);cell=read(cp);cid=cell['cell_id'];row=m['rows'][cid]
        require(token(cell)==signatures[ck],'Exact immutable vector/cache authority')
        inherited_receipt_checks=ck in prefix
        if not inherited_receipt_checks:
            values=cell['first_requested_vector'];key,meta=render(values,cid)
            ident=dict(cell_id=cid,effective_values=meta['physical_values_a_b_c_d_R_epsilon'],xyz=row['source_input']['xyz'],parameters=row['source_input']['parameters'],key_sha256=hashlib.sha256(key.encode()).hexdigest(),implementation=ec['implementation'],source_manifest=ec['source_manifest'],executables=ec['executables'],helpers={k:ec[k] for k in e.HELPERS})
            require(token(ident)==ck==cell['cache_key'] and cell['cache_identity']==ident and cell['render']==meta,'Remaining cache identity')
            require(token(values) in requests and cell['implementation']==ec['implementation'] and cell['root_release']==dc['root_release'],'Remaining actual producer/release')
            programs=['analyze']+(['testgrad'] if row['needs_gradient'] else [])
            require(len(cell['receipts'])<=len(programs) and set(cell['program_inputs'])<=set(programs),'Remaining finite program slots')
            for program,receipt in zip(programs,cell['receipts']):
                folder=path.parent/program;inp=cell['program_inputs'][program]
                for name,file in (('xyz','system.xyz'),('parameters','amoeba09.prm')):
                    require(verify(inp[name])==folder/file and inp[name]['sha256']==row['source_input'][name]['sha256'],'Remaining literal input')
                key_text,_=render(values,cid,row['needs_export'] if program=='analyze' else False)
                require(verify(inp['key'])==folder/'system.key' and verify(inp['key']).read_text()==key_text,'Remaining exact rendered key')
                argv=[ec['executables'][program]['path'],'system.xyz']+(['EM'] if program=='analyze' else ['Y','N'])
                require(receipt['argv']==argv and receipt['cwd']==str(folder),'Remaining native invocation')
                require(set(receipt['threads'])==set(THREADS) and all(receipt['threads'][k]=='1' for k in THREADS),'Remaining native threads')
                require(verify(receipt['output'])==folder/(program+'.out'),'Remaining output hash')
            if cell['status']=='complete':require(len(cell['receipts'])==len(programs) and all(q.get('returncode')==0 for q in cell['receipts']),'Remaining completed native denominator')
        native.append(dict(cache_key=ck,cell_id=cid,result=cp,status=cell['status'],error=cell.get('error'),programs=len(cell['receipts']),receipt_checks_inherited_from_failed_audit=inherited_receipt_checks))
        receipts.extend(cell['receipts'])
        if cell['status']!='complete':failed_cells.append(native[-1])
        if j%2000==0:print(f'Cache receipt inventory {j}/{len(files)}.',flush=True)
    require(len(receipts)==r['evaluator']['native_program_invocations']==ledger['native_program_invocations'],'Full native program accounting')
    require(not failed_cells and not ledger['stopped'] and not r['evaluator']['stopped'],'No hidden native failure')
    for rec in receipts:require(math.isfinite(rec['elapsed_seconds']) and rec['elapsed_seconds']>=0,'Finite receipt time')
    for v in r['cost'].values():require(isinstance(v,(float,int)) and math.isfinite(v) and v>=0,'Finite driver cost')
    initial=read(r['vectors'][r['initial_key']]['artifact'])
    trajectory=[dict(outer_step=it['outer_step'],origin_key=it['origin_key'],selected_key=it['selected_key'],objective_before=r['vectors'][it['origin_key']]['objective'],objective_after=r['vectors'][it['selected_key']]['objective'],actual_improvement=it['actual_improvement']) for it in r['iterations']]
    source_regression=read(c['regression_ready'])['implementation'];verify(source_regression)
    result=dict(status='passed_terminal_joint_fit_review_continuation',aggregate_audit_status='passed_using_preserved_completed_phases_and_narrow_continuation',
        original_audit_status=failed['status'],original_audit_passed=False,original_completed_check_count=failed['checks'],failed_review=c['failed_review'],
        original_configuration=c['original_configuration'],diagnosis=c['diagnosis'],all_selected_diagnosis=c['all_selected_diagnosis'],strict_continuation_failure=c['strict_continuation_failure'],fit_result=cfg['fit_result'],implementation=pin(__file__),
        driver_implementation=dc['implementation'],evaluator_implementation=ec['implementation'],selected=r['selected'],
        completed_phases_inherited=['terminal/config/source/release/domain/target joins','all673 vector identities, exact objective and response arithmetic, initial saved replay','all24 derivative/trial batches, BVLS/KKT saved arithmetic, strict trial/current selection, stopping reason','full evaluator ledger/vector/request joins','lexicographic cache receipt prefix through boundary (raw replay stopped at boundary)'],
        phases_not_repeated=['objective evaluation','EF checks','response assembly','Jacobian reconstruction','direction or optimizer calculations','trial selection'],
        boundary_cache_key=boundary,receipt_prefix_cells=len(prefix),receipt_continuation_cells=len(suffix),unique_native_cells=len(native),native_programs=len(receipts),native_cells=native,
        selected_raw_review=pin(out/'SELECTED_RAW_REPLAY.json'),selected_raw_cells=sum(x['status']=='passed_selected_raw_replay' for x in selected_rows.values()),
        all_vectors=len(r['vectors']),new_vectors=r['new_parameter_vectors'],completed_outer_steps=r['completed_outer_steps'],original_limits=r['limits'],
        initial_objective=initial['score']['objective'],initial_blocks=initial['score']['block_objectives'],selected_objective=r['selected_objective'],selected_blocks=r['selected_blocks'],
        original_EF_checks_passed=r['original_EF_checks_passed'],stopping_reason=r['stopping_reason'],bound_contacts=r['bound_contacts'],trajectory=trajectory,
        failed_cells=failed_cells,failed_vector_records=[x for x in r['vectors'].values() if x['status']!='complete'],
        regression_dependency=dict(source=source_regression,source_changed=False,all_selected_energy_moment_subsets_exact=True,derived_export_metadata_excluded_by_context=True,full_context_not_executed=True),
        cost=dict(driver=r['cost'],sum_native_elapsed_seconds=math.fsum(x['elapsed_seconds'] for x in receipts),failed_full_audit=failed['review_cost'],failed_strict_continuation=strict_failure['review_cost'],diagnostic=diag['cost'],all_selected_diagnostic=all_diag['cost'],scheduler_accounting_bound=False,
                  child_CPU_warning='Per-call concurrent child CPU deltas overlap; use driver child_CPU_seconds, not a receipt sum.'),
        new_native_calls=0,new_parameter_evaluations=0,new_optimization_calls=0,physical_accuracy_qualified=False,reserved_outcomes_opened=False,
        independence='Driver selection/objective normalization was independently reconstructed in the preserved audit. Evaluator/objective authored by this reviewer; their frozen pure parsers are reused, not presented as independent oracles.',
        limits=['Original unavailable audit remains failed; this receipt completes its unfinished phase.','Only replayed derived export diagnostics/chi/residual use the existing scout arithmetic allowance; stored objective chi and scores stay exact, as do energies, moments, u/v, native exports and all input tokens.','No physical tolerance, objective, fit protocol, source or candidate changed.','No fresh validation or consumed trimer-regression outcome was read.'])
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One local numerical thread')
    c=json.loads(a.config.read_text());require(c['implementation']==pin(__file__),'Continuation source pin');verify(c['declaration'])
    a.out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time()
    try:r=run(c,a.out)
    except Exception as exc:r=dict(status='unavailable_terminal_joint_fit_review_continuation',error=repr(exc),implementation=pin(__file__),failed_review=c['failed_review'],new_native_calls=0,new_parameter_evaluations=0)
    r.update(configuration=pin(a.config),checks=CHECKS,review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    save(a.out/'REVIEW.json',r)
    text='# Terminal fit audit continuation\n\n'+r['status']+'\n\n'
    if 'error' in r:text+=r['error']+'\n'
    else:text+=f"The original audit remains failed after {r['original_completed_check_count']} completed checks. The frozen terminal/config/source and its completed phases are reused. Only the unfinished cache receipts and selected raw rows are completed here.\n\n{r['all_vectors']} vectors; {r['completed_outer_steps']} outer steps; {r['native_programs']} native programs; {r['selected_raw_cells']} selected raw cells. Stop: {r['stopping_reason']}. Original25 E/F checks pass: {r['original_EF_checks_passed']}.\n\n"+'\n'.join('- '+x for x in r['limits'])+'\n'
    (a.out/'REPORT.md').write_text(text);print(json.dumps(dict(status=r['status'],review=pin(a.out/'REVIEW.json'))),flush=True)
    return 0 if r['status']=='passed_terminal_joint_fit_review_continuation' else 1

if __name__=='__main__':raise SystemExit(main())
