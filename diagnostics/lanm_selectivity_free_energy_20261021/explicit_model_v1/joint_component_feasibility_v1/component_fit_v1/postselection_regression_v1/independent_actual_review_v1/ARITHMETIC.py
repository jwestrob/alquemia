"""Direct saved-token audit of eight consumed cycles; no admission suite or native calls."""
import hashlib,json,math,os,re,sys,time,importlib.util,datetime
from types import SimpleNamespace
from decimal import Decimal,getcontext
from pathlib import Path
getcontext().prec=50
tick,cpu=time.monotonic(),time.process_time()
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[6]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/postselection_regression_v1'
PINS={}
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
    assert pin(p['path'])==p,p;PINS[p['path']]=p['sha256'];return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def txt(p):return verify(p).read_text()
def D(v):return Decimal(str(v).replace('D','E'))
def norm(v):return math.sqrt(math.fsum(float(x)**2 for x in v))
def arithmetic(*v):return 128*24*sys.float_info.epsilon*max([1.]+[abs(float(x)) for x in v])
for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:assert os.environ.get(key)=='1'
rp=pin(W/'run_v1/execution/RESULT.json');assert rp['sha256']=='d4685a09808f4a15ee32b29aea545ada882d8be3f8b53e1c85edb73cfd109db0'
r=read(rp);assert r['status']=='consumed_trimer_regression_complete' and r['native_program_invocations']==6 and r['error'] is None
verify(r['implementation']);m=read(r['manifest']);ec=read(r['configuration']['evaluator_config']);em=read(ec['manifest']);selected=read(r['selected'])
assert r['parameters']==selected['parameters'] and r['fit'] is False and r['reserved_outcomes_opened'] is False
assert len(r['rows'])==6 and len(r['cycles'])==8 and all(x['status']=='complete' and not x['inherited_native'] for x in r['rows'].values())
for cell in r['rows'].values():
    assert len(cell['receipts'])==1 and cell['receipts'][0]['returncode']==0
    assert cell['receipts'][0]['argv']==[ec['executables']['analyze']['path'],'system.xyz','EM']
    for p in cell['program_inputs']['analyze'].values():verify(p)
verify(ec['executables']['analyze'])
cfg=read(pin(W/'CONFIG_RUN.json'));assert cfg==r['configuration']
release=read(cfg['root_release']);source_review=read(cfg['source_review'])
assert release['status']=='released_component_postselection_regression' and release['maximum_native_calls']==6 and release['workers']==1 and release['reserved_outcomes_authorized'] is False
for key in ('implementation','manifest','fit_result','fit_review','evaluator_config','source_review','hostname','output_directory'):assert release[key]==cfg[key]
assert release['selected']==r['selected'] and r['implementation']==cfg['implementation']==m['implementation']==source_review['implementation']
assert source_review['status']=='passed_component_postselection_regression_source_review' and source_review['manifest']==r['manifest']
fit=read(r['fit_result']);fitreview=read(r['fit_review']);fitcfg=read(fit['configuration'])
assert fit['status']=='bounded_component_fit_complete' and fit['candidate_frozen'] and fit['error'] is None and fit['selected']==r['selected']
assert fitreview['status']=='passed_component_fit_actual_review' and fitreview['selected_native_admission_passed'] is True
assert fitreview['result']==r['fit_result'] and fitreview['selected']==r['selected'] and fitreview['configuration']==fit['configuration'] and fitreview['implementation']==fit['implementation']
assert fitcfg['parent_configuration']==cfg['evaluator_config'] and fitcfg['parent_evaluator']==ec['implementation'] and fitcfg['declaration']==m['campaign_declaration']
assert r['selected_metrics']==fit['selected_metrics'] and r['parameters']==fit['vectors'][fit['selected_key']]['parameters']
assert r['job_id'] is None and r['executor']['mode']=='login' and r['executor']['job_id'] is None and r['executor']['native_workers']==1
assert r['cost']['scheduler_allocation'] is None and r['cost']['scheduler_allocated_CPU_seconds']==0
previous=read(m['previous_regression']);previous_review=read(m['previous_regression_review'])
assert previous_review['status']=='passed_consumed_trimer_arithmetic_review' and previous_review['result']==m['previous_regression']
assert r['previous_selected_regression']==m['previous_regression'] and r['previous_selected_cycles']==previous['cycles']

# Qualified parser/render reuse is confined to six new cells; never invoke native.
def module(ref,name):
    spec=importlib.util.spec_from_file_location(name,verify(ref));obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj
def forbidden(*args,**kwargs):raise AssertionError('Saved review forbids molecular calls')
ev=module(ec['implementation'],'review_scientific_evaluator');h={k:module(ec[k],'regression_review_'+k) for k in ev.HELPERS}
h['native_helper'].invoke=forbidden
domain=read(ec['domain']);proxy=SimpleNamespace(h=h,m=em,config=ec)
qualification=read(ec['qualification_review']);build=read(ec['build_receipt'])
assert qualification['status']=='passed_local_response_native_qualification' and build['status']=='complete'
for key in ('source_manifest','build_receipt','executables'):assert qualification[key]==ec[key]
assert build['manifest']==ec['source_manifest'] and build['executables']==ec['executables']
for key in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed'):assert qualification[key] is True
last_end=None;native_replays={}
for cid,cell in r['rows'].items():
    row=m['rows'][cid];inp=cell['program_inputs']['analyze'];receipt=cell['receipts'][0]
    rendered,meta=ev.render(r['parameters'],row,domain,h['resolver'],exports=False)
    assert txt(inp['key'])==rendered and cell['render']==meta
    for key in ('xyz','parameters'):assert verify(inp[key]).read_bytes()==verify(row['source_input'][key]).read_bytes()
    assert receipt['executor']==cell['executor']==r['executor'] and receipt['timeout_seconds']==300 and all(v=='1' for v in receipt['threads'].values())
    directory=Path(receipt['cwd']);assert directory==Path(cfg['output_directory'])/'cells'/cid
    assert json.loads((directory/'RESULT.json').read_text())==cell
    start=datetime.datetime.fromisoformat(receipt['started_utc']);end=datetime.datetime.fromisoformat(receipt['ended_utc'])
    assert end>=start and (last_end is None or start>=last_end);last_end=end
    parsed=ev.Evaluator._parse(proxy,cell,row);assert all(cell[k]==v for k,v in parsed.items())
    native_replays[cid]=dict(status='passed_saved_native_input_output_replay',cell_result=pin(directory/'RESULT.json'),input=inp,output=receipt['output'])
targets=read(m['sources']['target_inventory']);old=read(targets['artifacts']['trimer_comparison'])
artifacts={'model_moments':read(targets['artifacts']['model_moments']),'trimer_model':read(old['inputs']['model_result'])}
baseline=read(m['sources']['baseline_result']);mp=read(m['sources']['moment_preparation'])
needed={k for c in m['cycles'].values() for k in c['model_terms']};values={};terms={};max_term_mu=0.

def raw_energy(mom):
    text=txt(mom['analyze_output']);v=re.findall(r'^\s*Total Potential Energy\s*:\s*([-+\d.EeDd]+)\s+Kcal/mole',text,re.M)
    assert len(v)==1 and float(D(v[0]))==mom['energy_kcal_mol'];return D(v[0]),text

for cid in sorted(needed):
    if cid in r['rows'] or cid in selected['cells']:
        cell=r['rows'].get(cid,selected['cells'].get(cid)); assert cell['status']=='complete' and not cell['inherited_native']
        row=m['rows'].get(cid,em['rows'].get(cid));mom=cell['moments'];energy,text=raw_energy(mom)
        n=row['atoms'];dual=[x.split()[2:] for x in text.splitlines() if x.startswith('LOCAL_RESPONSE_INDUCED')][-n:]
        assert len(dual)==n and all(len(x)==6 for x in dual)
        perm=[D(x) for x in row['permanent_lab_Debye']]
        mu=[perm[j]+D(em['units']['debye_per_eA'])*sum((D(x[j])+D(x[j+3]))/2 for x in dual) for j in range(3)]
        authority='six_new_trimers' if cid in r['rows'] else 'frozen_selected_artifact'
        stored=cell['mu_Debye'];assert float(energy)==cell['energy']
        source=mom['analyze_output']
    else:
        cell=em['static'][cid];src=cell['source']
        if 'artifact' in src:
            mom=artifacts[src['artifact']]
            for k in src['json_path']:mom=mom[k]
            perm=[D(x) for x in cell['permanent_lab_Debye']]
        else:
            assert src['result']==m['sources']['baseline_result'];mom=baseline['rows'][src['cell']]['moments']
            perm=[D(x) for x in mp['rows'][cid]['permanent_lab_Debye']]
        energy,text=raw_energy(mom);lines=[x.split() for x in txt(mom['induced_output']).splitlines() if x.strip()]
        n=int(lines[0][0]);assert len(lines)==n+1 and [int(x[0]) for x in lines[1:]]==list(range(1,n+1))
        mu=[perm[j]+sum(D(x[j+2]) for x in lines[1:]) for j in range(3)]
        assert float(energy)==cell['energy_kcal_mol'];stored=cell['mu_Debye'];authority='unchanged_saved_fragment';source=mom['analyze_output']
    discrepancy=max(abs(float(mu[j])-stored[j]) for j in range(3));max_term_mu=max(max_term_mu,discrepancy)
    assert discrepancy<=arithmetic(*mu,*stored)
    values[cid]=(energy,mu);terms[cid]=dict(authority=authority,raw_output=source,energy_decimal=str(energy),mu_Debye_decimal=[str(v) for v in mu],maximum_mu_arithmetic_difference=discrepancy)
cycles={};maxE=maxM=0.
for cid,c in m['cycles'].items():
    a=r['cycles'][cid];o=old['cycles'][cid]
    assert a['status']=='complete' and not a['missing'] and a['model_terms']==c['model_terms']==o['model_terms'] and len(c['model_terms'])==7 and sorted(c['model_terms'].values())==[-1,-1,-1,1,1,1,1]
    for key in ['reference_electronic_E3_kcal_mol','reference_mu3_Debye']:assert a[key]==c[key]==o[key]
    assert a['baseline_E3_kcal_mol']==c['model_E3_kcal_mol']==o['model_E3_kcal_mol'] and a['baseline_mu3_Debye']==c['model_mu3_Debye']==o['model_mu3_Debye']
    E=sum(D(co)*values[k][0] for k,co in c['model_terms'].items());M=[sum(D(co)*values[k][1][j] for k,co in c['model_terms'].items()) for j in range(3)]
    ed=abs(float(E)-a['candidate_E3_kcal_mol']);md=max(abs(float(M[j])-a['candidate_mu3_Debye'][j]) for j in range(3));maxE=max(maxE,ed);maxM=max(maxM,md)
    assert ed<=arithmetic(E,a['candidate_E3_kcal_mol']) and md<=arithmetic(*M,*a['candidate_mu3_Debye'])
    ref=D(c['reference_electronic_E3_kcal_mol']);bm=[D(x) for x in c['model_mu3_Debye']];qm=[D(x) for x in c['reference_mu3_Debye']]
    err=E-ref;be=D(c['model_E3_kcal_mol'])-ref;me=[M[j]-qm[j] for j in range(3)];bme=[bm[j]-qm[j] for j in range(3)]
    assert abs(float(err)-a['energy_error_kcal_mol'])<=arithmetic(err) and max(abs(float(me[j])-a['dipole_error_Debye'][j]) for j in range(3))<=arithmetic(*me)
    cycles[cid]=dict(candidate_E3_kcal_mol=float(E),candidate_E3_decimal=str(E),baseline_E3_kcal_mol=float(c['model_E3_kcal_mol']),reference_E3_kcal_mol=float(ref),candidate_energy_error=float(err),baseline_energy_error=float(be),absolute_energy_error_change=float(abs(err)-abs(be)),candidate_mu3_Debye=[float(x) for x in M],baseline_mu3_Debye=[float(x) for x in bm],reference_mu3_Debye=[float(x) for x in qm],candidate_mu_error_norm_Debye=norm(me),baseline_mu_error_norm_Debye=norm(bme),mu_error_norm_change_Debye=norm(me)-norm(bme),original_reference_print_bound_Debye_per_component=a['original_reference_print_bound_Debye_per_component'],original_model_print_bound_Debye_per_component=a['original_model_print_bound_Debye_per_component'],coefficients=c['model_terms'])
    prior=previous['cycles'][cid]
    for key in ('model_terms','reference_electronic_E3_kcal_mol','reference_mu3_Debye','baseline_E3_kcal_mol','baseline_mu3_Debye'):assert prior[key]==a[key]
    pe=D(prior['candidate_E3_kcal_mol'])-ref;pm=[D(x)-qm[j] for j,x in enumerate(prior['candidate_mu3_Debye'])]
    cycles[cid].update(previous_selected_E3_kcal_mol=prior['candidate_E3_kcal_mol'],previous_selected_mu3_Debye=prior['candidate_mu3_Debye'],previous_selected_energy_error=float(pe),previous_selected_mu_error_norm_Debye=norm(pm),absolute_energy_error_change_from_previous=float(abs(err)-abs(pe)),mu_error_norm_change_from_previous_Debye=norm(me)-norm(pm),candidate_mu_error_vector_Debye=[float(x) for x in me],previous_selected_mu_error_vector_Debye=[float(x) for x in pm],baseline_mu_error_vector_Debye=[float(x) for x in bme])
paired={}
for n in [8,9]:
    for context in ['near','far']:
        la=cycles[f'aquo{n}_central_La__{context}'];dy=cycles[f'aquo{n}_central_Dy__{context}'];p={}
        for kind in ['candidate','baseline','reference','previous_selected']:
            p[kind+'_E3_Dy_minus_La']=float(D(dy[kind+'_E3_kcal_mol'])-D(la[kind+'_E3_kcal_mol']))
            p[kind+'_mu3_Dy_minus_La_Debye']=[float(D(y)-D(x)) for x,y in zip(la[kind+'_mu3_Debye'],dy[kind+'_mu3_Debye'])]
        for kind in ['candidate','baseline','previous_selected']:
            p[kind+'_energy_error']=p[kind+'_E3_Dy_minus_La']-p['reference_E3_Dy_minus_La']
            p[kind+'_mu_error_norm_Debye']=norm([y-x for x,y in zip(p['reference_mu3_Dy_minus_La_Debye'],p[kind+'_mu3_Dy_minus_La_Debye'])])
        p['absolute_energy_error_change']=abs(p['candidate_energy_error'])-abs(p['baseline_energy_error']);p['mu_error_norm_change_Debye']=p['candidate_mu_error_norm_Debye']-p['baseline_mu_error_norm_Debye']
        p['absolute_energy_error_change_from_previous']=abs(p['candidate_energy_error'])-abs(p['previous_selected_energy_error'])
        p['mu_error_norm_change_from_previous_Debye']=p['candidate_mu_error_norm_Debye']-p['previous_selected_mu_error_norm_Debye']
        p['candidate_mu_error_vector_Debye']=[y-x for x,y in zip(p['reference_mu3_Dy_minus_La_Debye'],p['candidate_mu3_Dy_minus_La_Debye'])]
        paired[f'W{n}_{context}']=p
result=dict(status='passed_consumed_trimer_arithmetic_review',implementation=pin(__file__),result=rp,producer=r['implementation'],manifest=r['manifest'],selected=r['selected'],cycles=cycles,paired_Dy_minus_La=paired,raw_terms=terms,term_count=len(terms),cycle_count=8,paired_count=4,source_pins=[dict(path=k,sha256=v) for k,v in PINS.items()],maximum_energy_sum_replay_difference_kcal_mol=maxE,maximum_mu_sum_replay_difference_Debye=maxM,maximum_term_mu_arithmetic_difference_Debye=max_term_mu,energy_absolute_errors_worsened=sum(c['absolute_energy_error_change']>0 for c in cycles.values()),dipole_vector_norm_errors_improved=sum(c['mu_error_norm_change_Debye']<0 for c in cycles.values()),paired_energy_absolute_errors_worsened=sum(c['absolute_energy_error_change']>0 for c in paired.values()),new_native_calls=0,original_EF_checks_passed=r['original_EF_checks_passed'],scientific_threshold_added=False,physical_accuracy_qualified=False,limits=['Arithmetic review reuses admitted baseline/reference and permanent-frame values; it does not repeat native/QM admission or the full fit audit.','Eight consumed E3 cycles are outside the objective; their dependencies overlap objective cells, so this is not held-out validation.','No E3 acceptance threshold was declared; worsening is descriptive evidence for a progression decision.','PBE0 and frozen electronic-state limitations persist; no unique missing interaction or parameter mechanism is established.','Static fragments retain native six-decimal induced moments and the original printed-precision envelopes.'],cost=dict(saved_review_wall_seconds=time.monotonic()-tick,saved_review_CPU_seconds=time.process_time()-cpu,threads=1,producer_application=r['cost'],scheduler_allocation=None,scheduler_allocated_CPU_seconds=0))
out=W/'independent_actual_review_v1';out.mkdir(parents=True,exist_ok=False)
result.update(configuration=pin(W/'CONFIG_RUN.json'),root_release=cfg['root_release'],source_review=cfg['source_review'],fit_result=r['fit_result'],fit_review=r['fit_review'],previous_selected_regression=m['previous_regression'],previous_selected_review=m['previous_regression_review'],native_input_output_replay=native_replays,new_cell_native_admission_passed=True,reviewer_parent=pin(ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/regression_v1/independent_actual_review_v1/ARITHMETIC.py'),comparison_counts={
    'cycle_energy_abs_improved_from_previous':sum(v['absolute_energy_error_change_from_previous']<0 for v in cycles.values()),
    'cycle_mu_norm_improved_from_previous':sum(v['mu_error_norm_change_from_previous_Debye']<0 for v in cycles.values()),
    'paired_energy_abs_improved_from_previous':sum(v['absolute_energy_error_change_from_previous']<0 for v in paired.values()),
    'paired_mu_norm_improved_from_previous':sum(v['mu_error_norm_change_from_previous_Debye']<0 for v in paired.values()),
    'cycle_energy_abs_worse_than_original_baseline':sum(v['absolute_energy_error_change']>0 for v in cycles.values()),
    'cycle_mu_norm_worse_than_original_baseline':sum(v['mu_error_norm_change_Debye']>0 for v in cycles.values()),
    'paired_energy_abs_worse_than_original_baseline':sum(v['absolute_energy_error_change']>0 for v in paired.values()),
    'paired_mu_norm_worse_than_original_baseline':sum(v['mu_error_norm_change_Debye']>0 for v in paired.values())})
result['limits'][0]='Six new cells replay exact rendering/input bytes and the qualified native parser. Forty-two cycle terms replay raw tokens; admitted selected/reference/static authorities are reused without another full fit or quantum audit.'
lines=['Consumed development regression: eight E3/mu3 cycles; not holdout validation.',
       'Signed energy errors are model minus reference (kcal/mol); dipole errors below are Euclidean vector-error norms (Debye).',
       'Order: original baseline / prior selected joint fit / final component-fit selected.',
       '', 'Cycle | signed E3 error: baseline / prior / final | mu3 error norm: baseline / prior / final']
for cid,v in cycles.items():
    lines.append(f"{cid} | {v['baseline_energy_error']:.8f} / {v['previous_selected_energy_error']:.8f} / {v['candidate_energy_error']:.8f} | {v['baseline_mu_error_norm_Debye']:.8f} / {v['previous_selected_mu_error_norm_Debye']:.8f} / {v['candidate_mu_error_norm_Debye']:.8f}")
lines.extend(['','Paired Dy-minus-La | signed E3 error: baseline / prior / final | mu3 error norm: baseline / prior / final'])
for cid,v in paired.items():
    lines.append(f"{cid} | {v['baseline_energy_error']:.8f} / {v['previous_selected_energy_error']:.8f} / {v['candidate_energy_error']:.8f} | {v['baseline_mu_error_norm_Debye']:.8f} / {v['previous_selected_mu_error_norm_Debye']:.8f} / {v['candidate_mu_error_norm_Debye']:.8f}")
lines.extend(['','Comparison counts: '+json.dumps(result['comparison_counts'],sort_keys=True),
              'No physical E3/mu3 acceptance threshold was introduced. These cycles were consumed earlier and share objective dependencies.',
              'Original PBE0/frozen-electronic-state and print-precision limits remain; no unique mechanism or physical accuracy follows.',
              'Actual regression cost: '+json.dumps(r['cost'],sort_keys=True),
              'Saved review cost: '+json.dumps(result['cost'],sort_keys=True)])
(out/'REPORT.txt').write_text('\n'.join(lines)+'\n')
result['report']=pin(out/'REPORT.txt')
(out/'REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(review=pin(out/'REVIEW.json'),counts={k:result[k] for k in ['term_count','energy_absolute_errors_worsened','dipole_vector_norm_errors_improved','paired_energy_absolute_errors_worsened']},paired=paired,maximum_replay_differences=[maxE,maxM],cost=result['cost'])))
