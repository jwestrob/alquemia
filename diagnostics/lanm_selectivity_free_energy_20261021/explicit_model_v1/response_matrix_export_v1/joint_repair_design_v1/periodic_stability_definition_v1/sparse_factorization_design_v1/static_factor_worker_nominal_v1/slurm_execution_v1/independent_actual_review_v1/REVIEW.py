"""Saved-only review of job 20053; no numerical imports or backend calls."""
import datetime, hashlib, json, math, os, sys, time
from pathlib import Path

ROOT = Path.cwd().resolve()
REL = Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/periodic_stability_definition_v1/sparse_factorization_design_v1/static_factor_worker_nominal_v1/slurm_execution_v1')
D = ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL
W = ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL
OUT = W/'independent_actual_review_v1/actual_v1'
checks, verified = [], {}
def check(name, good):
    checks.append(dict(name=name, passed=bool(good)))
    if not good: raise AssertionError(name)
def pin(path):
    p = Path(path).resolve(); h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1<<20), b''): h.update(b)
    return dict(path=str(p), sha256=h.hexdigest())
def verify(p):
    if p['path'] not in verified:
        actual = pin(p['path']); check('pin '+p['path'], actual == p); verified[p['path']] = actual
    else: check('repeated pin '+p['path'], verified[p['path']] == p)
    return Path(p['path'])
def load(p): return json.loads(verify(p).read_text())
def direct_pins(obj):
    if isinstance(obj, dict):
        if set(obj) == {'path','sha256'}: verify(obj)
        else:
            for v in obj.values(): direct_pins(v)
    elif isinstance(obj, list):
        for v in obj: direct_pins(v)
def seconds(s):
    parts=s.split(':');return sum(float(v)*60**i for i,v in enumerate(reversed(parts)))

tick,cpu=time.monotonic(),time.process_time()
result=dict(status='failed_nominal_static_factor_actual_review', checks=checks, error=None,
            new_native_calls=0, matrix_assemblies=0, numeric_factorizations=0, solver_calls=0,
            numerical_backend_imports=0, nominal_only=True, e_app=None,
            application_error_bound_established=False, physical_accuracy_admitted=False, MD_admitted=False)
try:
    check('one thread', all(os.environ.get(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']))
    cfg=load(pin(D/'CONFIG.json')); rel=load(pin(W/'ROOT_RELEASE.json'))
    sup=load(pin(W/'actual_v1/RESULT.json')); wr=load(sup['worker_result'])
    start=load(pin(W/'actual_v1/START.json')); submission=load(pin(W/'submission_v1/SUBMISSION.json'))
    for obj in [cfg,rel,sup,wr,start,submission]: direct_pins(obj)
    inp=load(cfg['input_manifest']); ir=load(cfg['input_review']); sr=load(rel['source_review'])
    for obj in [inp,ir]: direct_pins(obj)
    assembly=load(ir['independent_assembly_review'])
    check('frozen source and config', cfg['implementation']['sha256']=='e73e9a788ca8ba310f3c55bc20d0c82faa6c3dc3f3f58e2da5ef15bea35d1139' and pin(D/'CONFIG.json')['sha256']=='4c4b741aa4b599cf91c814c8bd73ec3dee26af6e71d332c72ca2f1c9130d7fc3')
    check('release and source admission', rel['status']=='released_one_slurm_bound_static_factor_attempt' and sr['status']=='passed_slurm_nominal_static_factor_source_review' and sr['implementation']==cfg['implementation'] and sr['configuration']==pin(D/'CONFIG.json'))
    check('exact result/start/release joins', all(o['configuration']==pin(D/'CONFIG.json') and o['release']==pin(W/'ROOT_RELEASE.json') for o in [sup,wr,start]) and rel['implementation']==cfg['implementation']==sup['implementation']==wr['implementation'] and rel['resources']==cfg['resources'] and rel['output']==str(W/'actual_v1'))
    check('submission joins', submission['returncode']==0 and submission['job_id']=='20053' and submission['release']==sup['release'] and submission['argv'][-1]==rel['output'])
    check('terminal and no retry', sup['status']==wr['status']=='passed_nominal_static_sufficient_numerical_test' and wr['phase']=='finished' and sup['exit_code']==0 and not sup['wall_timeout'] and not sup['automatic_retry'] and wr['error'] is None)
    check('single symbolic numeric and ten solves', (wr['symbolic_attempts'],wr['numeric_attempts'],wr['solve_attempts'])==(1,1,10) and wr['new_native_calls']==sup['new_native_calls']==0)
    check('input lineage', wr['input_manifest']==cfg['input_manifest']==ir['input_manifest']==assembly['input_manifest'] and wr['input_review']==cfg['input_review'] and assembly['status']=='passed_candidate_sparse_assembly_actual_review')
    for k in ['candidate','composition_manifest','pme_export','nominal_contract','pme_result','pme_root_review','pme_bound_review']:
        check('candidate PME chain '+k, inp[k]==cfg[k]==ir[k])
    check('assembly source/result/arrays unchanged', inp['operator_producer']==ir['operator_producer']==assembly['producer_implementation'] and inp['operator_preparation_result']==ir['operator_preparation_result']==assembly['operator_preparation_result'] and inp['arrays']==assembly['arrays'])
    check('exact saved operator domain', inp['shape']==[35874,35874] and inp['nnz']==16371009 and inp['shift']==0 and inp['metal']=='La' and inp['cell_A']==[49.324]*3 and inp['cutoff_A']==9 and inp['ewald_inverse_A']==.4 and inp['all_alpha_positive'] and inp['mutual_scales_all_one'] and inp['metal_coordinate_indices']==[0,1,2])
    check('symmetry scope preserved', ir['normalized_symmetry_max_absolute']==0 and not ir['native_operator_parity_admitted'] and not assembly['native_operator_parity_admitted'])
    op=wr['operator']; unom=inp['U_B']+inp['U_R']; tau=1e-6*unom+inp['delta_odd']
    check('nominal bound and shift exact', all(op[k]==inp[k]==ir[k] for k in ['U_B','U_R','delta_odd']) and op['U_nominal']==unom and op['tau']==tau==.004903508115932018 and op['scope']=='nominal_saved_operator_only' and op['shape']==inp['shape'] and op['nnz']==inp['nnz'])
    old=load(pin(W.parent/'actual_v1/worker/RESULT.json'))
    check('same shifted matrix as preserved technical failure', old['operator']==op and old['numeric_attempts']==old['solve_attempts']==0 and old['status']=='failed_static_numeric_attempt')
    check('unknown application error retained', wr['nominal_only'] and wr['e_app'] is None and op['e_app'] is None and not wr['application_error_bound_established'] and not op['application_error_bound_established'] and not wr['physical_accuracy_admitted'] and not wr['MD_admitted'])
    check('clean backend', wr['backend_exception'] is None and wr['backend_warnings']==[] and wr['audit_warnings']==[])
    pre,post=wr['pre_numeric'],wr['post_numeric']; au=wr['audit']
    check('symbolic and numeric flags', pre['N']==35874 and not pre['is_numeric'] and pre['is_super'] and pre['is_ll'] and pre['is_lower'] and pre['order']=='amd' and post['N']==35874 and post['is_numeric'] and post['is_ll'] and not post['is_super'] and post['is_lower'] and post['order']=='amd')
    check('unexposed raw backend status retained', all(x['raw_minor'] is None and x['raw_status'] is None and x['unavailable_reason'] for x in [pre,post]))
    check('factor readout stable', au['before_copy']==au['after_copy']==post and au['original_factor_flags_preserved'] and au['finite_positive_diagonal'] and au['sparse_factor_shape']==[35874,35874] and au['sparse_factor_nnz']==399831681)
    check('complete permutation', len(au['permutation'])==35874 and sorted(au['permutation'])==list(range(35874)))
    bound=100*35874*sys.float_info.epsilon
    check('unchanged numerical bound', au['arithmetic_bound']==bound and math.isfinite(au['K_infinity_norm']) and au['K_infinity_norm']>0)
    names=[f'uniform_axis_{k}' for k in range(3)]+[f'metal_basis_{k}' for k in range(3)]+[f'rademacher_{k}' for k in range(4)]
    check('exact ten probe inventory', [r['name'] for r in au['probes']]==names)
    for p in au['probes']:
        good=all(math.isfinite(p[k]) and 0<=p[k]<=bound for k in ['factor_action_error','solve_backward_error']) and p['finite']
        check('saved probe '+p['name'], good and p['passed']==good)
    check('maxima and aggregate', au['maximum_factor_action_error']==max(p['factor_action_error'] for p in au['probes']) and au['maximum_solve_backward_error']==max(p['solve_backward_error'] for p in au['probes']) and au['passed'] and wr['sufficient_numerical_test_passed'])
    check('limits actual', wr['limits']['address_space_bytes']['actual']==[34359738368]*2 and wr['limits']['CPU_seconds']['actual']==[10800]*2 and wr['actual_threads']==dict(openblas=12,openmp=1) and 0<wr['cost']['wall_seconds']<cfg['resources']['wall_seconds'] and 0<wr['cost']['CPU_seconds']<10800 and wr['cost']['peak_RSS_KiB']*1024<34359738368 and sup['cost']['wall_seconds']<900)
    allocation=sup['allocation']
    check('exact allocation joins', start['allocation']==wr['allocation']==allocation and start['job_id']=='20053' and start['hostname']==allocation['hostname']=='node-64-768g-16')
    check('truthful scheduler resources', allocation['mode']=='slurm' and allocation['job_id']=='20053' and allocation['allocated_CPUs']==allocation['affinity_CPUs']==12 and allocation['nodes']==allocation['tasks']==1 and allocation['memory_grant_bytes']==48000*1024**2 and allocation['default_memory_MiB_per_CPU']==4000 and not allocation['requeue'] and not allocation['exclusive'])
    for k,v in [('account','banfield'),('partition','standard'),('qos','standard_expand')]: check('scheduler '+k, allocation[k]==v)
    fields='JobIDRaw,JobName,State,ExitCode,Start,End,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,UserCPU,SystemCPU,MaxRSS,MaxVMSize,ReqMem,AllocTRES,NodeList'.split(',')
    sacct=OUT/'SACCT.txt';rows=[dict(zip(fields,line.split('|'))) for line in sacct.read_text().splitlines()]
    check('terminal scheduler denominator', [r['JobIDRaw'] for r in rows]==['20053','20053.batch'] and all(r['State']=='COMPLETED' and r['ExitCode']=='0:0' for r in rows))
    job=rows[0];allocated=int(job['ElapsedRaw'])*int(job['AllocCPUS'])
    check('one allocation without batch double count', int(job['AllocCPUS'])==12 and allocated==int(job['CPUTimeRAW'])==948 and job['NodeList']==allocation['hostname'])
    check('scheduler measured CPU sum', abs(seconds(job['TotalCPU'])-seconds(job['UserCPU'])-seconds(job['SystemCPU']))<1e-10)
    result.update(status='passed_nominal_static_factor_actual_review', producer_implementation=cfg['implementation'], configuration=pin(D/'CONFIG.json'), result=pin(W/'actual_v1/RESULT.json'), worker_result=sup['worker_result'], root_release=sup['release'], input_manifest=cfg['input_manifest'], input_review=cfg['input_review'], candidate=cfg['candidate'], nominal_contract=cfg['nominal_contract'], source_review=rel['source_review'], scheduler_receipt=pin(sacct), submission=pin(W/'submission_v1/SUBMISSION.json'), operator=op, symbolic_attempts=1, numeric_attempts=1, solve_attempts=10, maximum_factor_action_error=au['maximum_factor_action_error'], maximum_solve_backward_error=au['maximum_solve_backward_error'], arithmetic_bound=bound, positive_finite_factor_diagonal=True, original_factor_flags_preserved=True, raw_backend_status_available=False, numerical_qualification_passed=True, native_operator_parity_admitted=False, normalized_symmetry_scope=ir['symmetry_scope'], worker_cost=wr['cost'], supervisor_cost=sup['cost'], scheduler=dict(rows=rows,allocation_elapsed_seconds=int(job['ElapsedRaw']),allocated_CPU_seconds=allocated,measured_CPU_seconds=seconds(job['TotalCPU']),maximum_batch_RSS_KiB=int(rows[1]['MaxRSS'].rstrip('K')),step_rows_not_added=True), preserved_original_failure=pin(W.parent/'actual_v1/worker/RESULT.json'), verified_pins=list(verified.values()), limitation='Independent saved metadata/hash/scalar-gate review; no repeated factor action or solve. Sampled action/backward residuals do not establish a whole-factor error norm or floating PME application bound. Static nominal pass is neither trajectory nor energy/force nor physical accuracy admission.')
except Exception as exc:
    result['error']=dict(type=type(exc).__name__,message=str(exc))
finally:
    result['implementation']=pin(__file__);result['check_count']=len(checks)
    result['review_cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu)
    result['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (OUT/'REVIEW.json').open('x') as f: json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:result[k] for k in ['status','check_count','error','review_cost']}))
    print(json.dumps(pin(OUT/'REVIEW.json')))
sys.exit(0 if result['error'] is None else 1)
