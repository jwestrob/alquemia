"""Saved arrays/metadata only; no assembly, native or factor API."""
from pathlib import Path
import hashlib,json,math,os,sys,time,traceback
import numpy as np

HERE=Path(__file__).resolve().parent
tick,cpu=time.monotonic(),time.process_time();checks=[];pins={}

def pin(path):
    path=Path(path).resolve();h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return dict(path=str(path),sha256=h.hexdigest())

def check(name,ok):
    checks.append(dict(name=name,passed=bool(ok)))
    if not ok:raise ValueError(name)

def read(p):
    check('pin '+p['path'],pin(p['path'])==p);pins[p['path']]=p
    return json.loads(Path(p['path']).read_text())

def close(a,b):return abs(a-b)<=1e-14*max(1.,abs(a),abs(b))

out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
result=dict(status='started',implementation=pin(__file__),checklist=pin(HERE/'CHECKLIST.md'),checks=checks,
            new_native_calls=0,matrix_assemblies=0,geometric_weight_evaluations=0,upper_bound_recomputations=0,
            numeric_factorizations=0,solver_calls=0,physical_accuracy_admitted=False,MD_admitted=False,error=None)
try:
    check('one numerical thread',all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')))
    ready=read(pin(HERE.parent/'READY.json'));cfg=read(ready['configuration']);inputs=read(ready['inputs'])
    prep=read(ready['preparation']);meta=read(ready['metadata']);source_review=read(pin(HERE.parent/'independent_source_review_v1/REVIEW.json'))
    actual_dir=Path(ready['output_manifest']).parents[1];rp=pin(actual_dir/'RESULT.json');sup=read(rp)
    worker=read(sup['worker_result']);manifest=read(sup['input_manifest']);release=read(sup['release'])
    start=read(pin(actual_dir/'START.json'))
    check('source/config/release joins',sup['implementation']==worker['implementation']==release['implementation']==ready['implementation'] and sup['configuration']==worker['configuration']==release['configuration']==ready['configuration'] and worker['release']==sup['release'])
    check('source review/released output',release['source_review']==pin(HERE.parent/'independent_source_review_v1/REVIEW.json') and source_review['status']=='passed_candidate_sparse_assembly_source_review' and release['status']=='released_one_terminal_candidate_La_real_assembly' and Path(release['output']).resolve()==actual_dir.resolve())
    check('terminal complete once',sup['status']==worker['status']=='completed_terminal_candidate_La_real_assembly' and sup['exit_code']==0 and sup['wall_timeout'] is False and sup['automatic_retry'] is False and worker['error'] is None and worker['assembly_attempts']==worker['actual_assembly_calls']==sup['actual_assembly_calls_reported']==1)
    check('no backend/native calls',sup['native_calls']==sup['numeric_factorizations']==worker['native_calls']==worker['numeric_factorizations']==worker['symbolic_analyses']==0)
    check('immutable worker manifest joins',manifest['operator_producer']==ready['implementation'] and manifest['operator_preparation_result']==sup['worker_result'] and manifest['arrays']==worker['arrays'] and manifest['source_preparation']==ready['preparation'])
    changed={'status','U_B','arrays','metal_coordinate_indices','operator_producer','operator_preparation_result'}
    check('all unchanged template fields',all(manifest[k]==v for k,v in inputs['template'].items() if k not in changed))
    check('correct manifest and matrix status',manifest['status']=='prepared_normalized_real_operator_for_single_static_factor' and manifest['matrix']=='normalized_B0_plus_H_unshifted_lower' and manifest['shift']==0.)
    h=worker['metal_stiffness'];check('stiffness manifest/preparation identity',manifest['metal_stiffness']==h and h['coefficients_a_b_c_d']==meta['response_coefficients'] and h['metal_site_index']==meta['metal_site_index'] and h['intrinsic_alpha']==meta['type_parameters']['900']['alpha'])
    a,b,c,d=h['coefficients_a_b_c_d'];s=h['s'];q=h['cross_sum']
    check('saved scalar environment relations',all(math.isfinite(v) for v in (s,q,h['eta'],h['g'],h['sigma'],h['normalized_beta'],h['native_beta'])) and s>=0 and q>=0 and close(h['eta'],q/(1+s*s)) and close(h['g'],2*s/(1+s)) and close(h['sigma'],h['g']*(c+d*h['eta'])))
    check('saved scalar normalized stiffness',close(h['normalized_beta'],math.exp(h['sigma'])-1.) and close(h['native_beta']*h['intrinsic_alpha'],h['normalized_beta']) and h['normalized_beta']>-1)
    wp=worker['water_weights'];check('weight file pin',pin(wp['path'])==wp);weights=np.load(wp['path'],mmap_mode='r',allow_pickle=False)
    check('saved weights finite/support count',weights.dtype==np.float64 and weights.shape==(3980,) and np.isfinite(weights).all() and np.all(weights>=0) and int(np.count_nonzero(weights))==h['supported_water_count'] and h['water_count']==3980)
    info=worker['assembly'];check('candidate U_B/graph/normalization metadata',manifest['U_B']==info['unshifted_site_upper'] and math.isfinite(manifest['U_B']) and info['sites']==11958 and info['dimension']==35874 and info['pairs']==1815015 and info['lower_entries']==16371009 and info['normalization']=='sqrt(alpha) D sqrt(alpha)' and info['triangle']=='lower' and info['shift']==0 and info['PME_included'] is False)
    arrays={}
    for k in ('data','indices','indptr'):
        p=manifest['arrays'][k];check('array file '+k,pin(p['path'])==p);arrays[k]=np.load(p['path'],mmap_mode='r',allow_pickle=False)
        digest=hashlib.sha256(memoryview(arrays[k]).cast('B')).hexdigest()
        check('array payload '+k,digest==worker['payload_hashes'][k] and (k=='data' or digest==manifest[k+'_payload_sha256']))
    data,idx,ptr=(arrays[k] for k in ('data','indices','indptr'));n=35874
    check('array dimensions/dtypes',data.dtype==np.float64 and data.shape==(16371009,) and idx.dtype==np.int32 and idx.shape==data.shape and ptr.dtype==np.int32 and ptr.shape==(n+1,))
    check('finite data and explicit diagonals',np.isfinite(data).all() and ptr[0]==0 and ptr[-1]==len(data) and np.all(np.diff(ptr)>0) and np.array_equal(idx[ptr[:-1]],np.arange(n)))
    # Exact graph hashes are the previously admitted full canonical/lower structure authority.
    graph=read(inputs['graph_result']);check('original qualified graph',all(worker['payload_hashes'][k]==graph['assembly'][k+'_sha256'] for k in ('indices','indptr')))
    xyz=Path(manifest['geometry']['path']);check('geometry pin',pin(xyz)==manifest['geometry']);types=[int(line.split()[5]) for line in xyz.read_text().splitlines()[2:]]
    alpha=np.array([meta['type_parameters'][str(t)]['alpha'] for t in types]);expected=1-alpha*(4*.4**3/(3*math.sqrt(math.pi)));expected[h['metal_site_index']]+=h['normalized_beta'];expected=np.repeat(expected,3)
    diagonal=data[ptr[:-1]];diag_error=float(np.max(np.abs(diagonal-expected)))
    check('normalized self plus candidate-H diagonals',diag_error<=1e-14 and manifest['metal_coordinate_indices']==worker['metal_coordinate_indices']==list(range(3*h['metal_site_index'],3*h['metal_site_index']+3)) and manifest['U_B']>=float(np.max(diagonal)))
    caps=cfg['resources'];cost=sup['cost'];wc=worker['cost']
    check('actual resource limits',release['resources']==caps and worker['actual_limits']=={'address_space_bytes':[caps['address_space_bytes']]*2,'CPU_seconds':[caps['CPU_seconds']]*2} and start['hostname']==caps['hostname'] and start['scheduler_job_id'] is None and start['configuration']==ready['configuration'] and start['release']==sup['release'])
    check('recorded time/memory within limits',cost['numerical_threads']==wc['numerical_threads']==1 and 0<wc['wall_seconds']<=cost['wall_seconds']<caps['wall_seconds'] and 0<wc['CPU_seconds']<=cost['child_user_CPU_seconds']+cost['child_system_CPU_seconds']<caps['CPU_seconds'] and cost['child_max_RSS_KiB']*1024<caps['address_space_bytes'] and info['estimated_assembly_bytes']<caps['assembly_bytes'])
    result.update(status='passed_candidate_sparse_assembly_actual_review',result=rp,input_manifest=sup['input_manifest'],operator_preparation_result=sup['worker_result'],producer_implementation=ready['implementation'],configuration=ready['configuration'],root_release=sup['release'],source_review=pin(HERE.parent/'independent_source_review_v1/REVIEW.json'),candidate=manifest['candidate'],composition_manifest=manifest['composition_manifest'],pme_export=manifest['pme_export'],pme_root_review=manifest['pme_root_review'],geometry=manifest['geometry'],arrays=manifest['arrays'],array_payload_hashes=worker['payload_hashes'],U_B=manifest['U_B'],U_R=manifest['U_R'],delta_odd=manifest['delta_odd'],e_app=None,application_error_bound_established=False,metal_stiffness=h,diagonal_maximum_replay_error=diag_error,normalized_symmetry_max_absolute=0.,normalized_symmetry_scope='Exact symmetric completion of the saved lower CSC by definition; no separately stored upper triangle or full native operator parity measured.',native_operator_parity_admitted=False,upper_bound_scope='Returned by exact qualified assembler; derivation source-reviewed, full bound not recomputed.',actual_cost=cost,worker_cost=wc)
except Exception as exc:
    result.update(status='failed_candidate_sparse_assembly_actual_review',error=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
result['review_cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1)
(out/'REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(status=result['status'],review=pin(out/'REVIEW.json'),checks=len(checks),error=result['error'],cost=result['review_cost'])))
raise SystemExit(0 if result['error'] is None else 1)
