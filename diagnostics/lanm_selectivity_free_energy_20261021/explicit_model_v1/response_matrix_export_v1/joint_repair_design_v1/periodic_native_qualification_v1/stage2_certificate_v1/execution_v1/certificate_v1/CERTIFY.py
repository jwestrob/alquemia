"""One changed-state assembly and nominal factor. No native program calls."""
import argparse
import ctypes
import math
import warnings
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from COMMON import *

N=35874
HERE=Path(__file__).resolve().parent
CAPS={'assembly':{'wall':180,'AS':8*1024**3,'CPU':180,'BLAS':1},
      'factor':{'wall':900,'AS':32*1024**3,'CPU':10800,'BLAS':12}}

def context(config,release,out,phase=None):
    cfg,m,rel,alloc=source_context(config,release,out,__file__,'certificate',12,12 if phase=='factor' else 1)
    require(cfg['state_id']!='central' and m['rows'][cfg['state_id']]['role']=='ACTIVE','only original32 changed states')
    require(rel['state_id']==cfg['state_id'] and rel['maximum_assemblies']==rel['maximum_numeric_factors']==1 and rel['maximum_solve_probes']==10,'exact one-state finite release')
    require(cfg['caps']==CAPS and cfg['nominal_contract']==m['nominal_contract'],'fixed resources/nominal contract')
    review=read(rel['source_review'])
    require(review['certificate_sources']==cfg['certificate_sources'],'separately reviewed producer helpers')
    for p in cfg['certificate_sources'].values():verify(p)
    require(cfg['certificate_sources']['CERTIFY.py']==pin(__file__),'exact certificate producer')
    for p in cfg['authorities'].values():verify(p)
    initializer=initialization(cfg,m)
    require(initializer['implementation']==m['execution_sources']['NATIVE.py'] and initializer['source_manifest']==m['native_sources'],'initializer exact execution source')
    require(initializer['source_inputs']==m['rows'][cfg['state_id']]['source_inputs'],'initializer original literal input')
    command=initializer['commands']['native'];require(command['returncode']==0 and initializer['native_attempts']==1,'exact successful initializer attempt')
    verify(command['stdout']);verify(command['stderr'])
    return cfg,m,rel,alloc,initializer

def phase_record(config,release,out,phase):
    directory=out/phase;directory.mkdir(exist_ok=False)
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),'phase':phase,'native_calls':0,'assembly_attempts':0,'symbolic_attempts':0,'numeric_attempts':0,'solve_attempts':0,'sufficient_numerical_test_passed':False,'error':None}
    return directory,result

def physical_metadata(row,m,native_reader):
    expected,cell=native_reader.sites(row['inputs']);meta=read(m['physics_metadata'])
    require(meta['all_alpha_positive'] is True and meta['metal_site_index']==0 and meta['sites']==11958 and meta['water_count']==3980,'frozen active map')
    require(meta['mutual_scales']=={'mutual-'+s+'-scale':1. for s in ('11','12','13','14')},'all-unity mutual scaling')
    coefficients={};scales={};pairs={};polar={}
    for key in ('parameters','key'):
        seen=set()
        for line in verify(row['inputs'][key]).read_text().splitlines():
            v=line.split('#',1)[0].split()
            if not v:continue
            tag=v[0].lower()
            if tag=='polar-local-response':
                types=tuple(map(int,v[1:3]));require((tag,types) not in seen,'duplicate response');seen.add((tag,types))
                require(types in ((900,402),(901,402)) and len(v)==7,'exact response roles');coefficients[str(types[0])]=list(map(float,v[3:]))
            elif tag=='polpair':
                types=tuple(sorted(map(int,v[1:3])));require((tag,types) not in seen,'duplicate pair');seen.add((tag,types));pairs[types]=float(v[3])
            elif tag.startswith('mutual-') and tag.endswith('-scale'):
                require(tag not in seen,'duplicate scale');seen.add(tag);scales[tag]=float(v[1])
            elif tag=='polarize' and int(v[1]) in (402,403,405,414,900,901):
                ty=v[1];require((tag,ty) not in seen,'duplicate polarize');seen.add((tag,ty));polar[ty]={'alpha':float(v[2]),'thole':float(v[3])}
            elif tag=='polar-donor-response':raise ValueError('donor chemistry outside fixed panel')
    require(polar==meta['type_parameters'] and scales==meta['mutual_scales'],'exact intrinsic alpha/damping/mutual scaling')
    require(pairs=={tuple(x['types']):x['gamma'] for x in meta['pair_overrides']},'exact pair damping')
    comp=read(m['composition'])
    require(set(coefficients)=={'900','901'},'both immutable coefficient slots')
    for metal,ty in [('La','900'),('Dy','901')]:require(coefficients[ty]==comp['effective_a_b_c_d_R_epsilon_POLPAIR'][metal][:4],'exact immutable candidate coefficients')
    require(expected[0]['type']==900 and sum(x['type']==900 for x in expected)==1 and sum(x['type']==402 for x in expected)==3980,'La/oxygen active map')
    return expected,cell,meta,coefficients,pairs

def assemble(cfg,m,initializer,directory,result,persist):
    import numpy as np
    native=load(m['execution_sources']['NATIVE.py'],'stage2_metadata_reader')
    expected,cell,meta,coefficients,pairs=physical_metadata(m['rows'][cfg['state_id']],m,native)
    snapshot_reader=load(m['helpers']['snapshot'],'stage2_snapshot_reader')
    with verify(initializer['pme_export']).open() as stream:snapshot=snapshot_reader.read(stream)
    state=load(cfg['certificate_sources']['STATE.py'],'stage2_state_correspondence')
    h,beta,direct,pdamp,thole,gamma=state.read_current(verify(initializer['state_output']),snapshot,expected,coefficients)
    alpha=np.asarray([x['alpha'] for x in expected]);types=np.asarray([x['type'] for x in expected]);xyz=np.asarray([x['xyz'] for x in expected])
    checks=[]
    state.compare(pdamp,alpha**(1./6.),'intrinsic damping radius',checks)
    require(np.array_equal(thole,np.asarray([meta['type_parameters'][str(t)]['thole'] for t in types])),'native Thole exact')
    require(np.all(gamma[types==402]==pairs[(402,900)]),'native pair override exact')
    expected_h,expected_beta=state.compact_response(expected,coefficients['900'],gamma,pdamp,cell[:3])
    state.compare(h,expected_h,'native shared h at literal cell',checks);state.compare(beta,expected_beta,'native H at literal cell',checks)
    require(np.all(beta[1:]==0.),'only metal H support')
    response=load(m['helpers']['response_environment'],'stage2_ordered_response');weights=[]
    for i,ty in enumerate(types):
        if ty!=402:continue
        vector=[]
        for axis,L in enumerate(cell[:3]):
            d=float(xyz[i,axis]-xyz[0,axis]);nearest=math.copysign(math.floor(abs(d/L)+.5),d/L);vector.append(d-L*nearest)
        weights.append(response.weight(math.sqrt(sum(v*v for v in vector))))
    _,_,c,d=coefficients['900'];stiffness=response.environment(weights,c,d)
    require(len(weights)==3980 and all(math.isfinite(x) for x in stiffness.values()),'full finite collective environment')
    state.compare(np.array([beta[0]*alpha[0]]),np.array([stiffness['normalized_beta']]),'normalized collective metal H',checks)
    result.update(state_checks=checks,metal_stiffness=stiffness,expected_state=initializer['state_output'],pme_export=initializer['pme_export'],cell_A=cell[:3],active_site_count=11958)
    # Snapshot arrays are no longer needed during the large pair query.
    del snapshot,h,beta,direct,expected_h,expected_beta
    asm=load(m['helpers']['assembly'],'stage2_qualified_assembly')
    result['assembly_attempts']=1;result['phase']='pair_query_and_assembly';persist()
    lower,info=asm.assemble(xyz,types,alpha,thole,cell[:3],cutoff=9.,ewald=.4,pair_overrides=pairs,stiffness={0:stiffness['normalized_beta']},shift=0.,max_assembly_bytes=4*1024**3,batch_pairs=4096,mutual_scales_all_one=True)
    require(lower.shape==(N,N) and lower.nnz>0 and str(lower.data.dtype)=='float64' and str(lower.indices.dtype)==str(lower.indptr.dtype)=='int32','canonical lower shape/dtypes')
    require(lower.has_canonical_format and np.all(np.isfinite(lower.data)) and math.isfinite(info['unshifted_site_upper']),'finite canonical assembly')
    arrays={};payload={}
    for name in ('data','indices','indptr'):
        value=getattr(lower,name);np.save(directory/(name+'.npy'),value,allow_pickle=False);arrays[name]=pin(directory/(name+'.npy'))
        payload[name]=hashlib.sha256(memoryview(np.ascontiguousarray(value)).cast('B')).hexdigest()
    require(np.array_equal(lower.indices[lower.indptr[:-1]],np.arange(N)),'explicit diagonal first')
    norms=initializer['PME_checks']['operator_norms']
    result.update(status='stage2_assembly_complete',phase='complete',arrays=arrays,payload_sha256=payload,shape=[N,N],nnz=int(lower.nnz),info=info,U_B=float(info['unshifted_site_upper']),U_R=float(norms['U_R_nominal']),delta_odd=float(norms['evenness_perturbation_nominal']),e_app=None,application_error_bound_established=False,metal_coordinate_indices=[0,1,2])

def factor(cfg,m,initializer,directory,result,persist,out):
    import numpy as np
    from scipy.sparse import csc_array
    assembly_pin=pin(out/'assembly/RESULT.json');a=read(assembly_pin)
    require(a['status']=='stage2_assembly_complete' and a['state_id']==cfg['state_id'] and a['initializer']==cfg['initializer'] and a['configuration']==pin(result['configuration']['path']),'matching completed assembly')
    require(a['shape']==[N,N] and a['expected_state']==initializer['state_output'] and a['pme_export']==initializer['pme_export'],'assembly exact native state')
    blas=ctypes.CDLL(str(Path(sys.prefix)/'lib/libopenblasp-r0.3.34.so'));omp=ctypes.CDLL(str(Path(sys.prefix)/'lib/libgomp.so.1'))
    blas.openblas_get_num_threads.restype=omp.omp_get_max_threads.restype=ctypes.c_int
    result['actual_threads']={'openblas':blas.openblas_get_num_threads(),'openmp':omp.omp_get_max_threads()}
    require(result['actual_threads']=={'openblas':12,'openmp':1},'actual numerical threads')
    backend=load(m['helpers']['backend'],'stage2_qualified_backend');nominal_helper=load(m['helpers']['slurm12'],'stage2_nominal_parameters')
    require(m['nominal_contract']['sha256']==nominal_helper.NOMINAL_SHA,'unchanged nominal contract')
    nominal=nominal_helper.nominal_parameters(a,a)
    require(nominal['U_R']==initializer['PME_checks']['operator_norms']['U_R_nominal'] and nominal['delta_odd']==initializer['PME_checks']['operator_norms']['evenness_perturbation_nominal'],'same initialized reciprocal bounds')
    arrays={name:np.load(verify(a['arrays'][name]),mmap_mode='r',allow_pickle=False) for name in ('data','indices','indptr')}
    for name,value in arrays.items():require(hashlib.sha256(memoryview(value).cast('B')).hexdigest()==a['payload_sha256'][name],'exact array payload')
    data,idx,ptr=(arrays[name] for name in ('data','indices','indptr'))
    require(data.dtype==np.float64 and idx.dtype==ptr.dtype==np.int32 and data.shape==idx.shape==(a['nnz'],) and ptr.shape==(N+1,),'exact array representation')
    k=csc_array((np.array(data,copy=True),np.array(idx,copy=True,order='C'),np.array(ptr,copy=True,order='C')),shape=(N,N),copy=False)
    backend.validate_lower(k);k.data[k.indptr[:-1]]-=nominal['tau'];backend.validate_lower(k)
    result.update(assembly=assembly_pin,operator={'shape':[N,N],'nnz':a['nnz'],**nominal,'shifted_data_sha256':hashlib.sha256(memoryview(k.data).cast('B')).hexdigest()},phase='symbolic')
    def on_phase(name):
        key={'symbolic_analyses':'symbolic_attempts','numeric_factorizations':'numeric_attempts'}[name];result[key]+=1
        require(result[key]<=1,'one factor attempt per phase');result['phase']=name;persist()
    def on_solve():result['solve_attempts']+=1;require(result['solve_attempts']<=10,'ten original probes maximum');persist()
    f,error,ww,pre=backend.factor_once(k,on_phase=on_phase)
    result.update(pre_numeric=pre,post_numeric=backend.metadata(f),backend_exception=error,backend_warnings=ww);persist()
    if error is not None or ww:
        result['status']='inconclusive_real_split_test' if error and error['type']=='CholmodNotPositiveDefiniteError' else ('inconclusive_numeric_warning' if not error else 'failed_static_numeric_attempt')
        result['phase']='discarded_factor_without_solve';del f;return
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter('always');audit,lower,p=backend.audit_factor(f,k,backend.ten_probes(N,(0,1,2)),on_solve=on_solve)
    result['audit']=audit;result['audit_warnings']=[{'type':type(x.message).__name__,'message':str(x.message)} for x in observed]
    passed=bool(audit['passed'] and not observed and result['solve_attempts']==10)
    result.update(sufficient_numerical_test_passed=passed,status='passed_nominal_static_sufficient_numerical_test' if passed else 'inconclusive_static_arithmetic_or_boundary',phase='complete')

def child(config,release,out,phase):
    directory,result=phase_record(config,release,out,phase);tick=time.monotonic();cpu=time.process_time()
    def persist():
        result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'openblas_threads':CAPS[phase]['BLAS'],'openmp_threads':1}
        save(directory/'RESULT.json',result)
    persist()
    try:
        parent=json.loads((out/'RESULT.json').read_text());intent=parent['commands'][phase]
        require(parent['implementation']==pin(__file__) and parent['configuration']==pin(config) and parent['release']==pin(release),'same supervisor intent')
        require(intent['status'] in ('intent_durable_before_launch','launched') and intent['argv']==[sys.executable,str(Path(__file__).resolve()),'--config',str(config),'--release',str(release),'--out',str(out),'--phase',phase],'exact child invocation')
        caps=CAPS[phase];set_limits(caps['AS'],caps['CPU'])
        cfg,m,rel,alloc,initializer=context(config,release,out,phase)
        result.update(state_id=cfg['state_id'],manifest=cfg['manifest'],inputs=m['rows'][cfg['state_id']]['inputs'],source_inputs=m['rows'][cfg['state_id']]['source_inputs'],source_manifest=m['native_sources'],build=cfg['build'],initializer=cfg['initializer'],allocation=alloc,limits=caps);persist()
        if phase=='assembly':assemble(cfg,m,initializer,directory,result,persist)
        else:factor(cfg,m,initializer,directory,result,persist,out)
    except BaseException as exc:result.update(status='stage2_'+phase+'_failed',error=repr(exc))
    persist();return 0 if result['status'] in ('stage2_assembly_complete','passed_nominal_static_sufficient_numerical_test') else 2

def supervise(config,release,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),'commands':{},'native_calls':0,'automatic_retry':False,'error':None}
    save(out/'RESULT.json',result)
    try:
        cfg,m,rel,alloc,initializer=context(config,release,out)
        result.update(state_id=cfg['state_id'],candidate=m['candidate'],manifest=cfg['manifest'],source_manifest=m['native_sources'],source_inputs=m['rows'][cfg['state_id']]['source_inputs'],inputs=m['rows'][cfg['state_id']]['inputs'],build=cfg['build'],initializer=cfg['initializer'],expected_state=initializer['state_output'],pme_export=initializer['pme_export'],allocation=alloc,e_app=None,application_error_bound_established=False,sufficient_numerical_test_passed=False)
        for phase in ('assembly','factor'):
            caps=CAPS[phase]
            argv=[sys.executable,str(Path(__file__).resolve()),'--config',str(config),'--release',str(release),'--out',str(out),'--phase',phase]
            invoke(result,out,phase,argv,HERE,caps['wall'],caps['AS'],caps['CPU'],caps['BLAS'])
            item=pin(out/phase/'RESULT.json');record=read(item);result[phase]=item
            wanted='stage2_assembly_complete' if phase=='assembly' else 'passed_nominal_static_sufficient_numerical_test'
            require(record['status']==wanted and record['state_id']==cfg['state_id'],'phase admission before continuation')
        result.update(status='stage2_state_nominal_certificate_passed',nominal=record['operator'],sufficient_numerical_test_passed=True)
    except BaseException as exc:result.update(status='stage2_state_certificate_failed',error=repr(exc))
    for phase in ('assembly','factor'):
        p=out/phase/'RESULT.json'
        if p.exists():result[phase]=pin(p)
    result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'controller_threads':1,'child_costs_in_commands':True,'scheduler_cost_separate':True}
    save(out/'RESULT.json',result);return 0 if result['status']=='stage2_state_nominal_certificate_passed' else 2

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--phase',choices=['assembly','factor']);a=p.parse_args()
    config,release,out=(x.resolve() for x in (a.config,a.release,a.out))
    raise SystemExit(child(config,release,out,a.phase) if a.phase else supervise(config,release,out))
