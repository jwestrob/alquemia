#!/usr/bin/env python3
"""One released saved-stencil preparation and one released convex diagnostic."""
import argparse,ast,copy,hashlib,json,math,os,socket,sys,time,types
from pathlib import Path
import clarabel,numpy as np,scipy
from scipy import sparse
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
NAMES=('E3/La','mu3/La','E3/Dy','mu3/Dy','E3/Dy_minus_La','mu3/Dy_minus_La')

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(isinstance(p,dict) and pin(p['path'])==p,'Changed or unbound pin');return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def token(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def save(p,v):
    p=Path(p);q=p.with_suffix('.tmp');q.write_text(json.dumps(plain(v),indent=2,allow_nan=False)+'\n');os.replace(q,p)
def plain(v):
    if isinstance(v,dict):return {str(k):plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):return [plain(x) for x in v]
    if isinstance(v,(float,np.floating)):return float(v) if math.isfinite(v) else repr(float(v))
    if isinstance(v,np.integer):return int(v)
    return v if v is None or isinstance(v,(str,int,bool)) else str(v)
def extract(p,names,ns):
    nodes=[n for n in ast.parse(verify(p).read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
    require({n.name for n in nodes}==set(names),'Exact helper functions');exec(compile(ast.Module(body=nodes,type_ignores=[]),p['path'],'exec'),ns)
def helpers(inv,base_pin):
    ns=dict(__file__=inv['sources']['math']['path'],math=math,np=np,ast=ast,copy=copy,Path=Path,hashlib=hashlib,json=json,require=require,token=token,FRACTIONS=(1.,.5,.25,.125),SLICES=('old_EF','pair_mu','C_mu','chi'))
    extract(inv['sources']['math'],['clone'],ns)
    module=types.SimpleNamespace(__file__=str(verify(inv['sources']['parent_core'])))
    edits=ns['clone'](module,ns,'certify_fraction',{13:15,12:14})
    extract(inv['sources']['original_cycle_helper'],['cycles'],ns)
    base=dict(np=np,math=math,clarabel=clarabel,read=read,require=require)
    extract(base_pin,['attributes','frozen_settings','cone_margins','diagnostics'],base)
    return ns,base,edits

def values_for(vector,inv,cycles):
    require(vector['status']=='complete','Complete admitted saved vector')
    vals={k:{'energy':v['energy_kcal_mol'],'mu_Debye':v['mu_Debye'],'status':'complete'} for k,v in inv['static_terms'].items()}
    vals.update({k:vector['cells'][k] for k in inv['required_dynamic_cell_ids']})
    rows=cycles(inv['cycle_dictionaries'],vals);require(all(v['status']=='complete' for v in rows.values()),'Both complete seven-term cycles')
    out=[]
    for metal in ('La','Dy'):
        v=rows['aquo9_central_'+metal+'__far'];out.extend([[v['candidate_E3_kcal_mol']],v['candidate_mu3_Debye']])
    return out

def build_groups(inv,cycles):
    center=read(inv['candidate']);require(token(center['parameters'])==inv['origin_key'],'Exact fixed origin')
    vals=values_for(center,inv,cycles);probes=[]
    for i,row in enumerate(inv['ordered_derivative_vectors']):
        v=read(row['artifact']);require(v['parameters']==row['parameters'] and token(v['parameters'])==row['key'] and row['index']==i,'Exact ordered saved stencil')
        probes.append(values_for(v,inv,cycles))
    require(len(probes)==28,'Exactly28 saved vectors')
    refs=[]
    for metal in ('La','Dy'):
        v=inv['cycle_dictionaries']['aquo9_central_'+metal+'__far'];refs.extend([[v['reference_electronic_E3_kcal_mol']],v['reference_mu3_Debye']])
    result=[]
    for i in range(4):
        r=[a-b for a,b in zip(vals[i],refs[i])]
        J=np.column_stack([(np.asarray(probes[2*j+1][i])-np.asarray(probes[2*j][i]))/span for j,span in enumerate(inv['actual_spans'])])
        result.append(dict(name=NAMES[i],r=r,J=J,model=vals[i],reference=refs[i]))
    for a,b in ((0,2),(1,3)):
        model=[x-y for x,y in zip(vals[b],vals[a])];reference=[x-y for x,y in zip(refs[b],refs[a])]
        result.append(dict(name=NAMES[len(result)],r=[x-y for x,y in zip(model,reference)],J=result[b]['J']-result[a]['J'],model=model,reference=reference))
    for g,authority in zip(result,inv['groups']):
        require(g['name']==authority['name'] and np.isfinite(g['J']).all(),'Finite exact group order')
        require(np.array_equal(np.asarray(g['r']),np.asarray(authority['r'])),'Literal admitted center residual')
        g['scale']=math.sqrt(math.fsum(v*v for v in g['r']));require(g['scale']>0 and math.isfinite(g['scale']),'Nonzero finite norm without epsilon floor')
    return result

def adapt(A,b,layout,groups,q0):
    require(A.shape[1]==15 and len(groups)==6 and [g['name'] for g in groups]==list(NAMES),'Fourteen d columns plus old s; exactly six groups')
    require(math.isfinite(q0) and q0>0 and np.isfinite(A.data).all() and np.isfinite(b).all(),'Finite saved conic matrix/cap')
    require([i for x in layout for i in range(x['start'],x['stop'])]==list(range(len(b))),'Ordered complete original row map')
    fixed=math.sqrt(q0);parts=[];rhs=[];new=[];cones=[];off=0;linear=[];removed=[]
    def add(name,kind,a,v,old=None):
        nonlocal off
        parts.append(sparse.csc_matrix(a));rhs.extend(v);new.append(dict(name=name,kind=kind,start=off,stop=off+len(v),original_rows=old));off+=len(v)
    for x in layout:
        ids=list(range(x['start'],x['stop']));a=A[ids,:14];v=b[ids]-np.asarray(A[ids,14].toarray()).ravel()*fixed
        a=sparse.hstack([a,sparse.csc_matrix((len(ids),1))],format='csc')
        if x['name'] in ('epigraph/lower','epigraph/upper'):removed.append(x['name']);continue
        if x['kind']=='SOC':add(x['name'],'SOC',a,v,ids);cones.append(clarabel.SecondOrderConeT(len(ids)))
        else:require(x['kind']=='nonnegative' and len(ids)==1,'Literal linear sides');linear.append((x['name'],a,v,ids))
    require(removed==['epigraph/lower','epigraph/upper'],'Remove only old epigraph bounds')
    for g in groups:
        r=np.asarray(g['r'],dtype=float);J=np.asarray(g['J'],dtype=float);s=float(g['scale'])
        require(J.shape==(len(r),14) and len(r) in (1,3) and s>0 and np.isfinite(r).all() and np.isfinite(J).all(),'Finite six group shape and positive scale')
        a=np.zeros((len(r)+1,15));a[0,14]=1.;a[1:,:14]=-J/s
        add('trimer/'+g['name'],'SOC',a,np.r_[1.,r/s]);cones.append(clarabel.SecondOrderConeT(len(r)+1))
    for name,a,v,ids in linear:add(name,'nonnegative',a,v,ids)
    lower=np.zeros((1,15));lower[0,14]=-1;add('fraction/lower','nonnegative',lower,[0.])
    add('fraction/upper','nonnegative',-lower,[1.]);cones.append(clarabel.NonnegativeConeT(len(linear)+2))
    return dict(P=sparse.csc_matrix((15,15)),A=sparse.vstack(parts,format='csc'),b=np.asarray(rhs),c=np.r_[np.zeros(14),-1.],cones=cones,layout=new,fixed_s=fixed,literal_q_ceiling=q0,removed=removed)

def certificate(p,groups,x,checker):
    x=np.asarray(x,dtype=float);require(x.shape==(15,) and np.isfinite(x).all(),'Finite14 displacement plus fraction')
    d=x[:14].tolist();t=float(x[14]);q0=p['raw_center_q']
    old=checker(p,d+[q0],1.,False)['certificate'];m={k:v for k,v in old['residual_constraints'].items() if k!='block/old_EF'}
    signed=p['signed_W9_z'];sv=signed['r_Debye']+math.fsum(float(a)*b for a,b in zip(signed['J_Debye'],d))
    m.update({'signed_W9/upper':signed['absolute_ceiling_Debye']-sv,'signed_W9/lower':signed['absolute_ceiling_Debye']+sv,'fraction/lower':t,'fraction/upper':1-t})
    params=[v+w*dv for v,w,dv in zip(p['center'],p['width'],d)]
    for j,dv in enumerate(d):
        for name,value in [('trust/lower',dv+.025),('trust/upper',.025-dv),('box/lower',params[j]-p['lower'][j]),('box/upper',p['upper'][j]-params[j])]:m[name+'/'+str(j)]=value
    details={}
    for g in groups:
        v=[float(r)+math.fsum(float(a)*b for a,b in zip(row,d)) for r,row in zip(g['r'],g['J'])];norm=math.sqrt(math.fsum(a*a for a in v));margin=(1-t)*g['scale']-norm
        m['trimer/'+g['name']]=margin;details[g['name']]=dict(signed_residual=v,reference=g['reference'],predicted_model=[a+b for a,b in zip(v,g['reference'])],norm=norm,center_norm=g['scale'],margin=margin,achieved_fraction=1-norm/g['scale'])
    q=max(old['blocks'][k]['ratio'] for k in ('P_energy','B_energy'));m['literal_q_cap']=q0-q
    return dict(parameters=params,normalized_displacement=d,fraction=t,q=q,literal_q_ceiling=q0,all25_predicted=old['all25_predicted'],original25=old['original25'],blocks=old['blocks'],signed_W9_residual=sv,groups=details,margins=m,minimum_margin=min(m.values()),exact_sign_feasible=all(v>=0 for v in m.values()) and old['all25_predicted'],native_admitted=False,global_family_claim=False)

def fingerprint(cfg):return token({k:cfg[k] for k in ('implementation','inventory','design_ready','conic_helpers','attempt_path','limits','runtime_declaration','design_admission')})
def release_gate(cfg,mode,fp):
    release=read(cfg['root_release'])
    require(release['status']=='released_trimer_local_span_'+mode and release['mode']==mode and release['configuration_fingerprint']==fp,'Exact released mode/fingerprint')
    require(all(release[k]==cfg[k] for k in ('implementation','output_directory','source_review','runtime_declaration')) and release['maximum_solver_calls']==(1 if mode=='solve' else 0),'Exact separate root release output/source/limits')

def gate(cfg,mode):
    require(mode in ('prepare','solve') and cfg['implementation']==pin(__file__),'Exact runtime source/mode')
    require(cfg['limits']=={'solver_calls':1,'native_calls':0,'new_probes':0,'threads':1,'solver_time_limit_seconds':60},'Exact finite limits')
    require(all(os.environ.get(k)=='1' for k in THREADS),'One-thread environment')
    require(not os.environ.get('SLURM_JOB_ID') and socket.gethostname()=='biotite','Truthful one-thread login source/runtime')
    fp=fingerprint(cfg);require(cfg['configuration_fingerprint']==fp,'Scientific configuration fingerprint')
    review=read(cfg['source_review']);require(review['status']=='passed_trimer_local_span_runtime_source_review' and review['implementation']==cfg['implementation'] and review['configuration_fingerprint']==fp,'Frozen source review')
    release_gate(cfg,mode,fp)
    inv=read(cfg['inventory']);ready=read(cfg['design_ready']);require(ready['inventory']==cfg['inventory'],'Frozen design inventory')
    declaration=read(cfg['runtime_declaration']);design=read(cfg['design_admission'])
    require(declaration['input_authority']==cfg['inventory'] and declaration['design_admission']==cfg['design_admission'] and design['ready']==cfg['design_ready'] and design['inventory']==cfg['inventory'],'Root design/runtime authority')
    require(declaration['maximum_real_solver_calls']==1 and declaration['maximum_numerical_threads']==1 and declaration['new_native_calls']==declaration['new_probe_calls']==declaration['new_reference_calls']==0,'Frozen root finite scope')
    require(read(inv['trimer_admission'])['result']==inv['trimer_result'] and read(inv['fit_root_admission'])['status'].startswith('admitted'),'Admitted saved terminal inputs')
    for source in inv['sources'].values():verify(source)
    env=read(inv['sources']['environment']);require(sys.prefix==env['prefix'] and {'clarabel':clarabel.__version__,'numpy':np.__version__,'scipy':scipy.__version__}==env['versions'],'Pinned existing environment')
    for file in env['files']:verify(file)
    return inv

def once(cfg,model,settings,out,result,base,factory=None):
    require(result['solver_calls']==0 and result['solver_constructions']==0,'One attempt only')
    lock=Path(cfg['attempt_path']);lock.parent.mkdir(parents=True,exist_ok=True)
    with lock.open('x') as f:json.dump({'configuration_fingerprint':cfg['configuration_fingerprint'],'implementation':cfg['implementation'],'PID':os.getpid(),'maximum_attempts':1},f);f.flush();os.fsync(f.fileno())
    result.update(status='attempt_reserved',attempt=pin(lock),solver_constructions=1);save(out/'RESULT.json',result)
    solver=(factory or clarabel.DefaultSolver)(model['P'],model['c'],model['A'],model['b'],model['cones'],settings)
    require(base['attributes'](solver.get_settings())==base['attributes'](settings),'All effective settings exact')
    solver.print_to_buffer();result.update(status='solve_entered',solver_calls=1);save(out/'RESULT.json',result)
    tick,cpu=time.monotonic(),time.process_time();sol=solver.solve();result['solver_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1};info=solver.get_info();(out/'SOLVER.log').write_text(solver.get_print_buffer())
    data=base['attributes'](info);data['linsolver']=base['attributes'](info.linsolver)
    result.update(solution=base['attributes'](sol),info=data,solver_status=str(sol.status),solver_success=str(sol.status)=='Solved',independent_residuals=base['diagnostics'](model,sol))
    return np.asarray(sol.x,dtype=float)

def run(path,mode):
    cfg=json.loads(Path(path).read_text());cfg['_config_path']=str(Path(path).resolve());out=Path(cfg['output_directory']);out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time()
    result=dict(status='started',mode=mode,implementation=pin(__file__),configuration=pin(path),solver_calls=0,solver_constructions=0,native_calls=0,new_probes=0,error=None)
    try:
        inv=gate(cfg,mode);h,base,edits=helpers(inv,cfg['conic_helpers']);p=read(inv['linearization']);q0=inv['literal_q_ceiling'];require(p['raw_center_q']==q0,'Literal original q cap')
        layout=read(inv['row_layout']);data=np.load(verify(inv['conic_data']),allow_pickle=False);A=sparse.csc_matrix((data['A_data'],data['A_indices'],data['A_indptr']),shape=tuple(data['A_shape']));b=data['b'];require(A.shape==(1437,15),'Saved step8 shape')
        if mode=='prepare':
            groups=build_groups(inv,h['cycles']);save(out/'GROUPS.json',groups)
        else:
            prep=read(cfg['preparation']);require(prep['status']=='trimer_local_span_prepared' and prep['implementation']==cfg['implementation'] and prep['configuration_fingerprint']==cfg['configuration_fingerprint'],'Exact successful preparation')
            groups=read(prep['groups'])
        model=adapt(A,b,layout,groups,q0);require(model['A'].shape==(1455,15) and len(model['cones'])==27,'Frozen dimensions')
        zero=certificate(p,groups,np.zeros(15),h['certify_fraction']);zm=base['cone_margins'](model,model['b']);result.update(configuration_fingerprint=cfg['configuration_fingerprint'],inventory=cfg['inventory'],zero_scalar=zero,zero_conic_margins=zm,dimension_edits=edits)
        save(out/'RESULT.json',result);require(zero['exact_sign_feasible'] and min(zm.values())>=0,'Zero witness fails; preserve margins without padding')
        settings=base['frozen_settings'](inv['sources']['settings']);result['effective_settings']=base['attributes'](settings)
        np.savez(out/'CONIC_DATA.npz',A_data=model['A'].data,A_indices=model['A'].indices,A_indptr=model['A'].indptr,A_shape=model['A'].shape,b=model['b'],c=model['c']);save(out/'ROW_LAYOUT.json',model['layout'])
        result.update(conic_data=pin(out/'CONIC_DATA.npz'),row_layout=pin(out/'ROW_LAYOUT.json'))
        if mode=='prepare':result.update(status='trimer_local_span_prepared',groups=pin(out/'GROUPS.json'))
        else:
            expected=np.load(verify(prep['conic_data']),allow_pickle=False)
            for k,v in [('A_data',model['A'].data),('A_indices',model['A'].indices),('A_indptr',model['A'].indptr),('A_shape',model['A'].shape),('b',model['b']),('c',model['c'])]:require(np.array_equal(expected[k],v),'Exact prepared matrix '+k)
            require(read(prep['row_layout'])==model['layout'],'Exact prepared row layout')
            x=once(cfg,model,settings,out,result,base)
            if x.shape==(15,) and np.isfinite(x).all():result['witness']=certificate(p,groups,x,h['certify_fraction'])
            else:result['witness_rejection']='Malformed/nonfinite solver point; no fallback'
            result.update(status='trimer_local_span_solve_complete',selected_candidate=None,native_admitted=False,global_infeasibility_claim=False)
    except Exception as exc:result.update(status='failed_trimer_local_span_no_retry',error=repr(exc));raise
    finally:
        result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'hostname':socket.gethostname(),'PID':os.getpid(),'job_id':None};save(out/'RESULT.json',result)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);parser.add_argument('--mode',choices=('prepare','solve'),required=True);args=parser.parse_args();run(args.config,args.mode)
