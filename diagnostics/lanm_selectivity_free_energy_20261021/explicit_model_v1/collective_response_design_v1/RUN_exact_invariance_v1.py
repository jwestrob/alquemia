#!/usr/bin/env python3
"""Fixed, unfit EXCHPOL response matrix; preparation never invokes a native program."""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import threading
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
Q=(.01,.05)
EXPONENT=4.0
R0=2.5
CELLS=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central') for m in ('La','Dy'))
CLUSTERS=tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))
CONTROLS=('zero_aquo9_La','zero_aquo9_Dy','shell_zero','shell_q05')
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
POINT=[3.919394849221829,.17662704152555941,3.640166322240097,.15457541153280463,.15863887317162664,.202479823314524]
MAX_PROGRAMS=116
MAX_FD=636

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    p=Path(path);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(t,p)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def fingerprint(cfg):return hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k not in ('root_release','preparation_result')},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def label(q):return 'q01' if q==.01 else 'q05' if q==.05 else 'zero'

def overlap(r):
    require(math.isfinite(r) and r>0,'Invalid distance');t=.5*EXPONENT*r;s=(1+t+t*t/3)*math.exp(-t)
    return s*s

def render(key,q):
    require(q in (0.,*Q),'Undeclared amplitude')
    require(not re.search(r'(?im)^\s*(?:exchpol|exchange-polar|polar-print)\b',key),'Existing EXCHPOL/print record')
    k=q/overlap(R0)
    extra=['# Fixed collective-response diagnostic; original bytes above unchanged','exchange-polar S2U','polar-print',
        f'exchpol 34 {k:.17g} 1 4 0','exchpol 35 0 0 4 0','exchpol 500 0 1 4 1','exchpol 501 0 1 4 1']
    return key+'\n'+'\n'.join(extra)+'\n'

def water_shell(text,pair):
    rows=pair.typed_rows(text);require(len(rows)==28 and int(rows[0][5])==500,'Expected W9 La source')
    out=[]
    for r in rows[1:]:
        require(int(r[5]) in (36,37) and all(int(v)>1 for v in r[6:]),'Water shell graph/type differs')
        out.append([str(int(r[0])-1),*r[1:6],*[str(int(v)-1) for v in r[6:]]])
    return '27  literal W9 metal-absent water control\n'+'\n'.join(' '.join(r) for r in out)+'\n'

def tensors(text,q,pair,resolver,prm,key):
    """Source reconstruction of native site matrices, not native tensor readback."""
    rows=pair.typed_rows(text);p=resolver.parse_parameters(prm,key);metal=[i for i,r in enumerate(rows) if int(r[5]) in (500,501)]
    require(len(metal)<=1,'Unexpected multiple metals');out=[]
    for r in rows:
        typ=int(r[5]);require(typ in (36,37,500,501),'Unqualified represented type')
        require(p['atom'][typ]['class']=={36:34,37:35,500:500,501:501}[typ],'EXCHPOL class mapping differs')
        require(p['polarize'][typ]['alpha_A3']>0,'Native dexpol requires positive alpha at every represented site')
    for i,r in enumerate(rows):
        matrix=[[float(j==k) for k in range(3)] for j in range(3)];increment=0.;distance=None
        if int(r[5])==36 and metal:
            a=[float(v) for v in r[2:5]];b=[float(v) for v in rows[metal[0]][2:5]];v=[x-y for x,y in zip(a,b)];distance=math.sqrt(sum(x*x for x in v))
            require(distance<5.4,'Diagnostic reached native repulsion taper; reconstruction would need its polynomial')
            increment=q*overlap(distance)/overlap(R0);u=[x/distance for x in v]
            matrix=[[float(j==k)+increment*u[j]*u[k] for k in range(3)] for j in range(3)]
        require(increment>=0 and all(math.isfinite(x) for rr in matrix for x in rr),'Invalid tensor')
        # Explicit Cholesky verifies the full off-diagonal matrix, not just positive diagonals.
        L=[[0.]*3 for _ in range(3)]
        for j in range(3):
            for k in range(j+1):
                x=matrix[j][k]-sum(L[j][h]*L[k][h] for h in range(k))
                if j==k:require(x>0,'Nonpositive tensor');L[j][k]=math.sqrt(x)
                else:L[j][k]=x/L[k][k]
        out.append({'atom_index':i+1,'type':int(r[5]),'matrix':matrix,'eigenvalues':[1.,1.,1.+increment],
            'maximum_absolute_off_diagonal':max(abs(matrix[j][k]) for j in range(3) for k in range(3) if j!=k),'M_O_distance_A':distance})
    return {'origin':'pinned native formula applied to exact prepared geometry; not executable tensor readback','positive_definite':True,'sites':out}

def default_config():
    e=json.loads((D/'explicit_model_v1/paired_water_coordination_joint_v1/EVALUATOR_CONFIG.json').read_text())
    c=json.loads((D/'explicit_model_v1/paired_water_coordination_joint_v1/CONFIG_PREPARE.json').read_text())
    src=W/'explicit_model_v1/tinker_source/source'
    return {'protocol':'collective_response_exchpol_v1','implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),
        'data':pin(W/'explicit_model_v1/paired_water_coordination_joint_v1/metal_alpha_response_v1/preparation_v2/DATA.json'),
        'pair_adapter':pin(D/'reference_chemistry_v1/matched_pair_decomposition_v1/model/RUN.py'),
        'pair_result':pin(W/'reference_chemistry_v1/matched_pair_decomposition_v1/model/run_v1/execution/RESULT.json'),
        'native_helper':e['native_helper'],'resolver':e['resolver'],'executables':e['executables'],'objective':c['objective'],'body_helper':c['body_helper'],
        'physics_sources':{n:pin(src/n) for n in ('kexpol.f','kpolar.f','prmkey.f','alterpol.f','induce.f','epolar.f','epolar1.f','epolar3.f','dexpol.f','damping.f','cutoffs.f','switch.f','analyze.f','analysis.f')},
        'q_values':list(Q),'native_overlap_exponent_A_inverse':EXPONENT,'normalization_distance_A':R0,
        'maximum_native_programs':MAX_PROGRAMS,'maximum_internal_FD_energies':MAX_FD,'workers':4,'root_release':None,'preparation_result':None}

def checked_config(cfg,threads=False):
    require(cfg['protocol']=='collective_response_exchpol_v1' and cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Source/protocol differs')
    require(cfg['q_values']==list(Q) and cfg['native_overlap_exponent_A_inverse']==EXPONENT and cfg['normalization_distance_A']==R0,'Frozen probe settings differ')
    require(cfg['maximum_native_programs']==MAX_PROGRAMS and cfg['maximum_internal_FD_energies']==MAX_FD and cfg['workers']==4,'Scope differs')
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
    for group in ('executables','physics_sources'):
        for v in cfg[group].values():verify(v)
    if threads:require(all(os.environ.get(k)=='1' for k in THREADS),'One-thread environment required')

def context(cfg):
    data=read(cfg['data']);pair=module(cfg['pair_adapter'],'exch_pair');native=module(cfg['native_helper'],'exch_native');resolver=module(cfg['resolver'],'exch_resolver')
    require(data['baseline_parameters']==POINT and set(data['rows'])==set(CELLS),'Frozen ten-cell point differs')
    review=read(data['provenance']['selected_review']);verify(data['provenance']['candidate_result'])
    require(data['provenance']['candidate_result']['sha256']=='5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832','Wrong selected result')
    require(review['status']=='selected_ten_raw_review_passed' and review['all_selected_raw_replays_exact'] and review['exact_selected_score_replayed'] and review['all_development_controls_passed'] is False and review['error'] is None,'Raw/science admission differs')
    for gid in CELLS:
        m=data['baseline_model'][gid];n=m['native'];raw=native.parse(verify(n['analyze_output']),verify(n['gradient_output']),{'atoms':data['rows'][gid]['geometry']['atoms'],'derivative_test':n['finite_difference'] is not None})
        require(raw['total_energy_kcal_mol']==m['energy'] and raw['force_kcal_mol_A']==m['forces'],'Frozen E/F replay differs')
        verify(m['provenance']['native_authority'])
    prior=read(cfg['pair_result']);pm=read(prior['manifest']);admitted=pair.admit_result(prior['config'],pm,prior)
    require(set(admitted['rows'])==set(pair.IDS),'Missing literal dimers')
    objective=module(cfg['objective'],'exch_objective');body=module(cfg['body_helper'],'exch_body')
    score=objective.evaluate(data,data['baseline_model'],body);require(score==data['baseline_score'],'Saved baseline score differs')
    return data,pair,native,resolver,pm,prior,objective,body

def prepare(cfg,out):
    checked_config(cfg,True);data,pair,native,resolver,pm,prior,objective,body=context(cfg);rows={}
    def add(cid,kind,q,source,xyz_text,atoms,stage,gid=None,pid=None,fd=False):
        base=verify(source['key']).read_text();prm=verify(source['parameters']).read_text();rendered=render(base,q)
        require(not re.search(r'(?im)^\s*(?:exchpol|exchange-polar)\b',prm),'Parameter EXCHPOL already active')
        folder=out/'inputs'/cid;folder.mkdir(parents=True);(folder/'system.xyz').write_text(xyz_text);(folder/'system.key').write_text(rendered);shutil.copyfile(verify(source['parameters']),folder/'amoeba09.prm')
        require(re.findall(r'(?im)^\s*parameters\s+(\S+)',rendered)==['amoeba09.prm'],'Parameter basename differs')
        rows[cid]={'id':cid,'kind':kind,'q':q,'gid':gid,'pair_id':pid,'atoms':atoms,'stage':stage,'derivative_test':fd,
            'input':{k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))},'source_input':source,
            'tensor_control':tensors(xyz_text,q,pair,resolver,prm,base),'native_programs':1 if kind=='pair' else 2,'internal_FD_energies':6*atoms if fd else 0}
    for m in ('La','Dy'):
        gid='aquo9_central_'+m;inp=data['baseline_model'][gid]['provenance']['effective_input'];add('zero_aquo9_'+m,'zero',0.,inp,verify(inp['xyz']).read_text(),28,'controls',gid=gid)
    inp=data['baseline_model']['aquo9_central_La']['provenance']['effective_input'];shell=water_shell(verify(inp['xyz']).read_text(),pair)
    for cid,q in (('shell_zero',0.),('shell_q05',.05)):add(cid,'shell',q,inp,shell,27,'controls')
    for q in Q:
        for gid in CELLS:
            inp=data['baseline_model'][gid]['provenance']['effective_input'];fd=q==.05 and gid in CLUSTERS
            add(label(q)+'__'+gid,'fixture',q,inp,verify(inp['xyz']).read_text(),data['rows'][gid]['geometry']['atoms'],'derivatives' if fd else 'matrix',gid=gid,fd=fd)
        for pid in pair.IDS:
            inp=pm['rows'][pid]['input'];add(label(q)+'__'+pid,'pair',q,inp,verify(inp['xyz']).read_text(),4,'matrix',gid=pm['rows'][pid]['parent_geometry_id'],pid=pid)
    require(len(rows)==92 and sum(r['native_programs'] for r in rows.values())==MAX_PROGRAMS and sum(r['internal_FD_energies'] for r in rows.values())==MAX_FD,'Inventory/cost differs')
    manifest={'status':'prepared','implementation':cfg['implementation'],'config_fingerprint':fingerprint(cfg),'data':cfg['data'],'pair_result':cfg['pair_result'],
        'candidate_result':data['provenance']['candidate_result'],'selected_review':data['provenance']['selected_review'],'rows':rows,
        'baseline_cycle':prior['model_cycle'],'analytic_bare_ions':prior['analytic_bare_ions'],
        'overlap':{'native_type':'S2U','native_exponent_A_inverse':EXPONENT,'t':'0.5*sqrt(b_i*b_j)*r','S2':'((1+t+t*t/3)*exp(-t))**2 * size_i*size_j','normalization_r_A':R0,'normalization_S2':overlap(R0),'O_springs':{label(q):q/overlap(R0) for q in Q}},
        'stages':{s:[cid for cid,r in rows.items() if r['stage']==s] for s in ('controls','derivatives','matrix')},
        'maximum_native_programs':MAX_PROGRAMS,'maximum_internal_FD_energies':MAX_FD,'new_native_calls':0,'candidate_selection':False}
    save(out/'MANIFEST.json',manifest);return manifest

def bound_inputs(cfg,manifest):
    checked_config(cfg);require(manifest['implementation']==cfg['implementation'] and manifest['config_fingerprint']==fingerprint(cfg),'Manifest source/config differs')
    require(manifest['data']==cfg['data'] and manifest['pair_result']==cfg['pair_result'] and len(manifest['rows'])==92,'Manifest lineage/inventory differs')
    pair=module(cfg['pair_adapter'],'exch_bound_pair')
    for row in manifest['rows'].values():
        inp=row['input'];src=row['source_input'];base=verify(src['key']).read_text()
        require(verify(inp['key']).read_text()==render(base,row['q']) and verify(inp['parameters']).read_bytes()==verify(src['parameters']).read_bytes(),'Rendered parameters differ')
        expected=water_shell(verify(src['xyz']).read_text(),pair) if row['kind']=='shell' else verify(src['xyz']).read_text()
        require(verify(inp['xyz']).read_text()==expected,'Prepared geometry differs')
    require(sum(r['native_programs'] for r in manifest['rows'].values())==MAX_PROGRAMS and sum(r['internal_FD_energies'] for r in manifest['rows'].values())==MAX_FD,'Manifest limits differ')
    return manifest

def induction(path):
    text=Path(path).read_text();matches=re.findall(r'Induced Dipoles\s*:\s*Iterations\s+(\d+)\s+RMS Residual\s+([-+0-9.Ee]+)',text)
    require(matches,'Native induction summary missing');values=[(int(i),float(r)) for i,r in matches]
    require(all(i>=1 and math.isfinite(r) and 0<=r<=1e-8 for i,r in values),'Induction residual not admitted')
    return {'solutions':len(values),'maximum_iterations':max(i for i,r in values),'maximum_printed_RMS_residual_Debye':max(r for i,r in values),'tolerance_Debye':1e-8,'printed_resolution_Debye':1e-10}

class Counter:
    def __init__(self):self.lock=threading.RLock();self.stop=threading.Event();self.count=0;self.fd=0
    def reserve(self,fd=0):
        with self.lock:
            require(not self.stop.is_set() and self.count<MAX_PROGRAMS and self.fd+fd<=MAX_FD,'Global stop or finite ceiling');self.count+=1;self.fd+=fd

def release_gate(cfg,prep,manifest,runtime=True):
    require(cfg['root_release'] is not None,'Root release required');r=read(cfg['root_release'])
    require(r['status']=='approved_collective_response_exchpol_diagnostic' and r['implementation']==cfg['implementation'] and r['config_fingerprint']==fingerprint(cfg) and r['preparation_result']==cfg['preparation_result'],'Release lineage differs')
    require(r['maximum_native_programs']==MAX_PROGRAMS and r['maximum_internal_FD_energies']==MAX_FD and r['workers']==4 and r['fit_authorized'] is False,'Release scope differs')
    require(prep['status']=='prepared' and prep['new_native_calls']==0 and prep['manifest']==pin(verify(prep['manifest'])),'Preparation not admitted')
    if runtime:require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4,'One four-CPU task required')

def prechecks(rows,data):
    require(all(rows[c]['status']=='complete' for c in CONTROLS),'Incomplete pre-gates')
    for m in ('La','Dy'):
        a=rows['zero_aquo9_'+m]['model'];b=data['baseline_model']['aquo9_central_'+m]
        require(a['energy']==b['energy'] and a['forces']==b['forces'],'Zero-q printed E/F replay failed')
    a=rows['shell_zero']['model'];b=rows['shell_q05']['model']
    require(a['energy']==b['energy'] and a['forces']==b['forces'],'Metal-absent water printed E/F invariance failed')
    return {'passed':True,'zero_q_exact_metals':['La','Dy'],'metal_absent_water_exact':True,'comparison':'native printed energy and Cartesian force arrays'}

def summaries(rows,manifest,data,pair,objective,body):
    result={}
    for q in Q:
        lab=label(q);model={gid:rows[lab+'__'+gid]['model'] for gid in CELLS};score=objective.evaluate(data,model,body)
        fragments={pid:{**rows[lab+'__'+pid],'parent_geometry_id':manifest['rows'][lab+'__'+pid]['gid']} for pid in pair.IDS}
        clusters={gid:{'energy_kcal_mol':model[gid]['energy']} for gid in CLUSTERS};cycle=pair.model_cycle(fragments,clusters,manifest['analytic_bare_ions'])
        components={}
        names=set().union(*(rows[lab+'__'+pid]['components'] for pid in pair.IDS),*(model[gid]['native']['components'] for gid in CLUSTERS))
        for name in sorted(names):
            signs={f'aquo{n}_central_{m}':1 if (n,m) in ((9,'Dy'),(8,'La')) else -1 for n in (8,9) for m in ('La','Dy')}
            C=math.fsum(signs[g]*model[g]['native']['components'].get(name,{}).get('energy_kcal_mol',0.) for g in CLUSTERS)
            P=math.fsum(signs[manifest['rows'][lab+'__'+pid]['gid']]*rows[lab+'__'+pid]['components'].get(name,{}).get('energy_kcal_mol',0.) for pid in pair.IDS)
            components[name]={'C':C,'P':P,'B':C-P}
        result[lab]={'q':q,'score':score,'cycle':cycle,'components':components,'cycle_change':{k:cycle[k]-manifest['baseline_cycle'][k] for k in ('C','P','B')},'selected':False}
    result['scale_dependence']={k:result['q05']['cycle_change'][k]/.05-result['q01']['cycle_change'][k]/.01 for k in ('C','P','B')}
    return result

def run(cfg,out):
    start=time.monotonic();cpu=time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN);checked_config(cfg,True);prep=read(cfg['preparation_result']);manifest=bound_inputs(cfg,read(prep['manifest']));release_gate(cfg,prep,manifest)
    data=read(cfg['data']);pair=module(cfg['pair_adapter'],'exch_run_pair');native=module(cfg['native_helper'],'exch_run_native');objective=module(cfg['objective'],'exch_run_obj');body=module(cfg['body_helper'],'exch_run_body');counter=Counter()
    result={'status':'running','implementation':cfg['implementation'],'config':cfg,'manifest':prep['manifest'],'job_id':os.environ['SLURM_JOB_ID'],
        'rows':{cid:{'id':cid,'status':'not_attempted'} for cid in manifest['rows']},'prechecks':None,'summaries':None,'native_program_invocations':0,'internal_FD_energies_requested':0,'internal_FD_energies_completed':0,'global_stop':False,'error':None,'candidate_selection':False,'fit_targets_changed':False}
    def persist():
        with counter.lock:
            result.update(native_program_invocations=counter.count,internal_FD_energies_requested=counter.fd,global_stop=counter.stop.is_set());save(out/'RESULT.json',result)
    def worker(cid):
        if counter.stop.is_set():return
        row=manifest['rows'][cid];folder=out/'cells'/cid;record={'id':cid,'status':'started','q':row['q'],'input':row['input'],'receipts':[]}
        try:
            folder.mkdir(parents=True)
            for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),folder/n)
            record['runtime_input']={k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            with counter.lock:result['rows'][cid]={'id':cid,'status':'started'};save(folder/'RESULT.json',record)
            def invoke(argv,name,fd=0):
                with counter.lock:counter.reserve(fd);save(out/'COUNTERS.json',{'native_program_invocations':counter.count,'internal_FD_energies_requested':counter.fd})
                try:native.invoke(argv,folder,folder/name,record['receipts']);return induction(folder/name)
                except Exception:counter.stop.set();raise
            checks={'analyze':invoke([cfg['executables']['analyze']['path'],'system.xyz','E'],'analyze.out')}
            if row['kind']=='pair':parsed=pair.parse_analyze(folder/'analyze.out')
            else:
                args=[cfg['executables']['testgrad']['path'],'system.xyz','Y']+(['Y','0.00001'] if row['derivative_test'] else ['N'])
                checks['testgrad']=invoke(args,'testgrad.out',row['internal_FD_energies']);parsed=native.parse(folder/'analyze.out',folder/'testgrad.out',row)
                if row['derivative_test']:
                    with counter.lock:result['internal_FD_energies_completed']+=row['internal_FD_energies']
                    record['numerical_check']=parsed['finite_difference'];require(parsed['finite_difference']['passed'],'Cartesian derivative gate failed')
                record['model']={'energy':parsed['total_energy_kcal_mol'],'forces':parsed['force_kcal_mol_A'],'native':parsed}
            record.update(status='complete',energy_kcal_mol=parsed['total_energy_kcal_mol'],components=parsed['components'],native=parsed,induction=checks)
        except Exception as exc:
            counter.stop.set();record.update(status='failed_or_unavailable',error=repr(exc))
            with counter.lock:result['error']=repr(exc)
        finally:
            with counter.lock:
                result['rows'][cid]=record
                if folder.exists():save(folder/'RESULT.json',record)
                persist()
    try:
        persist()
        for stage in ('controls','derivatives','matrix'):
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(worker,manifest['stages'][stage]))
            require(not counter.stop.is_set() and all(result['rows'][c]['status']=='complete' for c in manifest['stages'][stage]),'Incomplete stage: '+stage)
            if stage=='controls':result['prechecks']=prechecks(result['rows'],data)
            persist()
        require(counter.count==MAX_PROGRAMS and counter.fd==MAX_FD and result['internal_FD_energies_completed']==MAX_FD,'Final count differs')
        result['summaries']=summaries(result['rows'],manifest,data,pair,objective,body);result['status']='collective_response_diagnostic_complete'
    except Exception as exc:counter.stop.set();result.update(status='failed_or_unavailable',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,'aggregate_child_user_CPU_seconds':ch.ru_utime-child.ru_utime,'aggregate_child_system_CPU_seconds':ch.ru_stime-child.ru_stime,'allocated_CPUs':4,'per_receipt_child_CPU_deltas_not_summed':True};persist()
    return result

def admit_result(cfg,manifest,result):
    """Pure saved replay; no molecular invocation and no parameter selection."""
    bound_inputs(cfg,manifest);release_gate(cfg,read(cfg['preparation_result']),manifest,runtime=False)
    require(result['implementation']==cfg['implementation'] and result['config']==cfg and read(result['manifest'])==manifest,'Result lineage differs')
    require(result['status']=='collective_response_diagnostic_complete' and result['error'] is None and result['global_stop'] is False,'Incomplete/failed result')
    require(set(result['rows'])==set(manifest['rows']) and result['native_program_invocations']==MAX_PROGRAMS and result['internal_FD_energies_requested']==MAX_FD and result['internal_FD_energies_completed']==MAX_FD,'Full result denominator/count differs')
    data=read(cfg['data']);pair=module(cfg['pair_adapter'],'exch_admit_pair');native=module(cfg['native_helper'],'exch_admit_native');objective=module(cfg['objective'],'exch_admit_obj');body=module(cfg['body_helper'],'exch_admit_body')
    for cid,row in result['rows'].items():
        frozen=manifest['rows'][cid];require(row['status']=='complete' and row['id']==cid and row['q']==frozen['q'] and row['input']==frozen['input'],'Native row identity differs')
        for k in ('xyz','key','parameters'):require(verify(row['runtime_input'][k]).read_bytes()==verify(frozen['input'][k]).read_bytes(),'Native input bytes differ')
        a=verify(row['native']['analyze_output']);checks={'analyze':induction(a)}
        if frozen['kind']=='pair':parsed=pair.parse_analyze(a)
        else:
            g=verify(row['native']['gradient_output']);parsed=native.parse(a,g,frozen);checks['testgrad']=induction(g)
            require(row['model']=={'energy':parsed['total_energy_kcal_mol'],'forces':parsed['force_kcal_mol_A'],'native':parsed},'Raw E/F replay differs')
            if frozen['derivative_test']:require(parsed['finite_difference']['passed'] and row['numerical_check']==parsed['finite_difference'],'FD replay differs')
        require(parsed==row['native'] and checks==row['induction'] and row['energy_kcal_mol']==parsed['total_energy_kcal_mol'] and row['components']==parsed['components'],'Raw component/induction replay differs')
        expected=[[cfg['executables']['analyze']['path'],'system.xyz','E']]
        if frozen['kind']!='pair':expected.append([cfg['executables']['testgrad']['path'],'system.xyz','Y']+(['Y','0.00001'] if frozen['derivative_test'] else ['N']))
        require(len(row['receipts'])==len(expected),'Receipt count differs')
        for i,(receipt,argv) in enumerate(zip(row['receipts'],expected)):
            require(receipt['argv']==argv and receipt['returncode']==0 and receipt['cwd']==str(Path(row['runtime_input']['xyz']['path']).parent) and receipt['output']==row['native']['analyze_output' if i==0 else 'gradient_output'] and all(receipt['threads'][k]=='1' for k in THREADS),'Native receipt differs')
    require(prechecks(result['rows'],data)==result['prechecks'],'Precheck replay differs')
    summary=summaries(result['rows'],manifest,data,pair,objective,body);require(summary==result['summaries'],'Score/decomposition replay differs')
    return {'status':'saved_collective_response_replay_passed','complete_rows':92,'new_native_calls':0,'summaries':summary,'candidate_selection':False}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('config','prepare','run'));p.add_argument('--config');p.add_argument('--output',required=True);a=p.parse_args()
    if a.stage=='config':save(a.output,default_config());raise SystemExit(0)
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);cfg=json.loads(Path(a.config).read_text());start=time.monotonic();cpu=time.process_time()
    try:
        if a.stage=='prepare':
            require(cfg['root_release'] is None and cfg['preparation_result'] is None,'Preparation must be unreleased');prepare(cfg,out)
            r={'status':'prepared','implementation':pin(__file__),'config':pin(a.config),'manifest':pin(out/'MANIFEST.json'),'new_native_calls':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu};save(out/'RESULT.json',r)
        else:r=run(cfg,out)
    except Exception as exc:r={'status':'precall_or_preparation_failed','implementation':pin(__file__),'error':repr(exc)};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'output':str(out)}));raise SystemExit(0 if r['status'] in ('prepared','collective_response_diagnostic_complete') else 1)
