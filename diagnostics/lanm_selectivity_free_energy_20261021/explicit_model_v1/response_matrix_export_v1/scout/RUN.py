#!/usr/bin/env python3
"""Two native reporting calls; reusable saved matrix qualification only."""
import argparse
import copy
from decimal import Decimal
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import time
import numpy as np

HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
FLAG = 'RESPONSE-MATRIX-EXPORT'
NAMES = ('response_matrix_001.dat','response_matrix_002.dat')

def require(ok, message):
    if not ok: raise ValueError(message)

def pin(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def verify(p):
    require(pin(p['path'])==p,'Changed pinned file: '+p['path']);return Path(p['path'])

def read(p): return json.loads(verify(p).read_text())

def save(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)

def module(p):
    s=importlib.util.spec_from_file_location('frozen_scout_field_helper',verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def numeric(values):
    a=np.array([float(v.replace('D','E')) for v in values]);require(np.isfinite(a).all(),'Nonfinite exported number');return a

def arithmetic(dim,*values):
    scale=max([1.]+[float(np.max(np.abs(v))) for v in values if np.asarray(v).size])
    return 128*dim*np.finfo(float).eps*scale

def settings(row):
    alpha={};peek=1.
    for source in ('parameters','key'):
        for line in verify(row['input'][source]).read_text().splitlines():
            t=line.split('#',1)[0].split()
            if not t:continue
            if t[0].lower()=='polarize':alpha[int(t[1])]=float(t[2])
            if t[0].lower()=='pcg-peek':peek=float(t[1])
    require(math.isfinite(peek),'Nonfinite native peek')
    return alpha,peek

def parity(new_moments,baseline_moments,exact=False):
    """Mandatory energy/components; unavailable historical moment arrays explicit."""
    report={'passed':True,'exact_requested':exact,'checks':{},'untested_channels':[]}
    def compare(name,a,b,tol):
        x,y=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
        require(x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all(),'Parity shape/nonfinite: '+name)
        delta=float(np.max(np.abs(x-y))) if x.size else 0.
        limit=0. if exact else tol
        report['checks'][name]={'maximum_difference':delta,'tolerance':limit,'printed_exact':bool(np.array_equal(x,y)),'passed':delta<=limit}
        report['passed'] &= delta<=limit
    compare('total_energy',new_moments['energy_kcal_mol'],baseline_moments['energy_kcal_mol'],1e-8)
    require(set(new_moments['components'])==set(baseline_moments['components']),'Energy component inventory differs')
    for name,old in baseline_moments['components'].items():
        new=new_moments['components'][name]
        require(new['interactions']==old['interactions'],'Component interaction count differs')
        compare('component_'+name,new['energy_kcal_mol'],old['energy_kcal_mol'],1e-8)
    for name,tol in [('induced_Debye',1e-6),('total_COM_Debye',1e-3)]:
        if baseline_moments.get(name) is not None:compare(name,new_moments[name],baseline_moments[name],tol)
        else:report['untested_channels'].append(name)
    return report

def parse_export(path,row,units):
    """Read only the fixed V1 full denominator, preserving native row order."""
    path=verify(path) if isinstance(path,dict) else Path(path)
    lines=[x.split() for x in path.read_text().splitlines() if x.strip()]
    n=row['atoms'];dim=3*n
    require(1<=n<=28 and len(lines)==8+n+2*dim,'Export line denominator')
    require(lines[0]==['TINKER_RESPONSE_MATRIX_V1'] and lines[-1]==['END_TINKER_RESPONSE_MATRIX_V1'],'Export protocol markers')
    require(lines[1][0]=='COUNTS' and len(lines[1])==5,'Export counts grammar')
    seq,nn,npole,dd=map(int,lines[1][1:]);require(seq in (1,2) and (nn,npole,dd)==(n,n,dim),'Export scope/dimension')
    require(lines[2]==['UNITS','A=Angstrom^-3','RHS=e/Angstrom^2','U=eA'],'Export units')
    require(lines[3][0]=='CONSTANTS' and len(lines[3])==4,'Constants grammar')
    constants=numeric(lines[3][1:]);require(constants.tolist()==[units['electric'],1.,units['debye_per_eA']],'Native constants differ')
    require(lines[4][0]=='EXTERNAL' and len(lines[4])==5 and lines[4][1]=='F' and np.array_equal(numeric(lines[4][2:]),np.zeros(3)),'Zero-field export required')
    require(lines[5][0]=='CHECKS' and len(lines[5])==7,'Export checks grammar');checks=numeric(lines[5][1:]);require((checks>=0).all(),'Negative native check')
    require(lines[6]==['RESTORED','T'],'Exact native restoration required')
    xyz=[s.split() for s in verify(row['input']['xyz']).read_text().splitlines()[1:] if s.strip()]
    require(len(xyz)==n and [int(s[5]) for s in xyz]==row['types'],'Input site/type inventory')
    alpha,peek=settings(row);sites=[];seen=set()
    for pos,t in enumerate(lines[7:7+n]):
        require(t[0]=='SITE' and len(t)==9 and int(t[1])==pos+1,'Export site grammar/order')
        atom,typ=int(t[2]),int(t[3]);require(1<=atom<=n and atom not in seen,'Duplicate/missing atom');seen.add(atom)
        a=numeric(t[5:]);expected=numeric([xyz[atom-1][2],xyz[atom-1][3],xyz[atom-1][4]])
        require(typ==row['types'][atom-1] and t[4]==row['native_names'][atom-1],'Export atom type/name map')
        require(a[0]==alpha[typ] and a[0]>=1e-8 and np.max(abs(a[1:]-expected))<=arithmetic(dim,expected),'Export alpha/coordinate map')
        sites.append(dict(position=pos+1,atom=atom,type=typ,name=t[4],alpha_A3=float(a[0]),xyz_A=a[1:].tolist()))
    vectors=[];matrix=[]
    for i in range(dim):
        v,a=lines[7+n+2*i:9+n+2*i]
        require(v[0]=='VECTOR' and len(v)==6 and int(v[1])==i+1,'Vector denominator/order')
        require(a[0]=='A' and len(a)==dim+2 and int(a[1])==i+1,'Matrix row denominator/order')
        vectors.append(numeric(v[2:]).tolist());matrix.append(numeric(a[2:]).tolist())
    return dict(source=pin(path),sequence=seq,sites=sites,vectors=vectors,A=matrix,checks=checks.tolist(),restored=True,peek=peek)

def induction_values(output):
    values=re.findall(r'Induced Dipoles\s*:\s*Iterations\s+(\d+)\s+RMS Residual\s+(\S+)',verify(output).read_text())
    require(len(values)==2,'Exactly two EM induction summaries')
    q=[float(v) for _,v in values];require(all(math.isfinite(v) and 0<=v<=1e-8 for v in q),'Original induction gate');return q

def half_print(token):
    return float(Decimal('0.5')*Decimal(10)**Decimal(token.replace('D','E')).as_tuple().exponent)

def qualify_exports(export_paths,row,units,new_moments,field_targets=None,report=None):
    """Generic1..28site saved qualification; field_targets only for this scout."""
    report={} if report is None else report;report.update(checks={},exports=[],field_predictions={})
    def gate(name,error,limit):
        passed=math.isfinite(float(error)) and math.isfinite(float(limit)) and 0<=error<=limit
        report['checks'][name]={'error':float(error),'bound':float(limit),'passed':bool(passed)}
        require(passed,'Numerical export qualification: '+name)
    require(len(export_paths)==2,'Exactly two exports')
    parsed=[parse_export(p,row,units) for p in export_paths];report['exports']=parsed
    require([p['sequence'] for p in parsed]==[1,2] and parsed[0]['sites']==parsed[1]['sites'],'Two solve/site identities')
    q=induction_values(new_moments['analyze_output']);n=row['atoms'];dim=3*n;D=units['debye_per_eA'];K=units['electric'];B=units['bohr_A'];H=units['hartree_kcal_mol']
    for k,p in enumerate(parsed):
        A=np.array(p['A']);v=np.array(p['vectors']);ed,ep,u,w=v.T
        alpha=np.repeat([s['alpha_A3'] for s in p['sites']],3);sqrt=np.sqrt(alpha)
        ar=arithmetic(dim,A,ed,ep,u,w);an=math.sqrt(dim)*ar
        sym=float(np.max(abs(A-A.T)));gate(f'{k}_symmetry',sym,ar)
        lam=np.linalg.eigvalsh(A);scaled=np.linalg.eigvalsh(sqrt[:,None]*A*sqrt[None,:])
        require(np.isfinite(lam).all() and lam[0]>0 and scaled[0]>1e-6*scaled[-1]>0,'Full matrix stability failed')
        p['spectrum']={'A':lam.tolist(),'normalized_A':scaled.tolist(),'condition_A':float(lam[-1]/lam[0]),'condition_normalized':float(scaled[-1]/scaled[0])}
        ud=np.linalg.solve(A,ed);wd=np.linalg.solve(A,ep)
        r=A@u-ed;rp=A@w-ep;dr=A@ud-ed;dp=A@wd-ep
        amp=float(np.linalg.norm(np.eye(dim)-p['peek']*A*alpha[None,:],2))
        cg=amp*math.sqrt(n)*(q[k]+.5e-10)/D+an
        gate(f'{k}_Au_Ed',np.linalg.norm(r),cg);gate(f'{k}_Av_Ep',np.linalg.norm(rp),cg)
        bound_u=(np.linalg.norm(r)+np.linalg.norm(dr)+2*an)/lam[0]
        bound_v=(np.linalg.norm(rp)+np.linalg.norm(dp)+2*an)/lam[0]
        gate(f'{k}_native_vs_dense_u',np.linalg.norm(u-ud),bound_u);gate(f'{k}_native_vs_dense_v',np.linalg.norm(w-wd),bound_v)
        gate(f'{k}_equal_dual_RHS',np.max(abs(ed-ep)),ar)
        gate(f'{k}_equal_dual_solution',np.linalg.norm(u-w),bound_u+bound_v+an)
        for j,name in enumerate(('channel_operator','direct_RHS','native_symmetry','native_operator_action')):gate(f'{k}_{name}',p['checks'][j],ar)
        gate(f'{k}_native_residual_diagnostic_u',abs(p['checks'][4]-np.max(abs(r))),ar)
        gate(f'{k}_native_residual_diagnostic_v',abs(p['checks'][5]-np.max(abs(rp))),ar)
        p['numerical_bounds']=dict(arithmetic=ar,arithmetic_norm=an,post_peek_amplification=amp,CG_post_peek_residual_norm=cg,native_vs_dense_u=bound_u,native_vs_dense_v=bound_v)
        if k==0:
            text=verify(new_moments['analyze_output']).read_text()
            token=re.findall(r'^\s*Polarization\s+([-+\d.EeDd]+)\s+\d+\s*$',text,re.M)
            if token:require(len(token)==1,'Polarization component denominator');printed=float(token[0]);print_error=half_print(token[0])
            else:
                require('Polarization' not in new_moments['components'] and not np.any(ed) and not np.any(ep),'Missing nonzero polarization component')
                printed=0.;print_error=0.
            native_E=-.5*K*float(u@ep);dense_E=-.5*K*float(ud@ep);ea=arithmetic(dim,native_E,printed)
            gate('native_polarization_energy',abs(native_E-printed),print_error+ea)
            gate('dense_polarization_energy',abs(dense_E-printed),print_error+ea+.5*K*np.linalg.norm(ep)*bound_u)
            report['polarization_energy']=dict(native_contraction=native_E,dense_contraction=dense_E,printed=printed,print_half_unit=print_error)
        else:
            mapped=np.empty((n,3))
            for pos,s in enumerate(p['sites']):mapped[s['atom']-1]=D*u[3*pos:3*pos+3]
            gate('native_uind_site_print',np.max(abs(mapped-np.array(new_moments['induced_Debye']))),.5e-6+D*ar)
    a,b=parsed;gate('two_solve_matrix',np.max(abs(np.array(a['A'])-b['A'])),arithmetic(dim,a['A'],b['A']))
    gate('two_solve_RHS',np.max(abs(np.array(a['vectors'])[:,:2]-np.array(b['vectors'])[:,:2])),arithmetic(dim,a['vectors'],b['vectors']))
    report['two_solve_exact']={'matrix':a['A']==b['A'],'RHS':np.array_equal(np.array(a['vectors'])[:,:2],np.array(b['vectors'])[:,:2])}
    if field_targets:
        require(set(field_targets)=={'minus_h','minus_half','plus_half','plus_h'},'Exact four old fields')
        require({k:v['field_au'] for k,v in field_targets.items()}=={'minus_h':-.001,'minus_half':-.0005,'plus_half':.0005,'plus_h':.001},'Frozen field sizes/signs')
        A=np.array(b['A']);ed=np.array(b['vectors'])[:,0];lmin=b['spectrum']['A'][0];ar=b['numerical_bounds']['arithmetic'];an=math.sqrt(dim)*ar
        z=np.tile([0.,0.,1.],n);coupling=H/(K*B);column=np.linalg.solve(A,z);chi=column.reshape(n,3).sum(axis=0)*coupling/B
        field_errors={}
        for label,target in field_targets.items():
            f=target['field_au'];moment=target['moments'];predict=np.linalg.solve(A,ed+z*coupling*f)
            oldq=max(induction_values(moment['analyze_output']))
            residual_bound=b['numerical_bounds']['post_peek_amplification']*math.sqrt(n)*(oldq+.5e-10)/D+an
            dense_error=np.linalg.norm(A@predict-ed-z*coupling*f)+an
            bound=(residual_bound+dense_error)/lmin
            mapped=np.empty((n,3))
            for pos,s in enumerate(b['sites']):mapped[s['atom']-1]=D*predict[3*pos:3*pos+3]
            gate('old_field_sites_'+label,np.max(abs(mapped-np.array(moment['induced_Debye']))),.5e-6+D*bound)
            field_errors[label]=bound
            report['field_predictions'][label]={'field_au':f,'site_dipoles_Debye':mapped.tolist(),'solution_error_bound_eA':bound}
        report['susceptibility']={'predicted_column_a0_cubed':chi.tolist(),'steps':{}}
        for tag,h in [('h',.001),('half',.0005)]:
            plus=field_targets['plus_'+tag]['moments'];minus=field_targets['minus_'+tag]['moments']
            observed=(np.array(plus['induced_Debye'])-minus['induced_Debye']).sum(axis=0)/(2*h*D*B)
            bound=n*.5e-6/(h*D*B)+math.sqrt(n)*(field_errors['plus_'+tag]+field_errors['minus_'+tag])/(2*h*B)+arithmetic(dim,chi)
            gate('old_field_chi_'+tag,np.max(abs(chi-observed)),bound)
            report['susceptibility']['steps'][tag]={'old_site_column_a0_cubed':observed.tolist(),'error_bound_a0_cubed':bound}
    report['passed']=all(x['passed'] for x in report['checks'].values());return report

def old_context(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['declaration']==pin(HERE/'DECLARATION.md'),'Scout source binding')
    m=read(cfg['old_manifest']);r=read(cfg['old_model_result']);review=read(cfg['old_model_review'])
    require(review['status']=='passed_independent_actual_model_susceptibility_review' and review['input_result']==cfg['old_model_result'] and review['manifest']==cfg['old_manifest'],'Old actual admission')
    require(r['status']=='conditional_model_susceptibility_complete' and r['native_program_invocations']==19 and r['error'] is None and r['controls_passed'] and r['summary']['numerical_checks_passed'],'Old19 incomplete')
    old=module(m['implementation']);helpers=old.helpers(m['sources']);return m,r,old,helpers

def prepare(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();m,r,old,h=old_context(cfg);row=m['rows']['La_dimer']
    require(row['types']==[500,36,37,37] and row['atoms']==4,'Exact La water1 dimer')
    exporter=read(cfg['exporter_ready']);require(exporter['source_manifest']==cfg['source_manifest'],'Export source join');read(cfg['source_manifest']);verify(exporter['helper']);verify(exporter['patch'])
    zero=old.parse(verify(row['baseline']['analyze_output']),row,0.,h)
    require(parity(zero,row['baseline'],exact=True)['passed'],'Old zero raw replay')
    fields={}
    for tag in ('minus_h','minus_half','plus_half','plus_h'):
        rec=r['rows']['La_dimer__'+tag];q=old.parse(verify(rec['moments']['analyze_output']),row,rec['job']['field_au'][2],h)
        require(q==rec['moments'],'Old field raw replay differs');fields[tag]={'field_au':rec['job']['field_au'][2],'moments':q}
    key=verify(row['input']['key']).read_text();tokens=[s.split('#',1)[0].split() for s in key.splitlines()]
    require(not any(t and t[0].lower() in ('response-matrix-export','external-field','exfld-freq','polar-pair-response','exchpol','exchange-polar') for t in tokens),'Unexpected baseline keyword')
    jobs={}
    for state in ('off','on'):
        folder=out/'inputs'/state;folder.mkdir(parents=True)
        for k,name in [('xyz','system.xyz'),('parameters','amoeba09.prm')]:shutil.copyfile(verify(row['input'][k]),folder/name)
        (folder/'system.key').write_text(key if state=='off' else key+'\n'+FLAG+'\n')
        jobs[state]={'input':{k:pin(folder/name) for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]},'export_enabled':state=='on'}
    result={'status':'prepared_unbound_binary','implementation':pin(__file__),'configuration':cfg,'declaration':cfg['declaration'],'row':row,'units':m['conversion'],'old_zero':zero,'old_fields':fields,'jobs':jobs,'maximum_native_programs':2,'workers':1,'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu}}
    save(out/'MANIFEST.json',result)
    template={'implementation':pin(__file__),'declaration':cfg['declaration'],'manifest':pin(out/'MANIFEST.json'),'source_manifest':cfg['source_manifest'],'analyze':None,'build_receipt':None,'root_release':None,'maximum_native_programs':2,'workers':1}
    save(out/'CONFIG_RUN_TEMPLATE.json',template);return result

def context(cfg,runtime=False):
    require(cfg['implementation']==pin(__file__) and cfg['declaration']==pin(HERE/'DECLARATION.md'),'Runtime source')
    m=read(cfg['manifest']);require(m['implementation']==cfg['implementation'] and m['new_native_calls']==0 and m['configuration']['source_manifest']==cfg['source_manifest'],'Prepared source/input join')
    _,_,old,h=old_context(m['configuration']);read(cfg['source_manifest'])
    require(cfg['analyze'] is not None and cfg['build_receipt'] is not None and cfg['root_release'] is not None,'Root-bound reviewed executable required')
    verify(cfg['analyze']);build=read(cfg['build_receipt'])
    require(build['status']=='complete' and build['manifest']==cfg['source_manifest'] and build['executables']['analyze']==cfg['analyze'],'Completed build/source/executable join')
    rel=read(cfg['root_release'])
    require(rel['status']=='approved_response_matrix_scout' and all(rel[k]==cfg[k] for k in ('implementation','manifest','source_manifest','analyze','build_receipt')),'Root release join')
    require(rel['maximum_native_programs']==cfg['maximum_native_programs']==2 and rel['workers']==cfg['workers']==1 and rel['force_qualification'] is False,'Finite release scope')
    require(set(m['jobs'])=={'off','on'} and m['maximum_native_programs']==2,'Two prepared slots')
    key=verify(m['row']['input']['key']).read_text()
    for state,job in m['jobs'].items():
        require(verify(job['input']['key']).read_text()==(key if state=='off' else key+'\n'+FLAG+'\n'),'Only reporting keyword allowed')
        for k in ('xyz','parameters'):require(verify(job['input'][k]).read_bytes()==verify(m['row']['input'][k]).read_bytes(),'Changed geometry/parameters')
    if runtime:require(os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_NTASKS','1')=='1' and os.environ.get('SLURM_CPUS_PER_TASK')=='1' and all(os.environ.get(k)=='1' for k in THREADS),'One Slurm CPU/task/thread')
    return m,old,h,rel

def run(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'status':'running','implementation':pin(__file__),'configuration':cfg,'rows':{s:{'status':'not_attempted'} for s in ('off','on')},'native_program_invocations':0,'error':None,'qualification':{},'new_reference_calls':0}
    save(out/'RESULT.json',result)
    try:
        m,old,h,rel=context(cfg,True);require(str(out.resolve())==rel['execution_root'],'Exact root execution location')
        for state in ('off','on'):
            rec=result['rows'][state];rec.update(status='started',receipts=[]);folder=out/state;folder.mkdir(exist_ok=False)
            for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:shutil.copyfile(verify(m['jobs'][state]['input'][k]),folder/name)
            rec['runtime_input']={k:pin(folder/name) for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]}
            save(folder/'START.json',{'implementation':cfg['implementation'],'configuration':cfg,'state':state,'runtime_input':rec['runtime_input'],'job_id':os.environ['SLURM_JOB_ID']});rec['start']=pin(folder/'START.json')
            require(result['native_program_invocations']<2,'Two-call cap');result['native_program_invocations']+=1;save(out/'RESULT.json',result)
            h[2].invoke([cfg['analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',rec['receipts'])
            rec['moments']=old.parse(folder/'analyze.out',m['row'],0.,h)
            names=re.findall(r'Response Matrix Written To\s*:\s*(\S+)',(folder/'analyze.out').read_text())
            exports=sorted(p.name for p in folder.glob('response_matrix_*.dat'))
            if state=='off':
                require(not names and not exports,'Disabled export produced files')
                rec['parity']=parity(rec['moments'],m['old_zero']);require(rec['parity']['passed'],'Rebuilt baseline parity failed')
            else:
                require(names==list(NAMES) and exports==list(NAMES),'Two announced export files required')
                rec['exports']=[pin(folder/n) for n in NAMES]
                rec['parity']=parity(rec['moments'],result['rows']['off']['moments'],exact=True);require(rec['parity']['passed'],'Reporting changed native outputs')
                qualify_exports(rec['exports'],m['row'],m['units'],rec['moments'],m['old_fields'],result['qualification'])
            rec['status']='complete';save(folder/'RESULT.json',rec);save(out/'RESULT.json',result)
        result['status']='response_matrix_scout_complete'
    except Exception as exc:
        result['status']='failed_or_unqualified';result['error']=repr(exc)
        for rec in result['rows'].values():
            if rec['status']=='started':rec['status']='failed';rec['error']=repr(exc)
    now=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-tick,'parent_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':now.ru_utime-child.ru_utime,'child_system_CPU_seconds':now.ru_stime-child.ru_stime,'workers':1}
    save(out/'RESULT.json',result);return result

def admit_result(cfg,manifest,result):
    """Saved scout replay for the separately gated remaining fixture runner."""
    m,old,h,rel=context(cfg)
    require(m==manifest and result['configuration']==cfg and result['implementation']==cfg['implementation'],'Saved scout source/config')
    require(result['status']=='response_matrix_scout_complete' and result['error'] is None and result['native_program_invocations']==2 and set(result['rows'])=={'off','on'},'Complete two-call scout required')
    replay={}
    for state in ('off','on'):
        rec=result['rows'][state];require(rec['status']=='complete' and len(rec['receipts'])==1,'Saved native slot unavailable')
        folder=verify(rec['runtime_input']['xyz']).parent
        require(folder.resolve()==Path(rel['execution_root']).resolve()/state,'Saved execution path differs')
        for key in ('xyz','key','parameters'):require(verify(rec['runtime_input'][key]).read_bytes()==verify(m['jobs'][state]['input'][key]).read_bytes(),'Saved native physical input differs')
        st=read(rec['start']);require(st['configuration']==cfg and st['state']==state and st['runtime_input']==rec['runtime_input'],'Saved START join')
        rr=rec['receipts'][0]
        require(rr['argv']==[cfg['analyze']['path'],'system.xyz','EM'] and rr['cwd']==str(folder) and rr['returncode']==0 and rr['output']==pin(folder/'analyze.out') and set(rr['threads'])==set(THREADS) and all(rr['threads'][k]=='1' for k in THREADS),'Saved native receipt')
        q=old.parse(folder/'analyze.out',m['row'],0.,h);require(q==rec['moments'],'Saved native moments replay')
        baseline=m['old_zero'] if state=='off' else replay['off']
        require(parity(q,baseline,exact=state=='on')==rec['parity'] and rec['parity']['passed'],'Saved parity differs')
        names=re.findall(r'Response Matrix Written To\s*:\s*(\S+)',(folder/'analyze.out').read_text())
        files=sorted(p.name for p in folder.glob('response_matrix_*.dat'))
        require((not names and not files) if state=='off' else names==list(NAMES) and files==list(NAMES),'Saved export denominator')
        if state=='on':require(rec['exports']==[pin(folder/n) for n in NAMES],'Saved raw export pins')
        replay[state]=q
    numerical=qualify_exports(result['rows']['on']['exports'],m['row'],m['units'],replay['on'],m['old_fields'])
    require(numerical['passed'] and result['qualification']['passed'] and set(numerical['checks'])==set(result['qualification']['checks']),'Saved numerical qualification failed')
    for a,b in zip(numerical['exports'],result['qualification']['exports']):
        require(all(a[k]==b[k] for k in ('source','sequence','sites','vectors','A','checks','restored','peek')),'Stored raw matrix differs')
    return {'status':'admitted_response_matrix_scout','on_row':result['rows']['on'],'qualification':numerical,'new_native_calls':0}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=('prepare','run'));ap.add_argument('--config',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    try:
        cfg=json.loads(args.config.read_text());r=prepare(cfg,args.out) if args.action=='prepare' else run(cfg,args.out)
    except Exception as exc:
        saved=json.loads((args.out/'RESULT.json').read_text()) if (args.out/'RESULT.json').exists() else {}
        save(args.out/'FAILURE.json',{'error':repr(exc),'native_calls_attempted':saved.get('native_program_invocations',0),'configuration':pin(args.config) if args.config.exists() else None});raise
    print(json.dumps({'status':r['status'],'output':str(args.out.resolve())}));return 0 if r['status'] in ('prepared_unbound_binary','response_matrix_scout_complete') else 1

if __name__=='__main__':raise SystemExit(main())
