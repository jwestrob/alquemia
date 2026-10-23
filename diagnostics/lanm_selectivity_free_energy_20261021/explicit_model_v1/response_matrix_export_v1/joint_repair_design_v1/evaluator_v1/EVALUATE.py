#!/usr/bin/env python3
"""Finite joint evaluator: exact fixture/cache joins and existing native parsers."""
import argparse
import concurrent.futures
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import threading
import time
import numpy as np

HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
HELPERS = ('objective','old_objective','body_helper','moment_adapter','induction_helper','native_helper','resolver','scout_helper')
SOURCES = HELPERS + ('domain','fit_protocol','declaration','target_inventory','targets','development_data',
                    'fixture_inventory','baseline_result','baseline_review','baseline_manifest','moment_preparation')
MAX_VECTORS, MAX_PROGRAMS, WORKERS = 673, 24000, 12
SCOUT = 'aquo9_central_La__dimer_01'


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(q):
    require(pin(q['path']) == q, 'Changed pin: '+q['path'])
    return Path(q['path'])


def read(q): return json.loads(verify(q).read_text())


def save(path, value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)


def digest(value): return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def fingerprint(cfg): return digest({k:v for k,v in cfg.items() if k not in ('root_release','preparation_result')})


def module(q, name):
    spec=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def checked(values, domain):
    values=list(values);require(len(values)==12,'Exactly12 box parameters')
    require(all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) and p['lower']<=v<=p['upper'] for v,p in zip(values,domain['parameters'])),'Parameter outside declared finite box')
    return [0. if v==0 else float(v) for v in values]


def effective(values, metal, domain):
    v=checked(values,domain);i=0 if metal=='La' else 4;j=8 if metal=='La' else 10;S=domain['S']
    a,b,u,w=v[i:i+4]
    return [a,b,S*(u+w)/4,S*(u-w)/4,v[j],v[j+1]]


def render(values, row, domain, resolver, exports=None):
    """Change one active metal-O VDWPR and add one physical response record."""
    vals=effective(values,row['metal'],domain);mt=500 if row['metal']=='La' else 501
    original=verify(row['source_input']['key']).read_text();prm=verify(row['source_input']['parameters']).read_text()
    lines=original.splitlines(keepends=True);out=[];matches=0;names=[]
    for line in lines:
        t=line.split('#',1)[0].split()
        if not t:out.append(line);continue
        name=t[0].lower();names.append(name)
        require(name not in ('polar-local-response','polar-pair-response','external-field','exfld-freq','response-matrix-export','exchpol','exchange-polar'),'Unexpected response/field keyword in baseline')
        if name=='vdwpr' and sorted(map(int,t[1:3]))==[34,mt]:
            require(len(t)==5,'Ambiguous metal-O VDWPR');matches+=1
            out.append(f'vdwpr {t[1]} {t[2]} {vals[4]:.17g} {vals[5]:.17g}\n')
        else:out.append(line)
    require(matches==1,'Exactly one selected metal-O pair override required')
    for name in ('polar-print','save-uinduce'):
        require(names.count(name)<=1,'Duplicate reporting flag')
        if name not in names:out.append('\n'+name+'\n')
    out.append(f'\nPOLAR-LOCAL-RESPONSE {mt} 36 '+' '.join(f'{x:.17g}' for x in vals[:4])+'\n')
    enabled=row['needs_export'] if exports is None else exports
    if enabled:out.append('RESPONSE-MATRIX-EXPORT\n')
    key=''.join(out);before=resolver.parse_parameters(prm,original);after=resolver.parse_parameters(prm,key)
    for name in before:
        if name!='vdwpr':require(before[name]==after[name],'Non-VDW physical parameter changed: '+name)
    allowed=(34,mt)
    require({k:v for k,v in before['vdwpr'].items() if k!=allowed}=={k:v for k,v in after['vdwpr'].items() if k!=allowed},'Nonwater VDW pair changed')
    pair=resolver.resolve_pair(after,mt,36)
    require([pair['Rmin_pair_A'],pair['epsilon_kcal_mol'],pair['thole']]==[vals[4],vals[5],domain['fixed_native_POLPAIR'][row['metal']]],'Resolved pair differs')
    return key,{'metal':row['metal'],'physical_values_a_b_c_d_R_epsilon':vals,'resolved_pair':pair,'export_enabled':enabled,'POLPAIR_unchanged':True}


def local_readback(text, row, physical, units, solves=2):
    """Admit the finite native reporting grammar; do not recalculate the functional."""
    lines=[s.split() for s in text.splitlines() if s.strip().startswith('LOCAL_RESPONSE_')]
    if not any(physical[:4]):
        require(not lines,'Zero response must retain native bypass');return None
    names=('PARAMETERS','PAIR','ENVIRONMENT','SPECTRUM','INDUCED','DUAL')
    groups={name:[] for name in names}
    for tokens in lines:
        name=tokens[0].removeprefix('LOCAL_RESPONSE_');require(name in groups,'Unsupported local response readback')
        groups[name].append(tokens[1:])
    n=row['atoms'];mt=500 if row['metal']=='La' else 501
    metals=[i+1 for i,t in enumerate(row['types']) if t==mt];oxygen=[i+1 for i,t in enumerate(row['types']) if t==36]
    require(len(metals)==1,'Exactly one represented metal');metal=metals[0]
    expected={'PARAMETERS':solves,'PAIR':solves*len(oxygen),'ENVIRONMENT':solves,'SPECTRUM':solves,'INDUCED':solves*n,'DUAL':solves}
    require({k:len(v) for k,v in groups.items()}==expected,'Local response reporting denominator')
    def nums(t,count):
        require(len(t)==count,'Local response numeric grammar')
        a=np.asarray([float(x.replace('D','E')) for x in t]);require(np.isfinite(a).all(),'Nonfinite local response readback');return a
    def close(a,b):
        a,b=np.asarray(a),np.asarray(b)
        return bool(np.max(abs(a-b))<=128*3*n*np.finfo(float).eps*max(1.,float(np.max(abs(a))),float(np.max(abs(b)))))
    for t in groups['PARAMETERS']:
        require(len(t)==6 and list(map(int,t[:2]))==[mt,36] and close(nums(t[2:],4),physical[:4]),'Actual physical response coefficients differ')
    for k,t in enumerate(groups['PAIR']):
        oi=oxygen[k%len(oxygen)];require(len(t)==5 and list(map(int,t[:2]))==[metal,oi],'Literal native metal-O map')
        gamma,kappa,r=nums(t[2:],3);distance=np.linalg.norm(np.array(row['positions_A'][metal-1])-row['positions_A'][oi-1])
        expected_kappa=row['fixed_POLPAIR']/math.sqrt(row['alpha_A3'][metal-1]*row['alpha_A3'][oi-1])
        require(gamma>0 and kappa>0 and close([gamma,kappa,r],[row['fixed_POLPAIR'],expected_kappa,distance]),'Native pair gamma/damping length/geometry changed')
    for t in groups['ENVIRONMENT']:
        require(len(t)==5 and int(t[0])==metal,'Native environment metal');nums(t[1:],4)
    spectra=[]
    for t in groups['SPECTRUM']:
        lo,hi,sym,resid,orth=nums(t,5)
        require(hi>0 and lo>1e-6*hi and 0<=sym<=1e-10 and 0<=resid<=1e-10*max(1.,hi) and 0<=orth<=1e-10,'Native response spectrum/stopping gate')
        spectra.append([float(x) for x in (lo,hi,sym,resid,orth)])
    arrays=[]
    for solve in range(solves):
        uv=[]
        for i,t in enumerate(groups['INDUCED'][solve*n:(solve+1)*n]):
            require(len(t)==7 and int(t[0])==i+1,'Native full induced atom order');uv.append(nums(t[1:],6).tolist())
        uv=np.asarray(uv);dual=nums(groups['DUAL'][solve],2)
        require((dual>=0).all() and close([np.max(abs(uv[:,:3]-uv[:,3:]))],[dual[1]]),'Native dual readback closure')
        arrays.append(uv.tolist())
    return {'parameters':physical[:4],'native_rows':groups,'spectra':spectra,'induced_u_v_eA':arrays,
            'conjugate_induced_Debye':(units['debye_per_eA']*.5*(np.asarray(arrays[-1])[:,:3]+np.asarray(arrays[-1])[:,3:])).tolist(),
            'dual_RHS_and_solution_are_retained_without_assuming_exact_equality':True}


def prepare(cfg, out):
    tick,cpu=time.monotonic(),time.process_time()
    require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Evaluator source/spec')
    for k in SOURCES:verify(cfg[k])
    domain=read(cfg['domain']);inv=read(cfg['fixture_inventory']);data=read(cfg['development_data']);target=read(cfg['target_inventory'])
    full=read(cfg['baseline_result']);review=read(cfg['baseline_review']);rm=read(cfg['baseline_manifest']);mp=read(cfg['moment_preparation'])
    require(review['status']=='passed_baseline_matrix_inventory_review' and review['result']==cfg['baseline_result'] and review['manifest']==cfg['baseline_manifest']==full['manifest'],'Actual54 baseline admission')
    require(full['status']=='baseline_matrix_inventory_complete' and full['complete_fixtures']==54 and full['admitted_snapshots']==108 and full['error'] is None and full['global_stop'] is False,'Complete admitted54 baseline')
    require(rm['inventory']==mp['inventory']==cfg['fixture_inventory'] and mp['targets']==cfg['target_inventory'] and rm['units']==mp['units'],'Baseline inventory/units/frame lineage')
    require(rm['scout_helper']==cfg['scout_helper'] and all(rm['helpers'][k]==cfg[k] for k in ('moment_adapter','induction_helper','native_helper','resolver')),'Qualified native helper identity')
    require(domain['starting_candidate']==inv['candidate_result']==target['fixed_model_candidate'] and data['baseline_parameters']==target['fixed_six_parameters'],'Selected17676 lineage')
    require(domain['root_declaration']==cfg['declaration'] and len(domain['parameters'])==12,'Declared domain lineage')
    h={k:module(cfg[k],'joint_prepare_'+k) for k in HELPERS}
    require(h['objective'].targets_from_inventory(target)==read(cfg['targets']),'Prepared objective targets unchanged')
    initial=checked([p['initial'] for p in domain['parameters']],domain)
    aliases=inv['field_aliases'];pairs=target['groups']['pair_and_collective_dipoles']['pair_cycles']
    pair_ids=[next(cid for cid,c in r['model_terms'].items() if c==1) for r in pairs]
    require(len(pair_ids)==34 and all('__dimer_' in cid for cid in pair_ids),'Literal34 pair fixture map')
    ids=list(data['rows'])+pair_ids
    ids += [aliases[s]['fixture'] for s in ('La_trimer','Dy_trimer')]
    require(len(ids)==len(set(ids))==46 and len(data['rows'])==10,'Exact46 objective fixtures')
    def native_record(cid):
        r=full['rows'][cid];require(r['status']=='complete','Missing baseline fixture')
        return r['admission']['on_row'] if cid==SCOUT else r
    rows={};baselines={}
    for cid in ids:
        row=copy.deepcopy(rm['rows'][cid]);support=mp['rows'][cid]
        require(row['source_input']==support['input'] and row['types']==support['types'] and row['positions_A']==support['positions_A'],'Literal fixture/frame map')
        row.update(metal=inv['rows'][cid]['metal'],needs_gradient=cid in data['rows'],needs_export=cid in {x['fixture'] for x in aliases.values()},permanent_lab_Debye=support['permanent_lab_Debye'],energy_replay_kcal_mol=None)
        row['fixed_POLPAIR']=domain['fixed_native_POLPAIR'][row['metal']]
        old=native_record(cid);mom=copy.deepcopy(old['moments'])
        for k in ('analyze_output','induced_output'):verify(mom[k])
        mu=[row['permanent_lab_Debye'][j]+mom['induced_sum_Debye'][j] for j in range(3)]
        require(max(abs(mu[j]-mom['total_lab_Debye'][j]) for j in range(3))<=.0005+row['atoms']*.5e-6+2e-9,'Baseline laboratory moment check')
        base={'status':'complete','inherited_native':True,'energy':mom['energy_kcal_mol'],'moments':mom,'mu_Debye':mu,
              'producer_result':cfg['baseline_result'],'producer_cell_id':cid,'producer_input':row['source_input'],'producer_review':cfg['baseline_review']}
        if row['needs_gradient']:
            model=copy.deepcopy(data['baseline_model'][cid]);require(model['energy']==base['energy'],'Baseline EF/moment energy join')
            for k in ('analyze_output','gradient_output'):verify(model['native'][k])
            base['model']=model
        if row['needs_export']:
            state=next(s for s,v in aliases.items() if v['fixture']==cid);q=review['field_alias_consistency'][state]
            require(q['fixture']==cid and q['passed'],'Qualified four-state susceptibility alias')
            base['chi_a0_cubed']=q['susceptibility']['predicted_column_a0_cubed'][2];base['chi_authority']={'review':cfg['baseline_review'],'state':state}
        rows[cid]=row;baselines[cid]=base
        key,meta=render(initial,row,domain,h['resolver']);folder=out/'inputs'/cid;folder.mkdir(parents=True)
        for k,name in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['source_input'][k]),folder/name)
        (folder/'system.key').write_text(key);row['prepared_initial_input']={k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))};row['initial_render']=meta
    static=copy.deepcopy(mp['static'])
    for cid in ('ion_La','ion_Dy'):
        old=native_record(cid);mom=old['moments'];support=mp['rows'][cid]
        static[cid]={'energy_kcal_mol':mom['energy_kcal_mol'],'mu_Debye':[support['permanent_lab_Debye'][j]+mom['induced_sum_Debye'][j] for j in range(3)],'source':{'result':cfg['baseline_result'],'cell':cid},'input':rm['rows'][cid]['source_input']}
    require(len(static)==48 and sum(r['needs_gradient'] for r in rows.values())==10 and sum(r['needs_export'] for r in rows.values())==4,'46/10/4 plus48 exact cache')
    m={'status':'prepared_joint_evaluator_unbound_native','implementation':cfg['implementation'],'specification':cfg['specification'],'sources':{k:cfg[k] for k in SOURCES},
       'rows':rows,'baseline_cells':baselines,'static':static,'field_aliases':aliases,'units':rm['units'],'initial_vector':initial,
       'counts':{'objective_cells':46,'EF_cells':10,'chi_exports':4,'static_metal_free':46,'static_ions':2,'programs_per_full_changed_vector':56,'programs_per_changed_metal':28,'excluded_other_trimers':6},
       'limits':{'max_parameter_vectors':MAX_VECTORS,'max_native_programs':MAX_PROGRAMS,'workers':WORKERS},'new_native_calls':0,
       'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    save(out/'MANIFEST.json',m);return m


class Evaluator:
    def __init__(self, config, cache_dir, root_release=None):
        self.config=copy.deepcopy(config);self.fingerprint=fingerprint(config);self.root_release=root_release or config.get('root_release')
        require(config['implementation']==pin(__file__) and config['specification']==pin(HERE/'SPEC.md'),'Evaluator source changed')
        self.m=read(config['manifest']);require(self.m['implementation']==config['implementation'] and self.m['specification']==config['specification'],'Prepared evaluator source differs')
        require(self.m['sources']=={k:config[k] for k in SOURCES},'Prepared scientific source bindings differ')
        require(config['max_parameter_vectors']==MAX_VECTORS and config['max_native_programs']==MAX_PROGRAMS and config['workers']==WORKERS,'Declared finite caps differ')
        require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread per process')
        self.h={k:module(config[k],'joint_runtime_'+k) for k in HELPERS};self.domain=read(config['domain'])
        self.data=read(config['development_data']);self.targets=read(config['targets']);self.target_inventory=read(config['target_inventory'])
        self.initial_vector=list(self.m['initial_vector']);self.family={'parameters':copy.deepcopy(self.domain['parameters']),'initial_vector':self.initial_vector[:],'domain':config['domain']}
        self.cache_dir=Path(cache_dir).resolve();self.cache_dir.mkdir(parents=True,exist_ok=True)
        self.lock=threading.RLock();self.evaluation_lock=threading.Lock();self.stop_event=threading.Event();self.slots=threading.BoundedSemaphore(WORKERS)
        self.cell_locks={};self.cache={};self.rendered={};self._gated=False;self.ledger_path=self.cache_dir/'LEDGER.json'
        if self.ledger_path.exists():
            self.ledger=json.loads(self.ledger_path.read_text());require(self.ledger['config_sha256']==self.fingerprint,'Ledger configuration mismatch')
            if self.ledger['stopped']:self.stop_event.set()
        else:self.ledger={'config_sha256':self.fingerprint,'vector_keys':[],'native_program_invocations':0,'requests':[],'stopped':False};save(self.ledger_path,self.ledger)

    def checked_values(self, values): return checked(values,self.domain)

    def render(self, values, cid, exports=None):
        row=self.m['rows'][cid];enabled=row['needs_export'] if exports is None else exports
        # The46 geometries share four literal baseline keys. Resolve identical
        # rendered bytes once, retaining the same validation and returned data.
        for k in ('key','parameters'):verify(row['source_input'][k])
        token=(row['source_input']['key']['sha256'],row['source_input']['parameters']['sha256'],row['metal'],
               tuple(effective(values,row['metal'],self.domain)),enabled)
        with self.lock:
            if token not in self.rendered:self.rendered[token]=render(values,row,self.domain,self.h['resolver'],enabled)
            return copy.deepcopy(self.rendered[token])

    def cache_identity(self, values, cid):
        key,meta=self.render(values,cid);row=self.m['rows'][cid]
        q={'cell_id':cid,'effective_values':meta['physical_values_a_b_c_d_R_epsilon'],'xyz':row['source_input']['xyz'],'parameters':row['source_input']['parameters'],
           'key_sha256':hashlib.sha256(key.encode()).hexdigest(),'implementation':self.config['implementation'],
           'source_manifest':self.config.get('source_manifest'),'executables':self.config.get('executables'),'helpers':{k:self.config[k] for k in HELPERS}}
        return digest(q),q,key,meta

    def _stop(self, exc):
        self.stop_event.set()
        with self.lock:self.ledger.update(stopped=True,failure=repr(exc));save(self.ledger_path,self.ledger)

    def _gate(self):
        require(not self.stop_event.is_set(),'Global native stop active')
        require(os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_NTASKS','1')=='1' and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==12,'Exactly one12-CPU Slurm task')
        require(all(os.environ.get(k)=='1' for k in THREADS),'One native numerical thread')
        if self._gated:return
        cfg=self.config;require(self.root_release is not None and all(isinstance(cfg.get(k),dict) for k in ('source_manifest','build_receipt','executables','qualification_review')),'Actual reviewed native bindings unavailable')
        require(set(cfg['executables'])=={'analyze','testgrad'},'Exact two native executables')
        for q in cfg['executables'].values():verify(q)
        build=read(cfg['build_receipt']);qualification=read(cfg['qualification_review']);release=read(self.root_release)
        verify(qualification['implementation']);verify(qualification['result'])
        require(build['status']=='complete' and build['manifest']==cfg['source_manifest'] and all(build['executables'][k]==v for k,v in cfg['executables'].items()),'Completed build/source/executable joins')
        verify(cfg['source_manifest'])
        require(qualification['status']=='passed_local_response_native_qualification' and all(qualification[k]==cfg[k] for k in ('source_manifest','build_receipt','executables')),'Actual native qualification joins')
        require(all(qualification[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Native qualification incomplete')
        require(release['status']=='released_bounded_joint_response_fit' and release['evaluator']==cfg['implementation'] and release['evaluator_config_sha256']==self.fingerprint,'Root release/evaluator binding')
        require(all(release[k]==cfg[k] for k in ('manifest','domain','fit_protocol','source_manifest','build_receipt','executables','qualification_review')),'Root scientific/native bindings')
        require(release['limits']=={'outer_steps':24,'new_vectors':672,'native_programs':MAX_PROGRAMS,'workers':WORKERS} and release['cache_directory']==str(self.cache_dir),'Root finite allocation/cache binding')
        require(release['reserved_outcomes_authorized'] is False,'Reserved outcomes remain excluded')
        self._gated=True

    def _invoke(self, program, folder, record):
        with self.slots:
            with self.lock:
                self._gate();require(self.ledger['native_program_invocations']<MAX_PROGRAMS,'Native program hard ceiling')
                self.ledger['native_program_invocations']+=1;save(self.ledger_path,self.ledger)
            args=['system.xyz','EM'] if program=='analyze' else ['system.xyz','Y','N']
            try:self.h['native_helper'].invoke([self.config['executables'][program]['path']]+args,folder,folder/(program+'.out'),record['receipts'])
            except Exception as exc:self._stop(exc);raise

    def _parse(self, record, row):
        h=self.h;a=verify(record['program_inputs']['analyze']['xyz']).parent/'analyze.out';runtime=copy.deepcopy(row);runtime['input']=record['program_inputs']['analyze']
        mom=h['moment_adapter'].parse_moments(a,h['moment_adapter'].induced_file(a),runtime,h['induction_helper'].induction)
        local=local_readback(a.read_text(),row,record['render']['physical_values_a_b_c_d_R_epsilon'],self.m['units'])
        if local is not None:
            uv=np.asarray(local['induced_u_v_eA'][-1]);D=self.m['units']['debye_per_eA']
            require(np.max(abs(D*uv[:,:3]-mom['induced_Debye']))<=.5e-6+self.h['scout_helper'].arithmetic(3*row['atoms'],uv)*D,'Full induced readback/.uind closure')
            induced=np.sum(local['conjugate_induced_Debye'],axis=0)
        else:induced=mom['induced_sum_Debye']
        mu=[float(row['permanent_lab_Debye'][j]+induced[j]) for j in range(3)]
        require(max(abs(mu[j]-mom['total_lab_Debye'][j]) for j in range(3))<=.0005+row['atoms']*.5e-6+2e-9,'Conjugate native/common-origin moment consistency')
        parsed={'energy':mom['energy_kcal_mol'],'moments':mom,'mu_Debye':mu,'local_response_readback':local}
        if row['needs_export']:
            paths=sorted(a.parent.glob('response_matrix_*.dat'))
            require([p.name for p in paths]==['response_matrix_001.dat','response_matrix_002.dat'],'Exact two analyze exports')
            require(re.findall(r'Response Matrix Written To\s*:\s*(\S+)',a.read_text())==[p.name for p in paths],'Exact export announcements')
            q=h['scout_helper'].qualify_exports(paths,runtime,self.m['units'],mom)
            require(q['passed'],'Candidate matrix qualification')
            A=np.asarray(q['exports'][1]['A']);z=np.tile([0.,0.,1.],row['atoms']);x=np.linalg.solve(A,z)
            residual=float(np.max(abs(A@x-z)));bound=h['scout_helper'].arithmetic(3*row['atoms'],A,x)*max(1.,float(np.max(abs(A))),float(np.max(abs(x))))
            require(residual<=bound,'Susceptibility solve residual')
            units=self.m['units'];chi=units['hartree_kcal_mol']/(units['electric']*units['bohr_A']**2)*float(x.reshape(row['atoms'],3).sum(axis=0)[2])
            require(math.isfinite(chi),'Finite candidate susceptibility')
            parsed.update(chi_a0_cubed=chi,exports=[pin(p) for p in paths],export_qualification=q,chi_solve_residual=residual)
        if row['needs_gradient'] and 'testgrad' in record['program_inputs']:
            g=verify(record['program_inputs']['testgrad']['xyz']).parent/'testgrad.out'
            h['induction_helper'].induction(g)
            native=h['native_helper'].parse(a,g,{'atoms':row['atoms'],'derivative_test':False})
            parsed['model']={'energy':native['total_energy_kcal_mol'],'forces':native['force_kcal_mol_A'],'native':native}
        return parsed

    def _cell(self, values, cid):
        token,identity,key,meta=self.cache_identity(values,cid);row=self.m['rows'][cid]
        with self.lock:lock=self.cell_locks.setdefault(token,threading.Lock())
        with lock:
            if token in self.cache:record=self.cache[token]
            elif effective(values,row['metal'],self.domain)==effective(self.initial_vector,row['metal'],self.domain):
                record=copy.deepcopy(self.m['baseline_cells'][cid]);record.update(cache_identity=identity,request_render=meta,cache_key=token)
                self.cache[token]=record
            else:
                directory=self.cache_dir/'cells'/token;rp=directory/'RESULT.json'
                if rp.exists():
                    record=read(pin(rp));require(record['status']=='complete' and record['cache_identity']==identity,'Failed/incomplete/different cached cell; no retry')
                    require(len(record['receipts'])==1+int(row['needs_gradient']) and all(r['returncode']==0 for r in record['receipts']),'Cached receipt denominator')
                    self._saved_inputs(record,row,values,directory)
                    parsed=self._parse(record,row)
                    require(all(record[k]==v for k,v in parsed.items()),'Cached raw admission changed')
                    self.cache[token]=record
                else:
                    self._gate();directory.mkdir(parents=True,exist_ok=False)
                    record={'status':'started','cache_key':token,'cache_identity':identity,'cell_id':cid,'first_requested_vector':values[:],
                            'inherited_native':False,'implementation':self.config['implementation'],'render':meta,'receipts':[],
                            'program_inputs':{},'job_id':os.environ['SLURM_JOB_ID'],'root_release':self.root_release}
                    try:
                        for program in (('analyze','testgrad') if row['needs_gradient'] else ('analyze',)):
                            require(not self.stop_event.is_set(),'Global native stop active')
                            folder=directory/program;folder.mkdir();text=key if program=='analyze' else self.render(values,cid,exports=False)[0]
                            for k,n in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['source_input'][k]),folder/n)
                            (folder/'system.key').write_text(text);record['program_inputs'][program]={k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
                            save(rp,record);self._invoke(program,folder,record)
                            # Analyze admission happens before admitting a testgrad start.
                            try:parsed=self._parse(record,row)
                            except Exception as exc:self._stop(exc);raise
                            record.update(parsed);save(rp,record)
                        require(not row['needs_gradient'] or 'model' in record,'Ten native EF cells required')
                        record['status']='complete';self.cache[token]=record
                    except Exception as exc:self._stop(exc);record.update(status='failed',error=repr(exc));raise
                    finally:save(rp,record)
            result=copy.deepcopy(record);result['request_provenance']={'vector':values[:],'cache_key':token,'cache_identity':identity,'render':meta}
            return result

    def _saved_inputs(self, record, row, values, directory):
        programs=['analyze']+(['testgrad'] if row['needs_gradient'] else [])
        require(set(record['program_inputs'])==set(programs),'Saved program input denominator')
        require(record['implementation']==self.config['implementation'] and record['render']==self.render(values,record['cell_id'])[1],'Saved producer/render identity')
        for program,receipt in zip(programs,record['receipts']):
            folder=directory/program;inp=record['program_inputs'][program]
            expected_key=self.render(values,record['cell_id'],exports=row['needs_export'] if program=='analyze' else False)[0]
            require(verify(inp['key'])==folder/'system.key' and verify(inp['key']).read_text()==expected_key,'Saved exact rendered key')
            for k,name in (('xyz','system.xyz'),('parameters','amoeba09.prm')):
                require(verify(inp[k])==folder/name and inp[k]['sha256']==row['source_input'][k]['sha256'],'Saved literal geometry/parameters')
            args=['system.xyz','EM'] if program=='analyze' else ['system.xyz','Y','N']
            require(receipt['argv']==[self.config['executables'][program]['path']]+args and receipt['cwd']==str(folder) and receipt['returncode']==0 and receipt['output']==pin(folder/(program+'.out')),'Saved native argv/cwd/output')
            require(set(receipt['threads'])==set(THREADS) and all(receipt['threads'][k]=='1' for k in THREADS),'Saved native thread allocation')

    def _assemble(self, cells):
        require(set(cells)==set(self.m['rows']) and all(v['status']=='complete' for v in cells.values()),'Full46 objective denominator')
        values={cid:v['mu_Debye'] for cid,v in self.m['static'].items()};values.update({cid:v['mu_Debye'] for cid,v in cells.items()})
        group=self.target_inventory['groups']['pair_and_collective_dipoles']
        pair={r['id']:[math.fsum(c*values[cid][j] for cid,c in r['model_terms'].items()) for j in range(3)] for r in group['pair_cycles']}
        C={p:[values[p][j]-values['ion_'+p.rsplit('_',1)[1]][j]-values[p+'__shell'][j] for j in range(3)] for p in group['parents']}
        chi={state:cells[v['fixture']]['chi_a0_cubed'] for state,v in self.m['field_aliases'].items()}
        model={cid:copy.deepcopy(cells[cid]['model']) for cid in self.data['rows']}
        response={'pair_dipoles_Debye':pair,'C_dipoles_Debye':C,'chi_states_a0_cubed':chi}
        score=self.h['objective'].evaluate(self.data,model,response,self.targets,self.h['old_objective'],self.h['body_helper'])
        return {'model':model,'response':response,'score':score,'cells':cells}

    def evaluate(self, values):
        values=self.checked_values(values)
        with self.evaluation_lock:
            token=digest(values)
            with self.lock:
                require(not self.stop_event.is_set(),'Global evaluator stop active')
                if token not in self.ledger['vector_keys']:
                    require(len(self.ledger['vector_keys'])<MAX_VECTORS,'Parameter-vector hard ceiling');self.ledger['vector_keys'].append(token)
                request={'values':values[:],'vector_key':token,'status':'started','cells':{cid:{'status':'not_attempted'} for cid in self.m['rows']}}
                self.ledger['requests'].append(request);save(self.ledger_path,self.ledger)
            cells={}
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
                    futures={pool.submit(self._cell,values,cid):cid for cid in self.m['rows']}
                    first_error=None
                    for future in concurrent.futures.as_completed(futures):
                        cid=futures[future]
                        try:
                            value=future.result();cells[cid]=value
                            state={'status':'complete','cache_key':value['request_provenance']['cache_key'],'inherited_native':value['inherited_native']}
                        except Exception as exc:
                            self._stop(exc);state={'status':'failed_or_stopped','error':repr(exc)}
                            if first_error is None:first_error=exc
                        with self.lock:request['cells'][cid]=state;save(self.ledger_path,self.ledger)
                    if first_error is not None:raise first_error
                result=self._assemble(cells)
                with self.lock:request['status']='complete';request['objective']=result['score']['objective'];save(self.ledger_path,self.ledger)
                return result
            except Exception as exc:
                self._stop(exc)
                with self.lock:request.update(status='failed',error=repr(exc));save(self.ledger_path,self.ledger)
                raise

    def snapshot(self):
        with self.lock:
            return {'implementation':self.config['implementation'],'config_sha256':self.fingerprint,'family':copy.deepcopy(self.family),
                    'native_program_invocations':self.ledger['native_program_invocations'],'distinct_parameter_vectors':len(self.ledger['vector_keys']),
                    'stopped':self.stop_event.is_set(),'requests':copy.deepcopy(self.ledger['requests']),'ledger':pin(self.ledger_path),
                    'max_parameter_vectors':MAX_VECTORS,'max_native_programs':MAX_PROGRAMS,'workers':WORKERS,
                    'baseline_result':self.config['baseline_result'],'baseline_review':self.config['baseline_review'],
                    'candidate_forces_native_analytic_only':True,'reserved_outcomes_opened':False,
                    'native_child_CPU_receipts_may_overlap_between_concurrent_workers':True}


def main():
    tick,cpu=time.monotonic(),time.process_time()
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare',));p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One local numerical thread');a.out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads(a.config.read_text());m=prepare(cfg,a.out)
    runtime=dict(cfg,manifest=pin(a.out/'MANIFEST.json'));save(a.out/'CONFIG_EVALUATOR_TEMPLATE.json',runtime)
    ev=Evaluator(runtime,a.out/'initial_cache');initial=ev.evaluate(ev.initial_vector)
    require(ev.snapshot()['native_program_invocations']==0,'Initial point must reuse saved authorities')
    save(a.out/'INITIAL.json',initial);save(a.out/'RESULT.json',{'status':'prepared_joint_evaluator_zero_calls','manifest':pin(a.out/'MANIFEST.json'),
         'configuration':pin(a.out/'CONFIG_EVALUATOR_TEMPLATE.json'),'initial':pin(a.out/'INITIAL.json'),'snapshot':ev.snapshot(),'new_native_calls':0,'preparation_cost':m['cost'],
         'total_cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}})
    print(json.dumps({'result':pin(a.out/'RESULT.json'),'objective':initial['score']['objective'],'native_calls':0}))


if __name__=='__main__':main()
