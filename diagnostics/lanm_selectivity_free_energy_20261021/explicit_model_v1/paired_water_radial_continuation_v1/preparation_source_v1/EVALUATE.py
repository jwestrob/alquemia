#!/usr/bin/env python3
"""Eight consumed water cells; exact native cache, no objective/search."""
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

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];PROJECT=HERE.parents[1]
TRAINING=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_radial_repair_v1/training_v1/RESULT.json'
TRAINING_SHA='3aefeb74f40aaf51e327fdc5c0ffd85f7f00e6e769b07c080f3ece16c345192d'
BASELINE_SHA='702d97ababc232336e3b5405a1326cec4cafeccc9570ce9f514076cf7908e986'
ORIGINAL_SHA='72af9700048cc48e27ee1ad06bdd98cf9cefa974be083a286259e678d47d4003'
CELLS=('water_minus_La','water_minus_Dy','water_central_La','water_central_Dy','water_plus_La','water_plus_Dy','aquo8_central_La','aquo8_central_Dy')
ORIGINAL_INITIAL=(3.6985300908646654,.25,3.3221124657663292,.139)
INITIAL=(3.710148544089586,.2001607792689767,3.389554168967775,.16712500000000005)
EPSILON={'La':.24429636467571708,'Dy':.37332330078584813}
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(path,obj):
    p=Path(path);temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');os.replace(temp,p)
def config_fingerprint(cfg):return digest({k:v for k,v in cfg.items() if k!='root_release'})

def default_config():
    require(pin(TRAINING)['sha256']==TRAINING_SHA,'Training authority changed');data=json.loads(TRAINING.read_text())
    old=json.loads((HERE.parent/'matched_dry_v1/CONFIG_PILOT.template.json').read_text())
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'training_result':pin(TRAINING),
        'resolver':pin(PROJECT/'alchemical_path_v1/PATH.py'),'native_helper':old['native_helper'],'executables':old['executables'],
        'source_inputs':{gid:{k:data['rows'][gid]['input'][k] for k in ('xyz','key','parameters')} for gid in CELLS},
        'baseline_inputs':{gid:data['rows'][gid]['baseline_input'] for gid in CELLS},
        'parent_evaluator':pin(HERE.parent/'paired_water_repair_v1/EVALUATE.py'),
        'parameter_names':['La_R','La_damping','Dy_R','Dy_damping'],'initial_vector':list(INITIAL),'fixed_pair_epsilon':EPSILON,
        'prior_result':pin(ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_radial_repair_v1/run_v1/execution/RESULT.json'),
        'prior_evaluator':pin(HERE.parent/'paired_water_radial_repair_v1/EVALUATE.py'),'maximum_parameter_vectors':385,'maximum_native_program_invocations':4096,'maximum_native_concurrency':8,'root_release':None}

class Evaluator:
    """Thread-safe shared instance; each new cell uses serial analyze/testgrad."""
    def __init__(self,config,cache_dir,root_release=None):
        self.config=copy.deepcopy(json.loads(Path(config).read_text()) if isinstance(config,(str,Path)) else config)
        cfg=self.config;require(cfg['implementation']==pin(__file__),'Evaluator source differs');verify(cfg['specification'])
        require(cfg['training_result']['sha256']==TRAINING_SHA,'Wrong training source');data=read(cfg['training_result'])
        require(data['status']=='radial_training_complete' and data['complete']==8 and data['saved_controls_passed'] and set(data['rows'])==set(CELLS),'Training unavailable')
        require(cfg['parent_evaluator']['sha256']==ORIGINAL_SHA,'Original evaluator identity changed');verify(cfg['parent_evaluator'])
        self.data=copy.deepcopy(data);self.data['rows']={gid:self.data['rows'][gid] for gid in CELLS}
        self.resolver=module(cfg['resolver'],'paired_water_resolver');self.native=module(cfg['native_helper'],'paired_water_native')
        require(cfg['resolver']['sha256']=='69e1e4dc46cce86f145596180a4507d1e09b3e6473a83210efd99d19e16abba3','Resolver convention changed')
        for p in cfg['executables'].values():verify(p)
        require(cfg['initial_vector']==list(INITIAL) and cfg['fixed_pair_epsilon']==EPSILON,'Initial pair constants changed')
        require(cfg['parameter_names']==['La_R','La_damping','Dy_R','Dy_damping'],'Parameter order changed')
        require(cfg['maximum_parameter_vectors']==385 and cfg['maximum_native_program_invocations']==4096 and cfg['maximum_native_concurrency']==8,'Finite limits changed')
        self.root_release=root_release if root_release is not None else cfg.get('root_release')
        self.lock=threading.RLock();self.stop_event=threading.Event();self._cell_locks={};self._cache={};self._requests=[]
        self.cache_dir=Path(cache_dir).resolve();self.cache_dir.mkdir(parents=True,exist_ok=True)
        self.fingerprint=config_fingerprint(cfg);self._ledger_path=self.cache_dir/'LEDGER.json'
        if self._ledger_path.exists():
            self._ledger=json.loads(self._ledger_path.read_text());require(self._ledger['config_sha256']==self.fingerprint,'Cache context differs')
            if self._ledger['stopped']:self.stop_event.set()
        else:self._ledger={'config_sha256':self.fingerprint,'vector_keys':[],'native_program_invocations':0,'stopped':False,'requests':[]};save(self._ledger_path,self._ledger)
        allocated=int(os.environ.get('SLURM_CPUS_PER_TASK','1'))*int(os.environ.get('SLURM_NTASKS','1'))
        self._native_slots=threading.BoundedSemaphore(max(1,min(8,allocated)))
        self.family={'parameters':[],'fixed_pair_epsilon':copy.deepcopy(EPSILON),'initial_vector':list(INITIAL),'original_vector_for_bounds':list(ORIGINAL_INITIAL),'specification':cfg['specification']}
        for name,v in zip(cfg['parameter_names'],INITIAL):
            origin=ORIGINAL_INITIAL[len(self.family['parameters'])]
            lo,hi=(.05,.50) if name.endswith('damping') else (.85*origin,1.15*origin)
            self.family['parameters'].append({'name':name,'initial':v,'lower':lo,'upper':hi})
        self.baseline_review={}
        for gid,row in self.data['rows'].items():
            require(row['status']=='complete' and row['input']['id']==gid,'Training cell unavailable')
            for key in ('xyz','key','parameters'):
                require(row['input'][key]==cfg['source_inputs'][gid][key],'Input authority differs');verify(row['input'][key])
            require(row['baseline_result']['sha256']==BASELINE_SHA and row['baseline_parameters']==list(INITIAL) and row['baseline_parse_mode']=='analytic_only','Starting authority/vector/parse mode differs')
            baseline=read(row['baseline_result']);require(baseline['status']=='consumed_transfer_complete' and baseline['stage_A_passed'] and not baseline['global_stop'],'Starting transfer unavailable')
            actual=baseline['rows'][gid];require(actual['status']=='complete','Starting cell unavailable')
            require(all(actual['model'][k]==row['baseline'][k] for k in ('energy','forces','native')),'Starting model differs from exact transfer row')
            require(row['baseline_input']==cfg['baseline_inputs'][gid]==actual['input']['files'],'Baseline input authority differs')
            for k in ('xyz','key','parameters'):verify(row['baseline_input'][k])
            require(all(row['baseline_input'][k]['sha256']==row['input'][k]['sha256'] for k in ('xyz','parameters')),'Starting geometry/parameters differ')
            verify(row['original_baseline_result']);verify(row['reference_result']);verify(row['geometry']['xyz'])
            for key in ('output','engrad','input','panel'):verify(row['reference'][key])
            rr=read(row['reference_result']);require(rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values()),'Saved reference unavailable')
            require(rr['forces_kcal_mol_A']==row['reference']['forces'] and rr['energy_hartree']==row['reference']['energy_hartree'],'Reference training arrays differ')
            require(len(row['reference']['forces'])==row['input']['atoms'] and all(len(x)==3 and all(math.isfinite(v) for v in x) for x in row['reference']['forces']),'Reference force shape/values differ')
            native=row['baseline']['native'];replay=self.native.parse(verify(native['analyze_output']),verify(native['gradient_output']),dict(row['input'],derivative_test=False))
            require(replay==native and row['baseline']['energy']==native['total_energy_kcal_mol'] and row['baseline']['forces']==native['force_kcal_mol_A'],'Saved baseline replay differs')
            prm=verify(row['input']['parameters']).read_text();key=verify(row['input']['key']).read_text();p=self.resolver.parse_parameters(prm,key)
            metal=row['input']['metal'];mt=500 if metal=='La' else 501;offset=0 if metal=='La' else 2
            require(row['input']['types'][0]==mt and set(row['input']['types'][1:])=={36,37},'Only one ion and water admitted')
            require([t for t,a in p['atom'].items() if a['class']==34]==[36],'Water oxygen class has another type alias')
            pair=self.resolver.resolve_pair(p,mt,36)
            require(pair['Rmin_pair_A']==ORIGINAL_INITIAL[offset] and pair['thole']==ORIGINAL_INITIAL[offset+1] and pair['epsilon_kcal_mol']==EPSILON[metal],'Initial resolved pair differs')
            rendered,render=self.render(INITIAL,gid)
            require(rendered==verify(row['baseline_input']['key']).read_text() and render['baseline_key_reused'],'Starting rendered input differs')
            self.baseline_review[gid]={'original_resolved_pair':pair,'starting_resolved_pair':render['resolved_pair'],'baseline_input':row['baseline_input'],'native_authority':row['baseline_result'],'native_outputs':{k:native[k] for k in ('analyze_output','gradient_output')},'source_input':cfg['source_inputs'][gid]}

        self._inherited={};self._inherited_admitted=set()
        require(cfg['prior_result']['sha256']=='31cb5caa7f96a3f312cc3367964149923665faa16b41d8993e664ad22eeb259c','Wrong immutable continuation parent')
        require(cfg['prior_evaluator']['sha256']=='1cf6025e27291523ac7c7242a07245da22c8927d67f4323c87392974133eda56','Wrong prior evaluator');verify(cfg['prior_evaluator'])
        prior=read(cfg['prior_result']);pcfg=read(prior['config']);pecfg=read(pcfg['evaluator_config'])
        require(prior['status']=='bounded_radial_refinement_complete' and prior['error'] is None and prior['candidate_frozen'] and prior['original_results_unchanged'],'Prior result unavailable')
        require(pcfg['evaluator']==cfg['prior_evaluator'] and pcfg['training_result']==cfg['training_result'],'Prior evaluator/training differs')
        for k in ('resolver','native_helper','executables','source_inputs','baseline_inputs','fixed_pair_epsilon'):require(pecfg[k]==cfg[k],'Prior physical context differs: '+k)
        self.prior_ledger=prior['evaluator']['ledger'];verify(self.prior_ledger)
        self.continuation_start=prior['selected']['parameters']
        seen=set()
        for vector in prior['vectors'].values():
            require(vector['status']=='complete' and set(vector['model'])==set(CELLS),'Prior vector unavailable')
            for gid,value in vector['model'].items():
                pv=value['provenance']
                if pv['reused_baseline']:continue
                ap=pv['native_authority']
                if ap['path'] in seen:continue
                seen.add(ap['path']);record=read(ap);identity=record['cache_identity']
                require(record['status']=='complete' and record['cell_id']==gid and record['cache_key']==digest(identity)==pv['cache_key'],'Inherited native authority differs')
                require(Path(ap['path']).is_relative_to(Path(cfg['prior_result']['path']).parent/'cache/cells'),'Inherited native outside prior cache')
                require(identity['implementation']==cfg['prior_evaluator'],'Do not relabel prior native source')
                for k in ('resolver','native_helper','executables','training_result'):require(identity[k]==cfg[k],'Inherited native context differs: '+k)
                require(all(record[k]==value[k] for k in ('energy','forces','native')) and record['input']==pv['effective_input'],'Inherited saved values/input differ')
                require(len(record['receipts'])==2 and all(r['returncode']==0 for r in record['receipts']),'Inherited native receipts incomplete')
                require(record['source_input']==cfg['source_inputs'][gid] and identity['xyz']==cfg['source_inputs'][gid]['xyz'] and identity['parameters']==cfg['source_inputs'][gid]['parameters'],'Inherited source geometry/parameters differ')
                effective=dict(identity,implementation=cfg['implementation']);newkey=digest(effective)
                self._inherited[newkey]={'authority':ap,'record':record,'effective_identity':effective}
        require(len(self._inherited)*2==prior['new_native_program_invocations'],'Inherited complete native-cell count differs')

    def checked_values(self,values):
        values=list(values);require(len(values)==4,'Four parameters required')
        for value,p in zip(values,self.family['parameters']):require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and p['lower']<=value<=p['upper'],'Parameter outside frozen bounds')
        return [float(v) for v in values]

    def render(self,values,gid):
        values=self.checked_values(values);require(gid in CELLS,'Cell outside eight-water panel');row=self.data['rows'][gid];metal=row['input']['metal'];mt=500 if metal=='La' else 501;i=0 if metal=='La' else 2
        base=verify(row['input']['key']).read_text();prm=verify(row['input']['parameters']).read_text();before=self.resolver.parse_parameters(prm,base)
        key=base if values[i:i+2]==list(ORIGINAL_INITIAL[i:i+2]) else base.rstrip()+f'\n\n# Paired water development override: {metal} only\nvdwpr 34 {mt} {values[i]:.17g} {EPSILON[metal]:.17g}\npolpair 36 {mt} {values[i+1]:.17g}\n'
        after=self.resolver.parse_parameters(prm,key);pair=self.resolver.resolve_pair(after,mt,36)
        require(pair['Rmin_pair_A']==values[i] and pair['thole']==values[i+1] and pair['epsilon_kcal_mol']==EPSILON[metal],'Rendered pair differs')
        for category in ('settings','atom','vdw','polarize','metal_multipoles'):require(before[category]==after[category],'Non-pair parameter changed')
        for category,allowed in (('vdwpr',(34,mt)),('polpair',(36,mt))):
            require({k:v for k,v in before[category].items() if k!=allowed}=={k:v for k,v in after[category].items() if k!=allowed},'Non-water pair changed')
        require(self.resolver.resolve_pair(before,mt,37)==self.resolver.resolve_pair(after,mt,37),'Water-H pair changed')
        return key,{'metal':metal,'effective_values':values[i:i+2],'resolved_pair':pair,'original_key_reused':key==base,'baseline_key_reused':key==verify(self.config['baseline_inputs'][gid]['key']).read_text(),'nonwater_parameters_unchanged':True}

    def cache_identity(self,values,gid):
        key,render=self.render(values,gid);row=self.data['rows'][gid]
        authority={'cell_id':gid,'effective_values':render['effective_values'],'xyz':row['input']['xyz'],'parameters':row['input']['parameters'],
            'rendered_key_sha256':hashlib.sha256(key.encode()).hexdigest(),'implementation':self.config['implementation'],'resolver':self.config['resolver'],
            'native_helper':self.config['native_helper'],'executables':self.config['executables'],'training_result':self.config['training_result']}
        return digest(authority),authority,key,render

    def _gate(self):
        require(not self.stop_event.is_set(),'Global native stop is active')
        require(os.environ.get('SLURM_JOB_ID','').isdigit(),'New native calls require Slurm')
        require(int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'One Slurm task with exactly eight CPUs required')
        require(all(os.environ.get(key)=='1' for key in THREADS),'Every native numerical thread control must equal one')
        require(self.root_release is not None,'Explicit root release required');release=read(self.root_release)
        require(release['status']=='approved_paired_water_radial_continuation','Wrong root release')
        require(release['evaluator']==self.config['implementation'] and release['evaluator_config_sha256']==self.fingerprint,'Root release does not bind evaluator/config')
        require(release['max_parameter_vectors']==385 and release['max_native_program_invocations']==4096,'Root limits differ')
        require(release['validation_scoring_authorized'] is False and release['physical_accuracy_claimed'] is False,'Release scope differs')

    def _stopped(self,exc):
        self.stop_event.set()
        with self.lock:self._ledger['stopped']=True;self._ledger['failure']=repr(exc);save(self._ledger_path,self._ledger)

    def _program_admission(self,path,program):
        text=Path(path).read_text()
        require('Tinker' in text and 'Version 26.2' in text,'Native identity absent')
        require(not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b',text),'Native failure/nonfinite marker')
        energies=re.findall(r'Total Potential Energy\s*:\s*('+self.native.NUMBER+')',text)
        require(len(energies)==1 and math.isfinite(float(energies[0])),'Native finite energy unavailable')
        if program=='analyze':require('Energy Component Breakdown :' in text,'Native energy components unavailable')

    def _cell(self,values,gid):
        token,authority,key,render=self.cache_identity(values,gid)
        with self.lock:lock=self._cell_locks.setdefault(token,threading.Lock())
        with lock:
            row=self.data['rows'][gid];directory=self.cache_dir/'cells'/token;result_path=directory/'RESULT.json'
            if token in self._cache:record=self._cache[token]
            elif token in self._inherited:
                inherited=self._inherited[token];record=inherited['record'];verify(inherited['authority'])
                require(inherited['effective_identity']==authority and record['render']==render,'Inherited effective input differs')
                for p in record['input'].values():verify(p)
                require(record['input']['xyz']['sha256']==authority['xyz']['sha256'] and record['input']['parameters']['sha256']==authority['parameters']['sha256'] and record['input']['key']['sha256']==authority['rendered_key_sha256'],'Inherited physical input differs')
                for receipt in record['receipts']:verify(receipt['output'])
                native=record['native'];again=self.native.parse(verify(native['analyze_output']),verify(native['gradient_output']),dict(row['input'],derivative_test=False))
                require(again==native and record['energy']==native['total_energy_kcal_mol'] and record['forces']==native['force_kcal_mol_A'],'Inherited raw replay differs')
                record=copy.deepcopy(record);record['inherited_native_authority']=inherited['authority'];self._cache[token]=record;self._inherited_admitted.add(token)
            elif render['baseline_key_reused']:
                record={'status':'complete','cache_identity':authority,'cache_key':token,'cell_id':gid,'first_requested_vector':values[:],
                    'energy':row['baseline']['energy'],'forces':row['baseline']['forces'],'native':row['baseline']['native'],
                    'native_authority':row['baseline_result'],'baseline_input':row['baseline_input'],'source_input':self.config['source_inputs'][gid],'reused_baseline':True,'receipts':[],'render':render}
                self._cache[token]=record
            elif result_path.exists():
                record=json.loads(result_path.read_text());require(record['cache_identity']==authority and record['status']=='complete','Cached cell failed or differs; no automatic retry')
                for value in record['input'].values():verify(value)
                require(record['input']['xyz']['sha256']==authority['xyz']['sha256'] and record['input']['parameters']['sha256']==authority['parameters']['sha256'] and record['input']['key']['sha256']==authority['rendered_key_sha256'],'Cached effective input differs')
                require(record['render']==render and len(record['receipts'])==2 and all(x['returncode']==0 for x in record['receipts']),'Cached execution evidence differs')
                native=record['native'];again=self.native.parse(verify(native['analyze_output']),verify(native['gradient_output']),dict(row['input'],derivative_test=False))
                require(again==native and record['energy']==native['total_energy_kcal_mol'] and record['forces']==native['force_kcal_mol_A'],'Cached native replay differs');self._cache[token]=record
            else:
                self._gate();directory.mkdir(parents=True,exist_ok=False)
                record={'status':'started','cache_identity':authority,'cache_key':token,'cell_id':gid,'first_requested_vector':values[:],
                    'source_input':self.config['source_inputs'][gid],'render':render,'reused_baseline':False,'receipts':[],'job_id':os.environ['SLURM_JOB_ID']}
                try:
                    shutil.copyfile(verify(row['input']['xyz']),directory/'system.xyz');shutil.copyfile(verify(row['input']['parameters']),directory/'amoeba09.prm');(directory/'system.key').write_text(key)
                    record['input']={k:pin(directory/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
                    require(record['input']['xyz']['sha256']==row['input']['xyz']['sha256'] and record['input']['parameters']['sha256']==row['input']['parameters']['sha256'],'Copied geometry/parameter changed')
                    require(record['input']['key']['sha256']==authority['rendered_key_sha256'],'Rendered key changed');save(result_path,record)
                    for program,args in (('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])):
                        with self._native_slots:
                            with self.lock:
                                self._gate();require(self._ledger['native_program_invocations']<4096,'Native call ceiling reached');self._ledger['native_program_invocations']+=1;save(self._ledger_path,self._ledger)
                            try:
                                self.native.invoke([self.config['executables'][program]['path']]+args,directory,directory/(program+'.out'),record['receipts'])
                                self._program_admission(directory/(program+'.out'),program)
                                if program=='testgrad':raw=self.native.parse(directory/'analyze.out',directory/'testgrad.out',dict(row['input'],derivative_test=False))
                            except Exception as exc:self._stopped(exc);raise
                            save(result_path,record)
                    record.update(status='complete',energy=raw['total_energy_kcal_mol'],forces=raw['force_kcal_mol_A'],native=raw)
                    self._cache[token]=record
                except Exception as exc:
                    record.update(status='failed',error=repr(exc));self._stopped(exc);raise
                finally:save(result_path,record)
            value={key:copy.deepcopy(record[key]) for key in ('energy','forces','native')}
            value['provenance']={'requested_vector':values[:],'cache_key':token,'effective_parameters':render['effective_values'],
                'native_authority':record.get('native_authority') if record['reused_baseline'] else record.get('inherited_native_authority') or pin(result_path),'reused_baseline':record['reused_baseline'],
                'source_input':record['source_input'],'effective_input':record.get('baseline_input',record.get('input')),
                'render':render,'first_requested_vector':record['first_requested_vector'],
                'inherited_native':bool(record.get('inherited_native_authority')),'native_cache_key':record['cache_key'],
                'native_implementation':record['cache_identity']['implementation']}
            return value

    def evaluate(self,values):
        values=self.checked_values(values);token=digest(values)
        with self.lock:
            if token not in self._ledger['vector_keys']:
                require(len(self._ledger['vector_keys'])<385,'Parameter-vector ceiling reached');self._ledger['vector_keys'].append(token)
            request={'vector':values[:],'vector_key':token,'status':'started','cells':{gid:{'status':'not_attempted'} for gid in CELLS},'started_monotonic':time.monotonic()};self._ledger['requests'].append(request);save(self._ledger_path,self._ledger)
        model={}
        try:
            for gid in CELLS:
                model[gid]=self._cell(values,gid)
                with self.lock:request['cells'][gid]=model[gid]['provenance'];save(self._ledger_path,self._ledger)
            with self.lock:request['status']='complete'
            return model
        except Exception as exc:
            with self.lock:request.update(status='failed',error=repr(exc))
            self._stopped(exc);raise
        finally:
            with self.lock:request['elapsed_seconds']=time.monotonic()-request['started_monotonic'];save(self._ledger_path,self._ledger)

    def snapshot(self):
        with self.lock:
            return {'implementation':self.config['implementation'],'config_sha256':self.fingerprint,'family':copy.deepcopy(self.family),
                'native_program_invocations':self._ledger['native_program_invocations'],'distinct_parameter_vectors':len(self._ledger['vector_keys']),
                'stopped':self.stop_event.is_set(),'requests':copy.deepcopy(self._ledger['requests']),'baseline_review':copy.deepcopy(self.baseline_review),
                'inherited_native_cells_indexed':len(self._inherited),'inherited_native_cells_admitted':len(self._inherited_admitted),
                'prior_result':self.config['prior_result'],'prior_ledger':self.prior_ledger,
                'ledger':pin(self._ledger_path),'native_child_CPU_receipts_may_overlap_between_concurrent_workers':True,
                'validation_outcomes_opened':False,'physical_accuracy_qualified':False}
