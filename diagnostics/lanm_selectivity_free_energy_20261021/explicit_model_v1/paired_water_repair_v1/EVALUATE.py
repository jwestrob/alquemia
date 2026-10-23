#!/usr/bin/env python3
"""One paired-water parameter vector; exact native cache, no objective/search."""
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
TRAINING=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/training_v1/RESULT.json'
TRAINING_SHA='99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb'
CELLS=('water_central_La','water_central_Dy','aquo8_central_La','aquo8_central_Dy')
INITIAL=(3.6985300908646654,.25,3.3221124657663292,.139)
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
        'parameter_names':['La_R','La_damping','Dy_R','Dy_damping'],'initial_vector':list(INITIAL),'fixed_pair_epsilon':EPSILON,
        'maximum_parameter_vectors':241,'maximum_native_program_invocations':1920,'maximum_native_concurrency':8,'root_release':None}

class Evaluator:
    """Thread-safe shared instance; each new cell uses serial analyze/testgrad."""
    def __init__(self,config,cache_dir,root_release=None):
        self.config=copy.deepcopy(json.loads(Path(config).read_text()) if isinstance(config,(str,Path)) else config)
        cfg=self.config;require(cfg['implementation']==pin(__file__),'Evaluator source differs');verify(cfg['specification'])
        require(cfg['training_result']['sha256']==TRAINING_SHA,'Wrong training source');data=read(cfg['training_result'])
        require(data['status']=='training_complete' and data['complete']==28 and data['saved_controls_passed'],'Training unavailable')
        self.data=copy.deepcopy(data);self.data['rows']={gid:self.data['rows'][gid] for gid in CELLS}
        self.resolver=module(cfg['resolver'],'paired_water_resolver');self.native=module(cfg['native_helper'],'paired_water_native')
        require(cfg['resolver']['sha256']=='69e1e4dc46cce86f145596180a4507d1e09b3e6473a83210efd99d19e16abba3','Resolver convention changed')
        for p in cfg['executables'].values():verify(p)
        require(cfg['initial_vector']==list(INITIAL) and cfg['fixed_pair_epsilon']==EPSILON,'Initial pair constants changed')
        require(cfg['parameter_names']==['La_R','La_damping','Dy_R','Dy_damping'],'Parameter order changed')
        require(cfg['maximum_parameter_vectors']==241 and cfg['maximum_native_program_invocations']==1920 and cfg['maximum_native_concurrency']==8,'Finite limits changed')
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
        self.family={'parameters':[],'fixed_pair_epsilon':copy.deepcopy(EPSILON),'initial_vector':list(INITIAL),'specification':cfg['specification']}
        for name,v in zip(cfg['parameter_names'],INITIAL):
            lo,hi=(.05,.50) if name.endswith('damping') else (.85*v,1.15*v)
            self.family['parameters'].append({'name':name,'initial':v,'lower':lo,'upper':hi})
        self.baseline_review={}
        for gid,row in self.data['rows'].items():
            require(row['status']=='complete' and row['input']['id']==gid,'Training cell unavailable')
            for key in ('xyz','key','parameters'):
                require(row['input'][key]==cfg['source_inputs'][gid][key],'Input authority differs');verify(row['input'][key])
            verify(row['baseline_result']);verify(row['reference_result']);verify(row['geometry']['xyz'])
            for key in ('output','engrad','input','panel'):verify(row['reference'][key])
            rr=read(row['reference_result']);require(rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values()),'Saved reference unavailable')
            require(rr['forces_kcal_mol_A']==row['reference']['forces'] and rr['energy_hartree']==row['reference']['energy_hartree'],'Reference training arrays differ')
            require(len(row['reference']['forces'])==row['input']['atoms'] and all(len(x)==3 and all(math.isfinite(v) for v in x) for x in row['reference']['forces']),'Reference force shape/values differ')
            native=row['baseline']['native'];replay=self.native.parse(verify(native['analyze_output']),verify(native['gradient_output']),row['input'])
            require(replay==native and row['baseline']['energy']==native['total_energy_kcal_mol'] and row['baseline']['forces']==native['force_kcal_mol_A'],'Saved baseline replay differs')
            prm=verify(row['input']['parameters']).read_text();key=verify(row['input']['key']).read_text();p=self.resolver.parse_parameters(prm,key)
            metal=row['input']['metal'];mt=500 if metal=='La' else 501;offset=0 if metal=='La' else 2
            require(row['input']['types'][0]==mt and set(row['input']['types'][1:])=={36,37},'Only one ion and water admitted')
            require([t for t,a in p['atom'].items() if a['class']==34]==[36],'Water oxygen class has another type alias')
            pair=self.resolver.resolve_pair(p,mt,36)
            require(pair['Rmin_pair_A']==INITIAL[offset] and pair['thole']==INITIAL[offset+1] and pair['epsilon_kcal_mol']==EPSILON[metal],'Initial resolved pair differs')
            self.baseline_review[gid]={'resolved_pair':pair,'native_authority':row['baseline_result'],'native_outputs':{k:native[k] for k in ('analyze_output','gradient_output')},'source_input':cfg['source_inputs'][gid]}

    def checked_values(self,values):
        values=list(values);require(len(values)==4,'Four parameters required')
        for value,p in zip(values,self.family['parameters']):require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and p['lower']<=value<=p['upper'],'Parameter outside frozen bounds')
        return [float(v) for v in values]

    def render(self,values,gid):
        values=self.checked_values(values);require(gid in CELLS,'Cell outside four-water panel');row=self.data['rows'][gid];metal=row['input']['metal'];mt=500 if metal=='La' else 501;i=0 if metal=='La' else 2
        base=verify(row['input']['key']).read_text();prm=verify(row['input']['parameters']).read_text();before=self.resolver.parse_parameters(prm,base)
        key=base if values[i:i+2]==list(INITIAL[i:i+2]) else base.rstrip()+f'\n\n# Paired water development override: {metal} only\nvdwpr 34 {mt} {values[i]:.17g} {EPSILON[metal]:.17g}\npolpair 36 {mt} {values[i+1]:.17g}\n'
        after=self.resolver.parse_parameters(prm,key);pair=self.resolver.resolve_pair(after,mt,36)
        require(pair['Rmin_pair_A']==values[i] and pair['thole']==values[i+1] and pair['epsilon_kcal_mol']==EPSILON[metal],'Rendered pair differs')
        for category in ('settings','atom','vdw','polarize','metal_multipoles'):require(before[category]==after[category],'Non-pair parameter changed')
        for category,allowed in (('vdwpr',(34,mt)),('polpair',(36,mt))):
            require({k:v for k,v in before[category].items() if k!=allowed}=={k:v for k,v in after[category].items() if k!=allowed},'Non-water pair changed')
        require(self.resolver.resolve_pair(before,mt,37)==self.resolver.resolve_pair(after,mt,37),'Water-H pair changed')
        return key,{'metal':metal,'effective_values':values[i:i+2],'resolved_pair':pair,'original_key_reused':key==base,'nonwater_parameters_unchanged':True}

    def cache_identity(self,values,gid):
        key,render=self.render(values,gid);row=self.data['rows'][gid]
        authority={'cell_id':gid,'effective_values':render['effective_values'],'xyz':row['input']['xyz'],'parameters':row['input']['parameters'],
            'rendered_key_sha256':hashlib.sha256(key.encode()).hexdigest(),'implementation':self.config['implementation'],'resolver':self.config['resolver'],
            'native_helper':self.config['native_helper'],'executables':self.config['executables'],'training_result':self.config['training_result']}
        return digest(authority),authority,key,render

    def _gate(self):
        require(not self.stop_event.is_set(),'Global native stop is active')
        require(os.environ.get('SLURM_JOB_ID','').isdigit(),'New native calls require Slurm')
        require(int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=1,'One Slurm task with actual CPU allocation required')
        require(all(os.environ.get(key)=='1' for key in THREADS),'Every native numerical thread control must equal one')
        require(self.root_release is not None,'Explicit root release required');release=read(self.root_release)
        require(release['status']=='approved_paired_water_representability','Wrong root release')
        require(release['evaluator']==self.config['implementation'] and release['evaluator_config_sha256']==self.fingerprint,'Root release does not bind evaluator/config')
        require(release['max_parameter_vectors']==241 and release['max_native_program_invocations']==1920,'Root limits differ')
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
            elif render['original_key_reused']:
                record={'status':'complete','cache_identity':authority,'cache_key':token,'cell_id':gid,'first_requested_vector':values[:],
                    'energy':row['baseline']['energy'],'forces':row['baseline']['forces'],'native':row['baseline']['native'],
                    'native_authority':row['baseline_result'],'source_input':self.config['source_inputs'][gid],'reused_baseline':True,'receipts':[],'render':render}
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
                                self._gate();require(self._ledger['native_program_invocations']<1920,'Native call ceiling reached');self._ledger['native_program_invocations']+=1;save(self._ledger_path,self._ledger)
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
                'native_authority':record.get('native_authority') if record['reused_baseline'] else pin(result_path),'reused_baseline':record['reused_baseline'],
                'source_input':record['source_input'],'render':render,'first_requested_vector':record['first_requested_vector']}
            return value

    def evaluate(self,values):
        values=self.checked_values(values);token=digest(values)
        with self.lock:
            if token not in self._ledger['vector_keys']:
                require(len(self._ledger['vector_keys'])<241,'Parameter-vector ceiling reached');self._ledger['vector_keys'].append(token)
            request={'vector':values[:],'vector_key':token,'status':'started','cells':{gid:{'status':'not_attempted'} for gid in CELLS},'started_monotonic':time.monotonic()};self._ledger['requests'].append(request);save(self._ledger_path,self._ledger)
        model={}
        try:
            for gid in CELLS:
                model[gid]=self._cell(values,gid)
                with self.lock:request['cells'][gid]=model[gid]['provenance'];save(self._ledger_path,self._ledger)
            request['status']='complete';return model
        except Exception as exc:
            request.update(status='failed',error=repr(exc));self._stopped(exc);raise
        finally:
            with self.lock:request['elapsed_seconds']=time.monotonic()-request['started_monotonic'];save(self._ledger_path,self._ledger)

    def snapshot(self):
        with self.lock:
            return {'implementation':self.config['implementation'],'config_sha256':self.fingerprint,'family':copy.deepcopy(self.family),
                'native_program_invocations':self._ledger['native_program_invocations'],'distinct_parameter_vectors':len(self._ledger['vector_keys']),
                'stopped':self.stop_event.is_set(),'requests':copy.deepcopy(self._ledger['requests']),'baseline_review':copy.deepcopy(self.baseline_review),
                'ledger':pin(self._ledger_path),'native_child_CPU_receipts_may_overlap_between_concurrent_workers':True,
                'validation_outcomes_opened':False,'physical_accuracy_qualified':False}
