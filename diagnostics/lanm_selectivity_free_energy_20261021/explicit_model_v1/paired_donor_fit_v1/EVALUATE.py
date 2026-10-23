#!/usr/bin/env python3
"""Eight donor-pair parameters; frozen qualified water, saved28 baseline, no search."""
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
WORK=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
TRAINING=WORK/'matched_dry_v1/training_v1/RESULT.json'
TRAINING_SHA='99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb'
TRANSFER_SHA='e3eeb870aac22f408c8b38a5f050404cfee5a8f4c0468633c07122fbc5e07254'
WATER_SHA='6027e319a261eb85171fb06e30aa3e519411f6371a7905eeca0e8a79812c48fc'
FAMILIES=('acetate_mono','acetate_bi','acetamide')
GEOMETRIES=('water_central','aquo8_central')+tuple('hydrated_'+f+'_'+q for f in FAMILIES for q in ('minus','central','plus'))+tuple('matched_dry_'+f+'_central' for f in FAMILIES)
CELLS=tuple(q+'_'+m for q in GEOMETRIES for m in ('La','Dy'))
WATER_CENTERS=tuple(q+'_'+m for q in ('water_central','aquo8_central') for m in ('La','Dy'))
WATER_CELLS=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central') for m in ('La','Dy'))
INITIAL=(3.6318,.1518,3.3,.2987,3.347237658604727,.139,3.2657502966294456,.139)
WATER=(3.906136716492464,.28782327521600637,3.5865909437359376,.2837564269122712,.20765190997435956,.3173248056679709)
EPSILON={'La_carboxylate':.228,'La_amide':.3019,'Dy_carboxylate':.3795688850865599,'Dy_amide':.3795688850865599}
NAMES=[m+'_'+f+'_'+v for m in ('La','Dy') for f in ('carboxylate','amide') for v in ('R','damping')]
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
    transfer_pin=pin(WORK/'paired_water_six_transfer_v1/run_v1/execution/RESULT.json')
    require(transfer_pin['sha256']==TRANSFER_SHA,'Transfer authority changed');transfer=read(transfer_pin);bound=read(transfer['bound'])
    old=json.loads((HERE.parent/'matched_dry_v1/CONFIG_PILOT.template.json').read_text())
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'training_result':pin(TRAINING),
        'water_result':pin(WORK/'paired_water_six_parameter_v1/run_v1/execution/RESULT.json'),'transfer_result':transfer_pin,
        'resolver':pin(PROJECT/'alchemical_path_v1/PATH.py'),'native_helper':old['native_helper'],'executables':old['executables'],
        'source_inputs':{gid:{k:data['rows'][gid]['input'][k] for k in ('xyz','key','parameters')} for gid in CELLS},
        'baseline_inputs':{gid:bound['inputs'][gid]['files'] for gid in CELLS},
        'parameter_names':NAMES,'initial_vector':list(INITIAL),'fixed_pair_epsilon':EPSILON,'frozen_water_parameters':list(WATER),
        'maximum_parameter_vectors':481,'maximum_native_program_invocations':9216,'maximum_native_concurrency':8,'root_release':None}

class Evaluator:
    """Shared thread-safe native cache; one serial analyze/testgrad pair per new cell."""
    def __init__(self,config,cache_dir,root_release=None):
        self.config=copy.deepcopy(json.loads(Path(config).read_text()) if isinstance(config,(str,Path)) else config)
        cfg=self.config;require(cfg['implementation']==pin(__file__),'Evaluator source differs');verify(cfg['specification'])
        require(cfg['training_result']['sha256']==TRAINING_SHA,'Wrong training source');data=read(cfg['training_result'])
        require(data['status']=='training_complete' and data['complete']==28 and data['saved_controls_passed'] and set(data['rows'])==set(CELLS),'Training unavailable')
        self.data=copy.deepcopy(data)  # Original consumed baseline/reference metadata stays unchanged.
        self.resolver=module(cfg['resolver'],'paired_donor_resolver');self.native=module(cfg['native_helper'],'paired_donor_native')
        require(cfg['resolver']['sha256']=='69e1e4dc46cce86f145596180a4507d1e09b3e6473a83210efd99d19e16abba3','Resolver convention changed')
        require(cfg['native_helper']['sha256']=='4189602145853bc39d8dfaf30915b2a6b1d4632984e0c50fc0bf321644fa306c','Native helper changed')
        for p in cfg['executables'].values():verify(p)
        require(cfg['initial_vector']==list(INITIAL) and cfg['fixed_pair_epsilon']==EPSILON and cfg['frozen_water_parameters']==list(WATER),'Initial/fixed parameters changed')
        require(cfg['parameter_names']==NAMES,'Parameter order changed')
        require(cfg['maximum_parameter_vectors']==481 and cfg['maximum_native_program_invocations']==9216 and cfg['maximum_native_concurrency']==8,'Finite limits changed')
        self.initial_vector=list(INITIAL);self.continuation_start=self.initial_vector[:]
        self.family={'parameters':[],'fixed_pair_epsilon':copy.deepcopy(EPSILON),'initial_vector':list(INITIAL),
            'frozen_water_parameters':list(WATER),'specification':cfg['specification']}
        for name,v in zip(NAMES,INITIAL):
            lo,hi=(.05,.50) if name.endswith('damping') else (.85*v,1.15*v)
            self.family['parameters'].append({'name':name,'initial':v,'lower':lo,'upper':hi})
        self.root_release=root_release if root_release is not None else cfg.get('root_release')
        self.lock=threading.RLock();self.stop_event=threading.Event();self._cell_locks={};self._cache={}
        self.cache_dir=Path(cache_dir).resolve();self.cache_dir.mkdir(parents=True,exist_ok=True)
        self.fingerprint=config_fingerprint(cfg);self._ledger_path=self.cache_dir/'LEDGER.json'
        if self._ledger_path.exists():
            self._ledger=json.loads(self._ledger_path.read_text());require(self._ledger['config_sha256']==self.fingerprint,'Cache context differs')
            if self._ledger['stopped']:self.stop_event.set()
        else:self._ledger={'config_sha256':self.fingerprint,'vector_keys':[],'native_program_invocations':0,'stopped':False,'requests':[]};save(self._ledger_path,self._ledger)
        allocated=int(os.environ.get('SLURM_CPUS_PER_TASK','1'))*int(os.environ.get('SLURM_NTASKS','1'))
        self._native_slots=threading.BoundedSemaphore(max(1,min(8,allocated)))
        require(cfg['transfer_result']['sha256']==TRANSFER_SHA and cfg['water_result']['sha256']==WATER_SHA,'Wrong actual baseline authorities')
        transfer=read(cfg['transfer_result']);water=read(cfg['water_result']);bound=read(transfer['bound'])
        require(transfer['status']=='consumed_six_parameter_transfer_complete' and transfer['error'] is None and transfer['stage_A_passed'] and not transfer['global_stop'] and transfer['job_id']=='17612','Transfer unavailable')
        require(water['status']=='bounded_six_parameter_water_complete' and water['error'] is None and water['development_checks_passed'] and water['candidate_frozen'] and water['job_id']=='17605','Passing water unavailable')
        require(water['selected']['parameters']==list(WATER) and water['selected']['score']['development_checks_passed'] and set(water['selected']['model'])==set(WATER_CELLS),'Water selection differs')
        require(bound['candidate_result']==cfg['water_result'] and bound['parameters']==list(WATER),'Transfer/water join differs')
        all_data=read(bound['data']);self.baselines={};self.baseline_review={};self._base_text={};self._base_parameters={}
        for gid,row in self.data['rows'].items():
            require(row['status']=='complete' and row['input']['id']==gid,'Training cell unavailable')
            require(row==all_data['rows'][gid],'Original consumed row changed')
            for k in ('xyz','key','parameters'):
                require(row['input'][k]==cfg['source_inputs'][gid][k],'Original input differs');verify(row['input'][k])
            actual=transfer['rows'][gid];require(actual['status']=='complete','Transfer cell unavailable')
            inp=cfg['baseline_inputs'][gid]
            require(inp==bound['inputs'][gid]['files'],'Baseline input authority differs')
            if 'input' in actual:require(inp==actual['input']['files'],'Actual native input differs')
            else:
                require(gid.startswith('matched_dry_') and actual['reused_from']==bound['reuse'][gid]['authority'],'Unchanged dry authority differs');verify(actual['reused_from'])
                require(all(actual['model'][k]==bound['reuse'][gid]['value'][k] for k in ('energy','forces','native')),'Dry reused values differ')
            for k in ('xyz','key','parameters'):verify(inp[k])
            require(all(inp[k]['sha256']==row['input'][k]['sha256'] for k in ('xyz','parameters')),'Baseline geometry/parameter file differs')
            verify(row['baseline_result']);verify(row['reference_result']);verify(row['geometry']['xyz'])
            for k in ('output','engrad','input','panel'):verify(row['reference'][k])
            rr=read(row['reference_result']);require(rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values()),'Saved reference unavailable')
            require(rr['forces_kcal_mol_A']==row['reference']['forces'] and rr['energy_hartree']==row['reference']['energy_hartree'],'Reference arrays differ')
            native=actual['model']['native'];again=self.native.parse(verify(native['analyze_output']),verify(native['gradient_output']),dict(row['input'],derivative_test=native['finite_difference'] is not None))
            require(again==native and actual['model']['energy']==native['total_energy_kcal_mol'] and actual['model']['forces']==native['force_kcal_mol_A'],'Actual baseline raw replay differs')
            self.baselines[gid]=copy.deepcopy(actual['model'])
            prm=verify(inp['parameters']).read_text();key=verify(inp['key']).read_text();before=self.resolver.parse_parameters(prm,key)
            self._base_text[gid]=(prm,key);self._base_parameters[gid]=before
            metal=row['input']['metal'];mt=500 if metal=='La' else 501
            require(row['input']['types'][0]==mt and all(t not in (500,501) for t in row['input']['types'][1:]),'One declared metal first required')
            family,j,cl,typ=self.recipient(gid)
            require((403 in row['input']['types'])==(family=='carboxylate') and (105 in row['input']['types'])==(family=='amide'),'Donor family/type differs')
            pair=None
            if family:
                require(before['atom'][typ]['class']==cl,'Donor type/class differs')
                require(all(before['atom'][t]['class']!=cl or t==typ for t in set(row['input']['types'])),'Recipient class alias in actual fixture')
                pair=self.resolver.resolve_pair(before,mt,typ)
                require([pair['Rmin_pair_A'],pair['thole']]==list(INITIAL[j:j+2]) and pair['epsilon_kcal_mol']==EPSILON[metal+'_'+family],'Original donor pair differs')
            if 36 in row['input']['types']:
                wp=self.resolver.resolve_pair(before,mt,36);wi=0 if metal=='La' else 2
                require([wp['Rmin_pair_A'],wp['thole'],wp['epsilon_kcal_mol']]==[WATER[wi],WATER[wi+1],WATER[4 if metal=='La' else 5]],'Frozen water pair differs')
            if family=='amide' and metal=='La':
                cp=self.resolver.resolve_pair(before,mt,104);require(cp['nonmetal_class']==53 and cp['Rmin_pair_A']==3.9541 and cp['epsilon_kcal_mol']==.8295,'La amide-C override differs')
            self.baseline_review[gid]={'transfer_authority':cfg['transfer_result'],'original_source_input':cfg['source_inputs'][gid],
                'baseline_input':inp,'resolved_donor_pair':pair,
                'reused_transfer_origin':{k:v for k,v in bound['reuse'].get(gid,{}).items() if k!='value'},
                'native_outputs':{k:native[k] for k in ('analyze_output','gradient_output')}}
        self.water_invariants={'passed':True,'water_result':cfg['water_result'],'parameters':list(WATER),
            'original_controls':copy.deepcopy(water['selected']['score']),'records':{},'new_native_calls':0}
        for gid in WATER_CELLS:
            value=water['selected']['model'][gid];pv=value['provenance'];record=read(pv['native_authority']);row=all_data['rows'][gid]
            require(all(record[k]==value[k] for k in ('energy','forces','native')),'Water native authority differs')
            native=value['native'];again=self.native.parse(verify(native['analyze_output']),verify(native['gradient_output']),dict(row['input'],derivative_test=native['finite_difference'] is not None))
            require(again==native and value['energy']==native['total_energy_kcal_mol'] and value['forces']==native['force_kcal_mol_A'],'Water raw replay differs')
            for p in pv['effective_input'].values():verify(p)
            require(all(pv['effective_input'][k]['sha256']==bound['inputs'][gid]['files'][k]['sha256'] for k in ('xyz','key','parameters')),'Water physical input differs')
            require(all(value[k]==transfer['rows'][gid]['model'][k] for k in ('energy','forces','native')),'Transfer changed selected water')
            self.water_invariants['records'][gid]={'native_authority':pv['native_authority'],'native_provenance':pv,'energy':value['energy'],
                'forces':value['forces'],'native':native,'unchanged':True}
        self.water_baseline={gid:copy.deepcopy(self.baselines[gid]) for gid in WATER_CENTERS}

    def checked_values(self,values):
        values=list(values);require(len(values)==8,'Eight parameters required')
        for value,p in zip(values,self.family['parameters']):require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and p['lower']<=value<=p['upper'],'Parameter outside frozen bounds')
        return [float(v) for v in values]

    def recipient(self,gid):
        require(gid in CELLS,'Cell outside consumed28');metal=self.data['rows'][gid]['input']['metal'];offset=0 if metal=='La' else 4
        if 'acetate_' in gid:return 'carboxylate',offset,403,403
        if 'acetamide_' in gid:return 'amide',offset+2,55,105
        return None,None,None,None

    def render(self,values,gid):
        values=self.checked_values(values);family,i,cl,typ=self.recipient(gid);row=self.data['rows'][gid];metal=row['input']['metal'];mt=500 if metal=='La' else 501
        prm,base=self._base_text[gid];before=self._base_parameters[gid];effective=[] if family is None else values[i:i+2]
        key=base
        if family and effective!=list(INITIAL[i:i+2]):
            key=base.rstrip()+f'\n\n# Paired donor development override: {metal} {family} only\nvdwpr {cl} {mt} {values[i]:.17g} {EPSILON[metal+"_"+family]:.17g}\npolpair {typ} {mt} {values[i+1]:.17g}\n'
        after=before if key==base else self.resolver.parse_parameters(prm,key);pair=None
        if family:
            pair=self.resolver.resolve_pair(after,mt,typ)
            require([pair['Rmin_pair_A'],pair['thole']]==effective and pair['epsilon_kcal_mol']==EPSILON[metal+'_'+family],'Rendered donor pair differs')
        for category in ('settings','atom','vdw','polarize','metal_multipoles'):require(before[category]==after[category],'Non-pair parameter changed')
        for category,allowed in (('vdwpr',(cl,mt)),('polpair',(typ,mt))):
            require({k:v for k,v in before[category].items() if k!=allowed}=={k:v for k,v in after[category].items() if k!=allowed},'Unrelated pair changed')
        for t in sorted(set(row['input']['types'][1:])):
            if t!=typ:require(self.resolver.resolve_pair(before,mt,t)==self.resolver.resolve_pair(after,mt,t),'Other atom pair changed')
        return key,{'metal':metal,'donor_family':family,'effective_values':effective,'resolved_pair':pair,
            'baseline_key_reused':key==base,'fixed_water_and_other_parameters_unchanged':True}

    def cache_identity(self,values,gid):
        key,render=self.render(values,gid);row=self.data['rows'][gid]
        authority={'cell_id':gid,'effective_values':render['effective_values'],'donor_family':render['donor_family'],
            'xyz':row['input']['xyz'],'parameters':row['input']['parameters'],'baseline_input':self.config['baseline_inputs'][gid],
            'rendered_key_sha256':hashlib.sha256(key.encode()).hexdigest(),'implementation':self.config['implementation'],'resolver':self.config['resolver'],
            'native_helper':self.config['native_helper'],'executables':self.config['executables'],'training_result':self.config['training_result'],
            'water_result':self.config['water_result'],'transfer_result':self.config['transfer_result']}
        return digest(authority),authority,key,render

    def _gate(self):
        require(not self.stop_event.is_set(),'Global native stop is active')
        require(os.environ.get('SLURM_JOB_ID','').isdigit(),'New native calls require Slurm')
        require(int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'One Slurm task with exactly eight CPUs required')
        require(all(os.environ.get(key)=='1' for key in THREADS),'Every native numerical thread control must equal one')
        require(self.root_release is not None,'Explicit root release required');release=read(self.root_release)
        require(release['status']=='approved_paired_donor_development','Wrong root release')
        require(release['evaluator']==self.config['implementation'] and release['evaluator_config_sha256']==self.fingerprint,'Root release does not bind evaluator/config')
        require(release['max_parameter_vectors']==481 and release['max_native_program_invocations']==9216,'Root limits differ')
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
            elif render['baseline_key_reused']:
                record={'status':'complete','cache_identity':authority,'cache_key':token,'cell_id':gid,'first_requested_vector':values[:],
                    'energy':self.baselines[gid]['energy'],'forces':self.baselines[gid]['forces'],'native':self.baselines[gid]['native'],
                    'native_authority':self.config['transfer_result'],'baseline_input':self.config['baseline_inputs'][gid],'source_input':self.config['source_inputs'][gid],'reused_baseline':True,'receipts':[],'render':render,'baseline_origin':self.baseline_review[gid]}
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
                                self._gate();require(self._ledger['native_program_invocations']<9216,'Native call ceiling reached');self._ledger['native_program_invocations']+=1;save(self._ledger_path,self._ledger)
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
                'baseline_origin':record.get('baseline_origin'),'native_cache_key':None if record['reused_baseline'] else record['cache_key'],
                'native_implementation':None if record['reused_baseline'] else record['cache_identity']['implementation'],
                'request_cache_implementation':self.config['implementation']}
            return value

    def evaluate(self,values):
        values=self.checked_values(values);token=digest(values)
        with self.lock:
            if token not in self._ledger['vector_keys']:
                require(len(self._ledger['vector_keys'])<481,'Parameter-vector ceiling reached');self._ledger['vector_keys'].append(token)
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
                'training_result':self.config['training_result'],'water_result':self.config['water_result'],'transfer_result':self.config['transfer_result'],
                'water_invariants_passed':self.water_invariants['passed'],'ledger':pin(self._ledger_path),
                'native_child_CPU_receipts_may_overlap_between_concurrent_workers':True,
                'validation_outcomes_opened':False,'physical_accuracy_qualified':False}
