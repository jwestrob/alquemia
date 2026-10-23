#!/usr/bin/env python3
"""Ten consumed W8/W9 water cells; six pair parameters, no search or objective."""
import argparse
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

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'; D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'
OLD=HERE.parent/'paired_water_six_parameter_v1'
OLD_CELLS=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central') for m in ('La','Dy'))
CELLS=OLD_CELLS+('aquo9_central_La','aquo9_central_Dy')
START=[3.906136716492464,.28782327521600637,3.5865909437359376,.2837564269122712,.20765190997435956,.3173248056679709]
ORIGINAL_INITIAL=[3.6985300908646654,.25,3.3221124657663292,.139]
EPSILON={'La':.24429636467571708,'Dy':.37332330078584813}
NAMES=['La_R','La_damping','Dy_R','Dy_damping','La_epsilon','Dy_epsilon']
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
TRAINING=W/'explicit_model_v1/paired_water_coordination_joint_v1/training_v1/RESULT.json'
PARENT=W/'explicit_model_v1/paired_water_six_parameter_v1/run_v1/execution/RESULT.json'
COORD=W/'aqueous_qualification_design_v1/coordination_screen_v1'

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

def source_pins():
    return {'original_training_result':pin(W/'explicit_model_v1/paired_water_radial_repair_v1/training_v1/RESULT.json'),
        'prior_result':pin(PARENT),'prior_review':pin(W/'explicit_model_v1/paired_water_six_parameter_v1/result_review_v1/REVIEW.json'),
        'coordination_result':pin(COORD/'run_v1/execution/RESULT.json'),'coordination_bound':pin(COORD/'preparation_v2/BOUND.json'),
        'coordination_review':pin(COORD/'actual_review_v1/REVIEW.json'),'coordination_provenance_notice':pin(COORD/'run_v1/PROVENANCE_NOTICE.json'),
        'parent_evaluator':pin(OLD/'EVALUATE.py')}

def checked_sources(s):
    expected={'original_training_result':'3aefeb74f40aaf51e327fdc5c0ffd85f7f00e6e769b07c080f3ece16c345192d',
      'prior_result':'6027e319a261eb85171fb06e30aa3e519411f6371a7905eeca0e8a79812c48fc',
      'prior_review':'a3147adcd7448012437e3ac275872fa58fe65e98a114619a95ec01e17a0df187',
      'coordination_result':'d8de1f0ba01852762f6a012ca241b7072603140eecf68ba07947110d3f04f34d',
      'coordination_bound':'9310cdd2e6ad278c20aa77461ecb3d470a629c270e04602080267a886200aa9d',
      'coordination_review':'6e8e2b4e5f9bb761ad8bbb80e68283724addbe7e471c6351b82dd074f37e0a4f',
      'parent_evaluator':'c73e62c1b6eaba44c3a595c2f14f682c8311c721966d5e97e27d612c6b69c3dd'}
    for k,v in expected.items():require(s[k]['sha256']==v,'Wrong frozen source: '+k);verify(s[k])
    notice=read(s['coordination_provenance_notice'])
    require(notice['reconstruction_attempted'] is False and notice['science_rerun'] is False,'Historical provenance notice changed')
    old=read(s['original_training_result']);prior=read(s['prior_result']);coord=read(s['coordination_result']);bound=read(s['coordination_bound'])
    require(old['status']=='radial_training_complete' and old['complete']==8 and old['saved_controls_passed'],'Old training unavailable')
    require(prior['status']=='bounded_six_parameter_water_complete' and prior['error'] is None and prior['development_checks_passed'] and prior['selected']['parameters']==START,'Wrong exact17605 selection')
    require(coord['status']=='complete' and coord['error'] is None and coord['job_id']=='17670' and coord['native_program_invocations']==4 and coord['bound']==s['coordination_bound'],'W9 execution unavailable')
    require(bound['status']=='prepared_no_calls' and bound['candidate_result']==s['prior_result'] and bound['candidate_parameters']==START,'W9 frozen candidate differs')
    verify(coord['implementation']);verify(bound['implementation']);verify(bound['config']);verify(coord['release'])
    # The pre-submission root review named by that release was overwritten. The
    # intact scientific artifacts and explicit notice are used; no reconstruction.
    return old,prior,coord,bound

def assemble_training(out):
    s=source_pins();old,prior,coord,bound=checked_sources(s);rows=copy.deepcopy(old['rows'])
    original=read(read(bound['config'])['original_manifest']);original={r['id']:r for r in original['rows']}
    for gid in OLD_CELLS:
        r=rows[gid];value=prior['selected']['model'][gid];r['initial_model']=copy.deepcopy(value)
        r['initial_input']=value['provenance']['effective_input'];r['initial_result']=s['prior_result'];r['initial_source_kind']='17605_selected'
    for gid in CELLS[8:]:
        r=copy.deepcopy(bound['rows'][gid]);original_input=copy.deepcopy(original[gid]);r['input']=original_input
        # The original SI key is the rendering origin; the exact17670 candidate
        # key remains the admitted initial input. No candidate override stacks.
        r.update(status='complete',initial_model=copy.deepcopy(coord['rows'][gid]['model']),initial_input=bound['rows'][gid]['files'],
            initial_result=s['coordination_result'],initial_source_kind='17670_W9',reference_result=r['reference']['raw_result'])
        rows[gid]=r
    data={'status':'coordination_joint_training_complete','complete':10,'rows':rows,'sources':s,'initial_parameters':START,
        'old_training_result':s['original_training_result'],'family':copy.deepcopy(prior['family']),
        'orientation_in_training':False,'new_native_calls':0,'physical_accuracy_qualified':False,
        'historical_review_continuity_broken':True,'saved_controls_passed':True}
    out=Path(out);out.mkdir(parents=True,exist_ok=False);save(out/'RESULT.json',data);return pin(out/'RESULT.json')

def default_config():
    s=source_pins();data=read(pin(TRAINING));oldcfg=read(read(read(s['prior_result'])['config'])['evaluator_config'])
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_declaration':pin(HERE/'ROOT_DECLARATION.md'),'training_result':pin(TRAINING),**s,
        'resolver':oldcfg['resolver'],'native_helper':oldcfg['native_helper'],'executables':oldcfg['executables'],
        'source_inputs':{gid:{k:data['rows'][gid]['input'][k] for k in ('xyz','key','parameters')} for gid in CELLS},
        'baseline_inputs':{gid:data['rows'][gid]['initial_input'] for gid in CELLS},
        'initial_vector':START,'original_pair_epsilon':EPSILON,'parameter_names':NAMES,
        'maximum_parameter_vectors':385,'maximum_native_program_invocations':5120,'maximum_native_concurrency':8,'root_release':None}

class Evaluator:
    def __init__(self,config,cache_dir,root_release=None):
        cfg=copy.deepcopy(json.loads(Path(config).read_text()) if isinstance(config,(str,Path)) else config);self.config=cfg
        require(cfg['implementation']==pin(__file__),'Evaluator changed');verify(cfg['specification']);verify(cfg['root_declaration'])
        old,prior,coord,bound=checked_sources(cfg);data=read(cfg['training_result'])
        require(data['status']=='coordination_joint_training_complete' and data['complete']==10 and set(data['rows'])==set(CELLS) and data['sources']=={k:cfg[k] for k in data['sources']},'Ten-cell data binding differs')
        require(data['initial_parameters']==START and not data['orientation_in_training'],'Training start/domain differs')
        require(cfg['initial_vector']==START and cfg['parameter_names']==NAMES and cfg['original_pair_epsilon']==EPSILON,'Family changed')
        require([cfg[k] for k in ('maximum_parameter_vectors','maximum_native_program_invocations','maximum_native_concurrency')]==[385,5120,8],'Finite limits changed')
        self.data=copy.deepcopy(data);self.initial_model={g:copy.deepcopy(r['initial_model']) for g,r in data['rows'].items()};self.continuation_start=START[:]
        self.native=module(cfg['native_helper'],'joint_water_native');self.resolver=module(cfg['resolver'],'joint_water_resolver')
        require(cfg['resolver']['sha256']=='69e1e4dc46cce86f145596180a4507d1e09b3e6473a83210efd99d19e16abba3','Resolver changed')
        parentcfg=read(read(prior['config'])['evaluator_config'])
        for k in ('native_helper','resolver','executables'):require(cfg[k]==parentcfg[k]==bound[k],'Native physical context differs: '+k) if k!='resolver' else require(cfg[k]==parentcfg[k],'Resolver differs')
        for p in cfg['executables'].values():verify(p)
        self.family=copy.deepcopy(prior['family']);self.family.update(initial_vector=START[:],specification=cfg['specification'])
        require([p['name'] for p in self.family['parameters']]==NAMES,'Family order differs')
        for j,p in enumerate(self.family['parameters']):
            origin=ORIGINAL_INITIAL[j] if j<4 else EPSILON['La' if j==4 else 'Dy']
            lo,hi=(.5*origin,1.5*origin) if j>=4 else ((.05,.5) if j%2 else (.85*origin,1.15*origin))
            require((p['lower'],p['upper'])==(lo,hi),'Original bounds changed');p['initial']=START[j]
        self.root_release=root_release if root_release is not None else cfg.get('root_release')
        self.lock=threading.RLock();self.stop_event=threading.Event();self._cell_locks={};self._cache={};self._inherited={};self._admitted=set()
        self.cache_dir=Path(cache_dir).resolve();self.cache_dir.mkdir(parents=True,exist_ok=True);self.fingerprint=config_fingerprint(cfg)
        self._ledger_path=self.cache_dir/'LEDGER.json'
        if self._ledger_path.exists():
            self._ledger=json.loads(self._ledger_path.read_text());require(self._ledger['config_sha256']==self.fingerprint,'Cache context differs')
            if self._ledger['stopped']:self.stop_event.set()
        else:
            self._ledger={'config_sha256':self.fingerprint,'vector_keys':[],'native_program_invocations':0,'stopped':False,'requests':[]};save(self._ledger_path,self._ledger)
        allocated=int(os.environ.get('SLURM_CPUS_PER_TASK','1'))*int(os.environ.get('SLURM_NTASKS','1'))
        self._native_slots=threading.BoundedSemaphore(max(1,min(8,allocated)));self.baseline_review={}
        for gid,row in self.data['rows'].items():
            require(row['status']=='complete' and row['input']['id']==gid,'Training cell unavailable')
            if gid in OLD_CELLS:
                require(all(row[k]==old['rows'][gid][k] for k in ('geometry','input','bodies','reference','reference_result')),'Old8 reference/geometry changed')
                require(row['initial_model']==prior['selected']['model'][gid],'Old8 starting model differs')
            else:
                require(all(row[k]==bound['rows'][gid][k] for k in ('geometry','bodies','reference','reference_result') if k in bound['rows'][gid]),'W9 data differs')
                require(row['initial_model']==coord['rows'][gid]['model'] and row['initial_input']==bound['rows'][gid]['files'],'W9 starting model/input differs')
            for k in ('xyz','key','parameters'):
                require(row['input'][k]==cfg['source_inputs'][gid][k] and row['initial_input'][k]==cfg['baseline_inputs'][gid][k],'Input binding differs');verify(row['input'][k]);verify(row['initial_input'][k])
            require(all(row['input'][k]['sha256']==row['initial_input'][k]['sha256'] for k in ('xyz','parameters')),'Initial geometry/parameters differ')
            g=row['geometry'];lines=verify(row['input']['xyz']).read_text().splitlines();verify(g['xyz'])
            require(int(lines[0].split()[0])==g['atoms']==row['input']['atoms'] and len(lines)==g['atoms']+1,'Typed count differs')
            for i,line in enumerate(lines[1:]):
                t=line.split();require(int(t[0])==i+1 and t[1]==g['symbols'][i] and [float(v) for v in t[2:5]]==g['positions_A'][i] and int(t[5])==row['input']['types'][i] and [int(v)-1 for v in t[6:]]==row['input']['bonds_zero_based'][i],'Typed geometry/graph differs')
            raw=read(row['reference_result']);require(raw['status']=='complete' and raw['returncode']==0 and all(raw['checks'].values()),'Reference unavailable')
            require(raw['forces_kcal_mol_A']==row['reference']['forces'] and math.isclose(raw['energy_hartree']*627.5094740631,row['reference']['energy'],rel_tol=0,abs_tol=1e-7),'Reference E/F changed')
            for k in ('output','engrad','input','panel'):
                if k in row['reference']:verify(row['reference'][k])
            metal=row['input']['metal'];mt=500 if metal=='La' else 501;i=0 if metal=='La' else 2
            p=self.resolver.parse_parameters(verify(row['input']['parameters']).read_text(),verify(row['input']['key']).read_text());pair=self.resolver.resolve_pair(p,mt,36)
            require(row['input']['types'][0]==mt and set(row['input']['types'][1:])=={36,37} and [t for t,a in p['atom'].items() if a['class']==34]==[36],'Water typing differs')
            require([pair['Rmin_pair_A'],pair['thole'],pair['epsilon_kcal_mol']]==ORIGINAL_INITIAL[i:i+2]+[EPSILON[metal]],'Original water pair differs')
            key,render=self.render(START,gid);require(key==verify(row['initial_input']['key']).read_text(),'Starting key is not exact executed input')
            self._replay(row['initial_model'],row)
            self.baseline_review[gid]={'source_result':row['initial_result'],'source_kind':row['initial_source_kind'],'input':row['initial_input'],'resolved_pair':render['resolved_pair'],'raw_replayed':True}
        # Index exact already computed physical inputs once. No old ledger writes,
        # no whole-history reads per call, no orientation outcomes.
        for vector in prior['vectors'].values():
            require(vector['status']=='complete' and set(vector['model'])==set(OLD_CELLS),'Prior vector incomplete')
            for gid,value in vector['model'].items():
                pv=value['provenance'];inp=pv['effective_input'];self._register(gid,vector['parameters'],value,inp,pv['native_authority'],sourceresult=cfg['prior_result'],kind='17605_inventory',producer=pv)
        for gid in CELLS[8:]:
            row=self.data['rows'][gid];actual=coord['rows'][gid];require(actual['status']=='complete' and len(actual['receipts'])==2 and all(r['returncode']==0 for r in actual['receipts']),'W9 native receipts differ')
            runtime={k:pin(Path(cfg['coordination_result']['path']).parent/gid/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            require(all(runtime[k]['sha256']==row['initial_input'][k]['sha256'] for k in runtime),'W9 runtime input differs')
            for receipt in actual['receipts']:verify(receipt['output'])
            self._register(gid,START,row['initial_model'],runtime,cfg['coordination_result'],sourceresult=cfg['coordination_result'],kind='17670_W9',producer={'native_implementation':coord['implementation'],'native_authority':cfg['coordination_result']})

    def checked_values(self,values):
        values=list(values);require(len(values)==6,'Six parameters required')
        for value,p in zip(values,self.family['parameters']):require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and p['lower']<=value<=p['upper'],'Parameter outside original bounds')
        return [float(v) for v in values]

    def render(self,values,gid):
        values=self.checked_values(values);require(gid in CELLS,'Cell outside ten-water panel');row=self.data['rows'][gid]
        metal=row['input']['metal'];mt=500 if metal=='La' else 501;i=0 if metal=='La' else 2;epsilon=values[4 if metal=='La' else 5]
        base=verify(row['input']['key']).read_text();prm=verify(row['input']['parameters']).read_text();before=self.resolver.parse_parameters(prm,base)
        key=base if values[i:i+2]==ORIGINAL_INITIAL[i:i+2] and epsilon==EPSILON[metal] else base.rstrip()+f'\n\n# Paired water development override: {metal} only\nvdwpr 34 {mt} {values[i]:.17g} {epsilon:.17g}\npolpair 36 {mt} {values[i+1]:.17g}\n'
        after=self.resolver.parse_parameters(prm,key);pair=self.resolver.resolve_pair(after,mt,36)
        require([pair['Rmin_pair_A'],pair['thole'],pair['epsilon_kcal_mol']]==values[i:i+2]+[epsilon],'Rendered pair differs')
        for k in ('settings','atom','vdw','polarize','metal_multipoles'):require(before[k]==after[k],'Nonpair term changed')
        for k,allowed in (('vdwpr',(34,mt)),('polpair',(36,mt))):require({a:b for a,b in before[k].items() if a!=allowed}=={a:b for a,b in after[k].items() if a!=allowed},'Nonwater pair changed')
        require(self.resolver.resolve_pair(before,mt,37)==self.resolver.resolve_pair(after,mt,37),'Water-H pair changed')
        return key,{'metal':metal,'effective_values':values[i:i+2]+[epsilon],'resolved_pair':pair,'original_key_reused':key==base,'baseline_key_reused':key==verify(row['initial_input']['key']).read_text(),'nonwater_parameters_unchanged':True}

    def cache_identity(self,values,gid):
        key,render=self.render(values,gid);row=self.data['rows'][gid]
        authority={'cell_id':gid,'effective_values':render['effective_values'],'xyz':row['input']['xyz'],'parameters':row['input']['parameters'],
            'rendered_key_sha256':hashlib.sha256(key.encode()).hexdigest(),'implementation':self.config['implementation'],'resolver':self.config['resolver'],
            'native_helper':self.config['native_helper'],'executables':self.config['executables'],'training_result':self.config['training_result']}
        return digest(authority),authority,key,render

    def _register(self,gid,values,value,inp,authority,sourceresult,kind,producer):
        # Native file hashes and raw data are admitted on first use. All index
        # members are from the exact immutable result inventory verified above.
        i=0 if self.data['rows'][gid]['input']['metal']=='La' else 2;effective=list(values[i:i+2])+[values[4 if i==0 else 5]]
        identity=(gid,tuple(effective),inp['xyz']['sha256'],inp['key']['sha256'],inp['parameters']['sha256'])
        if identity in self._inherited:
            other=self._inherited[identity]['model'];require(all(other[k]==value[k] for k in ('energy','forces','native')),'Conflicting inherited physical result');return
        self._inherited[identity]={'model':copy.deepcopy(value),'input':inp,'authority':authority,'source_result':sourceresult,'kind':kind,'producer':copy.deepcopy(producer)}

    def _replay(self,value,row):
        native=value['native'];parsed=self.native.parse(verify(native['analyze_output']),verify(native['gradient_output']),dict(row['input'],derivative_test=False))
        require(parsed==native and value['energy']==parsed['total_energy_kcal_mol'] and value['forces']==parsed['force_kcal_mol_A'],'Saved native E/F replay differs')

    def _gate(self):
        require(not self.stop_event.is_set(),'Global native stop active')
        require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'Exactly one eight-CPU Slurm task required')
        require(all(os.environ.get(k)=='1' for k in THREADS),'Native threads must equal one')
        require(self.root_release is not None,'Root release required');r=read(self.root_release)
        require(r['status']=='approved_paired_water_coordination_joint' and r['evaluator']==self.config['implementation'] and r['evaluator_config_sha256']==self.fingerprint,'Release binding differs')
        require(r['max_parameter_vectors']==385 and r['max_native_program_invocations']==5120 and not r['validation_scoring_authorized'] and not r['physical_accuracy_claimed'],'Release scope differs')

    def _stopped(self,exc):
        self.stop_event.set()
        with self.lock:self._ledger.update(stopped=True,failure=repr(exc));save(self._ledger_path,self._ledger)

    def _program_admission(self,path,program):
        text=Path(path).read_text();require('Tinker' in text and 'Version 26.2' in text,'Native identity absent')
        require(not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b',text),'Native failure/nonfinite marker')
        energies=re.findall(r'Total Potential Energy\s*:\s*('+self.native.NUMBER+')',text)
        require(len(energies)==1 and math.isfinite(float(energies[0])),'Native finite energy absent')
        if program=='analyze':require('Energy Component Breakdown :' in text,'Native components absent')

    def _cell(self,values,gid):
        token,authority,key,render=self.cache_identity(values,gid)
        with self.lock:lock=self._cell_locks.setdefault(token,threading.Lock())
        with lock:
            row=self.data['rows'][gid];directory=self.cache_dir/'cells'/token;rp=directory/'RESULT.json'
            identity=(gid,tuple(render['effective_values']),authority['xyz']['sha256'],authority['rendered_key_sha256'],authority['parameters']['sha256'])
            if token in self._cache:record=self._cache[token]
            elif identity in self._inherited:
                old=self._inherited[identity];verify(old['authority'])
                for p in old['input'].values():verify(p)
                self._replay(old['model'],row)
                record={'status':'complete',**{k:copy.deepcopy(old['model'][k]) for k in ('energy','forces','native')},'input':old['input'],
                    'native_authority':old['authority'],'native_source_result':old['source_result'],'native_source_kind':old['kind'],'producer_provenance':old['producer'],
                    'native_implementation':old['producer'].get('native_implementation'),'inherited_native':True,'render':render,'first_requested_vector':None}
                self._cache[token]=record;self._admitted.add(identity)
            elif rp.exists():
                record=read(pin(rp));require(record['status']=='complete' and record['cache_identity']==authority,'Cached cell failed/differs; no retry')
                for p in record['input'].values():verify(p)
                require(all(record['input'][k]['sha256']==authority[k]['sha256'] for k in ('xyz','parameters')) and record['input']['key']['sha256']==authority['rendered_key_sha256'],'Cached input changed')
                require(len(record['receipts'])==2 and all(r['returncode']==0 for r in record['receipts']),'Native receipts unavailable');self._replay(record,row);self._cache[token]=record
            else:
                self._gate();directory.mkdir(parents=True,exist_ok=False)
                record={'status':'started','cell_id':gid,'cache_key':token,'cache_identity':authority,'first_requested_vector':values[:],
                    'source_input':self.config['source_inputs'][gid],'render':render,'inherited_native':False,'native_implementation':self.config['implementation'],'receipts':[],'job_id':os.environ['SLURM_JOB_ID']}
                try:
                    shutil.copyfile(verify(row['input']['xyz']),directory/'system.xyz');shutil.copyfile(verify(row['input']['parameters']),directory/'amoeba09.prm');(directory/'system.key').write_text(key)
                    record['input']={k:pin(directory/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
                    require(all(record['input'][k]['sha256']==authority[k]['sha256'] for k in ('xyz','parameters')) and record['input']['key']['sha256']==authority['rendered_key_sha256'],'Copied physical input changed');save(rp,record)
                    for program,args in (('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])):
                        with self._native_slots:
                            with self.lock:
                                self._gate();require(self._ledger['native_program_invocations']<5120,'Native program ceiling');self._ledger['native_program_invocations']+=1;save(self._ledger_path,self._ledger)
                            try:
                                self.native.invoke([self.config['executables'][program]['path']]+args,directory,directory/(program+'.out'),record['receipts']);self._program_admission(directory/(program+'.out'),program)
                                if program=='testgrad':raw=self.native.parse(directory/'analyze.out',directory/'testgrad.out',dict(row['input'],derivative_test=False))
                            except Exception as exc:self._stopped(exc);raise
                            save(rp,record)
                    record.update(status='complete',energy=raw['total_energy_kcal_mol'],forces=raw['force_kcal_mol_A'],native=raw);self._cache[token]=record
                except Exception as exc:record.update(status='failed',error=repr(exc));self._stopped(exc);raise
                finally:save(rp,record)
            value={k:copy.deepcopy(record[k]) for k in ('energy','forces','native')}
            value['provenance']={'requested_vector':values[:],'cache_key':token,'effective_parameters':render['effective_values'],
                'native_authority':record.get('native_authority') or pin(rp),'native_implementation':record.get('native_implementation'),
                'native_source_result':record.get('native_source_result'),'native_source_kind':record.get('native_source_kind'),
                'producer_provenance':record.get('producer_provenance'),'inherited_native':record['inherited_native'],
                'effective_input':record['input'],'source_input':self.config['source_inputs'][gid],'render':render,
                'first_requested_vector':record['first_requested_vector'],'request_cache_implementation':self.config['implementation']}
            return value

    def evaluate(self,values):
        values=self.checked_values(values);token=digest(values)
        with self.lock:
            if token not in self._ledger['vector_keys']:
                require(len(self._ledger['vector_keys'])<385,'Parameter-vector ceiling');self._ledger['vector_keys'].append(token)
            request={'vector':values[:],'vector_key':token,'status':'started','cells':{g:{'status':'not_attempted'} for g in CELLS},'started_monotonic':time.monotonic()}
            self._ledger['requests'].append(request);save(self._ledger_path,self._ledger)
        model={}
        try:
            for gid in CELLS:
                model[gid]=self._cell(values,gid)
                # Durable compact request progress; full E/F/render provenance is
                # already in native cell receipts and returned driver models.
                p=model[gid]['provenance']
                with self.lock:request['cells'][gid]={'status':'complete','cache_key':p['cache_key'],'native_authority':p['native_authority']};save(self._ledger_path,self._ledger)
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
                'stopped':self.stop_event.is_set(),'requests':copy.deepcopy(self._ledger['requests']),'baseline_review':self.baseline_review,
                'inherited_native_cells_indexed':len(self._inherited),'inherited_native_cells_admitted':len(self._admitted),'ledger':pin(self._ledger_path),
                'prior_result':self.config['prior_result'],'coordination_result':self.config['coordination_result'],
                'validation_outcomes_opened':False,'physical_accuracy_qualified':False,'native_child_CPU_receipts_may_overlap_between_concurrent_workers':True}

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['training','config','prepare']);p.add_argument('--out',required=True);p.add_argument('--config');a=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'Preparation needs one numerical thread');start=time.monotonic();cpu=time.process_time()
    if a.action=='training':print(json.dumps(assemble_training(a.out)));return
    if a.action=='config':save(a.out,default_config());print(json.dumps(pin(a.out)));return
    out=Path(a.out);out.mkdir(parents=True,exist_ok=False);cfg=json.loads(Path(a.config).read_text());ev=Evaluator(cfg,out/'cache');model=ev.evaluate(START)
    require(all(all(model[g][k]==ev.initial_model[g][k] for k in ('energy','forces','native')) for g in CELLS),'Initial ten-cell replay differs')
    result={'status':'prepared_no_calls','implementation':pin(__file__),'config':pin(a.config),'training_result':cfg['training_result'],
        'initial_parameters':START,'initial_model':model,'family':ev.family,'evaluator':ev.snapshot(),'new_native_program_invocations':0,
        'initial_ten_raw_replay_passed':True,'orientation_outcomes_opened':False,'historical_review_continuity_broken':True,
        'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1}
    require(ev.snapshot()['native_program_invocations']==0,'Preparation made calls');save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}))
if __name__=='__main__':main()
