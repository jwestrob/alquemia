"""Four new radial probes; immutable zero-metal views and admitted native journal."""
import copy,hashlib,json,math,os,re,socket,types
from pathlib import Path
import numpy as np

def require(ok,msg):
    if not ok:raise ValueError(msg)
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def radial(values,metal):return values[12 if metal=='La' else 13] if len(values)==14 else 0.

def evaluator_class(old,base,checks,release_gate):
    facade=types.SimpleNamespace(**vars(base));facade.effective=lambda v,m,d:base.effective(v[:12],m,d)+[radial(v,m)]
    Parent=old.evaluator_class(facade)
    class RadialEvaluator(Parent):
        def __init__(self,cfg,path,mode,inv,center,parent_config):
            require(mode in ('prepare','run') and not Path(path).exists(),'Fresh four-vector cache')
            self.batch=copy.deepcopy(cfg);self.mode=mode;self.batch_fingerprint=cfg['fingerprint'];self.allowed={r['key']:{'parameters':r['parameters'],'role':'radial_probe','ordinal':r['ordinal']} for r in inv['new_probes']}
            self.complete_vectors={};self.imported_vector_keys=set();self.registry={};self.loaded_saved_vectors={};self.saved_metal_index={};self.old_config=copy.deepcopy(parent_config)
            self.executor_identity={'mode':'login','job_id':None,'hostname':socket.gethostname(),'controller_pid':os.getpid(),'run_id':cfg['run_id'],'implementation':cfg['implementation'],'adapter':cfg['adapter'],'parent_implementation':cfg['parent_evaluator']}
            runtime=copy.deepcopy(parent_config);runtime.update(source_manifest=cfg['native_source_manifest'],executables=cfg['executables'],build_receipt=cfg['build'],qualification_review=cfg.get('qualification_review'))
            base.Evaluator.__init__(self,runtime,path,root_release=cfg.get('root_release'))
            self.family['radial_extension']={'append_order':['La_radial_drive','Dy_radial_drive'],'bounds':[[-1.,1.],[-1.,1.]],'initial':[0.,0.],'original12_family_preserved':True}
            self.initial_vector+= [0.,0.];self.program_attempts=[];self.cell_attempts={};self.started_cells=set();self.imported_cell_keys=set();self.center=copy.deepcopy(center)
            for cid,record in center['cells'].items():
                row=self.m['rows'][cid];oldkey,oldmeta=base.render(center['parameters'],row,self.domain,self.h['resolver'])
                q={'cell_id':cid,'effective_values':oldmeta['physical_values_a_b_c_d_R_epsilon'],'xyz':row['source_input']['xyz'],'parameters':row['source_input']['parameters'],'key_sha256':hashlib.sha256(oldkey.encode()).hexdigest(),'implementation':parent_config['implementation'],'source_manifest':parent_config['source_manifest'],'executables':parent_config['executables'],'helpers':{k:parent_config[k] for k in base.HELPERS}}
                require(record['status']=='complete' and record['cache_identity']==q and record['cache_key']==digest(q) and record['request_provenance']['render']==oldmeta,'Original immutable center cell identity')
                for program,inp in record['program_inputs'].items():
                    expected=base.render(center['parameters'],row,self.domain,self.h['resolver'],row['needs_export'] if program=='analyze' else False)[0]
                    require(base.verify(inp['key']).read_text()==expected,'Original saved rendered key')
                    for k in ('xyz','parameters'):require(base.verify(inp[k]).read_bytes()==base.verify(row['source_input'][k]).read_bytes(),'Original literal geometry/parameters')
                for receipt in record['receipts']:base.verify(receipt['output']);require(receipt['returncode']==0,'Original successful native receipt')
                token,identity,key,meta=self.cache_identity(center['parameters']+[0.,0.],cid)
                require(key==oldkey+'POLAR-RADIAL-DRIVE '+str(500 if row['metal']=='La' else 501)+' 36 0\n','Only zero record added to requested key')
                self.cache[token]=copy.deepcopy(record);self.imported_cell_keys.add(token)
            require(len(self.imported_cell_keys)==46,'Exactly46 immutable zero-center cells')
            self.ledger.update(limits=copy.deepcopy(cfg['limits']),new_policy=cfg['adapter'],configuration_fingerprint=self.batch_fingerprint,native_returncodes_observed=0,initial_imported_unique_cells=46)
            self.persistence=old.J.Journal(self.cache_dir/'EVENTS.jsonl',self.ledger,self.lock,{'implementation':cfg['adapter'],'configuration_fingerprint':self.batch_fingerprint,'root_release':cfg.get('root_release')});self.persistence.materialize(self.ledger_path,expected=self.ledger)
            for r in inv['new_probes']:require(self.required_programs(r['parameters'])==28,'Exactly28 programs after verified unchanged-metal reuse')

        def checked_values(self,v):
            require(len(v)==14 and all(isinstance(e,(int,float)) and not isinstance(e,bool) and math.isfinite(e) and -1<=e<=1 for e in v[12:]),'Fourteen finite declared parameters')
            return base.checked(v[:12],self.domain)+[float(e) for e in v[12:]]

        def render(self,values,cid,exports=None):
            values=self.checked_values(values);row=self.m['rows'][cid];text,meta=base.Evaluator.render(self,values[:12],cid,exports)
            e=radial(values,row['metal']);mt=500 if row['metal']=='La' else 501;text+='POLAR-RADIAL-DRIVE '+str(mt)+' 36 '+format(e,'.17g')+'\n'
            meta.update(radial_drive=e,radial_adapter=self.batch['adapter'],original12_render_unchanged=True);return text,meta

        def cache_identity(self,values,cid):
            _,q,key,meta=base.Evaluator.cache_identity(self,values,cid);q.update(radial_adapter=self.batch['adapter'],radial_drive=meta['radial_drive']);return digest(q),q,key,meta

        def _saved_metal(self,*args):return None
        def _import_saved(self,*args):raise ValueError('Only explicitly verified center zero views may be imported')
        def authorize(self,*args):raise ValueError('Four fixed radial probes only')

        def _cell(self,values,cid):
            token,identity,key,meta=self.cache_identity(values,cid)
            if token in self.imported_cell_keys:
                require(radial(values,self.m['rows'][cid]['metal'])==0.,'No nonzero cross-binary reuse')
                answer=copy.deepcopy(self.cache[token])
                answer['radial_zero_reuse']={'old_artifact':self.batch['center'],'root_cache_bridge':self.batch.get('root_cache_bridge'),'requested_identity':identity,'requested_cache_key':token,'requested_vector':list(values),'requested_render':meta,'original_producer_and_raw_cell_retained':True,'new_native_calls':0}
                return answer
            require(facade.effective(values,self.m['rows'][cid]['metal'],self.domain)!=facade.effective(self.initial_vector,self.m['rows'][cid]['metal'],self.domain),'No implicit inherited baseline fallback')
            return super()._cell(values,cid)

        def _parse(self,record,row):
            parsed=base.Evaluator._parse(self,record,row)
            e=record['render']['radial_drive'];physical=record['render']['physical_values_a_b_c_d_R_epsilon'];rr=copy.deepcopy(row);rr.update(metal_type=500 if row['metal']=='La' else 501,physical=physical,physical_POLPAIR=row['fixed_POLPAIR'],e=e)
            h={'matrix':self.h['scout_helper']};m={'units':self.m['units']};readbacks={}
            for program,inp in record['program_inputs'].items():
                output=base.verify(inp['xyz']).parent/(program+'.out');job={'program':program,'stage':'active','profile':'active','cell':record['cell_id'],'fd':False}
                readbacks[program]=checks.readback(output,job,rr,m,h,base)
            parsed['radial_response_readback']=readbacks;parsed['radial_adapter']=self.batch['adapter'];return parsed

        def _gate(self):
            require(self.mode=='run' and not self.stop_event.is_set(),'No native in preparation or after failure')
            require(self.ledger['native_program_invocations']<112,'Hard112 invocation-boundary ceiling')
            require(all(os.environ.get(k)=='1' for k in self.batch['threads']) and socket.gethostname()==self.batch['hostname'] and not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'Truthful12-worker login/one-thread native')
            if not self._gated:release_gate(self.batch);self._gated=True

        def evaluate(self,values):
            values=self.checked_values(values);key=digest(values);require(key in self.allowed and self.allowed[key]['parameters']==values,'Only four exact radial vectors')
            require(key not in self.ledger['vector_keys'] and len(self.ledger['vector_keys'])<4 and self.required_programs(values)==28 and self.ledger['native_program_invocations']+28<=112,'No retry/whole-vector28 budget')
            before=self.ledger['native_program_invocations'];answer=super().evaluate(values)
            require(self.ledger['native_program_invocations']-before==28 and len(answer['cells'])==46,'Exact changed-metal native denominator');return answer

        def snapshot(self):
            with self.lock:self.persistence.materialize(self.ledger_path,expected=self.ledger)
            s=base.Evaluator.snapshot(self);s.update(active_limits=self.batch['limits'],new_parameter_vectors=len(self.ledger['vector_keys']),native_returncodes_observed=self.ledger['native_returncodes_observed'],imported_unique_cells=46,old_mutable_ledger_read=False,executor=copy.deepcopy(self.executor_identity));return s
    return RadialEvaluator
