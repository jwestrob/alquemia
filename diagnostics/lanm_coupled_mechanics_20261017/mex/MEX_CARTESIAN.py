#!/usr/bin/env python3
"""Full-coordinate Mex composite adapter; inherited calculators, new coordinate contract.

Metadata-only --stage bind writes a pinned config. All other stages require Slurm.
The old seven-mode evaluator supplies calculators only: new q=[] records are
wrapped with explicit Cartesian provenance and never admitted through its q7 gate.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
HERE = Path(__file__).resolve().parent
POOL = ROOT/'workspaces/lanm_coupled_response_20261016/mex_continuation_pool_v1/execution/RESULT.json'
PROTOCOL = 'mex_all_physical_cartesian_same_energy_v1'
COMPONENTS = ('MACE_vacuum','GFN2_field_difference','MM_ledger','fixed_FF_charge_ALPB')
SETTINGS = {'FD_steps_A':[0.0002,0.0001], 'FD_absolute_kcal_mol_A':0.005,
            'FD_relative':0.0001, 'replay_energy_kcal_mol':0.000001,
            'replay_gradient_kcal_mol_A':0.000001, 'maximum_new_composites':34,
            'directions':['core_backbone','boundary_anchor','external_water','spectator'],
            'physical_atoms':2049, 'electronic_atoms':209, 'external_charges':1840}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'], 'Changed pinned artifact: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write('\n')


def load(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def bind(path):
    """Read/hash existing metadata only; no molecular or geometric evaluation."""
    pool=json.loads(POOL.read_text()); pcfg=read(pool['config'])
    starts={row['metal']:row['candidate'] for row in pool['terminal_slots']}
    require(pool['status']=='complete' and len(starts)==2,'Old common5 unavailable')
    cfg={'protocol':PROTOCOL,'implementation':pin(__file__),
         'plan':pin(HERE.parent/'PLAN.md'),'specification':pin(HERE/'SPEC.md'),
         'pool_result':pin(POOL),'pool_implementation':pool['implementation'],
         'base_implementation':pcfg['frozen_scout_implementation'],
         'starts':starts,'settings':SETTINGS,
         'scope':'Same Mex model, all 2049 physical coordinates; conditional EF3 exchange with three Nd spectators',
         'DFT_calls':0,'production_change':False}
    write(path,cfg)
    return cfg


class MexEvaluator:
    """Reusable total-energy/gradient evaluation at arbitrary physical coordinates.

    Numeric use requires the pinned MACE interpreter in a Slurm allocation.
    Constructor replays old artifacts only. evaluate(..., fresh=True) deliberately
    bypasses this object's exact cache for qualification, never overwriting it.
    """
    def __init__(self, config, output):
        require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
        self.cfg=read(config) if isinstance(config,dict) and set(config)=={'path','sha256'} else config
        cfg=self.cfg
        require(cfg['protocol']==PROTOCOL and cfg['settings']==SETTINGS,'Cartesian protocol differs')
        require(pin(__file__)==cfg['implementation'],'Cartesian implementation differs')
        verify(cfg['plan']); verify(cfg['specification'])
        self.out=Path(output); self.out.mkdir(parents=True,exist_ok=False)
        self.pool=read(cfg['pool_result']); require(self.pool['status']=='complete','Parent pool unavailable')
        self.pcfg=read(self.pool['config'])
        require(cfg['pool_implementation']==self.pool['implementation'],'Pool source differs')
        poolmodule=load('mex_cartesian_old_pool',verify(cfg['pool_implementation']))
        self.base=load('mex_cartesian_base',verify(cfg['base_implementation']))
        require(cfg['base_implementation']==self.pcfg['frozen_scout_implementation'],'Base lineage differs')
        self.loaded=self.base.inputs(self.pcfg,self.out)
        response,fcfg,world,mech,pilot,points=self.loaded
        _,_,terminal_rows,terminal_slots=poolmodule.parents(self.base,self.pcfg,pilot,points,mech)
        require(len(terminal_rows)==2 and all(s['status']=='available' for s in terminal_slots),'Terminal starts unavailable')
        require({s['metal']:s['candidate'] for s in terminal_slots}==cfg['starts'],'Frozen terminal pin differs')
        self.np=world[6]; self.ids=mech.ids; self.index={s:i for i,s in enumerate(self.ids)}
        self.counters={k:0 for k in ('new_composite_evaluations_attempted','new_MACE_calls_attempted',
            'new_GFN2_calls_attempted','new_MM_configurations_attempted','new_MM_group_queries_attempted',
            'new_ALPB_calls_attempted','new_native_model_loads')}
        self.counters['observed_GFN2_child_CPU_seconds']=0.
        oldcfg=dict(self.pcfg,implementation=cfg['base_implementation'])
        self.old=self.base.Evaluator(oldcfg,self.out,self.counters,self.loaded)
        # This gate applies only to the finite archived q7 points; no Cartesian
        # result ever passes through old.admit_saved().
        for row in self.pool['matrix']:
            for cell in row['cells'].values():
                self.old.admit_saved(cell['point'])
        self.starts={m:self.old.coordinates(read(cfg['starts'][m])) for m in ('La','Dy')}
        self.state_metadata={}
        for metal in ('La','Dy'):
            task=world[11][metal]['GFN2']['embedded']['task']
            self.state_metadata[metal]={
                'metal':metal,'oxidation_state':3,'QM_charge':0,'whole_system_charge':2,
                'physical_multiplicity':1 if metal=='La' else 6,
                'MACE_backend_spin_convention':'multiplicity 2S+1 in atoms.info.spin; inherited audited native adapter',
                'MACE_spin_value':1 if metal=='La' else 6,
                'GFN2_state':{k:task[k] for k in ('charge','multiplicity','uhf','physical_multiplicity','electronic_state','valence_electrons')},
                'GFN2_task':world[1][metal]['field_rows']['embedded'],
                'spectators':[{'id':a['id'],'element':'Nd','oxidation_state':3,'charge_e':3.,
                    'site':a['site'],'position_policy':'all physical coordinates may move',
                    'representation':'formal +3 charge plus explicit pure 12-6 TIP3P IOD LJ from the pinned MM ledger'}
                    for a in world[9]['spectator_occupancy']],
                'fixed_charge_response':'No environment-charge polarization; fixed charges move with their atoms'}
        self.identity={'protocol':PROTOCOL,'wrapper':cfg['implementation'],'parent_pool':cfg['pool_result'],
            'base':cfg['base_implementation'],'source_inputs':world[1]['La']['mapping']['inputs'],
            'MACE_config':world[3],'GFN_runtime':world[8][4],
            'MM_ledgers':{m:world[1][m]['classical']['ledger'] for m in ('La','Dy')},
            'solvent_config':world[2],'solvent_binary':world[12]['build']['binary'],
            'force_assembly':pin(world[5].__file__), 'physical_source_ids':self.ids,
            'state_metadata':self.state_metadata}
        self.identity_hash=hashlib.sha256(canonical(self.identity)).hexdigest()
        self._frozen_states={m:canonical(s) for m,s in self.state_metadata.items()}
        # Replace only the native cache address, retaining all exact archived cells.
        points_to_import=list(self.old.cache.values()); self.old.cache={}
        self.old.key=self.key
        for record in points_to_import:
            point=read(record); self.old.cache[self.key(point['metal'],self.old.coordinates(point))]=record
        self.requests=[]; self.new_points=[]
        write(self.out/'IDENTITY.json',self.identity)
        self.fixtures={m:{name:vector.tolist() for name,vector in directions(self,m).items()} for m in ('La','Dy')}
        write(self.out/'FD_DIRECTIONS.json',{'rule':SETTINGS['directions'],'starts':cfg['starts'],'vectors':self.fixtures,
            'source_ids':self.ids,'selected_before_fresh_calls':True,'new_molecular_calls':0})
        key_checks=self.metadata_cache_checks()
        write(self.out/'CARTESIAN_PREFLIGHT.json',{
            'status':'passed','starts':cfg['starts'],'physical_atoms':len(self.ids),
            'identity_sha256':self.identity_hash,'new_molecular_calls':0,
            'archived_cache_points':len(self.old.cache),'settings':SETTINGS,
            'directions':pin(self.out/'FD_DIRECTIONS.json'),'metadata_cache_checks':key_checks,
            'old_mode_metadata_policy':'Inherited modes/q7 describe the archived source only. New q=[] has no geometric constraint; physical coordinates are authoritative.',
            'state_metadata':self.state_metadata})

    @property
    def new_evaluations(self):
        return self.counters['new_composite_evaluations_attempted']

    def key(self,metal,positions):
        require(metal in ('La','Dy'),'Unsupported target state')
        y=self.np.asarray(positions,dtype='<f8')
        require(y.shape==(2049,3) and self.np.isfinite(y).all(),'Invalid physical coordinates')
        return self._key_from_identity(self.identity_hash,self.state_metadata[metal],metal,y)

    def _key_from_identity(self,identity,state,metal,y):
        state_hash=hashlib.sha256(canonical(state)).hexdigest()
        return hashlib.sha256((identity+':'+state_hash+':'+metal).encode()+y.tobytes(order='C')).hexdigest()

    def metadata_cache_checks(self):
        checks=[]
        for metal in ('La','Dy'):
            y=self.starts[metal]; original=self.key(metal,y)
            require(original==self.key(metal,y.copy()),'Unchanged coordinate cache differs')
            for name,field in (('charge','QM_charge'),('physical_spin','physical_multiplicity'),('backend_spin','MACE_spin_value')):
                state=json.loads(self._frozen_states[metal]);state[field]+=1
                changed=self._key_from_identity(self.identity_hash,state,metal,y)
                checks.append({'metal':metal,'change':name,'invalidated':changed!=original})
            altered=dict(self.identity,GFN_runtime={'deliberate_metadata_only_change':self.identity['GFN_runtime']})
            altered_hash=hashlib.sha256(canonical(altered)).hexdigest()
            changed=self._key_from_identity(altered_hash,self.state_metadata[metal],metal,y)
            checks.append({'metal':metal,'change':'parameter_or_protocol_identity','invalidated':changed!=original})
        require(all(c['invalidated'] for c in checks),'Incomplete cache state/protocol identity')
        return checks

    def structural_input_check(self,positions):
        """Evaluator-domain check only; does not replace optimizer physical guards."""
        np=self.np; y=np.asarray(positions,dtype=float)
        require(y.shape==(2049,3) and np.isfinite(y).all(),'Invalid full physical array')
        distances=[]
        for cap in self.old.cells['La']['mapping']['caps']:
            r,o=cap['retained_index'],cap['omitted_index']
            distance=float(np.linalg.norm(y[r]-y[o])); distances.append(distance)
            require(distance>0.25,'Singular/near-singular boundary anchors')
        return {'passed':True,'scope':'finite coordinates and nonsingular boundary only; optimizer must supply chemistry/state guard',
                'cap_anchor_distances_A':distances,'coordinate_mode':'all_physical_cartesian'}

    def components(self,native):
        old=self.old.old
        return {'MACE_vacuum':native['MACE_eV']*old.EV,
                'GFN2_field_difference':(native['GFN2_embedded_Eh']-native['GFN2_isolated_Eh'])*old.HA,
                'MM_ledger':math.fsum(native['MM_groups_kcal_mol'].values()),
                'fixed_FF_charge_ALPB':native['ALPB_Eh']*old.HA}

    def package(self,point,record,positions,reused):
        gradient=read(point['assembled_gradient'])
        require(gradient['source_ids']==self.ids,'Gradient physical identity differs')
        g=self.np.asarray(gradient['gradient_kcal_mol_A']); components=self.components(point['native_energies'])
        require(g.shape==(2049,3) and self.np.isfinite(g).all(),'Gradient unavailable')
        require(set(gradient['component_gradients_kcal_mol_A'])==set(components),'Component gradient ledger differs')
        # Numeric assembly is replayed rather than inferred from endpoint success.
        closure=sum(self.np.asarray(gradient['component_gradients_kcal_mol_A'][k]) for k in COMPONENTS)
        require(self.np.array_equal(g,closure),'Total/component gradient closure differs')
        return {'total_energy_kcal_mol':math.fsum(components.values()),'gradient_kcal_mol_A':g,
                'component_energies_kcal_mol':components,'native_energies':point['native_energies'],
                'component_gradients_kcal_mol_A':gradient['component_gradients_kcal_mol_A'],
                'provenance':{'protocol':PROTOCOL,'point':record,'identity_sha256':self.identity_hash,
                    'cache_key':self.key(point['metal'],positions),'reused':reused,
                    'state':self.state_metadata[point['metal']]}}

    def evaluate(self,positions,metal,guard=None,fresh=False):
        np=self.np; y=np.asarray(positions,dtype=float); basic=self.structural_input_check(y)
        require(metal in self._frozen_states and canonical(self.state_metadata[metal])==self._frozen_states[metal],
                'State mutation requires a separately bound supported protocol')
        if guard is not None:
            require(guard.get('pass',guard.get('passed',False)),'Caller geometry guard failed')
        key=self.key(metal,y); previous=self.old.cache.pop(key,None) if fresh else None
        before=self.new_evaluations
        try:
            point,record=self.old.evaluate(metal,np.empty(0),y,dict(basic,caller_guard=guard,
                 coordinate_protocol=PROTOCOL,identity_sha256=self.identity_hash))
            require(np.array_equal(self.old.coordinates(point),y),'Returned coordinates differ')
            value=self.package(point,record,y,self.new_evaluations==before)
            event={'protocol':PROTOCOL,'point':record,'state':self.state_metadata[metal],
                   'identity_sha256':self.identity_hash,'cache_key':key,'fresh_requested':fresh,
                   'all_physical_coordinates_A':y.tolist(),'gradient':point['assembled_gradient'],
                   'source_ids':self.ids,'native_energies':point['native_energies'],
                   'reused':value['provenance']['reused'],'coordinate_dimensions':6147,
                   'old_modes_are_archival_only':True}
            path=self.out/('CARTESIAN_REQUEST_%04d.json'%len(self.requests));write(path,event)
            self.requests.append(pin(path))
            if self.new_evaluations>before:self.new_points.append(pin(path))
            return value
        finally:
            if fresh and previous is not None:
                self.old.cache[key]=previous

    def difference(self,new,old):
        parts=self.old.response.work(new['native_energies'],old['native_energies'],self.old.old)
        return dict(parts,composite_kcal_mol=math.fsum(parts.values()))

    def restore_cartesian(self,record):
        """Explicit checkpoint import, without the inherited q7/cache guard."""
        saved=read(record)
        require(saved['protocol']==PROTOCOL and saved['identity_sha256']==self.identity_hash,
                'Checkpoint model/protocol differs')
        require(saved['source_ids']==self.ids,'Checkpoint physical inventory differs')
        metal=saved['state']['metal'];require(saved['state']==self.state_metadata[metal],'Checkpoint state differs')
        point=read(saved['point']);positions=self.np.asarray(saved['all_physical_coordinates_A'])
        require(point['status']=='complete' and point['metal']==metal and point['native_energies']==saved['native_energies'],
                'Checkpoint native cell differs')
        require(self.np.array_equal(self.old.coordinates(point),positions),'Checkpoint physical coordinates differ')
        for name in ('cell','trial','MACE','GFN2_embedded','GFN2_isolated','ALPB','MM_result','assembled_gradient'):
            verify(point[name])
        require(saved['cache_key']==self.key(metal,positions),'Checkpoint environment cache identity differs')
        value=self.package(point,saved['point'],positions,True)
        self.old.cache[saved['cache_key']]=saved['point']
        return value

    def bound(self,metal):
        parent=self
        class Bound:
            @property
            def new_evaluations(self): return parent.new_evaluations
            def evaluate(self,positions): return parent.evaluate(positions,metal)
            def difference(self,new,old): return parent.difference(new,old)
        return Bound()


def directions(evaluator,metal):
    """Source/identity-defined fixtures fixed before evaluating displaced states."""
    np=evaluator.np; y=evaluator.starts[metal]; ix=evaluator.index
    def unit(v):
        norm=float(np.linalg.norm(v));require(norm>1e-8,'Degenerate fixture direction');return v/norm
    def single(sid,v):
        result=np.zeros_like(y);result[ix[sid]]=unit(v);return result
    result={}
    # Peptide C motion is no longer a rigid sidechain or axis endpoint motion.
    result['core_backbone']=single('A/84//C',y[ix['A/84//O']]-y[ix['A/84//C']])
    u=unit(y[ix['A/83//CA']]-y[ix['A/83//C']]);v=y[ix['A/83//N']]-y[ix['A/83//CA']]
    result['boundary_anchor']=single('A/83//CA',v-u*(u@v))
    physical=evaluator.old.cells[metal]['physical'];field={a['id'] for a in evaluator.old.cells[metal]['mapping']['field']}
    native=read(evaluator.old.export['atoms']); groups={}
    for atom in native:
        if atom['id'] in field: groups.setdefault(atom['id'].rsplit('/',1)[0],[]).append(atom)
    waters=[]
    for rid,atoms in groups.items():
        if len(atoms)==3 and sorted(a['element'] for a in atoms)==['H','H','O']:
            waters.append((rid,atoms))
    require(waters,'No wholly external physical water')
    rid,atoms=sorted(waters)[0];oxygen=next(a['id'] for a in atoms if a['element']=='O')
    target=evaluator.old.cells[metal]['target_id']; vector=unit(y[ix[oxygen]]-y[ix[target]])
    water=np.zeros_like(y)
    # Norm one in full Cartesian space; h is full displacement path length.
    for atom in atoms:water[ix[atom['id']]]=vector/math.sqrt(3.)
    result['external_water']=water
    spectator=sorted(a['id'] for a in physical if a['element']=='Nd')[0]
    result['spectator']=single(spectator,y[ix[spectator]]-y[ix[target]])
    return result


def qualification(ev,result):
    np=ev.np; slots=[];result['checks']=slots; result['all_declared_passed']=False
    result['declared_checks']=[{'metal':m,'check':name,'status':'not_evaluated'}
        for m in ('La','Dy') for name in ['fresh_central_replay',*SETTINGS['directions']]]
    for metal in ('La','Dy'):
        y=ev.starts[metal]; saved=ev.evaluate(y,metal); repeat=ev.evaluate(y.copy(),metal)
        require(saved['provenance']['cache_key']==repeat['provenance']['cache_key'] and
                saved['native_energies']==repeat['native_energies'] and
                np.array_equal(saved['gradient_kcal_mol_A'],repeat['gradient_kcal_mol_A']), 'Unchanged cache replay differs')
        fresh=ev.evaluate(y,metal,fresh=True); work=ev.difference(fresh,saved)
        errors={k:abs(v) for k,v in work.items()};gerr=float(np.max(np.abs(fresh['gradient_kcal_mol_A']-saved['gradient_kcal_mol_A'])))
        slots.append({'metal':metal,'check':'fresh_central_replay','work_errors_kcal_mol':errors,
                      'gradient_max_error_kcal_mol_A':gerr,'passed':max(errors.values())<=SETTINGS['replay_energy_kcal_mol'] and gerr<=SETTINGS['replay_gradient_kcal_mol_A']})
        next(c for c in result['declared_checks'] if c['metal']==metal and c['check']=='fresh_central_replay')['status']='passed' if slots[-1]['passed'] else 'failed'
        fixtures={name:np.asarray(vector) for name,vector in ev.fixtures[metal].items()}
        for name,v in fixtures.items():
            desc={'metal':metal,'check':name,'moving_source_ids':[ev.ids[i] for i in np.flatnonzero(np.any(v!=0,axis=1))],
                  'direction_A_per_A':v.tolist(),'steps':[]}
            slots.append(desc)
            predicted={k:float(np.sum(np.asarray(g)*v)) for k,g in fresh['component_gradients_kcal_mol_A'].items()}
            predicted['composite_kcal_mol']=float(np.sum(fresh['gradient_kcal_mol_A']*v))
            for step in SETTINGS['FD_steps_A']:
                minus=ev.evaluate(y-step*v,metal); plus=ev.evaluate(y+step*v,metal)
                require(minus['provenance']['cache_key']!=plus['provenance']['cache_key']!=saved['provenance']['cache_key'],'Coordinate cache did not invalidate')
                numeric={k:d/(2*step) for k,d in ev.difference(plus,minus).items()}
                parts={}
                for k,analytic in predicted.items():
                    residual=numeric[k]-analytic;tolerance=SETTINGS['FD_absolute_kcal_mol_A']+SETTINGS['FD_relative']*max(abs(numeric[k]),abs(analytic))
                    parts[k]={'analytic':analytic,'finite_difference':numeric[k],'residual':residual,'tolerance':tolerance,'passed':abs(residual)<=tolerance}
                desc['steps'].append({'h_A':step,'components_kcal_mol_A':parts,'passed':all(v['passed'] for v in parts.values())})
            desc['passed']=all(s['passed'] for s in desc['steps'])
            next(c for c in result['declared_checks'] if c['metal']==metal and c['check']==name)['status']='passed' if desc['passed'] else 'failed'
    require(ev.new_evaluations<=SETTINGS['maximum_new_composites'],'Qualification call ceiling exceeded')
    result['all_declared_passed']=len(slots)==10 and all(s['passed'] for s in slots)
    result['status']='passed' if result['all_declared_passed'] else 'failed'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',required=True,choices=('bind','preflight','qualification'))
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.stage=='bind': bind(args.config);return
    require(args.output is not None,'--output required')
    start=time.monotonic(); cpu=time.process_time()
    result={'status':'running','stage':args.stage,'config':pin(args.config),'implementation':pin(__file__),
            'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'host':socket.gethostname(),'new_DFT_calls':0}
    ev=None
    try:
        ev=MexEvaluator(json.loads(args.config.read_text()),args.output)
        result['preflight']=pin(args.output/'CARTESIAN_PREFLIGHT.json')
        if args.stage=='qualification':
            require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),'Allocated GPU required')
            require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=2,'At least two allocated CPUs required')
            qualification(ev,result)
        else:result['status']='passed'
    except Exception as exc:
        result.update(status='failed',failure=repr(exc),traceback=traceback.format_exc())
    finally:
        if ev is not None:result.update(counters=ev.counters,requests=ev.requests,new_points=ev.new_points,identity_sha256=ev.identity_hash)
        result.update(elapsed_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
        args.output.mkdir(parents=True,exist_ok=True)
        write(args.output/'RESULT.json',result)
    raise SystemExit(0 if result['status']=='passed' else 1)


if __name__=='__main__':main()
