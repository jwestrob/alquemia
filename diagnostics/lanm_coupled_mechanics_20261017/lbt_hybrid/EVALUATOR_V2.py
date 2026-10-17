#!/usr/bin/env python3
"""LBT local electronic / MM model, using existing pinned native calculators.

V2 changes only derived cap-Jacobian replay; original physics and gates remain.
Only --stage bind is metadata-only. Every other stage requires Slurm.
No optimization is implemented. Root owns execution and stage admission.
"""
from concurrent.futures import ThreadPoolExecutor
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import sys
import time
import traceback

ROOT=Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
HERE=Path(__file__).resolve().parent
OLD=ROOT/'workspaces/lanm_coupled_response_20261016'
PROTOCOL='lbt_local_electronic_MM_physical_v1'
EV=23.06054783061903; HA=627.509474; BOHR=0.529177210903
COMPONENTS=('MACE_vacuum','GFN2_field_difference','MM_ledger','fixed_FF_charge_ALPB')
STATES=[(s,r,m) for s in ('7CCO','7CCN') for r in ('baseline','expanded') for m in ('La','Dy')]
FD={'steps_A':[0.0002,0.0001], 'absolute_kcal_mol_A':0.005, 'relative':0.0001,
    'directions':['boundary_anchor','external_charge'], 'scout_composites':9}
CONTROL='$embedding\n interface=orca\n input=external.pc\n gradient=pcgrad\n$end\n'


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve(); h=hashlib.sha256()
    with p.open('rb') as stream:
        for b in iter(lambda:stream.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
    return pin(path)


def module(name,p):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def verify_tree(value):
    if isinstance(value,dict):
        if set(value)=={'path','sha256'}:verify(value)
        else:
            for v in value.values():verify_tree(v)
    elif isinstance(value,list):
        for v in value:verify_tree(v)


def environment(settings):
    e=os.environ.copy()
    for k,v in settings.items():
        if v is None:e.pop(k,None)
        else:e[k]=v
    for k in ('OMP_PLACES','OMP_PROC_BIND'):e.pop(k,None)
    return e


def bind(path):
    # Reuse the exact failed configuration; change only versioned replay metadata.
    old_config=pin(ROOT/'workspaces/lanm_coupled_mechanics_20261017/lbt_hybrid_scout_v1/CONFIG.json')
    cfg=read(old_config)
    diagnostic=pin(ROOT/'workspaces/lanm_coupled_mechanics_20261017/lbt_hybrid_mapping_diagnosis_v1/execution/RESULT.json')
    dc=read(read(diagnostic)['config'])
    require(dc['hybrid_config']['sha256']==old_config['sha256'],'Actual diagnosis configuration differs')
    cfg['implementation']=pin(__file__)
    cfg['mapping_replay_helper']=pin(HERE/'MAPPING_REPLAY_V2.py')
    cfg['mapping_recovery_spec']=pin(HERE/'RECOVERY_V2.md')
    cfg['mapping_recovery']={'original_config':old_config,'diagnostic':diagnostic,
        'failed_result':dc['failed_result'],'passed_preflight':dc['passed_preflight']}
    write(path,cfg);return cfg


def xyz_text(symbols,x):
    return str(len(symbols))+'\nLBT fixed physical coordinate mapping\n'+''.join(
        '%s %.17g %.17g %.17g\n'%(s,*r) for s,r in zip(symbols,x))


def mm_worker(cfg,out):
    """Exact prepared native XML, evaluated in the separate pinned MM environment."""
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'],'MM interpreter differs')
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'MM adapter differs')
    import numpy as np
    import openmm as mm
    from openmm import unit
    require(mm.__version__==cfg['version'],'MM version differs')
    system=mm.XmlSerializer.deserialize(verify(cfg['system']).read_text())
    require(system.getNumParticles()==256 and system.getNumConstraints()==0,'MM physical inventory/constraints differ')
    require(all(0<=f.getForceGroup()<4 for f in system.getForces()),'Unexpected MM group')
    x=np.asarray(read(cfg['positions'])['positions_A'],float)
    require(x.shape==(256,3) and np.isfinite(x).all(),'MM coordinate shape')
    # Group names are checked against the frozen native builder, never inferred from energy.
    sys.path.insert(0,str(verify(cfg['builder']).parent))
    helper=module('hybrid_mm_builder',cfg['builder'])
    names=tuple(helper.GROUPS);require(len(names)==4,'Four native groups required')
    integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
    context.setPositions(x*.1*unit.nanometer);energies={};gradients={}
    for i,n in enumerate(names):
        s=context.getState(getEnergy=True,getForces=True,groups=1<<i)
        energies[n]=s.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        gradients[n]=(-np.asarray(s.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom))).tolist()
    require(all(math.isfinite(e) for e in energies.values()) and np.isfinite(list(gradients.values())).all(),'Nonfinite MM output')
    return {'status':'complete','system':cfg['system'],'positions':cfg['positions'],'group_order':list(names),
            'energy_kcal_mol':math.fsum(energies.values()),'component_energies_kcal_mol':energies,
            'gradient_kcal_mol_A':sum(np.asarray(gradients[k]) for k in names).tolist(),
            'component_gradients_kcal_mol_A':gradients,'platform':'Reference','openmm_version':mm.__version__,
            'configurations_attempted':1,'group_queries_attempted':4}


class HybridEvaluator:
    def __init__(self,cfg,output):
        require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
        self.cfg=cfg;require(cfg['protocol']==PROTOCOL and cfg['FD']==FD and cfg['declared_states']==[list(s) for s in STATES],'Protocol differs')
        require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Wrapper source differs')
        verify(cfg['spec']);verify(cfg['plan']);verify(cfg['parameter_loading_source'])
        self.out=Path(output);self.out.mkdir(parents=True,exist_ok=False)
        import numpy as np
        self.np=np; self.native=None;self.requests=[];self.points=[];self.cache={};self.solvent_cache={}
        self.cost={k:0 for k in ('composites_attempted','MACE_calls_attempted','MACE_model_loads',
            'GFN_calls_attempted','MM_configurations_attempted','MM_group_queries_attempted','ALPB_calls_attempted','ALPB_cache_hits')}
        self.cost['observed_native_child_CPU_seconds']=0.
        self.prep=read(cfg['preparation_result']);self.pcfg=read(cfg['preparation_config'])
        require(self.prep['status']=='complete' and self.prep['slurm_job_id']=='12953' and len(self.prep['rows'])==4
            and self.prep['new_molecular_calls']==0 and self.prep['new_MM_energy_calls']==0,'Partition preparation unavailable')
        require(self.prep['config']==cfg['preparation_config'] and self.prep['implementation']==cfg['preparation_implementation'],'Partition lineage differs')
        verify_tree(self.pcfg)
        self.mapping=module('lbt_local_mapping',cfg['preparation_implementation'])
        sys.path.insert(0,str(verify(cfg['common_helper']).parent))
        self.gradient=module('lbt_hybrid_gradient',cfg['gradient_helper'])
        # These two pure exported functions require only these dependencies; do not
        # call initialize(), which also imports OpenMM into the separate MACE venv.
        self.mapping.np=np; self.mapping.cap_jacobians=self.gradient.cap_jacobians
        self.admit=module('lbt_hybrid_gfn_admission',cfg['admission'])
        self.backend=module('lbt_hybrid_native_backend',cfg['GFN_backend'])
        self.warm=module('lbt_hybrid_warm_MACE',cfg['native_MACE_helper'])
        self.mcfg=read(cfg['MACE_config']);mr=read(cfg['MACE_evidence'])
        require(mr['status']=='complete' and mr['slurm_job_id']=='12500' and mr['model_parameters_unchanged']
            and mr['model']==cfg['model']==self.mcfg['model'] and mr['versions']==cfg['expected_MACE_versions']
            and mr['config']==cfg['MACE_config'],'MACE precedent differs')
        verify_tree(self.mcfg['helper_pins']);verify(cfg['model']['checkpoint'])
        require(pin(sys.executable)['sha256']==self.mcfg['python']['sha256'],'MACE interpreter differs')
        self.scfg=read(cfg['GFN_parent_config']);self.runtime=cfg['GFN_runtime'];verify_tree(self.runtime)
        require(self.scfg['runtime']==self.runtime and self.scfg['tolerances']==self.admit.TOL,'Standalone runtime/tolerances differ')
        common=module('lbt_hybrid_units',cfg['common_helper'])
        require(common.HA_TO_KCAL==HA and common.BOHR_TO_A==BOHR,'Inherited gradient units differ')
        self.qualification=read(cfg['GFN_qualification']);qcfg=read(self.qualification['config'])
        manifest=read(self.qualification['manifest']);self.profile=read(cfg['GFN_profile'])
        require(self.qualification['status']=='passed' and self.qualification['slurm_job_id']=='12488'
            and self.qualification['accepted']==self.qualification['attempted']==8
            and all(p['passed'] for p in self.qualification['parity']+self.qualification['solvent_component_parity'])
            and manifest['profile']==cfg['GFN_profile'],'La/Dy runtime qualification differs')
        require(qcfg['staging']['standalone_settings']=={'accuracy':.0001,'maxiter':500,'electronic_temperature_K':300,'alpb_reference':'gsolv'},'GFN settings differ')
        for row in self.qualification['rows']:
            if row['task']['backend']=='standalone':
                require(row['status']=='accepted' and all(row['task'][k]==v for k,v in self.runtime.items()),'Qualified runtime binding differs')
        require(cfg['GFN_profile']['sha256']=='635e111bac31e7c5755bf612c8e148ce73b401c110ff31887fd5e31b590fa85f','Full corrected La+Dy profile differs')
        self.env=environment(self.runtime['environment'])
        self.gam={}
        for b in re.split(r'(?m)^\$Z=\s*',verify(self.runtime['parameters']).read_text())[1:]:
            m=re.search(r'(?m)^\s*GAM=\s*([-+\d.Ee]+)\s*$',b)
            if m:self.gam[int(b.split()[0])]=float(m.group(1))
        require(all(self.gam[self.admit.Z[s]]==self.profile['element'][s]['gam'] for s in ('H','C','N','O','La','Dy')),'Parameter/profile hardness mismatch')
        self.alpb=read(cfg['ALPB_evidence']);self.alpb_helper=module('lbt_hybrid_solvent',cfg['ALPB_helper'])
        require(self.alpb['status']=='complete' and self.alpb['slurm_job_id']=='12715'
            and self.alpb['implementation']==cfg['ALPB_helper'] and self.alpb['build']['binary']==cfg['ALPB_binary'],'ALPB lineage differs')
        verify_tree(self.alpb['build']);self.alpb_env=environment(self.alpb['native_environment'])
        self.mapping_replay=module('lbt_hybrid_mapping_replay_v2',cfg['mapping_replay_helper'])
        verify(cfg['mapping_recovery_spec']);recovery=cfg['mapping_recovery']
        original=read(recovery['original_config']);failed=read(recovery['failed_result']);passed=read(recovery['passed_preflight'])
        require(original['implementation']['sha256']=='a482559ca4febd240a40ec7d6042f4fa24ec4358e4f5ba05770bcfaacea4b784'
            and failed['slurm_job_id']=='12975' and failed['status']=='execution_failure'
            and failed['failure']=="ValueError('Exact origin mapping replay differs')"
            and passed['slurm_job_id']=='12972' and passed['status']=='complete'
            and failed['config']['sha256']==passed['config']['sha256']==recovery['original_config']['sha256'],
            'Original failed/preflight lineage differs')
        verify(original['implementation'])
        require(all(cfg[k]==v for k,v in original.items() if k!='implementation'),'V2 changed original scientific configuration')
        self.mapping_diagnostic_checks=self.mapping_replay.validate_diagnostic(recovery['diagnostic'],read)
        self.mapping_replays=[]
        self.partitions={};self.tasks={};self.origins={};self.states={}
        for row in self.prep['rows']:
            require(row['status']=='prepared','Preparation row unavailable')
            p=read(row['preparation']);partition=read(p['partition']);ledger=read(p['ledger'])
            verify_tree(p);verify_tree(ledger)
            sid,region=row['source_id'],row['region']; nq,nf=(172,84) if region=='baseline' else (210,46)
            require(p['status']=='prepared' and len(p['tasks'])==2 and p['real_QM_atoms']+p['cap_atoms']==nq
                and p['field_atoms']==nf and p['cap_atoms']==2 and all(p['common_MM_component_definition_checks'].values())
                and not ledger['C4_induction_included'] and not ledger['QM_MM_Coulomb_included'],'Partition/ledger contract differs')
            self.partitions[(sid,region)]=(p,partition,ledger)
            particles=read(ledger['artifacts']['particles.json'])
            require([a['id'] for a in particles]==[a['id'] for a in partition['physical_atoms']]
                and len(particles)==256 and particles[-1]['id']=='metal','Native MM/physical particle order differs')
            x=np.asarray([a['xyz_A'] for a in partition['physical_atoms']],float)
            mapping=self.mapping.update_geometry(partition,x)
            map_check=self.mapping_replay.compare(read(p['mapping']),mapping)
            self.mapping_replays.append({'kind':'origin','source_id':sid,'region':region,
                'archived_mapping':p['mapping'],'check':map_check})
            self.origins[(sid,region)]=x
            for t in p['tasks']:
                state=(sid,region,t['metal']);self.tasks[state]=t
                sy,y=self.admit.xyz(t['xyz']);mult=1 if t['metal']=='La' else 6
                ne=sum(sum(self.profile['element'][s]['refocc']) for s in sy)-t['charge']
                require(len(sy)==nq and sy[0]==t['metal'] and t['charge']==-2 and t['multiplicity']==mult
                    and t['physical_spin_2S']==mult-1 and t['metal_index']==0 and t['gfn2_uhf']==0
                    and ne==int(ne) and int(ne)%2==0 and np.array_equal(y,[a['xyz_A'] for a in mapping['core']]),'Electronic state or coordinates differ')
                require(sum(self.admit.Z[s] for s in sy)-t['charge']==t['all_electron_count']
                    and (t['all_electron_count']-mult+1)%2==0,'Physical electron count differs')
                charges=read(t['global_solvent_charges'])
                require(charges['physical_ids']==[a['id'] for a in partition['physical_atoms']]
                    and len(charges['charges_e'])==256 and abs(math.fsum(charges['charges_e'])+2)<=1e-6
                    and charges['charges_e'][-1]==3.,'Whole fixed-charge state differs')
                self.states[state]={'source':sid,'region':region,'metal':t['metal'],'oxidation_state':3,
                    'QM_charge':-2,'whole_charge':-2,'field_charge_e':partition['field_charge'],
                    'MACE_multiplicity':mult,'MACE_spin_convention':'multiplicity 2S+1',
                    'GFN_uhf':0,'GFN_effective_multiplicity':1,'GFN_valence_electrons':int(ne),
                    'physical_all_electrons':t['all_electron_count'],'physical_state_qualified':False,
                    'GFN_spin_independent_f_in_core':True,'DFT_ECP_not_used':True}
        require(set(self.tasks)==set(STATES),'Declared eight-state denominator differs')
        self.identity={'protocol':PROTOCOL,'implementation':cfg['implementation'],'admission':cfg['admission'],
            'spec':cfg['spec'],'plan':cfg['plan'],
            'preparation':cfg['preparation_result'],'mapping':cfg['preparation_implementation'],
            'mapping_replay_helper':cfg['mapping_replay_helper'],'mapping_recovery_spec':cfg['mapping_recovery_spec'],
            'mapping_recovery':cfg['mapping_recovery'],'mapping_replay_policy':self.mapping_replay.POLICY,
            'MACE':cfg['MACE_config'],'MACE_model':cfg['model'],'MACE_helper':cfg['native_MACE_helper'],
            'GFN_runtime':self.runtime,'GFN_profile':cfg['GFN_profile'],'GFN_qualification':cfg['GFN_qualification'],
            'GFN_backend':cfg['GFN_backend'],'gradient_helper':cfg['gradient_helper'],
            'MM_python':cfg['MM_python'],'MM_preparation':cfg['preparation_config'],
            'ALPB_binary':cfg['ALPB_binary'],'ALPB_build':self.alpb['build'],'ALPB_helper':cfg['ALPB_helper'],
            'states':{'/'.join(k):v for k,v in self.states.items()}}
        self.identity_hash=digest(self.identity);self.frozen_states=copy.deepcopy(self.states)
        self.replays=self.replay_saved()
        self.directions=self.freeze_directions()
        self.preflight={'status':'passed','identity_sha256':self.identity_hash,'states':self.identity['states'],
            'declared_origin_cells':8,'scout_composites':9,'new_molecular_calls':0,'new_MM_calls':0,
            'replays':self.replays,'FD_directions':write(self.out/'FD_DIRECTIONS.json',self.directions),
            'cache_checks':self.cache_checks(),'optimization_qualified':False,
            'mapping_replay_checks':self.mapping_replays,'mapping_diagnostic_checks':self.mapping_diagnostic_checks}
        write(self.out/'IDENTITY.json',self.identity);write(self.out/'PREFLIGHT.json',self.preflight)

    def key(self,sid,region,metal,x):
        state=(sid,region,metal);require(self.states[state]==self.frozen_states[state],'State metadata changed after construction')
        a=self.np.asarray(x,dtype='<f8');require(a.shape==(256,3) and self.np.isfinite(a).all(),'Physical coordinates invalid')
        return hashlib.sha256((self.identity_hash+':'+':'.join(state)).encode()+a.tobytes()).hexdigest()

    def solvent_key(self,sid,metal,x,q):
        identity={'binary':self.cfg['ALPB_binary'],'build':self.alpb['build'],'charges':q,
                  'symbols':[a['element'] for a in self.partitions[(sid,'baseline')][1]['physical_atoms'][:-1]]+[metal],
                  'physical_charge':-2,'metal':metal,'physical_multiplicity':1 if metal=='La' else 6}
        return hashlib.sha256(digest(identity).encode()+self.np.asarray(x,dtype='<f8').tobytes()).hexdigest()

    def cache_checks(self):
        s=STATES[0];x=self.origins[s[:2]]; k=self.key(*s,x); y=x.copy();y[-1,0]+=FD['steps_A'][0]
        require(self.key(*s,x.copy())==k and self.key(*s,y)!=k,'Coordinate cache invalidation failed')
        altered=copy.deepcopy(self.identity);altered['GFN_profile']['sha256']='changed_for_key_test'
        require(digest(altered)!=self.identity_hash,'Parameter identity omitted')
        require(all(digest(dict(self.states[s],**{field:self.states[s][field]+1}))!=digest(self.states[s])
                    for field in ('QM_charge','MACE_multiplicity','GFN_uhf')),'State identity omitted')
        return {'repeat_equal':True,'coordinate_change':True,'parameter_change':True,'charge_spin_state_change':True,'new_calls':0}

    def replay_saved(self):
        np=self.np;old=read(self.cfg['GFN_evidence']);checks=[]
        require(old['status']=='complete' and old['accepted']==old['attempted']==4,'Mex GFN precedent unavailable')
        for r in old['rows']:
            t=dict(r['task'],field_charge_e=2.); d=verify(r['execution']).parent
            got=self.admit.admit(d,t,read(r['execution']),backend=self.backend,gradients=self.gradient,np=np)
            require(got['energy_hartree']==r['energy_hartree'] and got['charges_e']==r['charges_e']
                and got['qm_gradient']['hartree_per_bohr']==r['qm_gradient']['hartree_per_bohr'],'Generalized native GFN replay differs')
            if t['kind']=='embedded':require(got['pointcharge_gradient']['hartree_per_bohr']==r['pointcharge_gradient']['hartree_per_bohr'],'Generalized PC derivative replay differs')
            checks.append({'kind':'GFN_native_replay','task':t['task_id'],'result':r['result'],'pass':True})
        for r in self.alpb['rows']:
            sid,metal=r['source_id'],r['metal'];require(r['status']=='complete','Old ALPB unavailable')
            t=next(t for t in self.mcfg['tasks'] if t['task_id']==r['MACE_task_id'])
            require(t['charge']==-2 and t['multiplicity']==(1 if metal=='La' else 6),'Old solvent state differs')
            sy,x=self.admit.xyz(t['xyz']);p,part,_=self.partitions[(sid,'baseline')]
            q=read(self.tasks[(sid,'baseline',metal)]['global_solvent_charges'])['charges_e']
            require(sy==[a['element'] for a in part['physical_atoms'][:-1]]+[metal]
                and np.array_equal(x,self.origins[(sid,'baseline')]),'Old solvent geometry mismatch')
            raw=read(r['ALPB']);verify_tree(raw['artifacts'])
            oldq=[float(v) for v in verify(raw['artifacts']['charges']).read_text().split()]
            require(oldq==q and raw['rounding_bounds_hartree']==[0.,0.] and abs(raw['charge_e']+2)<=1e-5,'Old fixed solvent charge mismatch')
            grad=np.loadtxt(verify(raw['artifacts']['gradient.dat']))
            require(grad.shape==(256,3) and np.isfinite(grad).all(),'Old solvent derivative missing')
            key=self.solvent_key(sid,metal,x,q)
            self.solvent_cache[key]={'energy_hartree':raw['total_hartree'],'gradient_hartree_per_bohr':grad.tolist(),
                'components_hartree':raw['components_hartree'],'record':r['ALPB'],'source':'exact_saved_origin'}
            require(r['native_energies']['ALPB_Eh']==raw['total_hartree'],'Old solvent energy replay')
            checks.append({'kind':'ALPB_origin_replay','source':sid,'metal':metal,'result':r['ALPB'],'pass':True})
        require(len(self.solvent_cache)==4 and len(checks)==8,'Replay denominator differs');return checks

    def freeze_directions(self):
        np=self.np;p,part,_=self.partitions[('7CCO','baseline')];x=self.origins[('7CCO','baseline')]
        ids=[a['id'] for a in part['physical_atoms']];ix={s:i for i,s in enumerate(ids)};c=part['caps'][0]
        a,b=ix[c['retained_source_id']],ix[c['omitted_source_id']];v=x[b]-x[a];v/=np.linalg.norm(v)
        axis=np.eye(3)[int(np.argmin(np.abs(v)))];perp=np.cross(v,axis);perp/=np.linalg.norm(perp)
        field=sorted((a for a in part['field'] if a['id'].rsplit('/',1)[-1] in ('N','O') and abs(a['charge_e'])>1e-12),key=lambda a:a['id'])
        require(field,'External backbone charged atom missing');f=field[0];j=f['source_index'];rad=x[255]-x[j];rad/=np.linalg.norm(rad)
        vectors={}
        for name,i,d in [('boundary_anchor',b,perp),('external_charge',j,rad)]:
            y=np.zeros((256,3));y[i]=d;vectors[name]={'physical_id':ids[i],'index':i,'vector':y.tolist()}
        return {'source':'7CCO','region':'baseline','metal':'La','physical_ids':ids,'directions':vectors,
            'origin':self.tasks[('7CCO','baseline','La')]['physical_source'],'settings':FD,'frozen_before_calls':True}

    def mace(self,task,symbols,x):
        if self.native is None:
            tick=time.monotonic();self.native=self.warm.NativeMACE(self.mcfg);self.cost['MACE_model_loads']+=1
            self.cost['MACE_model_load_seconds']=time.monotonic()-tick
        n=len(symbols);v=self.native;mult=task['multiplicity'];q=task['charge']
        require(n==task['atom_count'] in (172,210) and symbols[0]==task['metal'] and symbols.count(task['metal'])==1
            and q==-2 and mult==(1 if task['metal']=='La' else 6),'Prepared MACE state differs')
        require(set(symbols)<={'H','C','N','O',task['metal']} and {v.Z[s] for s in symbols}<=set(v.model.atomic_numbers.tolist()),'MACE element unsupported')
        require((sum(v.Z[s] for s in symbols)-q-mult+1)%2==0,'Physical spin parity')
        atoms=v.Atoms(symbols,positions=x,pbc=False);atoms.info.update(charge=q,spin=mult)
        batch=v.calc._atoms_to_batch(atoms);state=v.input_batch(v.calc,atoms,q,mult,batch,0)
        self.alpb_helper.state(state,n,q,mult,0)
        tick=time.monotonic();self.cost['MACE_calls_attempted']+=1
        with v.torch.enable_grad():
            v.calls+=1
            value=v.model(batch.to_dict(),training=False,compute_force=True,compute_virials=False,compute_stress=False,
                compute_displacement=False,compute_hessian=False,compute_edge_forces=False,compute_atomic_stresses=False)
        e=float(value['energy'].detach().cpu().item());f=value['forces'].detach().cpu()
        require(math.isfinite(e) and tuple(f.shape)==(n,3) and bool(v.torch.isfinite(f).all()),'Native MACE nonfinite')
        require([p._version for p in v.model.parameters()]==v.parameter_versions,'MACE parameters changed')
        return {'energy_eV':e,'forces_eV_A':f.tolist(),'input_state_check':state,'force_convention':'negative vacuum energy derivative',
                'wall_seconds':time.monotonic()-tick,'model':self.cfg['model']}

    def native_gfn(self,task,mapping,d,kind):
        d.mkdir();symbols=[task['metal']]+[a['element'] for a in mapping['core'][1:]]
        (d/'core.xyz').write_text(xyz_text(symbols,[a['xyz_A'] for a in mapping['core']]))
        # The exact complete qualified parameter text is staged for every call.
        (d/'param_gfn2-xtb.txt').write_bytes(verify(self.runtime['parameters']).read_bytes())
        staged={'core.xyz':pin(d/'core.xyz'),'parameters':pin(d/'param_gfn2-xtb.txt')}
        require(staged['parameters']['sha256']==self.runtime['parameters']['sha256'],'Parameter staging changed bytes')
        command=[str(verify(self.runtime['executable'])),'core.xyz','--gfn','2','--chrg',str(task['charge']),
            '--uhf','0','--etemp','300','--acc','0.0001','--iterations','500','--norestart','--json','--pop','--grad']
        field={'rows':mapping['field'],'parameter_source':self.runtime['parameters'],
               'potential_reference':'finite hardness-damped potential zero at infinity; no PC self-energy'}
        fp=write(d/'FIELD.json',field)
        nf=len(mapping['field']) if kind=='embedded' else 0
        if kind=='embedded':
            lines=[str(nf)]+[' '.join(format(v,'.17g') for v in (a['charge_e'],*a['xyz_A'],self.gam[self.admit.Z[a['element']]])) for a in mapping['field']]
            (d/'external.pc').write_text('\n'.join(lines)+'\n');(d/'control.inp').write_text(CONTROL)
            command+=['--input','control.inp'];staged.update({n:pin(d/n) for n in ('external.pc','control.inp')})
        state=self.states[(task['source_id'],task['region'],task['metal'])]
        t=dict(task,kind=kind,xyz=staged['core.xyz'],multiplicity=1,uhf=0,physical_multiplicity=task['multiplicity'],
            valence_electrons=state['GFN_valence_electrons'],field=fp,pointcharge_count=nf,
            field_charge_e=state['field_charge_e'] if kind=='embedded' else 0.,command=command,staged_files=staged,
            runtime=self.runtime,full_corrected_profile=self.cfg['GFN_profile'],parameter_loading_source=self.cfg['parameter_loading_source'])
        t['prepared_DFT_metadata_not_executed']={k:t.pop(k) for k in ('frozen_f_ecp_core_electrons','frozen_f_explicit_electrons')}
        t['backend_spin_convention']='GFN2 UHF=alpha-minus-beta=0, effective multiplicity1; physical Dy multiplicity6 remains metadata'
        t['electronic_state']='original published spin-independent GFN2 valence model, f-in-core; no DFT ECP calculation'
        write(d/'TASK.json',t);e=dict(self.env,XTBPATH=str(d.resolve()))
        result={'status':'unavailable','task':pin(d/'TASK.json'),'attempted':False}
        try:
            verify_tree(self.runtime);verify_tree(staged);result['attempted']=True
            execution=self.backend.call(command,d,e);result['execution']=pin(d/'EXECUTION.json')
            verify_tree(staged)
            result.update(self.admit.admit(d,t,execution,backend=self.backend,gradients=self.gradient,np=self.np),status='accepted')
        except Exception as exc:result.update(failure=repr(exc),traceback=traceback.format_exc())
        result['artifacts']={p.name:pin(p) for p in d.iterdir() if p.is_file()}
        result['record']=write(d/'RESULT.json',result);return result

    def solvent(self,sid,metal,x,q,task,d):
        key=self.solvent_key(sid,metal,x,q)
        if key in self.solvent_cache:self.cost['ALPB_cache_hits']+=1;return self.solvent_cache[key]
        verify_tree(self.alpb['build']);self.cost['ALPB_calls_attempted']+=1
        sy=[a['element'] for a in self.partitions[(sid,'baseline')][1]['physical_atoms'][:-1]]+[metal]
        raw=self.alpb_helper.solvent(self.cfg['ALPB_binary'],{'symbols':sy,'x':x.tolist(),'q':q,'task':task},d,self.alpb_env)
        grad=self.np.loadtxt(verify(raw['artifacts']['gradient.dat']))
        require(grad.shape==(256,3) and self.np.isfinite(grad).all(),'Native solvent gradient missing')
        value={'energy_hartree':raw['total_hartree'],'gradient_hartree_per_bohr':grad.tolist(),
            'components_hartree':raw['components_hartree'],'record':pin(d/'RESULT.json'),'source':'new_native_call'}
        self.solvent_cache[key]=value;return value

    def evaluate(self,sid,region,metal,x):
        np=self.np;x=np.asarray(x,float);key=self.key(sid,region,metal,x)
        if key in self.cache:
            require(self.cache[key] is not None,'Exact previous state failed; no automatic retry')
            return self.cache[key]
        self.cost['composites_attempted']+=1;self.cache[key]=None;tick=time.monotonic()
        d=self.out/('point_%03d_%s_%s_%s'%(self.cost['composites_attempted'],sid,region,metal));d.mkdir()
        p,part,ledger=self.partitions[(sid,region)];mapping=self.mapping.update_geometry(part,x)
        task=dict(self.tasks[(sid,region,metal)]);q=read(task['global_solvent_charges'])['charges_e']
        xp=write(d/'POSITIONS.json',{'source_ids':[a['id'] for a in part['physical_atoms']],'positions_A':x.tolist()})
        source_task=copy.deepcopy(task);task['physical_source']=write(d/'PHYSICAL.json',mapping['physical_atoms'])
        sy=[metal]+[a['element'] for a in mapping['core'][1:]]
        (d/'core.xyz').write_text(xyz_text(sy,[a['xyz_A'] for a in mapping['core']]))
        task.update(xyz=pin(d/'core.xyz'),origin_task=source_task,task_id=source_task['task_id']+'_'+d.name)
        tp=write(d/'TASK.json',task)
        mp=write(d/'MAPPING.json',mapping);receipt={'status':'unavailable','cache_key':key,'identity_sha256':self.identity_hash,
            'state':self.states[(sid,region,metal)],'positions':xp,'mapping':mp,'source_task':source_task,'task':tp,'component_receipts':{}}
        try:
            mmcfg={'implementation':self.cfg['implementation'],'python':self.cfg['MM_python'],'version':self.cfg['openmm_version'],
                'system':p['systems'][metal]['system'],'positions':xp,'builder':self.pcfg['helpers']['metal_environment_component_checks.py']}
            cp=write(d/'MM_CONFIG.json',mmcfg);driver=d/'mm_driver';driver.mkdir()
            command=[str(verify(mmcfg['python'])),str(verify(self.cfg['implementation'])),'--stage','mm','--config',str(verify(cp)),'--output',str(d/'MM')]
            self.cost['MM_configurations_attempted']+=1
            ex=self.backend.call(command,driver,self.env);self.cost['observed_native_child_CPU_seconds']+=ex['cpu_user_seconds']+ex['cpu_system_seconds']
            receipt['component_receipts']['MM_execution']=pin(driver/'EXECUTION.json')
            require(ex['returncode']==0 and ex['owned_process_group_empty'],'MM subprocess failed')
            mm=read(pin(d/'MM/RESULT.json'));self.cost['MM_group_queries_attempted']+=mm.get('group_queries_attempted',0)
            require(mm['status']=='complete' and mm['system']==mmcfg['system'] and mm['positions']==xp,'MM result identity differs')
            receipt['component_receipts']['MM']=pin(d/'MM/RESULT.json')
            sy=[metal]+[a['element'] for a in mapping['core'][1:]];qx=np.asarray([a['xyz_A'] for a in mapping['core']])
            ml=self.mace(task,sy,qx);ml.update(task=tp,xyz=task['xyz'],physical_source=task['physical_source'])
            receipt['component_receipts']['MACE']=write(d/'MACE.json',ml)
            def call_gfn(kind):
                try:return self.native_gfn(task,mapping,d/kind,kind)
                except Exception as exc:
                    # Preparation/launch failure still has an explicit denominator
                    # entry, while the other already-started native call is retained.
                    gd=d/kind;gd.mkdir(exist_ok=True)
                    failure={'status':'unavailable','attempted':(gd/'START.json').exists(),
                             'failure':repr(exc),'traceback':traceback.format_exc(),'kind':kind}
                    if (gd/'EXECUTION.json').is_file():failure['execution']=pin(gd/'EXECUTION.json')
                    failure['record']=write(gd/'STAGING_FAILURE.json',failure);return failure
            with ThreadPoolExecutor(max_workers=2) as pool:
                gf=list(pool.map(call_gfn,('embedded','isolated')))
            g=dict(zip(('embedded','isolated'),gf))
            for kind,r in g.items():
                self.cost['GFN_calls_attempted']+=int(r['attempted']);receipt['component_receipts']['GFN_'+kind]=r['record']
                if 'execution' in r:
                    ex=read(r['execution']);self.cost['observed_native_child_CPU_seconds']+=ex['cpu_user_seconds']+ex['cpu_system_seconds']
            require(all(r['status']=='accepted' for r in gf),'Native GFN component unavailable')
            sol=self.solvent(sid,metal,x,q,task,d/'ALPB');receipt['component_receipts']['ALPB']=sol['record']
            receipt['ALPB_reuse']=sol['source']
            nf=len(mapping['field']);zero=np.zeros((nf,3))
            gg={k:self.mapping.scatter(mapping,np.asarray(r['qm_gradient']['hartree_per_bohr']),
                np.asarray(r['pointcharge_gradient']['hartree_per_bohr']) if k=='embedded' else zero)*HA/BOHR for k,r in g.items()}
            gradients={'MACE_vacuum':self.mapping.scatter(mapping,-np.asarray(ml['forces_eV_A'])*EV,zero),
                'GFN2_field_difference':gg['embedded']-gg['isolated'],'MM_ledger':np.asarray(mm['gradient_kcal_mol_A']),
                'fixed_FF_charge_ALPB':np.asarray(sol['gradient_hartree_per_bohr'])*HA/BOHR}
            native={'MACE_eV':ml['energy_eV'],'GFN_embedded_Eh':g['embedded']['energy_hartree'],
                'GFN_isolated_Eh':g['isolated']['energy_hartree'],'MM_groups_kcal_mol':mm['component_energies_kcal_mol'],
                'ALPB_Eh':sol['energy_hartree']}
            components={'MACE_vacuum':native['MACE_eV']*EV,'GFN2_field_difference':(native['GFN_embedded_Eh']-native['GFN_isolated_Eh'])*HA,
                'MM_ledger':mm['energy_kcal_mol'],'fixed_FF_charge_ALPB':native['ALPB_Eh']*HA}
            gradient=sum(gradients[k] for k in COMPONENTS);require(gradient.shape==(256,3) and np.isfinite(gradient).all(),'Composite derivative invalid')
            point={'positions_A':x.tolist(),'gradient_kcal_mol_A':gradient.tolist(),'total_energy_kcal_mol':math.fsum(components.values()),
                'component_energies_kcal_mol':components,'component_gradients_kcal_mol_A':{k:gradients[k].tolist() for k in COMPONENTS},
                'native_energies':native,'source_id':sid,'region':region,'metal':metal,'state':receipt['state']}
            receipt.update(status='complete',point=point)
        except Exception as exc:receipt.update(failure=repr(exc),traceback=traceback.format_exc())
        receipt['elapsed_seconds']=time.monotonic()-tick;rp=write(d/'POINT.json',receipt);self.requests.append(rp)
        require(receipt['status']=='complete','Composite unavailable: '+receipt.get('failure','unknown'))
        point['provenance']={'point':rp,'identity_sha256':self.identity_hash,'cache_key':key,'component_receipts':receipt['component_receipts']}
        self.cache[key]=point;self.points.append(rp);return point

    def restore(self,record):
        r=read(record);require(r['status']=='complete' and r['identity_sha256']==self.identity_hash,'Restored model identity differs')
        verify_tree(r['component_receipts']);verify(r['mapping']);x=read(r['positions'])['positions_A'];p=copy.deepcopy(r['point'])
        require(p['positions_A']==x and r['cache_key']==self.key(p['source_id'],p['region'],p['metal'],x),'Restored coordinates differ')
        require(self.np.array_equal(p['gradient_kcal_mol_A'],sum(self.np.asarray(p['component_gradients_kcal_mol_A'][k]) for k in COMPONENTS))
            and p['total_energy_kcal_mol']==math.fsum(p['component_energies_kcal_mol'].values()),'Restored component closure differs')
        require(r['state']==p['state']==self.states[(p['source_id'],p['region'],p['metal'])],'Restored physical/backend state differs')
        receipts=r['component_receipts'];ml=read(receipts['MACE']);mm=read(receipts['MM'])
        emb=read(receipts['GFN_embedded']);iso=read(receipts['GFN_isolated']);sol=read(receipts['ALPB'])
        native={'MACE_eV':ml['energy_eV'],'GFN_embedded_Eh':emb['energy_hartree'],
                'GFN_isolated_Eh':iso['energy_hartree'],'MM_groups_kcal_mol':mm['component_energies_kcal_mol'],
                'ALPB_Eh':sol['total_hartree']}
        require(native==p['native_energies'] and emb['status']==iso['status']=='accepted' and mm['status']=='complete',
                'Restored native component evidence differs')
        mapping=read(r['mapping']);np=self.np;zero=np.zeros((len(mapping['field']),3))
        mapped=self.mapping.update_geometry(self.partitions[(p['source_id'],p['region'])][1],x)
        map_check=self.mapping_replay.compare(mapping,mapped)
        self.mapping_replays.append({'kind':'restored_point','point':record,'archived_mapping':r['mapping'],'check':map_check})
        # Raw derivative closure below deliberately continues using archived mapping.
        eg=self.mapping.scatter(mapping,np.asarray(emb['qm_gradient']['hartree_per_bohr']),np.asarray(emb['pointcharge_gradient']['hartree_per_bohr']))*HA/BOHR
        ig=self.mapping.scatter(mapping,np.asarray(iso['qm_gradient']['hartree_per_bohr']),zero)*HA/BOHR
        raw_gradients={'MACE_vacuum':self.mapping.scatter(mapping,-np.asarray(ml['forces_eV_A'])*EV,zero),
            'GFN2_field_difference':eg-ig,'MM_ledger':np.asarray(mm['gradient_kcal_mol_A']),
            'fixed_FF_charge_ALPB':np.loadtxt(verify(sol['artifacts']['gradient.dat']))*HA/BOHR}
        require(all(np.array_equal(raw_gradients[k],p['component_gradients_kcal_mol_A'][k]) for k in COMPONENTS),
                'Restored full component derivative evidence differs')
        p['provenance']={'point':record,'identity_sha256':self.identity_hash,'cache_key':r['cache_key'],'component_receipts':r['component_receipts']}
        self.cache[r['cache_key']]=p;return p


def difference(a,b):
    a,b=a['native_energies'],b['native_energies']
    parts={'MACE_vacuum':(a['MACE_eV']-b['MACE_eV'])*EV,
        'GFN2_field_difference':math.fsum((a['GFN_embedded_Eh']-b['GFN_embedded_Eh'],-a['GFN_isolated_Eh']+b['GFN_isolated_Eh']))*HA,
        'MM_ledger':math.fsum(a['MM_groups_kcal_mol'][k]-b['MM_groups_kcal_mol'][k] for k in a['MM_groups_kcal_mol']),
        'fixed_FF_charge_ALPB':(a['ALPB_Eh']-b['ALPB_Eh'])*HA}
    return dict(parts,composite_kcal_mol=math.fsum(parts.values()))


def scout(ev,result):
    np=ev.np;x=ev.origins[('7CCO','baseline')]
    slots=[{'cell':'center','status':'unavailable_not_attempted'}]
    slots += [{'cell':name+':'+format(h,'.17g')+':'+sign,'direction':name,'h_A':h,'sign':sign,
               'status':'unavailable_not_attempted'} for name in FD['directions'] for h in FD['steps_A'] for sign in ('plus','minus')]
    result['declared_cells']=slots
    def evaluate(slot,y):
        slot['status']='unavailable_attempted'
        try:
            value=ev.evaluate('7CCO','baseline','La',y)
            slot.update(status='complete',point=value['provenance']['point']);return value
        except Exception as exc:slot['failure']=repr(exc);raise
    center=evaluate(slots[0],x)
    result['center']=center['provenance']['point'];require(ev.evaluate('7CCO','baseline','La',x)==center,'Cached repeat differs')
    result['checks']=[]
    for name,data in ev.directions['directions'].items():
        v=np.asarray(data['vector']); checks=[]
        for h in FD['steps_A']:
            plus_slot=next(r for r in slots if r.get('direction')==name and r.get('h_A')==h and r.get('sign')=='plus')
            minus_slot=next(r for r in slots if r.get('direction')==name and r.get('h_A')==h and r.get('sign')=='minus')
            plus=evaluate(plus_slot,x+h*v);minus=evaluate(minus_slot,x-h*v)
            delta=difference(plus,minus);rows=[]
            for k in (*COMPONENTS,'composite_kcal_mol'):
                g=center['gradient_kcal_mol_A'] if k=='composite_kcal_mol' else center['component_gradients_kcal_mol_A'][k]
                analytic=float(np.sum(np.asarray(g)*v));finite=delta[k]/(2*h);residual=finite-analytic
                tol=FD['absolute_kcal_mol_A']+FD['relative']*max(abs(finite),abs(analytic))
                rows.append({'component':k,'analytic_kcal_mol_A':analytic,'FD_kcal_mol_A':finite,
                    'residual_kcal_mol_A':residual,'tolerance_kcal_mol_A':tol,'pass':abs(residual)<=tol})
            checks.append({'h_A':h,'plus':plus['provenance']['point'],'minus':minus['provenance']['point'],'components':rows,'pass':all(r['pass'] for r in rows)})
        result['checks'].append({'direction':name,'physical_id':data['physical_id'],'steps':checks,
            'step_halving_change_kcal_mol_A':{a['component']:b['FD_kcal_mol_A']-a['FD_kcal_mol_A'] for a,b in zip(checks[0]['components'],checks[1]['components'])},
            'pass':all(r['pass'] for r in checks)})
    require(ev.cost['composites_attempted']==9,'Scout finite denominator differs')
    result['qualification_passed']=all(r['pass'] for r in result['checks'])
    result['status']='complete' if result['qualification_passed'] else 'failed_derivative_gate'


def origins(ev,cfg,result):
    prior=read(cfg['scout_result']);require(prior['stage']=='scout' and prior['status']=='complete' and prior['qualification_passed']
        and prior['identity_sha256']==ev.identity_hash and prior['actual_cost']['composites_attempted']==9,'Actual scout admission required')
    ev.restore(prior['center']);result['rows']=[]
    for sid,region,metal in STATES:
        row={'source_id':sid,'region':region,'metal':metal,'status':'unavailable'};result['rows'].append(row)
        try:
            p=ev.evaluate(sid,region,metal,ev.origins[(sid,region)])
            row.update(status='complete',point=p['provenance']['point'],native_energies=p['native_energies'],
                component_energies_kcal_mol=p['component_energies_kcal_mol'],total_energy_kcal_mol=p['total_energy_kcal_mol'])
        except Exception as exc:row['failure']=repr(exc)
    require(ev.cost['composites_attempted']<=7,'Origin-stage finite denominator exceeded')
    result['status']='complete' if all(r['status']=='complete' for r in result['rows']) else 'unavailable'
    result['source_response_by_region']=[]
    cells={(r['source_id'],r['region'],r['metal']):r for r in result['rows']}
    for region in ('baseline','expanded'):
        r={'region':region,'status':'unavailable','definition':'[Dy(7CCN)-Dy(7CCO)]-[La(7CCN)-La(7CCO)]'}
        result['source_response_by_region'].append(r)
        if all(cells[(s,region,m)]['status']=='complete' for s in ('7CCO','7CCN') for m in ('La','Dy')):
            works={m:difference(cells[('7CCN',region,m)],cells[('7CCO',region,m)]) for m in ('La','Dy')}
            differential={k:works['Dy'][k]-works['La'][k] for k in (*COMPONENTS,'composite_kcal_mol')}
            r.update(status='available',source_response_kcal_mol=works,metal_dependent_source_response_kcal_mol=differential)
    result['interpretation']='Same-origin Hamiltonian comparison only; no absolute comparison across region partitions; no preference/affinity or optimization qualification.'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=('bind','preflight','scout','origins','mm'),required=True)
    p.add_argument('--config',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.stage=='bind':bind(a.output);return 0
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    cfg=json.loads(a.config.read_text());out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic()
    initial_usage=resource.getrusage(resource.RUSAGE_SELF)
    result={'status':'execution_failure','stage':a.stage,'protocol':PROTOCOL,'implementation':pin(__file__),'config':pin(a.config),
        'slurm_job_id':os.environ['SLURM_JOB_ID'],'declared_origins':8,'new_DFT_calls':0,'optimizer_steps':0};ev=None
    try:
        if a.stage=='mm':result.update(mm_worker(cfg,out))
        else:
            ev=HybridEvaluator(cfg,out/'evaluator');result.update(identity_sha256=ev.identity_hash,preflight=pin(ev.out/'PREFLIGHT.json'))
            if a.stage=='preflight':result['status']='complete'
            elif a.stage=='scout':scout(ev,result)
            else:origins(ev,cfg,result)
    except Exception as exc:result.update(failure=repr(exc),traceback=traceback.format_exc())
    if ev is not None:
        result.update(actual_cost=ev.cost,requests=ev.requests,points=ev.points,mapping_replay_checks=ev.mapping_replays)
        if ev.native is not None:
            result.update(native_MACE_calls=ev.native.calls,native_MACE_device=ev.native.device,native_MACE_versions=ev.native.versions,
                peak_CUDA_allocated_bytes=ev.native.torch.cuda.max_memory_allocated())
    result['elapsed_seconds']=time.monotonic()-tick
    usage=resource.getrusage(resource.RUSAGE_SELF)
    result['peak_host_RSS_KiB']=usage.ru_maxrss
    result['observed_parent_CPU_seconds']=usage.ru_utime+usage.ru_stime-initial_usage.ru_utime-initial_usage.ru_stime
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_PER_TASK','1'))
    write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}));return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
