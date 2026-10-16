#!/usr/bin/env python3
"""Finite same-basis Mex continuation. All scientific work requires Slurm.
CLI --stage preflight|scout|mm --config CONFIG --output FRESH.
Config: implementation,plan,pilot_result,search_pattern. Root alone submits.
Remaining search and paired pool use separately frozen adapters after the scout.
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
import socket
import sys
import time
import traceback

ROOT=Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
SETTINGS={'source_extent_A':.20,'target_extent_A':.20,'total_torsion_radian':.20,
 'trial_extent_A':.02,'trial_torsion_degrees':.5,'maximum_accepted_steps':12,'maximum_new_own_evaluations':24,
 'history':8,'armijo':1e-4,'trial_multipliers':[1.,.5,.25,.125],'metric_relative_cutoff':1e-10,
 'active_max_atom_kcal_mol_A':.1,'active_RMS_cartesian_kcal_mol_A':.03,
 'bond_tolerance_A':1e-10,'angle_cosine_tolerance':1e-10,'protected_distance_tolerance_A':1e-9,
 'fixed_atom_tolerance_A':1e-12,'displacement_roundoff_A':1e-12}
PILOT_SHA='2c5c50010f502ee0a09262cefa8934b5237716dc7290fe6a8cbc187a103eef32'
PATTERN_SHA='ee1db81bdf3db90c303f824d7a55b52b00c2bb64f2422903aed60c25931d8ca1'


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path']);return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def write(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


class Mechanics:
    """Pilot's complete physical guards, with only the declared search envelopes."""
    def __init__(self,world,modes,geometry):
        from scipy.spatial import cKDTree
        self.tree=cKDTree;self.np=np=world[6];c=world[1]['La'];export=world[10]
        self.c=c;self.x=c['x'];self.ids=c['ids'];self.ix={s:i for i,s in enumerate(self.ids)};self.target=self.ix[c['target_id']]
        self.modes=modes;self.axes=modes['axes'];x=self.x;ix=self.ix
        require(modes['physical_source_ids']==self.ids and modes['target_source_id']==c['target_id'],'Mode physical source differs')
        require(modes['mode_order']==['Leu83_N_CA','Asp84_CA_CB','Asp86_CA_CB','Glu95_CA_CB','target_x','target_y','target_z'],'Mode order differs')
        native=read(export['atoms']);metadata={a['id']:a for a in native};bonds=read(export['bonds'])
        require(modes['bonds']==export['bonds'],'Mode graph differs');protein=[a['id'] for a in native if a['kind']=='protein_source'];pset=set(protein)
        graph={sid:set() for sid in self.ids}
        for a,b in bonds:graph[a].add(b);graph[b].add(a)
        for m in self.axes:
            a,b=m['axis'];require(b in graph[a],'Axis not bonded');seen={b};pending=[b]
            while pending:
                now=pending.pop()
                for nxt in graph[now]-seen:
                    if {now,nxt}=={a,b}:continue
                    seen.add(nxt);pending.append(nxt)
            require(a not in seen and seen<=pset and sorted(seen)==m['moving_source_ids']
                and [ix[s] for s in sorted(seen)]==m['moving_indices'],'Frozen complete hinge membership differs')
        self.moving={self.target}|{i for m in self.axes for i in m['moving_indices']};self.fixed=sorted(set(range(2049))-self.moving)
        require(sorted(self.ids[i] for i in self.moving)==sorted(modes['moving_source_ids']),'Moving union differs')
        require(all(self.ids[i] in pset or i==self.target for i in self.moving),'Water/spectator motion')
        self.groups=modes['protected_groups']
        for m in self.axes:
            moving=set(m['moving_source_ids']);axis=set(m['axis'])
            for group in self.groups:
                require([ix[s] for s in group['source_ids']]==group['indices'],'Protected group mapping differs')
                rest=set(group['source_ids'])-axis;require(not(rest&moving and rest-moving),'Protected group split')
        pidx={s:i for i,s in enumerate(protein)};self.porder=[ix[s] for s in protein];self.other=[i for i,s in enumerate(self.ids) if s not in pidx]
        adapted=[{'source_id':s,'residue_id':s.rsplit('/',1)[0],'name':s.rsplit('/',1)[1],
                  'element':metadata[s]['element'],'resname':metadata[s]['resname']} for s in protein]
        pbonds=[(pidx[a],pidx[b]) for a,b in bonds if a in pidx and b in pidx]
        self.guard=geometry.Geometry(adapted,pbonds,x[self.porder],x[self.porder],
            [pidx[self.ids[i]] for i in sorted(self.moving) if self.ids[i] in pidx],
            [dict(c['physical'][i],xyz_A=x[i].tolist()) for i in self.other])
        self.bi=np.asarray([(ix[a],ix[b]) for a,b in bonds]);self.angles=np.asarray([(ix[a],ix[b],ix[d]) for b,ns in graph.items() for a in sorted(ns) for d in sorted(ns) if a<d])
        self.excluded=set()
        for s,ns in graph.items():
            near=set(ns)
            for n in ns:near.update(graph[n])
            self.excluded.update(tuple(sorted((ix[s],ix[n]))) for n in near if n!=s)
        self.lengths=np.linalg.norm(x[self.bi[:,0]]-x[self.bi[:,1]],axis=1);self.cos=self.cosines(x)
        self.j0=self.jacobian(x);matrix=self.j0.reshape(-1,7);self.norm=np.linalg.norm(matrix,axis=0)
        require(np.max(np.abs(self.j0-np.asarray(modes['origin_J'])))<=1e-12 and np.allclose(self.norm,modes['column_norms'],rtol=0,atol=1e-12),'Saved source metric differs')
        U=matrix/self.norm;ev,V=np.linalg.eigh(U.T@U);require(int((ev>1e-10*ev[-1]).sum())==7,'Original metric rank differs')
        self.T=((V/np.sqrt(ev))@V.T)/self.norm[:,None]
        require(np.max(np.abs((matrix@self.T).T@(matrix@self.T)-np.eye(7)))<=1e-10,'Whitening failed')
        self.eigenvalues=ev.tolist();self.geometry_module=geometry
        require(self.check(np.zeros(7),x,np.zeros(7),x)['passed'],'Original source guard failed')

    def cosines(self,y):
        np=self.np;a,b,d=self.angles.T;u,v=y[a]-y[b],y[d]-y[b]
        return np.einsum('ij,ij->i',u,v)/np.linalg.norm(u,axis=1)/np.linalg.norm(v,axis=1)

    def coordinates(self,q):
        np=self.np;y=self.x.copy()
        for k,m in enumerate(self.axes):
            a,b=[self.ix[s] for s in m['axis']];u=y[b]-y[a];u/=np.linalg.norm(u);mi=m['moving_indices'];pivot=y[b].copy();v=y[mi]-pivot
            y[mi]=pivot+v*np.cos(q[k])+np.cross(u,v)*np.sin(q[k])+np.outer(v@u,u)*(1-np.cos(q[k]))
        y[self.target]+=q[4:];return y

    def jacobian(self,y):
        np=self.np;j=np.zeros((2049,3,7))
        for k,m in enumerate(self.axes):
            a,b=[self.ix[s] for s in m['axis']];u=y[b]-y[a];u/=np.linalg.norm(u);mi=m['moving_indices']
            j[mi,:,k]=np.cross(u,y[mi]-y[b])
        j[self.target,:,4:]=np.eye(3);return j

    def check(self,q,y,previous_q,previous):
        np=self.np;x=self.x;require(y.shape==(2049,3) and q.shape==(7,) and np.isfinite(y).all() and np.isfinite(q).all(),'Nonfinite geometry')
        self.guard.other=y[self.other].copy();inherited=self.guard.check(y[self.porder]);contacts=[]
        for a,b in sorted(self.tree(y).query_pairs(1.)):
            if (a not in self.moving and b not in self.moving) or (a,b) in self.excluded:continue
            cutoff=.55 if 'H' in (self.c['physical'][a]['element'],self.c['physical'][b]['element']) else 1.
            d=float(np.linalg.norm(y[a]-y[b]));d0=float(np.linalg.norm(x[a]-x[b]))
            if d<min(cutoff,d0)-1e-7:contacts.append({'source_ids':[self.ids[a],self.ids[b]],'source_A':d0,'candidate_A':d})
        bond=float(np.max(np.abs(np.linalg.norm(y[self.bi[:,0]]-y[self.bi[:,1]],axis=1)-self.lengths)))
        angle=float(np.max(np.abs(self.cosines(y)-self.cos)));protected=0.
        for group in self.groups:
            ii=group['indices'];protected=max(protected,float(np.max(np.abs(np.linalg.norm(y[ii][:,None]-y[ii][None,:],axis=2)-np.linalg.norm(x[ii][:,None]-x[ii][None,:],axis=2)))))
        extent=float(np.linalg.norm(y-x,axis=1).max());target=float(np.linalg.norm(y[self.target]-x[self.target]));torsion=float(np.abs(q[:4]).max())
        trial=float(np.linalg.norm(y-previous,axis=1).max());dq=float(np.rad2deg(np.abs(q[:4]-previous_q[:4])).max())
        fixed=float(np.max(np.abs(y[self.fixed]-x[self.fixed]),initial=0));bounds={'source_extent':extent<=.2+1e-12,
            'target_extent':target<=.2+1e-12,'total_torsion':torsion<=.2+1e-12,'trial_extent':trial<=.02+1e-12,'trial_torsion':dq<=.5+1e-12}
        return {'passed':bool(inherited['pass'] and not contacts and bond<=1e-10 and angle<=1e-10 and protected<=1e-9 and fixed<=1e-12 and all(bounds.values())),
            'inherited':inherited,'all_physical_worsened_severe_contacts':contacts,'bond_max_error_A':bond,'angle_cosine_max_error':angle,
            'protected_distance_max_error_A':protected,'fixed_max_error_A':fixed,'max_source_displacement_A':extent,
            'target_source_displacement_A':target,'max_total_torsion_radian':torsion,'max_trial_displacement_A':trial,
            'max_trial_torsion_degrees':dq,'bounds':bounds,'immutable_source_reference':True}

    def diagnostic(self,y,g):
        np=self.np;matrix=self.jacobian(y).reshape(-1,7);norm=np.linalg.norm(matrix,axis=0);U=matrix/norm
        ev,V=np.linalg.eigh(U.T@U);require(int((ev>1e-10*ev[-1]).sum())==7,'Current metric rank changed')
        inv=(V/ev)@V.T;projected=(U@inv@(U.T@g.reshape(-1))).reshape(2049,3);residual=g-projected
        maximum=float(np.linalg.norm(projected,axis=1).max());rms=float(np.sqrt(np.mean(projected**2)))
        return {'active_span_max_atom_kcal_mol_A':maximum,'active_span_RMS_cartesian_kcal_mol_A':rms,
            'outside_active_span_L2_kcal_mol_A':float(np.linalg.norm(residual)),
            'current_axis_energy_loads_in_original_column_metric_kcal_mol_A':((matrix.T@g.reshape(-1))/self.norm).tolist(),
            'normalized_Gram_eigenvalues':ev.tolist(),'rank':7,'low_active_span_force':maximum<=.1 and rms<=.03,
            'stationary_minimum_claimed':False}


def inputs(cfg,out):
    require(cfg['pilot_result']['sha256']==PILOT_SHA and cfg['search_pattern']['sha256']==PATTERN_SHA,'Frozen parent differs')
    pilot=read(cfg['pilot_result']);verify(cfg['search_pattern'])
    require(pilot['status']=='passed' and pilot['slurm_job_id']=='12717' and pilot['comparison']['finite_pilot_passed']
        and len(pilot['rows'])==4 and all(r['status']=='complete' for r in pilot['rows']),'Pilot gate failed')
    pcfg=read(pilot['config']);response=module('mex_continuation_response',verify(pilot['implementation']))
    fd,fcfg,world,geometry=response.inputs(pcfg,out);np=world[6];modes=read(pilot['modes']);props=read(pilot['proposals'])
    require(pilot['modes']['sha256']=='b2f227d576fc9b1fb2f855408f5591d228230449f8cc67e0a4e70cceed1a2d1d'
        and pilot['proposals']['sha256']=='b7dc055a2853c0170c15339fb3e6f0c0b0e1bab44c7a8c4c65ac59194dc76c84','Basis/start pins differ')
    mech=Mechanics(world,modes,geometry);points={};checks=[]
    for row in pilot['rows']:
        require(read(row['cell'])=={k:v for k,v in row.items() if k!='cell'},'Pilot CELL differs')
        trial=read(row['trial']);mp=read(trial['mapping']);p=props[row['proposer']];physical=read(trial['physical'])
        y=np.asarray([a['xyz_A'] for a in physical]);q=np.asarray(p['q']);u=np.linalg.solve(mech.T,q)
        require(p['status']=='prepared' and [a['id'] for a in physical]==mech.ids and mp['physical_atoms']==physical,'Pilot physical map differs')
        require(np.array_equal(y,np.asarray(read(p['physical'])['xyz_A'])) and np.max(np.abs(mech.coordinates(q)-y))<=1e-12
            and np.max(np.abs(mech.coordinates(mech.T@u)-y))<=1e-12,'Pilot start replay differs')
        grad=read(row['assembled_gradient']);g=np.asarray(grad['gradient_kcal_mol_A'])
        require(grad['source_ids']==mech.ids and grad['units']=='kcal/mol/angstrom' and grad['quantity']=='gradient=dE/dR'
            and grad['physical_mapping']==trial['mapping'] and g.shape==(2049,3) and np.isfinite(g).all(),'Pilot gradient state/map differs')
        ml=read(row['MACE']);sol=read(row['ALPB']);mm=read(row['MM_result']);mmrow=next(r for r in mm['rows'] if r['task_id']==row['task_id'])
        gf={k:read(row['GFN2_'+k]) for k in ('embedded','isolated')};v=row['native_energies']
        require(ml['energy_eV']==v['MACE_eV'] and sol['total_hartree']==v['ALPB_Eh']
            and mmrow['components_kcal_mol']==v['MM_groups_kcal_mol'],'Pilot component energy linkage differs')
        for kind in gf:
            require(gf[kind]['status']=='accepted' and gf[kind]['energy_hartree']==v['GFN2_'+kind+'_Eh']
                and gf[kind]['task']['mapping']==trial['mapping'],'Pilot GFN state/map differs')
        gc=mech.check(q,y,q,y);require(gc['passed'],'Pilot start violates continuation source guard')
        point={'status':'complete','metal':row['metal'],'q':q.tolist(),'physical':trial['physical'],'trial':row['trial'],
            'native_energies':v,'assembled_gradient':row['assembled_gradient'],'cell':row['cell'],'origin':'admitted_pilot_12717'}
        points[(row['proposer'],row['metal'])]=point;checks.append({'proposer':row['proposer'],'metal':row['metal'],'passed':True,'geometry':gc})
    require(set(points)=={(p,m) for p in ('La','Dy') for m in ('La','Dy')},'Pilot paired denominator differs')
    write(out/'PREFLIGHT.json',{'status':'passed','pilot_result':cfg['pilot_result'],'modes':pilot['modes'],
        'starts':checks,'coordinate_scaling_T':mech.T.tolist(),'source_metric_eigenvalues':mech.eigenvalues,'new_molecular_calls':0})
    return response,fcfg,world,mech,pilot,points


# Same qualified Reference MM primitive, one evaluated state per request.
def mm_stage(cfg,out,result):
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'],'MM interpreter binary differs')
    verify(cfg['common_helper']);sys.path.insert(0,str(ROOT/'scripts'))
    builder=module('mex_fd_classical',verify(cfg['classical_helper']))
    import numpy as np
    import openmm as mm
    from openmm import unit
    require(mm.__version__==cfg['openmm_version'],'Original OpenMM version differs')
    platform=mm.Platform.getPlatformByName('Reference');prepared={};trial_data={}
    require(len(cfg['trials'])==1,'MM continuation cell denominator differs')
    # All imports, systems and coordinates are admitted before the first energy.
    for metal in ('La','Dy'):
        ledger=read(cfg['ledgers'][metal]);require(ledger['reused_builder']==cfg['classical_helper'],'MM builder lineage')
        artifacts={k:read(p) for k,p in ledger['artifacts'].items()}
        system,counts=builder.build_system(ledger,artifacts,metal)
        require(system.getNumParticles()==2049 and system.getNumConstraints()==0,'MM physical inventory/constraints')
        ids=[cfg['target_source_id'] if a['id']=='target_EF3' else a['id'] for a in artifacts['particles.json']]
        require(len(ids)==len(set(ids))==2049,'MM particle identity differs')
        prepared[metal]=(system,ids,counts)
    for p in cfg['trials']:
        t=read(p);physical=read(t['physical']);by={a['id']:a for a in physical};system,ids,counts=prepared[t['metal']]
        require(len(by)==2049 and set(by)==set(ids),'MM trial/source ID join differs')
        x=np.asarray([by[s]['xyz_A'] for s in ids]);require(x.shape==(2049,3) and np.isfinite(x).all(),'MM trial coordinate shape')
        trial_data[t['task_id']]=(t,x)
    result.update(runtime={'python':sys.executable,'numpy':np.__version__,'openmm':mm.__version__,
                          'openmm_full_version':mm.version.full_version,'platform':platform.getName()},rows=[],
                  configurations_attempted=0,group_queries_attempted=0)
    contexts={}
    for metal,(system,ids,counts) in prepared.items():
        integrator=mm.VerletIntegrator(.001);contexts[metal]=(mm.Context(system,integrator,platform),integrator)
    for tid,(t,x) in trial_data.items():
        row={'task_id':tid,'metal':t['metal'],'status':'unavailable','trial':next(p for p in cfg['trials'] if read(p)['task_id']==tid)}
        try:
            context=contexts[t['metal']][0];context.setPositions(x*.1*unit.nanometer);energies={};gradients={}
            result['configurations_attempted']+=1
            for i,name in enumerate(builder.GROUPS):
                result['group_queries_attempted']+=1
                state=context.getState(getEnergy=True,getForces=True,groups=1<<i)
                energies[name]=state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
                gradients[name]=(-np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom))).tolist()
            require(all(math.isfinite(v) for v in energies.values()) and all(np.isfinite(g).all() for g in gradients.values()),'Nonfinite MM energy/gradient')
            row.update(status='complete',components_kcal_mol=energies,total_kcal_mol=math.fsum(energies.values()),
                component_gradients_kcal_mol_A=gradients,source_ids=prepared[t['metal']][1])
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
        result['rows'].append(row)
    result['status']='complete' if all(r['status']=='complete' for r in result['rows']) else 'unavailable'


def relative(response,old,value,reference):
    parts=response.work(value,reference,old)
    return dict(parts,composite_kcal_mol=math.fsum(parts.values()))


class Evaluator:
    def __init__(self,cfg,out,result,loaded):
        self.cfg=cfg;self.out=out;self.result=result
        self.response,self.fcfg,self.world,self.mech,self.pilot,self.pilot_points=loaded
        self.old,self.cells,self.solcfg,self.mcfg,self.helper,self.mapper,self.np,self.local,self.fw,self.source,self.export,self.origins,self.previous=self.world
        self.backend,self.env=self.fw[6],self.fw[9];self.native=None;self.cache={};self.points=[];self.requests=[];self.new=0
        self.binary=verify(self.previous['build']['binary'])
        self.alpb_env=dict(self.env,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=self.solcfg['library_directory'])
        saved=out/'cached';saved.mkdir()
        for (proposer,metal),point in self.pilot_points.items():
            path=saved/(proposer+'_'+metal+'.json');write(path,point);self.admit_saved(pin(path))
        for metal in ('La','Dy'):
            write(saved/(metal+'_origin_physical.json'),self.cells[metal]['physical'])
            point={'status':'complete','metal':metal,'q':[0.]*7,'physical':pin(saved/(metal+'_origin_physical.json')),
                'native_energies':self.response.origin_native(metal,self.world),'assembled_gradient':self.origins[metal]['pin'],
                'origin':'admitted_static_12685','mapping':self.cells[metal]['mapping_pin']}
            path=saved/(metal+'_origin.json');write(path,point);self.admit_saved(pin(path))

    def key(self,metal,y):
        return metal+':'+hashlib.sha256(self.np.asarray(y,dtype='<f8').tobytes(order='C')).hexdigest()

    def coordinates(self,point):return self.np.asarray([a['xyz_A'] for a in read(point['physical'])])

    def admit_saved(self,record):
        p=read(record);np=self.np;require(p['status']=='complete' and p['metal'] in ('La','Dy'),'Unavailable cached point')
        ph=read(p['physical']);y=np.asarray([a['xyz_A'] for a in ph]);q=np.asarray(p['q']);g=read(p['assembled_gradient'])
        require([a['id'] for a in ph]==self.mech.ids and y.shape==(2049,3) and np.isfinite(y).all()
            and np.max(np.abs(self.mech.coordinates(q)-y))<=1e-12,'Cached physical source differs')
        require(g['source_ids']==self.mech.ids and np.asarray(g['gradient_kcal_mol_A']).shape==(2049,3)
            and np.isfinite(g['gradient_kcal_mol_A']).all(),'Cached gradient differs')
        require(self.mech.check(q,y,q,y)['passed'],'Cached source-relative guard failed')
        if 'cell' in p:
            cell=read(p['cell']);require(cell['status']=='complete' and cell['native_energies']==p['native_energies']
                and cell['assembled_gradient']==p['assembled_gradient'] and cell['trial']==p['trial'],'Cached CELL differs')
            trial=read(p['trial']);require(trial['metal']==p['metal'] and trial['physical']==p['physical'],'Cached state differs')
            for k in ('MACE','GFN2_embedded','GFN2_isolated','ALPB','MM_result'):verify(cell[k])
        key=self.key(p['metal'],y)
        if key in self.cache:
            prior=read(self.cache[key]);require(prior['native_energies']==p['native_energies'],'Exact cache energy disagreement')
        else:self.cache[key]=record
        self.points.append(record);return p

    def gradient(self,point):return self.np.asarray(read(point['assembled_gradient'])['gradient_kcal_mol_A'])

    def difference(self,value,reference):return relative(self.response,self.old,value['native_energies'],reference['native_energies'])

    def evaluate(self,metal,q,y,guard):
        np=self.np;key=self.key(metal,y)
        if key in self.cache:
            record=self.cache[key];p=read(record);self.requests.append({'metal':metal,'point':record,'reused':True});return p,record
        require(guard['passed'],'Guard unavailable before evaluation')
        self.new+=1;self.result['new_composite_evaluations_attempted']+=1
        d=self.out/('point_%03d_%s'%(self.new,metal));d.mkdir();tick=time.monotonic()
        point={'status':'unavailable','metal':metal,'q':q.tolist(),'geometry':guard,'origin':'new_continuation_cell'}
        try:
            # The frozen pilot preparer emits both representations; only the requested state is evaluated.
            write(d/'MODES.json',read(self.pilot['modes']));write(d/'PHYSICAL.json',{'source_ids':self.mech.ids,'xyz_A':y.tolist(),'origin':self.cells['La']['mapping']['inputs']})
            write(d/'GEOMETRY.json',guard)
            proposals={m:{'status':'unavailable_not_requested'} for m in ('La','Dy')}
            proposals[metal]={'status':'prepared','q':q.tolist(),'physical':pin(d/'PHYSICAL.json'),'geometry':pin(d/'GEOMETRY.json'),'modes':pin(d/'MODES.json')}
            write(d/'PROPOSALS.json',proposals);trials,tasks=self.response.prepare_cells(d,proposals,self.world)
            t=next(t for t in trials if t['metal']==metal);td=Path(t['record']['path']).parent
            own_tasks=[(read(t['GFN2'][k]),Path(t['GFN2'][k]['path']).parent) for k in ('embedded','isolated')]
            point.update(physical=t['physical'],trial=t['record'],mapping=t['mapping'],unused_paired_trial=next(v['record'] for v in trials if v['metal']!=metal))
            mmcfg={'implementation':self.cfg['implementation'],'plan':self.cfg['plan'],'python':self.fw[1]['python'],
                'common_helper':self.fw[1]['common_helper'],'classical_helper':self.fcfg['classical_helper'],
                'openmm_version':self.export['openmm_version'],'ledgers':{m:self.cells[m]['classical']['ledger'] for m in ('La','Dy')},
                'target_source_id':self.mech.c['target_id'],'trials':[t['record']]}
            write(d/'MM_CONFIG.json',mmcfg);driver=d/'mm_driver';driver.mkdir()
            command=[str(verify(mmcfg['python'])),str(verify(self.cfg['implementation'])),'--stage','mm','--config',str(d/'MM_CONFIG.json'),'--output',str(d/'MM')]
            execution=self.backend.call(command,driver,self.env);point['MM_execution']=pin(driver/'EXECUTION.json')
            mmresult=read(pin(d/'MM/RESULT.json'));point['MM_result']=pin(d/'MM/RESULT.json')
            self.result['new_MM_configurations_attempted']+=mmresult.get('configurations_attempted',0)
            self.result['new_MM_group_queries_attempted']+=mmresult.get('group_queries_attempted',0)
            require(execution['returncode']==0 and execution['owned_process_group_empty'] and mmresult['status']=='complete'
                and mmresult['config']==pin(d/'MM_CONFIG.json') and len(mmresult['rows'])==1,'MM cell unavailable')
            mmrow=mmresult['rows'][0];require(mmrow['task_id']==t['task_id'] and mmrow['metal']==metal,'MM cell state differs')
            if self.native is None:self.native=self.helper.NativeMACE(self.mcfg);self.result['new_native_model_loads']+=1
            c=dict(self.cells[metal]);c.update(qm_positions=np.asarray(self.local.xyz(t['xyz'])[1]),xyz=t['xyz'])
            self.result['new_MACE_calls_attempted']+=1;ml=self.old.mace_cell(self.native,c,metal);write(td/'MACE.json',ml)
            point['MACE']=pin(td/'MACE.json')
            def gfn(item):
                task,gd=item;row={'task':task,'status':'unavailable','attempted':False}
                try:
                    for p in [self.fw[4]['executable'],self.fw[4]['parameters'],*self.fw[4]['runtime_libraries'],*task['staged_files'].values()]:verify(p)
                    row['attempted']=True;ex=self.backend.call(task['command'],gd,self.env);row['execution']=pin(gd/'EXECUTION.json')
                    for p in task['staged_files'].values():verify(p)
                    row.update(self.local.admit(gd,task,ex,self.fw),status='accepted')
                except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
                row['artifacts']={p.name:pin(p) for p in gd.iterdir() if p.is_file()};write(gd/'RESULT.json',row);row['result']=pin(gd/'RESULT.json');return row
            with ThreadPoolExecutor(max_workers=2) as pool:gf=list(pool.map(gfn,own_tasks))
            self.result['new_GFN2_calls_attempted']+=sum(r['attempted'] for r in gf)
            for r in gf:
                if 'execution' in r:
                    ex=read(r['execution']);self.result['observed_GFN2_child_CPU_seconds']+=ex['cpu_user_seconds']+ex['cpu_system_seconds']
            g={k:r for k,r in zip(('embedded','isolated'),gf)}
            for k,r in g.items():point['GFN2_'+k]=r['result'];require(r['status']=='accepted','GFN2 cell unavailable')
            ph=read(t['physical']);c.update(physical=ph,x=np.asarray([a['xyz_A'] for a in ph]),mapping_pin=t['mapping'])
            self.result['new_ALPB_calls_attempted']+=1;sol=self.old.solvent_cell(self.binary,c,td/'solvent',self.alpb_env,np);point['ALPB']=pin(td/'solvent/RESULT.json')
            values={'MACE_eV':ml['energy_eV'],'GFN2_embedded_Eh':g['embedded']['energy_hartree'],
                'GFN2_isolated_Eh':g['isolated']['energy_hartree'],'MM_groups_kcal_mol':mmrow['components_kcal_mol'],'ALPB_Eh':sol['total_hartree']}
            mapped=self.mapper.scatter(read(t['mapping']),-np.asarray(ml['forces_eV_A'])*self.old.EV*self.old.BOHR/self.old.HA,np.zeros((1840,3)))
            mi={s:i for i,s in enumerate(mmrow['source_ids'])};order=[mi[s] for s in self.mech.ids]
            grads={'MACE_vacuum':np.asarray(mapped['gradient_kcal_mol_A']),
                'GFN2_field_difference':np.asarray(g['embedded']['mapped_physical_gradient']['gradient_kcal_mol_A'])-np.asarray(g['isolated']['mapped_physical_gradient']['gradient_kcal_mol_A']),
                'MM_ledger':sum(np.asarray(v)[order] for v in mmrow['component_gradients_kcal_mol_A'].values()),
                'fixed_FF_charge_ALPB':np.asarray(sol['gradient_hartree_per_bohr'])*self.old.HA/self.old.BOHR}
            total=sum(grads.values());require(total.shape==(2049,3) and np.isfinite(total).all(),'Nonfinite assembled gradient')
            write(td/'ASSEMBLED_GRADIENT.json',{'source_ids':self.mech.ids,'physical_mapping':t['mapping'],
                'quantity':'gradient=dE/dR','units':'kcal/mol/angstrom','gradient_kcal_mol_A':total.tolist(),
                'component_gradients_kcal_mol_A':{k:v.tolist() for k,v in grads.items()},'modes':self.pilot['modes'],
                'source_component_receipts':{k:point[k] for k in ('MACE','GFN2_embedded','GFN2_isolated','ALPB','MM_result')},'new_gradient_only_calls':0})
            point.update(status='complete',native_energies=values,assembled_gradient=pin(td/'ASSEMBLED_GRADIENT.json'))
            cell={k:point[k] for k in ('status','metal','trial','native_energies','assembled_gradient','MACE','GFN2_embedded','GFN2_isolated','ALPB','MM_result')}
            write(td/'CELL.json',cell);point['cell']=pin(td/'CELL.json')
        except Exception as exc:point.update(failure=repr(exc),traceback=traceback.format_exc())
        point['elapsed_seconds']=time.monotonic()-tick;write(d/'POINT.json',point);record=pin(d/'POINT.json')
        self.points.append(record);self.requests.append({'metal':metal,'point':record,'reused':False})
        if point['status']=='complete':self.cache[key]=record
        require(point['status']=='complete','Composite cell failed: '+point.get('failure','unknown'))
        return point,record

def search(evaluator,metal,out):
    """Original LBT finite L-BFGS loop, warm-started in the original source metric."""
    ev=evaluator;np=ev.np;mech=ev.mech;T=mech.T
    start=ev.pilot_points[(metal,metal)];q=np.asarray(start['q']);y=ev.coordinates(start);u=np.linalg.solve(T,q)
    record=ev.cache[ev.key(metal,y)];point=read(record);warm=point;gradient=ev.gradient(point)
    jac=mech.jacobian(y);gu=T.T@jac.reshape(-1,7).T@gradient.reshape(-1);history=[]
    current_guard=mech.check(q,y,q,y);require(current_guard['passed'],'Warm-start guard failed')
    original=ev.response.origin_native(metal,ev.world)
    result={'metal':metal,'status':'unavailable','candidate_admitted':False,'mode_order':mech.modes['mode_order'],
        'coordinate_scaling_T':T.tolist(),'source_column_norms':mech.norm.tolist(),'initial_history_length':0,
        'continuation_start':record,'authoritative_origin':ev.cells[metal]['mapping']['inputs'],
        'pilot_result':ev.cfg['pilot_result'],'settings':SETTINGS,'accepted_steps':0,'trace':[],
        'origin_policy':'immutable original source; exact admitted pilot endpoint and gradient are the warm start',
        'stationary_minimum_claimed':False}
    try:
        while True:
            diagnostic=mech.diagnostic(y,gradient)
            row={'accepted_step':result['accepted_steps'],'u_A':u.tolist(),'full_q':q.tolist(),'point':record,
                'work_from_original_source':relative(ev.response,ev.old,point['native_energies'],original),
                'work_from_continuation_start':ev.difference(point,warm),'geometry':current_guard,'force_diagnostic':diagnostic,'trials':[]}
            result['trace'].append(row)
            if diagnostic['low_active_span_force']:termination='low_active_span_force';break
            if result['accepted_steps']>=12:termination='accepted_step_limit';break
            if ev.new>=24:termination='evaluation_limit';break
            direction=mech.geometry_module.lbfgs_direction(gu,history)
            if float(gu@direction)>=-1e-14:history=[];direction=-gu
            dq=T@direction
            # Sum of all per-mode atom arc bounds; actual incremental geometry is also gated.
            speed=float(np.sum(np.linalg.norm(jac,axis=1)*np.abs(dq),axis=1).max())
            if not math.isfinite(speed) or speed<=1e-14:termination='zero_active_direction';break
            scale=min(1.,.02/speed);angle=float(np.abs(dq[:4]).max())
            if angle>0:scale=min(scale,float(np.deg2rad(.5))/angle)
            direction*=scale;slope=float(gu@direction);require(slope<0,'Own-energy direction is not descent')
            accepted=None
            for multiplier in SETTINGS['trial_multipliers']:
                trial_u=u+multiplier*direction;trial_q=T@trial_u;trial_y=mech.coordinates(trial_q)
                gc=mech.check(trial_q,trial_y,q,y)
                trial={'multiplier':multiplier,'u_A':trial_u.tolist(),'full_q':trial_q.tolist(),'geometry':gc,
                    'molecular_evaluation_requested':False,'accepted':False};row['trials'].append(trial)
                if not gc['passed']:trial['reason']='geometry_or_bound';continue
                if ev.new>=24 and ev.key(metal,trial_y) not in ev.cache:trial['reason']='evaluation_limit';break
                trial['molecular_evaluation_requested']=True;next_point,next_record=ev.evaluate(metal,trial_q,trial_y,gc)
                change=ev.difference(next_point,point);trial.update(point=next_record,work_from_current=change,armijo_rhs_kcal_mol=1e-4*multiplier*slope)
                if change['composite_kcal_mol']<=trial['armijo_rhs_kcal_mol']:
                    trial['accepted']=True;accepted=(trial_u,trial_q,trial_y,next_point,next_record,gc);break
                trial['reason']='own_composite_Armijo'
            write(out/('ITERATION_%02d.json'%result['accepted_steps']),row)
            if accepted is None:
                termination='evaluation_limit' if ev.new>=24 else 'no_admissible_decreasing_trial';break
            new_u,q,y,point,record,current_guard=accepted;gradient=ev.gradient(point);jac=mech.jacobian(y)
            new_gu=T.T@jac.reshape(-1,7).T@gradient.reshape(-1);s=new_u-u;difference=new_gu-gu;sy=float(s@difference)
            if sy>1e-10*np.linalg.norm(s)*np.linalg.norm(difference) and sy>1e-16:
                history.append((s,difference,1./sy));history=history[-8:]
            u,gu=new_u,new_gu;result['accepted_steps']+=1
        final=mech.check(q,y,q,y);require(final['passed'],'Final source guard failed')
        work=ev.difference(point,warm);require(work['composite_kcal_mol']<=1e-8,'Own continuation energy increased')
        result.update(status='finite_proposal_available',candidate_admitted=True,termination=termination,candidate=record,
            full_q=q.tolist(),u_A=u.tolist(),final_geometry=final,
            work_from_original_source=relative(ev.response,ev.old,point['native_energies'],original),
            work_from_continuation_start=work,final_force_diagnostic=diagnostic)
    except Exception as exc:
        result.update(termination='backend_or_numerical_failure',failure=repr(exc),traceback=traceback.format_exc(),last_archived_valid_point=record)
    result.update(new_own_composite_evaluations=ev.new,points=ev.points,requests=ev.requests)
    write(out/'SEARCH.json',result);return result


def run(cfg,out,result,stage):
    require(int(os.environ['SLURM_CPUS_PER_TASK'])>=2,'Two one-thread GFN workers require allocated CPUs')
    loaded=inputs(cfg,out);response,fcfg,world,mech,pilot,points=loaded
    result.update(preflight=pin(out/'PREFLIGHT.json'),modes=pilot['modes'],pilot_result=cfg['pilot_result'],
        original_common3=pilot['comparison']['common3'],force_check_result=pilot['force_check_result'],
        physical_model='MACE_local_vacuum + GFN2_embedded_minus_isolated + existing_MM + fixed_FF_global_ALPB')
    if stage=='preflight':result['status']='complete';return
    require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),'Allocated GPU required')
    evaluator=Evaluator(cfg,out,result,loaded);searchdir=out/'La_search';searchdir.mkdir()
    sr=search(evaluator,'La',searchdir)
    result.update(rows=[{'metal':'La','status':sr['status'],'candidate_admitted':sr['candidate_admitted'],
        'search':pin(searchdir/'SEARCH.json'),'points':evaluator.points}],points=evaluator.points,
        requests=evaluator.requests,accepted_steps=sr['accepted_steps'],
        status='complete' if sr['candidate_admitted'] else 'unavailable',
        continuation_denominator=[{'metal':'La','status':sr['status']},{'metal':'Dy','status':'not_launched'}])
    if evaluator.native is not None:
        result.update(native_model_calls=evaluator.native.calls,native_device=evaluator.native.device,
            model=evaluator.mcfg['model'],versions=evaluator.native.versions)
    require(result['new_composite_evaluations_attempted']<=24 and result['new_MACE_calls_attempted']<=24
        and result['new_GFN2_calls_attempted']<=48 and result['new_MM_configurations_attempted']<=24
        and result['new_MM_group_queries_attempted']<=96 and result['new_ALPB_calls_attempted']<=24,'Scout call ceiling exceeded')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',choices=('preflight','scout','mm'),default='scout')
    parser.add_argument('--config',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    cfg=json.loads(args.config.read_text());require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Implementation differs');verify(cfg['plan'])
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'running','stage':args.stage,'implementation':pin(__file__),'config':pin(args.config),'plan':cfg['plan'],
        'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'settings':SETTINGS,
        'new_composite_evaluations_attempted':0,'new_MACE_calls_attempted':0,'new_GFN2_calls_attempted':0,
        'new_MM_configurations_attempted':0,'new_MM_group_queries_attempted':0,'new_ALPB_calls_attempted':0,
        'new_native_model_loads':0,'observed_GFN2_child_CPU_seconds':0.,'new_DFT_calls':0,'new_wholeprotein_SCC_calls':0}
    try:
        if args.stage=='mm':mm_stage(cfg,out,result)
        else:run(cfg,out,result,args.stage)
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-tick;result['process_CPU_seconds']=time.process_time()-cpu
    if args.stage!='mm':
        result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
        if args.stage=='scout':result['application_allocated_GPU_seconds']=result['elapsed_seconds']
    write(out/'RESULT.json',result);return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
