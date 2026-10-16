#!/usr/bin/env python3
"""One own-gradient Mex proposal per metal; paired common3 scoring on Slurm.
CLI --stage run|mm --config CONFIG --output FRESH. No retries or optimization.
Config pins implementation,plan,force_check_result,geometry_helper,preparation_helper.
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
AXES=(('Leu83_N_CA','A/83//N','A/83//CA'),('Asp84_CA_CB','A/84//CA','A/84//CB'),
      ('Asp86_CA_CB','A/86//CA','A/86//CB'),('Glu95_CA_CB','A/95//CA','A/95//CB'))
SETTINGS={'max_atom_displacement_A':.02,'max_torsion_degrees':.5,'gram_relative_cutoff':1e-10,
          'bond_tolerance_A':1e-10,'angle_cosine_tolerance':1e-10,'protected_distance_tolerance_A':1e-9}
COMPONENTS=('MACE_vacuum','GFN2_field_difference','MM_ledger','fixed_FF_charge_ALPB')


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


def inputs(cfg,out):
    check=read(cfg['force_check_result'])
    require(check['status']=='passed' and check['slurm_job_id']=='12709' and check['comparison']['all_passed']
        and len(check['rows'])==16 and all(r['status']=='complete' for r in check['rows'])
        and check['new_MACE_calls_attempted']==8 and check['new_GFN2_calls_attempted']==24
        and check['new_MM_configurations_attempted']==16 and check['new_ALPB_calls_attempted']==16,'Derivative prerequisite differs')
    fcfg=read(check['config']);fd=module('mex_response_fd',verify(check['implementation']))
    world=fd.inputs(fcfg,out)
    old,cells,solcfg,mcfg,helper,mapper,np,local,fw,source,export,origins,previous=world
    prep=module('mex_response_preparation_pattern',verify(cfg['preparation_helper']))
    require(prep.SETTINGS['max_atom_displacement_A']==.02 and prep.SETTINGS['max_torsion_degrees']==.5
        and prep.SETTINGS['gram_relative_cutoff']==1e-10,'Inherited proposal metric/bounds differ')
    geometry=module('mex_response_geometry',verify(cfg['geometry_helper']))
    require({k:geometry.SETTINGS[k] for k in ('new_overlap_heavy_A','new_overlap_H_A','overlap_roundoff_A')}
        =={'new_overlap_heavy_A':1.,'new_overlap_H_A':.55,'overlap_roundoff_A':1e-7},'Original severe-contact rules differ')
    return fd,fcfg,world,geometry


def proposals(cfg,out,world,geometry):
    old,cells,solcfg,mcfg,helper,mapper,np,local,fw,source,export,origins,previous=world
    from scipy.spatial import cKDTree
    c=cells['La'];x=c['x'];ids=c['ids'];ix={sid:i for i,sid in enumerate(ids)};target=ix[c['target_id']]
    native=read(export['atoms']);metadata={a['id']:a for a in native};bonds=read(export['bonds'])
    protein=[a['id'] for a in native if a['kind']=='protein_source'];pset=set(protein);graph={sid:set() for sid in ids}
    for a,b in bonds:
        require(a in graph and b in graph,'Bond lacks physical source ID');graph[a].add(b);graph[b].add(a)
    for rid,resname in ((83,'LEU'),(84,'ASP'),(86,'ASP'),(95,'GLU')):
        require(metadata[f'A/{rid}//CA']['resname']==resname,'Declared residue identity differs')
    modes=[]
    for name,a,b in AXES:
        require(b in graph[a],'Axis is not an actual source bond');seen={b};pending=[b]
        while pending:
            now=pending.pop()
            for nxt in graph[now]-seen:
                if {now,nxt}=={a,b}:continue
                seen.add(nxt);pending.append(nxt)
        require(a not in seen and seen<=pset,'Axis does not cut a complete protein-only graph bridge')
        modes.append({'name':name,'axis':[a,b],'moving_source_ids':sorted(seen),'moving_indices':[ix[sid] for sid in sorted(seen)]})
    for k,parent in enumerate(modes):
        for child in modes[k+1:]:
            pset0=set(parent['moving_source_ids']);cset=set(child['moving_source_ids'])
            require(not pset0&cset or cset<=pset0 and set(child['axis'])<=pset0,'Mode ancestry/order differs')
    for k,rid,names in ((1,84,('CG','OD1','OD2')),(2,86,('CG','OD1','OD2')),(3,95,('CD','OE1','OE2'))):
        require({f'A/{rid}//{n}' for n in names}<=set(modes[k]['moving_source_ids']),'Incomplete donor carboxylate motion')
    moving={target}|set(i for m in modes for i in m['moving_indices']);fixed=sorted(set(range(2049))-moving)
    require(all(ids[i] in pset or i==target for i in moving),'Water/spectator would move')
    groups=[];residues={}
    for a in native:
        if a['id'] in pset:residues.setdefault(a['id'].rsplit('/',1)[0],{})[a['id'].rsplit('/',1)[1]]=a
    for rid,atoms in sorted(residues.items()):
        rn=next(iter(atoms.values()))['resname'];names={'ASP':('CG','OD1','OD2'),'GLU':('CD','OE1','OE2'),
             'ASN':('CG','OD1','ND2'),'GLN':('CD','OE1','NE2')}.get(rn)
        if names:
            require(all(n in atoms for n in names),'Incomplete protected source group')
            group={atoms[n]['id'] for n in names}
            group|={n for sid in tuple(group) for n in graph[sid] if metadata[n]['element']=='H'}
            groups.append({'residue':rid,'resname':rn,'source_ids':sorted(group),'indices':[ix[sid] for sid in sorted(group)]})
        ring={'PHE':('CG','CD1','CD2','CE1','CE2','CZ'),'TYR':('CG','CD1','CD2','CE1','CE2','CZ'),
              'TRP':('CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),
              'HIS':('CG','ND1','CD2','CE1','NE2'),'HID':('CG','ND1','CD2','CE1','NE2'),
              'HIE':('CG','ND1','CD2','CE1','NE2'),'HIP':('CG','ND1','CD2','CE1','NE2')}.get(rn)
        if ring:
            require(all(n in atoms for n in ring),'Incomplete protected ring');group={atoms[n]['id'] for n in ring}
            group|={n for sid in tuple(group) for n in graph[sid] if metadata[n]['element']=='H'}
            groups.append({'residue':rid,'resname':rn,'kind':'ring','source_ids':sorted(group),'indices':[ix[sid] for sid in sorted(group)]})
        if 'OXT' in atoms:
            group={atoms[n]['id'] for n in ('CA','C','O','OXT')}
            groups.append({'residue':rid,'kind':'terminal_carboxylate','source_ids':sorted(group),'indices':[ix[sid] for sid in sorted(group)]})
    for a,b in bonds:
        if a not in pset or b not in pset:continue
        if a.rsplit('/',1)[1]=='N' and b.rsplit('/',1)[1]=='C':a,b=b,a
        if a.rsplit('/',1)[1]!='C' or b.rsplit('/',1)[1]!='N':continue
        ra,rb=a.rsplit('/',1)[0],b.rsplit('/',1)[0]
        if ra==rb:continue
        group={residues[ra][n]['id'] for n in ('CA','C','O')}|{residues[rb][n]['id'] for n in ('N','CA')}
        group|={n for n in graph[b] if metadata[n]['element']=='H'}
        groups.append({'residues':[ra,rb],'kind':'peptide','source_ids':sorted(group),'indices':[ix[sid] for sid in sorted(group)]})
    for m in modes:
        moving_ids=set(m['moving_source_ids']);axis=set(m['axis'])
        for group in groups:
            rest=set(group['source_ids'])-axis
            require(not(rest&moving_ids and rest-moving_ids),'Axis splits a protected conjugated group')
    j0=np.zeros((2049,3,7))
    for k,m in enumerate(modes):
        a,b=[ix[s] for s in m['axis']];u=x[b]-x[a];u/=np.linalg.norm(u);mi=m['moving_indices'];j0[mi,:,k]=np.cross(u,x[mi]-x[b])
    j0[target,:,4:]=np.eye(3);matrix=j0.reshape(-1,7);norm=np.linalg.norm(matrix,axis=0)
    require(np.isfinite(matrix).all() and np.all(norm>0),'Invalid origin Jacobian')
    U=matrix/norm;gram=U.T@U;values,vectors=np.linalg.eigh(gram);keep=values>1e-10*values[-1]
    require(np.any(keep),'No retained source metric span');inverse=(vectors[:,keep]/values[keep])@vectors[:,keep].T
    mode_record={'source_inputs':c['mapping']['inputs'],'bonds':export['bonds'],'physical_source_ids':ids,'axes':modes,
        'target_source_id':c['target_id'],'mode_order':[m['name'] for m in modes]+['target_x','target_y','target_z'],
        'mode_units':['radian']*4+['angstrom']*3,'origin_J':j0.tolist(),'column_norms':norm.tolist(),
        'normalized_Gram':gram.tolist(),'eigenvalues':values.tolist(),'retained_rank':int(keep.sum()),
        'relative_cutoff':1e-10,'protected_groups':groups,'moving_source_ids':[ids[i] for i in sorted(moving)],'settings':SETTINGS}
    write(out/'MODES.json',mode_record)
    pidx={s:i for i,s in enumerate(protein)};porder=[ix[s] for s in protein]
    adapted=[{'source_id':s,'residue_id':s.rsplit('/',1)[0],'name':s.rsplit('/',1)[1],
              'element':metadata[s]['element'],'resname':metadata[s]['resname']} for s in protein]
    pbonds=[(pidx[a],pidx[b]) for a,b in bonds if a in pidx and b in pidx]
    other=[i for i in range(2049) if ids[i] not in pidx]
    guard=geometry.Geometry(adapted,pbonds,x[porder],x[porder],[pidx[ids[i]] for i in sorted(moving) if ids[i] in pidx],
        [dict(c['physical'][i],xyz_A=x[i].tolist()) for i in other])
    require(guard.check(x[porder])['pass'],'Origin covalent/geometry guards unavailable')
    bi=np.asarray([(ix[a],ix[b]) for a,b in bonds]);angles=np.asarray([(ix[a],ix[b],ix[d]) for b,ns in graph.items() for a in sorted(ns) for d in sorted(ns) if a<d])
    excluded=set()
    for s,ns in graph.items():
        near=set(ns)
        for n in ns:near.update(graph[n])
        excluded.update(tuple(sorted((ix[s],ix[n]))) for n in near if n!=s)
    def cosines(y):
        a,b,d=angles.T;u,v=y[a]-y[b],y[d]-y[b];return np.einsum('ij,ij->i',u,v)/np.linalg.norm(u,axis=1)/np.linalg.norm(v,axis=1)
    original_cos=cosines(x);original_lengths=np.linalg.norm(x[bi[:,0]]-x[bi[:,1]],axis=1)
    def coordinates(q):
        y=x.copy()
        for k,m in enumerate(modes):
            a,b=[ix[s] for s in m['axis']];u=y[b]-y[a];u/=np.linalg.norm(u);mi=m['moving_indices'];pivot=y[b].copy();v=y[mi]-pivot
            y[mi]=pivot+v*np.cos(q[k])+np.cross(u,v)*np.sin(q[k])+np.outer(v@u,u)*(1-np.cos(q[k]))
        y[target]+=q[4:];return y
    def check(y,q):
        guard.other=y[other].copy();inherited=guard.check(y[porder]);contacts=[]
        for a,b in sorted(cKDTree(y).query_pairs(1.)):
            if (a not in moving and b not in moving) or (a,b) in excluded:continue
            cutoff=.55 if 'H' in (c['physical'][a]['element'],c['physical'][b]['element']) else 1.
            d=float(np.linalg.norm(y[a]-y[b]));d0=float(np.linalg.norm(x[a]-x[b]))
            if d<min(cutoff,d0)-1e-7:contacts.append({'source_ids':[ids[a],ids[b]],'source_A':d0,'candidate_A':d})
        bond=float(np.max(np.abs(np.linalg.norm(y[bi[:,0]]-y[bi[:,1]],axis=1)-original_lengths)))
        angle=float(np.max(np.abs(cosines(y)-original_cos)));protected=0.
        for group in groups:
            ii=group['indices'];protected=max(protected,float(np.max(np.abs(np.linalg.norm(y[ii][:,None]-y[ii][None,:],axis=2)-np.linalg.norm(x[ii][:,None]-x[ii][None,:],axis=2)))))
        extent=float(np.linalg.norm(y-x,axis=1).max());torsion=float(np.rad2deg(np.abs(q[:4])).max());fixed_error=float(np.max(np.abs(y[fixed]-x[fixed]),initial=0))
        return {'passed':bool(inherited['pass'] and not contacts and bond<=1e-10 and angle<=1e-10 and protected<=1e-9
                and extent<=.02+1e-12 and torsion<=.5+1e-12 and fixed_error<=1e-12),
            'inherited':inherited,'all_physical_worsened_severe_contacts':contacts,'bond_max_error_A':bond,
            'angle_cosine_max_error':angle,'protected_distance_max_error_A':protected,'max_displacement_A':extent,
            'max_torsion_degrees':torsion,'fixed_max_error_A':fixed_error,'target_displacement_A':float(np.linalg.norm(y[target]-x[target]))}
    proposals={};gradients={m:np.asarray(origins[m]['assembled']['gradient_kcal_mol_A']) for m in ('La','Dy')}
    gradients.update(common=(gradients['La']+gradients['Dy'])/2,Dy_minus_La=gradients['Dy']-gradients['La'])
    for metal in ('La','Dy'):
        g=gradients[metal];coeff=-inverse@(U.T@g.reshape(-1));qdir=coeff/norm
        upper=np.sum(np.linalg.norm(j0,axis=1)*np.abs(qdir),axis=1);maximum=float(upper.max())
        row={'selector_metal':metal,'status':'unavailable','modes':pin(out/'MODES.json')};proposals[metal]=row
        if not math.isfinite(maximum) or maximum<=0:row['failure']='zero/nonfinite retained own direction';continue
        scale=.02/maximum;torsion=float(np.abs(qdir[:4]).max())
        if torsion>0:scale=min(scale,float(np.deg2rad(.5))/torsion)
        q=scale*qdir;tangent=(matrix@q).reshape(x.shape);pred={m:float(np.sum(gg*tangent)) for m,gg in gradients.items()}
        row.update(q=q.tolist(),scale=scale,predicted_work_kcal_mol=pred,own_component_predicted_work_kcal_mol={k:float(np.sum(np.asarray(origins[metal]['assembled']['component_gradients_kcal_mol_A'][k])*tangent)) for k in COMPONENTS})
        if pred[metal]>=0:row['failure']='own tangent is not descent';continue
        y=coordinates(q);gc=check(y,q);write(out/(metal+'_GEOMETRY.json'),gc);row['geometry']=pin(out/(metal+'_GEOMETRY.json'))
        write(out/(metal+'_PHYSICAL.json'),{'source_ids':ids,'xyz_A':y.tolist(),'origin':c['mapping']['inputs']});row['physical']=pin(out/(metal+'_PHYSICAL.json'))
        row['status']='prepared' if gc['passed'] else 'unavailable_geometry'
    write(out/'PROPOSALS.json',proposals)
    return proposals


def prepare_cells(out,proposals,world):
    old,cells,solcfg,mcfg,helper,mapper,np,local,fw,source,export,origins,previous=world
    ids=cells['La']['ids'];ix={s:i for i,s in enumerate(ids)};x=cells['La']['x'];native=read(export['atoms']);boundary=read(source['boundary_mapping'])
    root=out/'cells';root.mkdir();trials=[];tasks=[]
    for proposer,p in proposals.items():
        if p['status']!='prepared':continue
        y=np.asarray(read(p['physical'])['xyz_A']);require(y.shape==(2049,3) and np.isfinite(y).all(),'Proposal physical shape')
        changed=set(np.flatnonzero(np.any(y!=x,axis=1)).tolist())
        for metal in ('La','Dy'):
            c=cells[metal];tid=proposer+'_proposal_'+metal;d=root/tid;d.mkdir()
            physical=copy.deepcopy(c['physical'])
            for a in physical:a['xyz_A']=y[ix[a['id']]].tolist()
            write(d/'physical.json',physical)
            atoms=copy.deepcopy(native)
            for a in atoms:a['xyz_A']=y[ix[a['id']]].tolist()
            write(d/'atoms_charges.json',atoms);ex=copy.deepcopy(export);ex.update(atoms=pin(d/'atoms_charges.json'),origin_export=source['export']);write(d/'EXPORT.json',ex)
            env=copy.deepcopy(c['mapping']['field']);pc=verify(c['mapping']['pointcharges']).read_text().splitlines(keepends=True);field_changed=[]
            for k,a in enumerate(env):
                if ix[a['id']] not in changed:continue
                field_changed.append(k);a['xyz_A']=y[ix[a['id']]].tolist()
                if 'source_xyz_A' in a:a['source_xyz_A']=a['xyz_A']
                if 'source' in a:a['source']['xyz_A']=a['xyz_A']
                tokens=pc[k+1].split();require(len(tokens)==4,'Four-column source field layout')
                pc[k+1]=tokens[0]+' '+' '.join('%.17g'%v for v in a['xyz_A'])+'\n'
            write(d/'environment_atoms.json',env);(d/'environment.pc').write_text(''.join(pc))
            bd=copy.deepcopy(boundary);bd['original_boundary_mapping']=source['boundary_mapping'];bd['dipole_diagnostics_scope']='historical origin only; fixed shift membership'
            for sid in bd['retained_physical_CA_coordinates']:bd['retained_physical_CA_coordinates'][sid]=y[ix[sid]].tolist()
            write(d/'boundary_mapping.json',bd)
            core=copy.deepcopy(c['mapping']['core']);qm=c['qm_positions'].copy();changed_q=[]
            for a in core:
                qi=a['qm_index']
                if a['kind']!='cap':
                    sid=mapper.target_source_id(a) if a['kind']=='metal' else a['source_id']
                    if ix[sid] not in changed:continue
                    qm[qi]=y[ix[sid]];a['xyz_A']=qm[qi].tolist()
                    if a['kind']=='metal':a['source_site']['xyz_A']=a['xyz_A']
                else:
                    ri,oi=ix[a['retained_source_id']],ix[a['omitted_source_id']]
                    if ri not in changed and oi not in changed:continue
                    r,o=y[ri],y[oi];v=o-r;dist=float(np.linalg.norm(v));require(dist>0,'Coincident cap anchors');u=v/dist
                    jo=a['length_A']/dist*(np.eye(3)-np.outer(u,u));qm[qi]=r+a['length_A']*u
                    a.update(xyz_A=qm[qi].tolist(),retained_xyz_A=r.tolist(),omitted_xyz_A=o.tolist(),jacobian_retained=(np.eye(3)-jo).tolist(),jacobian_omitted=jo.tolist())
                changed_q.append(qi)
            write(d/'core_mapping.json',core);lines=verify(c['xyz']).read_text().splitlines(keepends=True)
            for qi in changed_q:lines[qi+2]=c['symbols'][qi]+' '+' '.join('%.17g'%v for v in qm[qi])+'\n'
            (d/'core.xyz').write_text(''.join(lines));symbols,roundtrip=local.xyz(pin(d/'core.xyz'))
            require(symbols==c['symbols'] and np.array_equal(roundtrip,qm),'QM round-trip serialization differs')
            endpoint=copy.deepcopy(source['configurations']['A']['endpoints'][metal]);endpoint['xyz']=pin(d/'core.xyz')
            conf={'endpoints':{metal:endpoint},'core_mapping':pin(d/'core_mapping.json'),'pointcharges':pin(d/'environment.pc')}
            inp=copy.deepcopy(source);inp.update(origin_inputs=c['mapping']['inputs'],coordinate_preparation='single_own_gradient_Mex_response_v1',
                export=pin(d/'EXPORT.json'),environment_atoms=pin(d/'environment_atoms.json'),boundary_mapping=pin(d/'boundary_mapping.json'),
                configurations={'A':conf},physical_source=pin(d/'atoms_charges.json'),physical_sources={'A':pin(d/'atoms_charges.json')},core_mapping={'A':pin(d/'core_mapping.json')})
            write(d/'INPUTS.json',inp);mp=mapper.prepare_mapping(d/'INPUTS.json','A',metal)
            require([a['id'] for a in mp['physical_atoms']]==ids and np.array_equal([a['xyz_A'] for a in mp['physical_atoms']],y),'Physical/QM/field/cap map disagreement')
            require([a['charge_e'] for a in mp['field']]==[a['charge_e'] for a in c['mapping']['field']],'Shifted field charges changed')
            write(d/'MAPPING.json',mp)
            desc=copy.deepcopy(read(origins[metal]['GFN2']['embedded']['task']['field']))
            desc.update(origin_field=origins[metal]['GFN2']['embedded']['task']['field'],source=pin(d/'environment.pc'),source_mapping=pin(d/'environment_atoms.json'),boundary_mapping=pin(d/'boundary_mapping.json'))
            for row in desc['rows']:row['xyz_A']=y[ix[row['source_id']]].tolist()
            write(d/'FIELD.json',desc)
            trial={'task_id':tid,'proposer':proposer,'metal':metal,'physical':pin(d/'physical.json'),'mapping':pin(d/'MAPPING.json'),
                'xyz':pin(d/'core.xyz'),'field':pin(d/'FIELD.json'),'source_inputs':pin(d/'INPUTS.json'),'proposal':pin(out/'PROPOSALS.json'),'GFN2':{}}
            for kind in ('embedded','isolated'):
                origin=origins[metal]['GFN2'][kind];task=copy.deepcopy(origin['task']);gd=d/('GFN2_'+kind);gd.mkdir();staged={}
                for fname,pin0 in task['staged_files'].items():
                    content=verify(pin0).read_bytes()
                    if fname=='core.xyz':content=(d/'core.xyz').read_bytes()
                    if fname=='external.pc':
                        raw=content.decode().splitlines(keepends=True)
                        for k in field_changed:
                            tok=raw[k+1].split();require(len(tok)==5,'Explicit GAM field layout')
                            raw[k+1]=' '.join([tok[0],*('%.17g'%v for v in env[k]['xyz_A']),tok[4]])+'\n'
                        content=''.join(raw).encode()
                    (gd/fname).write_bytes(content);staged[fname]=pin(gd/fname)
                task.update(task_id=tid+'_'+kind,xyz=staged['core.xyz'],source_xyz=trial['xyz'],mapping=trial['mapping'],field=trial['field'],
                    source_pointcharges=pin(d/'environment.pc'),staged_files=staged,origin_task=c['field_rows'][kind])
                task.pop('task_record',None);task.pop('preparation_task',None);write(gd/'TASK.json',task)
                trial['GFN2'][kind]=pin(gd/'TASK.json');tasks.append((task,gd))
            write(d/'TRIAL.json',trial);trial['record']=pin(d/'TRIAL.json');trials.append(trial)
    for proposer in ('La','Dy'):
        pair=[t for t in trials if t['proposer']==proposer]
        if pair:
            require(len(pair)==2 and np.array_equal([a['xyz_A'] for a in read(pair[0]['physical'])],[a['xyz_A'] for a in read(pair[1]['physical'])]),'Candidate La/Dy physical pairing differs')
            q0=local.xyz(pair[0]['xyz'])[1];q1=local.xyz(pair[1]['xyz'])[1];require(np.array_equal(q0,q1),'Candidate La/Dy local pairing differs')
    require(len(trials)<=4 and len(tasks)==2*len(trials),'Declared cell count differs')
    write(out/'PRECHECK.json',{'status':'passed','modes':pin(out/'MODES.json'),'proposals':pin(out/'PROPOSALS.json'),'trials':[t['record'] for t in trials],
        'declared_slots':[{'proposer':p,'metal':m,'status':'prepared' if proposals[p]['status']=='prepared' else 'unavailable_geometry'} for p in ('La','Dy') for m in ('La','Dy')],
        'new_molecular_calls':0,'maximum_new_calls':{'MACE':4,'GFN2':8,'MM_configurations':4,'MM_groups':16,'ALPB':4}})
    return trials,tasks


# Same qualified MM primitive; only the finite trial denominator is changed.
def mm_stage(cfg,out,result):
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'],'MM interpreter binary differs')
    verify(cfg['common_helper']);sys.path.insert(0,str(ROOT/'scripts'))
    builder=module('mex_fd_classical',verify(cfg['classical_helper']))
    import numpy as np
    import openmm as mm
    from openmm import unit
    require(mm.__version__==cfg['openmm_version'],'Original OpenMM version differs')
    platform=mm.Platform.getPlatformByName('Reference');prepared={};trial_data={}
    require(len(cfg['trials']) in (2,4),'MM response denominator differs')
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



def origin_native(metal,world):
    old,cells,solcfg,mcfg,helper,mapper,np,local,fw,source,export,origins,previous=world
    a=origins[metal];sol=read(a['assembled']['ALPB'])
    return {'MACE_eV':a['MACE']['energy_eV'],'GFN2_embedded_Eh':a['GFN2']['embedded']['energy_hartree'],
        'GFN2_isolated_Eh':a['GFN2']['isolated']['energy_hartree'],'MM_groups_kcal_mol':cells[metal]['classical']['components_kcal_mol'],
        'ALPB_Eh':sol['total_hartree']}


def work(value,origin,old):
    return {'MACE_vacuum':(value['MACE_eV']-origin['MACE_eV'])*old.EV,
        'GFN2_field_difference':math.fsum([value['GFN2_embedded_Eh']-origin['GFN2_embedded_Eh'],
                                         -(value['GFN2_isolated_Eh']-origin['GFN2_isolated_Eh'])])*old.HA,
        'MM_ledger':math.fsum(value['MM_groups_kcal_mol'][k]-origin['MM_groups_kcal_mol'][k] for k in origin['MM_groups_kcal_mol']),
        'fixed_FF_charge_ALPB':(value['ALPB_Eh']-origin['ALPB_Eh'])*old.HA}


def summarize(proposals,rows,world):
    old=world[0];by={(r['proposer'],r['metal']):r for r in rows};geometries=[]
    for proposer in ('La','Dy'):
        row={'proposer':proposer,'proposal_status':proposals[proposer]['status'],'status':'unavailable',
             'predicted_work_kcal_mol':proposals[proposer].get('predicted_work_kcal_mol')}
        if all(by[(proposer,m)]['status']=='complete' for m in ('La','Dy')):
            components={m:work(by[(proposer,m)]['native_energies'],origin_native(m,world),old) for m in ('La','Dy')}
            total={m:math.fsum(components[m].values()) for m in ('La','Dy')}
            row.update(status='available',component_work_kcal_mol=components,total_work_kcal_mol=total,
                common_work_kcal_mol=(total['La']+total['Dy'])/2,Dy_minus_La_work_kcal_mol=total['Dy']-total['La'],
                own_work_kcal_mol=total[proposer],own_descent=total[proposer]<=0,
                own_actual_minus_linear_kcal_mol=total[proposer]-proposals[proposer]['predicted_work_kcal_mol'][proposer])
        geometries.append(row)
    complete=all(r['status']=='available' for r in geometries);pool={'status':'unavailable','declared_menu':['origin','La_proposal','Dy_proposal']}
    if complete:
        selected={}
        for metal in ('La','Dy'):
            menu=[('origin',0.)]+[(r['proposer']+'_proposal',r['total_work_kcal_mol'][metal]) for r in geometries]
            name,value=min(menu,key=lambda r:r[1]);selected[metal]={'geometry':name,'work_from_origin_kcal_mol':value}
        pool.update(status='available',selections=selected,
            conditional_Dy_minus_La_accommodation_change_kcal_mol=selected['Dy']['work_from_origin_kcal_mol']-selected['La']['work_from_origin_kcal_mol'],
            common_selected_work_kcal_mol=(selected['La']['work_from_origin_kcal_mol']+selected['Dy']['work_from_origin_kcal_mol'])/2)
    passed=complete and all(r['own_descent'] for r in geometries)
    return {'geometries':geometries,'common3':pool,'all_four_cells_available':complete,'both_own_descent':passed,
        'finite_pilot_passed':passed,'automatic_expansion':False,
        'scope':'Single feasible response at fixed three-Nd spectators; no stationary minimum, protein affinity or cross-protein preference claim'}


def run(cfg,out,result):
    require(int(os.environ['SLURM_CPUS_PER_TASK'])>=4,'Four one-thread workers need allocated CPUs')
    fd,fcfg,world,geometry=inputs(cfg,out)
    old,cells,solcfg,mcfg,helper,mapper,np,local,fw,source,export,origins,previous=world
    proposal=proposals(cfg,out,world,geometry);trials,tasks=prepare_cells(out,proposal,world)
    result.update(modes=pin(out/'MODES.json'),proposals=pin(out/'PROPOSALS.json'),precheck=pin(out/'PRECHECK.json'),
        static_result=fcfg['static_result'],origin_replay=pin(out/'origin_replay/PRECHECK.json'),rows=[])
    for proposer in ('La','Dy'):
        for metal in ('La','Dy'):
            result['rows'].append({'task_id':proposer+'_proposal_'+metal,'proposer':proposer,'metal':metal,'status':'unavailable_geometry'})
    if not trials:
        result.update(status='unavailable_geometry',comparison=summarize(proposal,result['rows'],world));return
    backend,env=fw[6],fw[9]
    mmcfg={'implementation':cfg['implementation'],'plan':cfg['plan'],'force_check_result':cfg['force_check_result'],
        'python':fw[1]['python'],'common_helper':fw[1]['common_helper'],'classical_helper':fcfg['classical_helper'],
        'openmm_version':export['openmm_version'],'ledgers':{m:cells[m]['classical']['ledger'] for m in cells},
        'target_source_id':cells['La']['target_id'],'trials':[t['record'] for t in trials]}
    write(out/'MM_CONFIG.json',mmcfg);mmdir=out/'mm_driver';mmdir.mkdir()
    command=[str(verify(mmcfg['python'])),str(verify(cfg['implementation'])),'--stage','mm','--config',str(out/'MM_CONFIG.json'),'--output',str(out/'MM')]
    execution=backend.call(command,mmdir,env);result['MM_execution']=pin(mmdir/'EXECUTION.json')
    mmresult=read(pin(out/'MM/RESULT.json'));result['MM_result']=pin(out/'MM/RESULT.json')
    result['new_MM_configurations_attempted']=mmresult.get('configurations_attempted',0);result['new_MM_group_queries_attempted']=mmresult.get('group_queries_attempted',0)
    require(execution['returncode']==0 and execution['owned_process_group_empty'] and mmresult['status']=='complete'
        and mmresult['config']==pin(out/'MM_CONFIG.json'),'MM stage unavailable')
    mmby={r['task_id']:r for r in mmresult['rows']};require(set(mmby)=={t['task_id'] for t in trials},'MM result denominator')
    native=helper.NativeMACE(mcfg);ml={}
    for t in trials:
        row={'status':'unavailable'};ml[t['task_id']]=row;d=Path(t['record']['path']).parent
        try:
            c=dict(cells[t['metal']]);c.update(qm_positions=np.asarray(local.xyz(t['xyz'])[1]),xyz=t['xyz'])
            result['new_MACE_calls_attempted']+=1;value=old.mace_cell(native,c,t['metal']);write(d/'MACE.json',value)
            row.update(status='complete',value=value,receipt=pin(d/'MACE.json'))
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
    result.update(native_model_calls=native.calls,native_device=native.device,model=mcfg['model'],versions=native.versions)
    def gfn_cell(item):
        task,d=item;row={'task':task,'status':'unavailable','attempted':False}
        try:
            for p in [fw[4]['executable'],fw[4]['parameters'],*fw[4]['runtime_libraries']]:verify(p)
            for p in task['staged_files'].values():verify(p)
            row['attempted']=True;execution=backend.call(task['command'],d,env);row['execution']=pin(d/'EXECUTION.json')
            for p in task['staged_files'].values():verify(p)
            row.update(local.admit(d,task,execution,fw),status='accepted')
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
        row['artifacts']={p.name:pin(p) for p in d.iterdir() if p.is_file()};write(d/'RESULT.json',row);row['result']=pin(d/'RESULT.json');return row
    with ThreadPoolExecutor(max_workers=4) as pool:gfns=list(pool.map(gfn_cell,tasks))
    gby={r['task']['task_id']:r for r in gfns};result['GFN2_results']=[r['result'] for r in gfns]
    result['new_GFN2_calls_attempted']=sum(r['attempted'] for r in gfns)
    executions=[read(r['execution']) for r in gfns if 'execution' in r]
    result['observed_GFN2_child_CPU_seconds']=math.fsum(r['cpu_user_seconds']+r['cpu_system_seconds'] for r in executions)
    alpb_env=dict(env,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=solcfg['library_directory'])
    binary=verify(previous['build']['binary']);by={r['task_id']:r for r in result['rows']};modes=read(result['modes'])
    for t in trials:
        d=Path(t['record']['path']).parent;row=by[t['task_id']];row.update(status='unavailable',trial=t['record'])
        try:
            c=dict(cells[t['metal']]);physical=read(t['physical']);c.update(physical=physical,x=np.asarray([a['xyz_A'] for a in physical]),mapping_pin=t['mapping'])
            result['new_ALPB_calls_attempted']+=1;sol=old.solvent_cell(binary,c,d/'solvent',alpb_env,np);row['ALPB']=pin(d/'solvent/RESULT.json')
            m=ml[t['task_id']];require(m['status']=='complete','MACE unavailable');row['MACE']=m['receipt'];g={}
            for kind in ('embedded','isolated'):
                g[kind]=gby[read(t['GFN2'][kind])['task_id']];row['GFN2_'+kind]=g[kind]['result'];require(g[kind]['status']=='accepted','GFN2 unavailable')
            mmrow=mmby[t['task_id']];require(mmrow['status']=='complete','MM unavailable')
            values={'MACE_eV':m['value']['energy_eV'],'GFN2_embedded_Eh':g['embedded']['energy_hartree'],
                'GFN2_isolated_Eh':g['isolated']['energy_hartree'],'MM_groups_kcal_mol':mmrow['components_kcal_mol'],'ALPB_Eh':sol['total_hartree']}
            mp=read(t['mapping']);mlg=mapper.scatter(mp,-np.asarray(m['value']['forces_eV_A'])*old.EV*old.BOHR/old.HA,np.zeros((1840,3)))
            mmi={sid:i for i,sid in enumerate(mmrow['source_ids'])};order=[mmi[sid] for sid in c['ids']]
            grads={'MACE_vacuum':np.asarray(mlg['gradient_kcal_mol_A']),
                'GFN2_field_difference':np.asarray(g['embedded']['mapped_physical_gradient']['gradient_kcal_mol_A'])-np.asarray(g['isolated']['mapped_physical_gradient']['gradient_kcal_mol_A']),
                'MM_ledger':sum(np.asarray(v)[order] for v in mmrow['component_gradients_kcal_mol_A'].values()),
                'fixed_FF_charge_ALPB':np.asarray(sol['gradient_hartree_per_bohr'])*old.HA/old.BOHR}
            total=sum(grads.values());require(total.shape==(2049,3) and np.isfinite(total).all(),'Candidate total gradient invalid')
            y=c['x'];ix={sid:i for i,sid in enumerate(c['ids'])};j=np.zeros((2049,3,7))
            for k,mode in enumerate(modes['axes']):
                a,b=[ix[sid] for sid in mode['axis']];axis=y[b]-y[a];axis/=np.linalg.norm(axis);mi=mode['moving_indices']
                j[mi,:,k]=np.cross(axis,y[mi]-y[b])
            j[ix[c['target_id']],:,4:]=np.eye(3)
            loads=(j.reshape(-1,7).T@total.reshape(-1))/np.asarray(modes['column_norms'])
            write(d/'ASSEMBLED_GRADIENT.json',{'source_ids':c['ids'],'physical_mapping':t['mapping'],'quantity':'gradient=dE/dR',
                'units':'kcal/mol/angstrom','gradient_kcal_mol_A':total.tolist(),
                'component_gradients_kcal_mol_A':{k:v.tolist() for k,v in grads.items()},
                'current_axis_loads_in_original_column_metric_kcal_mol_A':loads.tolist(),'modes':result['modes'],
                'new_gradient_only_calls':0,'source_component_receipts':{'MACE':row['MACE'],'embedded':row['GFN2_embedded'],
                    'isolated':row['GFN2_isolated'],'ALPB':row['ALPB'],'MM':result['MM_result']}})
            row.update(status='complete',native_energies=values,MM_result=result['MM_result'])
            row.update(assembled_gradient=pin(d/'ASSEMBLED_GRADIENT.json'),current_axis_loads_in_original_column_metric_kcal_mol_A=loads.tolist())
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
        write(d/'CELL.json',row);row['cell']=pin(d/'CELL.json')
    result['comparison']=summarize(proposal,result['rows'],world)
    require(result['new_MACE_calls_attempted']<=4 and result['new_GFN2_calls_attempted']<=8
        and result['new_MM_configurations_attempted']<=4 and result['new_MM_group_queries_attempted']<=16
        and result['new_ALPB_calls_attempted']<=4,'Finite call ceiling exceeded')
    result['status']='passed' if result['comparison']['finite_pilot_passed'] else 'failed'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=('run','mm'),default='run')
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);cfg=json.loads(a.config.read_text());start=time.monotonic();cpu=time.process_time()
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Implementation differs');verify(cfg['plan'])
    result={'status':'running','stage':a.stage,'implementation':pin(__file__),'config':pin(a.config),'plan':cfg['plan'],
        'force_check_result':cfg['force_check_result'],'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),
        'new_MACE_calls_attempted':0,'new_GFN2_calls_attempted':0,'new_ALPB_calls_attempted':0,
        'new_MM_configurations_attempted':0,'new_MM_group_queries_attempted':0,'optimizer_steps':0,'own_proposals_declared':2,
        'new_DFT_calls':0,'new_wholeprotein_SCC_calls':0,'settings':SETTINGS}
    try:
        if a.stage=='mm':mm_stage(cfg,out,result)
        else:run(cfg,out,result)
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-start;result['process_CPU_seconds']=time.process_time()-cpu
    if a.stage=='run':
        result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
        result['application_allocated_GPU_seconds']=result['elapsed_seconds']
    write(out/'RESULT.json',result);return 0 if result['status'] in ('passed','complete') else 1


if __name__=='__main__':raise SystemExit(main())
