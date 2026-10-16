#!/usr/bin/env python3
"""Prepare a graph-defined whole-LBT3 torsion catalogue, without energy calls."""
from __future__ import annotations
import argparse
from collections import deque
import hashlib
import json
import os
from pathlib import Path
import sys
import time

SOURCES=('7CCO','7CCN')
RESNAMES='PHE ILE ASP THR ASN ASN ASP GLY TRP ILE GLU GLY ASP GLU LEU LEU ALA'.split()
RELATIVE_CUTOFF=1e-10


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(p):
    path=Path(p['path']).resolve(strict=True)
    require(path.is_file() and pin(path)['sha256']==p['sha256'],'missing/changed input: '+p['path'])
    return path


def read(p): return json.loads(verify(p).read_text())


def write(path,value):
    with Path(path).open('x') as handle:
        json.dump(value,handle,indent=2,sort_keys=True,allow_nan=False); handle.write('\n')


def component(graph,start,cut=None):
    found={start}; pending=[start]
    while pending:
        a=pending.pop()
        for b in graph[a]-found:
            if cut is not None and frozenset((a,b))==cut: continue
            found.add(b); pending.append(b)
    return found


def prepare_source(source,source_pin,out,Kinematics):
    import numpy as np
    p=read(source_pin)
    require(p['source_id']==source and p['status']=='prepared' and p['sequence']=='FIDTNNDGWIEGDELLA',
            'wrong actual source preparation')
    require(p['atoms']==256 and p['protein_atoms']==255 and p['charge']==-2 and p['protein_charge']==-5,
            'whole-source inventory/charge differs')
    require(p['waters']==p['spectators']==p['caps']==[],'unexpected waters/spectators/caps')
    atoms=read(p['physical_atoms']); native=read(p['native_peptide_atoms']); bonds=read(p['covalent_bonds'])
    verify(p['native_peptide_system']); verify(p['forcefield'])
    ids=[a['id'] for a in atoms]; ix={sid:i for i,sid in enumerate(ids)}
    require(len(ids)==256 and len(ix)==256 and ids[-1]=='metal' and atoms[-1]['element']=='M','physical map differs')
    require([a['id'] for a in native]==ids[:-1] and len(native)==255,'native particle map differs')
    require(all(a['element']==b['element'] and a['xyz_A']==b['xyz_A'] for a,b in zip(atoms,native)),
            'native and physical source coordinates differ')
    x=np.asarray([a['xyz_A'] for a in atoms],dtype=float)
    require(x.shape==(256,3) and np.isfinite(x).all(),'invalid actual coordinates')
    for metal,mult in (('La',1),('Dy',6)):
        task=p['endpoints']['whole'][metal]
        require(task['charge']==-2 and task['multiplicity']==mult and task['metal_index']==255,
                'paired physical state differs')
        text=verify(task['xyz']).read_text().splitlines()
        rows=[line.split() for line in text[2:] if line.strip()]
        require(int(text[0])==len(rows)==256,'whole XYZ atom count differs')
        require([(r[0],*[float(z) for z in r[1:]]) for r in rows]==
                [(metal if a['element']=='M' else a['element'],*a['xyz_A']) for a in atoms],
                'paired full-system XYZ differs from physical source')
    residues={i:{} for i in range(1,18)}
    for i,a in enumerate(atoms[:-1]):
        chain,resid,icode,name=a['id'].split('/'); resid=int(resid)
        require(chain=='A' and not icode and resid in residues and a['resname']==RESNAMES[resid-1],
                'source construct/numbering differs')
        require(name not in residues[resid],'duplicate residue atom')
        residues[resid][name]=i
    graph=[set() for _ in atoms]; pairs=[]
    for bond in bonds:
        a,b=[ix[bond[k]] for k in ('atom_a_id','atom_b_id')]
        require(a!=b and a<255 and b<255 and b not in graph[a],'invalid/duplicate native covalent edge')
        require(bond['native_particle_indices']==[a,b],'native/physical bond index mapping differs')
        graph[a].add(b); graph[b].add(a); pairs.append((a,b))
    require(component(graph,0)==set(range(255)) and not graph[255],'peptide connectivity or target edges differ')
    require(all(len(graph[i])==1 for i,a in enumerate(atoms) if a['element']=='H'),'hydrogen parent graph differs')
    for i in range(1,17):
        require(residues[i+1]['N'] in graph[residues[i]['C']],'missing true consecutive peptide bond')

    # Rigid conjugated groups include their attached H and peptide flanking CA.
    # An eligible axis may lie inside a group only if all remaining group atoms
    # stay on one side of the cut; no internal amide/carboxylate twist is allowed.
    protected=[]
    def protect(name,indices,hydrogen_parents=None):
        group=set(indices)
        parents=indices if hydrogen_parents is None else hydrogen_parents
        group.update(h for i in parents for h in graph[i] if atoms[h]['element']=='H')
        protected.append({'name':name,'indices':sorted(group),'source_ids':[ids[i] for i in sorted(group)]})
    for i in range(1,17):
        protect(f'peptide_{i}_{i+1}',[residues[i][n] for n in ('CA','C','O')]+
                [residues[i+1][n] for n in ('N','CA')],hydrogen_parents=[residues[i+1]['N']])
    protect('terminal_carboxylate_17',[residues[17][n] for n in ('CA','C','O','OXT')],hydrogen_parents=[])
    for i,r in residues.items():
        name=RESNAMES[i-1]
        if name=='ASP': protect(f'Asp{i}_carboxylate',[r[n] for n in ('CG','OD1','OD2')])
        if name=='GLU': protect(f'Glu{i}_carboxylate',[r[n] for n in ('CD','OE1','OE2')])
        if name=='ASN': protect(f'Asn{i}_amide',[r[n] for n in ('CG','OD1','ND2')])
        if name=='PHE': protect(f'Phe{i}_ring',[r[n] for n in ('CG','CD1','CD2','CE1','CE2','CZ')])
        if name=='TRP': protect(f'Trp{i}_ring',[r[n] for n in ('CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2')])

    modes=[]; eligibility=[]
    def candidate(name,family,resid,a,b,initial_reason=None,depth=0):
        row={'name':name,'family':family,'residue':resid,'axis_source_ids':[ids[a],ids[b]],
             'axis_indices':[a,b],'eligible':False}
        reason=initial_reason
        moving=component(graph,b,frozenset((a,b)))
        if reason is None and a in moving: reason='ring_or_nonbridge_edge'
        if reason is None and not any(atoms[j]['element']!='H' for j in moving-{a,b}):
            reason='only_hydrogens_move_off_the_axis'
        if reason is None:
            for group in protected:
                remaining=set(group['indices'])-{a,b}
                if remaining & moving and remaining-moving:
                    reason='would_split_protected_'+group['name']; break
        if reason is not None:
            row['reason']=reason; eligibility.append(row); return
        require(moving<=set(range(255)) and b in graph[a],'axis lacks actual peptide bridge')
        row.update(eligible=True,moving_source_ids=[ids[j] for j in sorted(moving)])
        eligibility.append(row)
        # Sort ancestors before descendants, with sidechains between phi/psi.
        order=(resid,{'phi':0,'sidechain':1,'psi':2}[family],depth,ids[a],ids[b])
        modes.append({'name':name,'family':family,'residue':resid,'unit':'radian',
                      'axis_indices':[a,b],'axis_source_ids':[ids[a],ids[b]],
                      'moving_indices':sorted(moving),'moving_source_ids':row['moving_source_ids'],
                      '_order':order})
    for resid,r in residues.items():
        candidate(f'A{resid}:phi:N-CA','phi',resid,r['N'],r['CA'],
                  'undefined_terminal_phi_no_preceding_peptide_C' if resid==1 else None)
        candidate(f'A{resid}:psi:CA-C','psi',resid,r['CA'],r['C'],
                  'undefined_terminal_psi_no_following_peptide_N' if resid==17 else None)
        local=set(r.values()); distance={r['CA']:0}; queue=deque([r['CA']])
        while queue:
            a=queue.popleft()
            for b in sorted(graph[a] & local):
                if atoms[b]['element']!='H' and b not in distance:
                    distance[b]=distance[a]+1; queue.append(b)
        for left,right in sorted(pairs):
            if left not in local or right not in local or atoms[left]['element']=='H' or atoms[right]['element']=='H': continue
            names={ids[left].split('/')[-1],ids[right].split('/')[-1]}
            if names in ({'N','CA'},{'CA','C'}): continue
            a,b=sorted((left,right),key=lambda j:(distance[j],ids[j]))
            # For this exact CHNO sequence every eligible sidechain rotor is C-C.
            # C-O/N conjugated bonds and backbone carbonyl bonds stay excluded.
            reason=None
            if atoms[a]['element']!='C' or atoms[b]['element']!='C':
                reason='heteroatom_or_conjugated_bond_not_a_sidechain_C_C_rotor'
            elif names <= {'N','CA','C','O','OXT'}:
                reason='backbone_bond_not_phi_or_psi'
            candidate(f'A{resid}:sidechain:{ids[a].split("/")[-1]}-{ids[b].split("/")[-1]}',
                      'sidechain',resid,a,b,reason,depth=distance[a])
    modes.sort(key=lambda m:m.pop('_order'))
    for k,mode in enumerate(modes):
        moving=set(mode['moving_indices'])
        for child in modes[k+1:]:
            other=set(child['moving_indices'])
            require(not moving & other or other<=moving and set(child['axis_indices'])<=moving,
                    'catalogue is not an ancestor-before-descendant/disjoint tree')
    for axis in range(3):
        vector=[0.,0.,0.]; vector[axis]=1.
        modes.append({'name':'metal_'+('x','y','z')[axis],'family':'metal_translation','unit':'angstrom',
                      'axis':vector,'moving_indices':[255],'moving_source_ids':['metal']})
    counts={family:sum(m['family']==family for m in modes) for family in ('phi','psi','sidechain','metal_translation')}
    require(counts=={'phi':16,'psi':16,'sidechain':29,'metal_translation':3},
            'eligible catalogue differs from declared source topology: '+str(counts))
    # Kinematics internally places target first in its auxiliary all-atom array;
    # every entry is a real physical atom, and every link is a source mapping.
    data={'positions_A':x.tolist(),'physical_ids':ids,'bonds':[list(b) for b in pairs],
          'modes':modes,'metal_index':255,'heavy_indices':[i for i,a in enumerate(atoms) if a['element']!='H'],
          'core_positions_A':x[[255,*range(255)]].tolist(),
          'core_links':[[i+1,'source',i] for i in range(255)],
          'auxiliary_atom_order':[255,*range(255)],
          'auxiliary_array_scope':'all 256 real atoms permuted target-first; no carved core or synthetic caps'}
    kin=Kinematics(data); y,aux,jac,aux_jac=kin.evaluate(np.zeros(len(modes)))
    require(np.allclose(y,x,rtol=0,atol=1e-12) and np.allclose(aux,x[data['auxiliary_atom_order']],rtol=0,atol=1e-12),
            'q0 kinematics does not reproduce exact source within floating roundoff')
    require(jac.shape==(64,256,3) and np.isfinite(jac).all() and
            np.allclose(aux_jac,jac[:,data['auxiliary_atom_order']],rtol=0,atol=1e-12),
            'analytic Jacobian/permutation invalid')
    matrix=jac.reshape(64,-1).T; scales=np.linalg.norm(matrix,axis=0)
    require(np.all(scales>1e-12),'zero physical displacement column')
    normalized=matrix/scales; gram=normalized.T@normalized
    eigenvalues,eigenvectors=np.linalg.eigh(gram)
    require(eigenvalues[-1]>0 and eigenvalues[0]>=-1e-10*eigenvalues[-1],'invalid physical displacement Gram')
    retained=eigenvalues>RELATIVE_CUTOFF*eigenvalues[-1]
    inverse=(eigenvectors[:,retained]/eigenvalues[retained])@eigenvectors[:,retained].T
    projector=eigenvectors[:,retained]@eigenvectors[:,retained].T
    out.mkdir(parents=True,exist_ok=False)
    write(out/'KINEMATICS.json',data)
    np.save(out/'JACOBIAN.npy',jac,allow_pickle=False)
    result={'source_id':source,'source_preparation':source_pin,'physical_source':p['physical_atoms'],
        'native_peptide_system':p['native_peptide_system'],'covalent_bonds':p['covalent_bonds'],
        'physical_source_ids':ids,'metal_index':255,'physical_state':{'charge':-2,'La_multiplicity':1,'Dy_multiplicity':6},
        'mode_order':[m['name'] for m in modes],'mode_counts':counts,'modes':modes,'eligibility':eligibility,
        'protected_groups':protected,'kinematics':pin(out/'KINEMATICS.json'),'origin_jacobian':pin(out/'JACOBIAN.npy'),
        'Jacobian_shape':['mode','physical_atom','Cartesian_component'],
        'Jacobian_units':'angstrom/radian for torsions; angstrom/angstrom for target translations',
        'column_norms':scales.tolist(),'physical_displacement_metric':'equal-weight Euclidean displacement of ALL256 physical atoms, including H',
        'normalized_Gram':gram.tolist(),'Gram_eigenvalues':eigenvalues.tolist(),
        'Gram_relative_cutoff':RELATIVE_CUTOFF,'retained_eigenvalue_mask':retained.tolist(),
        'retained_span_rank':int(retained.sum()),
        'retained_Gram_condition':float(eigenvalues[retained].max()/eigenvalues[retained].min()),
        'normalized_Gram_pseudoinverse':inverse.tolist(),'retained_coordinate_projector':projector.tolist(),
        'metric_interpretation':'kinematic overlap, not energy curvature or stiffness; no modes chosen using energy or force',
        'helper_derivative':'exact composed rotation/translation chain rule from Kinematics.evaluate; native source graph cuts',
        'finite_displaced_candidate_generated':False,'forces_loaded':False,'energy_calls':0,
        'constraints':'all graph bond lengths/angles, peptide and sidechain conjugated units preserved by rigid subtree rotations',
        'source_H_policy':'every deposited/completed H follows its covalent parent component; no bond normalization or independent H rotor'}
    write(out/'MODES.json',result)
    return {'source_id':source,'status':'prepared','modes':pin(out/'MODES.json'),
            'mode_order':result['mode_order'],'mode_counts':counts,'retained_span_rank':result['retained_span_rank']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('config','output','repo-root'): parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'all numerical geometry preparation requires Slurm')
    start=time.monotonic(); cfg=json.loads(args.config.read_text()); verify(cfg['plan'])
    prep=read(cfg['preparation']); require(prep['status']=='complete','completed whole-source preparation required')
    source_rows={r['source_id']:r for r in prep['rows']}
    require(set(source_rows)==set(SOURCES) and len(prep['rows'])==2,'exact two source denominator required')
    helpers={name:pin(args.repo_root/'scripts'/name) for name in
             ('mace_site_kinematics.py','affordable_response.py','affordable_common.py')}
    require(set(cfg['helpers'])==set(helpers),'explicit existing helper pins required')
    for name,p in cfg['helpers'].items():
        verify(p); require(p['sha256']==helpers[name]['sha256'],'helper differs from pinned source')
    sys.path.insert(0,str(args.repo_root.resolve()/'scripts'))
    from mace_site_kinematics import Kinematics
    out=args.output.resolve(); out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for source in SOURCES:
        try:
            require(source_rows[source]['status']=='prepared','unavailable source preparation')
            rows.append(prepare_source(source,source_rows[source]['preparation'],out/source,Kinematics))
        except Exception as exc:
            rows.append({'source_id':source,'status':'unavailable','reason':str(exc),'exception':type(exc).__name__})
        print(json.dumps(rows[-1]),flush=True)
    complete=all(r['status']=='prepared' for r in rows)
    same=complete and rows[0]['mode_order']==rows[1]['mode_order']
    result={'protocol_id':'lbt3_whole_graph_torsion_catalogue_v1','config':pin(args.config),
        'plan':cfg['plan'],'preparation':cfg['preparation'],'implementation':pin(__file__),'helpers':helpers,
        'status':'complete' if complete and same else 'incomplete','rows':rows,
        'identical_named_candidate_basis_across_sources':same,'same_modes_for_La_and_Dy':True,
        'selection':'complete eligible graph catalogue; small active basis selection deferred to actual forces',
        'source_specific_Gram_eigenvectors_are_not_shared_chemical_modes':True,
        'new_molecular_calls':0,'forces_loaded':False,'optimizer_steps':0,'finite_candidates_generated':0,
        'elapsed_seconds':time.monotonic()-start,'slurm_job_id':os.environ['SLURM_JOB_ID']}
    write(out/'RESULT.json',result)
    if result['status']!='complete': raise SystemExit(1)


if __name__=='__main__': main()
