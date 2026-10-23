"""Exploratory native-algebra Ewald real-block bounds; never an SPD admission."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time

for _k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[_k] = '1'

import numpy as np
from scipy.spatial import cKDTree
from scipy.special import erfc


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def checked(p):
    if pin(p['path']) != p:
        raise ValueError('pin mismatch: ' + p['path'])
    return Path(p['path'])


def tensor(displacement, ai, aj, gamma, ewald):
    """Vectorized umutual2a normalized pairT, all mutual scales=1."""
    r2 = np.sum(displacement**2, axis=1)
    if np.any(r2 <= 0):
        raise ValueError('coincident polarization sites')
    r = np.sqrt(r2)
    ex = np.exp(-(ewald*r)**2)
    bn0 = erfc(ewald*r)/r
    bn1 = (bn0 + 2*ewald/math.sqrt(math.pi)*ex)/r2
    bn2 = (3*bn1 + 4*ewald**3/math.sqrt(math.pi)*ex)/r2
    damp = gamma*r**3/np.sqrt(ai*aj)
    # Native dampthole sets both damping values to1 for damp>=50 or gamma=0.
    active = (gamma != 0) & (damp < 50)
    deficit3 = np.where(active, np.exp(-damp), 0.)
    deficit5 = np.where(active, (1+damp)*np.exp(-damp), 0.)
    rr3 = -bn1 + deficit3/r**3
    rr5 = bn2 - 3*deficit5/r**5
    result = rr3[:,None,None]*np.eye(3)[None,:,:]
    result += rr5[:,None,None]*displacement[:,:,None]*displacement[:,None,:]
    result *= np.sqrt(ai*aj)[:,None,None]
    if not np.isfinite(result).all():
        raise ValueError('nonfinite tensor')
    return result


def bounds(n, pairs, tensors, alpha, ewald, groups):
    """Only fixed site/water blocks; no comparison-matrix reweighting."""
    block = np.full(n, -1, dtype=int)
    slot = np.full(n, -1, dtype=int)
    for b, atoms in enumerate(groups):
        for s, atom in enumerate(atoms):
            if block[atom] != -1:
                raise ValueError('overlapping blocks')
            block[atom], slot[atom] = b, s
    if np.any(block < 0):
        raise ValueError('incomplete block partition')
    tself = 4*ewald**3/(3*math.sqrt(math.pi))
    diag = [np.diag(np.repeat(1-alpha[g]*tself, 3)) for g in groups]
    cross = {}
    for (i,j), tij in zip(pairs, tensors):
        a,b = int(block[i]),int(block[j])
        si,sj = 3*int(slot[i]),3*int(slot[j])
        if a == b:
            diag[a][si:si+3,sj:sj+3] -= tij
            diag[a][sj:sj+3,si:si+3] -= tij.T
        else:
            if a > b:
                a,b,si,sj,tij = b,a,sj,si,tij.T
            key = (a,b)
            if key not in cross:
                cross[key] = np.zeros((3*len(groups[a]),3*len(groups[b])))
            cross[key][si:si+3,sj:sj+3] -= tij
    row = np.zeros(len(groups))
    # Batch by dimensions to avoid one Python LAPACK dispatch per pair.
    shape_groups = {}
    for ab, mat in cross.items():
        shape_groups.setdefault(mat.shape, []).append(ab)
    for ab_list in shape_groups.values():
        for start in range(0,len(ab_list),4096):
            keys = ab_list[start:start+4096]
            mats = np.stack([cross[k] for k in keys])
            norms = np.linalg.svd(mats, compute_uv=False)[:,0]
            for (a,b), value in zip(keys,norms):
                row[a] += value
                row[b] += value
    lows = np.array([np.linalg.eigvalsh(m)[0] for m in diag])-row
    highs = np.array([np.linalg.eigvalsh(m)[-1] for m in diag])+row
    lo,hi = int(np.argmin(lows)),int(np.argmax(highs))
    return {'lower':float(lows[lo]), 'upper':float(highs[hi]),
            'lower_atoms_one_based':[i+1 for i in groups[lo]],
            'upper_atoms_one_based':[i+1 for i in groups[hi]],
            'lower_diagonal_min':float(lows[lo]+row[lo]),
            'lower_offblock_norm_sum':float(row[lo]),
            'blocks':len(groups),'cross_blocks':len(cross),
            'positive_real_lower_bound':bool(lows[lo]>0),
            'rigorous_rounding_certificate':False,
            'complete_PME_stability_admitted':False}


def read_fixture(fixture, cfg):
    lines = checked(fixture['xyz']).read_text().splitlines()
    n = int(lines[0].split()[0])
    cell = np.array([float(x) for x in lines[1].split()])
    if n != 11958 or len(lines) != n+2 or cell.tolist() != [49.324]*3+[90.]*3:
        raise ValueError('fixture size/cell mismatch')
    xyz, types, graph = [],[],[]
    for i,line in enumerate(lines[2:]):
        t = line.split()
        if int(t[0]) != i+1:
            raise ValueError('nonsequential atom indices')
        xyz.append([float(v) for v in t[2:5]])
        types.append(int(t[5]))
        graph.append([int(j)-1 for j in t[6:]])
    if sorted(types).count(fixture['metal_type']) != 1:
        raise ValueError('metal count')
    expected = {fixture['metal_type']:1,402:3980,403:7960,405:7,414:10}
    if {k:types.count(k) for k in set(types)} != expected:
        raise ValueError('type inventory')
    groups=[]
    for i,t in enumerate(types):
        if t == 402:
            hs=sorted(graph[i])
            if len(hs)!=2 or any(types[h]!=403 or graph[h]!=[i] for h in hs):
                raise ValueError('nonintact water')
            groups.append([i]+hs)
        elif t != 403:
            if graph[i]:
                raise ValueError('bonded singleton')
            groups.append([i])
    groups.sort(key=min)
    model=cfg['type_parameters']
    alpha=np.array([model[str(t)]['alpha'] for t in types])
    thole=np.array([model[str(t)]['thole'] for t in types])
    return np.array(xyz),np.array(types),alpha,thole,groups


def evaluate_fixture(fixture,cfg):
    xyz,types,alpha,thole,groups=read_fixture(fixture,cfg)
    cutoff,ewald=cfg['cutoff_A'],cfg['ewald_inverse_A']
    L=49.324
    if cutoff!=9. or ewald!=.4 or cutoff>=L/2:
        raise ValueError('undeclared Ewald settings')
    pairs=cKDTree(xyz % L,boxsize=L).query_pairs(cutoff,output_type='ndarray')
    disp=xyz[pairs[:,1]]-xyz[pairs[:,0]]
    disp-=L*np.rint(disp/L)
    gamma=np.minimum(thole[pairs[:,0]],thole[pairs[:,1]])
    gamma=np.where(gamma==0,np.maximum(thole[pairs[:,0]],thole[pairs[:,1]]),gamma)
    for metal,value in cfg['metal_oxygen_polpair'].items():
        left,right=types[pairs[:,0]],types[pairs[:,1]]
        mask=((left==int(metal))&(right==402))|((right==int(metal))&(left==402))
        gamma[mask]=value
    ts=tensor(disp,alpha[pairs[:,0]],alpha[pairs[:,1]],gamma,ewald)
    # Site diagonal is scalar; radial tensor has two equal transverse eigenvalues.
    norms=np.linalg.eigvalsh(ts)
    norms=np.max(np.abs(norms),axis=1)
    rows=np.bincount(pairs[:,0],weights=norms,minlength=len(xyz))
    rows+=np.bincount(pairs[:,1],weights=norms,minlength=len(xyz))
    diagonal=1-alpha*4*ewald**3/(3*math.sqrt(math.pi))
    atomlo,atomhi=diagonal-rows,diagonal+rows
    imin,imax=int(np.argmin(atomlo)),int(np.argmax(atomhi))
    site={'lower':float(atomlo[imin]),'upper':float(atomhi[imax]),
          'lower_atom_one_based':imin+1,'upper_atom_one_based':imax+1,
          'lower_type':int(types[imin]),'lower_diagonal':float(diagonal[imin]),
          'lower_row_norm_sum':float(rows[imin]),
          'positive_real_lower_bound':bool(atomlo[imin]>0)}
    return {'fixture':fixture,'pairs':len(pairs),'minimum_pair_distance_A':float(np.linalg.norm(disp,axis=1).min()),
            'atom_blocks':site,'water_ion_blocks':bounds(len(xyz),pairs,ts,alpha,ewald,groups)}


def controls():
    r=np.array([[2.,0.,0.]])
    a=np.array([1.])
    bare=tensor(r,a,a,np.array([0.]),0.)[0]
    if not np.allclose(bare,np.diag([.25,-.125,-.125]),rtol=0,atol=1e-15):
        raise AssertionError('vacuum dipole limit')
    rot=np.array([[0.,0.,2.]])
    other=tensor(rot,a,a,np.array([0.]),0.)[0]
    if not np.allclose(np.linalg.eigvalsh(bare),np.linalg.eigvalsh(other),rtol=0,atol=1e-15):
        raise AssertionError('rotation')
    coords=np.array([[0.,0.,0.],[2.,0.,0.],[0.,2.,0.]])
    pairs=np.array([[0,1],[0,2],[1,2]])
    aa=np.array([.4,.6,.3]);ewald=.4
    tt=tensor(coords[pairs[:,1]]-coords[pairs[:,0]],aa[pairs[:,0]],aa[pairs[:,1]],np.full(3,.39),ewald)
    full=np.diag(np.repeat(1-aa*4*ewald**3/(3*math.sqrt(math.pi)),3))
    for (i,j),t in zip(pairs,tt):
        full[3*i:3*i+3,3*j:3*j+3]-=t
        full[3*j:3*j+3,3*i:3*i+3]-=t.T
    eig=np.linalg.eigvalsh(full)
    for groups in ([[0],[1],[2]],[[0,1],[2]]):
        b=bounds(3,pairs,tt,aa,ewald,groups)
        if b['lower']>eig[0]+1e-14 or b['upper']<eig[-1]-1e-14:
            raise AssertionError('full matrix enclosure')
    return {'status':'passed_analytic_and_small_matrix_controls','controls':4,'new_native_calls':0}


def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['controls','run'])
    p.add_argument('--config');p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out)
    if out.exists():
        raise ValueError('output already exists')
    start,cpu=time.monotonic(),time.process_time()
    if args.mode=='controls': result=controls()
    else:
        cfg=json.loads(Path(args.config).read_text())
        for value in cfg['authorities'].values(): checked(value)
        checked(cfg['implementation'])
        if cfg['implementation'] != pin(__file__): raise ValueError('implementation mismatch')
        review=json.loads(checked(cfg['source_review']).read_text())
        if review['status']!='passed_real_bound_source_review' or review['implementation']!=pin(__file__):
            raise ValueError('source review not bound')
        result={'status':'complete_exploratory_real_bounds','configuration':pin(args.config),
                'baseline_only_H_zero':True,'rigorous_certificate':False,
                'complete_PME_stability_admitted':False,'new_native_calls':0,
                'fixtures':[evaluate_fixture(f,cfg) for f in cfg['fixtures']]}
    result.update(wall_seconds=time.monotonic()-start,cpu_seconds=time.process_time()-cpu,threads=1)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__': main()
