"""Fixed graph features and organic-only MMFF dummy-center descriptor.

The center is optimized as three extra coordinates, never an MMFF atom. No
restraint energy is included in s. All candidate pools are shared between sizes.
"""
import argparse
import itertools
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski
from scipy.optimize import minimize

from curate import HERE, W, donor_patterns, read, write

PAT = {k: Chem.MolFromSmarts(v) for k,v in donor_patterns().items()}
TOPOLOGY = ['chelate5','chelate6','chelate7plus','donor_excess8',
            'donor_path_rotors','donor_path_ring_bonds','macrocycle_size']
M1 = list(PAT) + ['macrocycle','ligand_charge']


def potential_donors(m):
    ds = set()
    for k, loc in [('carboxyl',2),('amine_N',0),('pyridine_N',0),('ether_O',0),
                   ('alcohol',0),('phenol',0),('amide_O',0),('carbonyl_O',0)]:
        for match in m.GetSubstructMatches(PAT[k]):
            ds.add(match[loc])
    return sorted(ds)


def topology(m, ds):
    distance = Chem.GetDistanceMatrix(m)
    parent = {d:d for d in ds}
    def find(d):
        while parent[d]!=d:
            d=parent[d]
        return d
    edges = sorted((int(distance[a,b]),a,b) for a,b in itertools.combinations(ds,2))
    features = {c:0 for c in TOPOLOGY}
    rotatable = {tuple(sorted(x)) for x in m.GetSubstructMatches(Lipinski.RotatableBondSmarts)}
    paths = set()
    for dist,a,b in edges:
        if dist>5 or find(a)==find(b):
            continue
        parent[find(a)] = find(b)
        size = dist+2
        if size>=5:
            features['chelate5' if size==5 else 'chelate6' if size==6 else 'chelate7plus']+=1
        path=Chem.GetShortestPath(m,a,b)
        paths.update(tuple(sorted(x)) for x in zip(path,path[1:]))
    features['donor_path_rotors']=len(paths & rotatable)
    features['donor_path_ring_bonds']=sum(m.GetBondBetweenAtoms(*x).IsInRing() for x in paths)
    features['donor_excess8']=max(0,len(ds)-8)
    features['macrocycle_size']=max([len(r) for r in m.GetRingInfo().AtomRings() if len(r)>=9] or [0])
    return features


def organic_species(row, neutral):
    # Acidic O assignments are unambiguous only when their number equals HnL.
    acid = set()
    for smarts,loc in [('[CX3](=O)[OX2H1]',2),('S(=O)(=O)[OX2H1]',3),
                       ('[P](=[O])-[OX2H1]',2),('[OX2H1]-[c]',0)]:
        for match in neutral.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
            acid.add(match[loc])
    if len(acid) != -int(row['ligand_charge']):
        raise ValueError('ambiguous_deprotonation_or_tautomer')
    if neutral.HasSubstructMatch(Chem.MolFromSmarts('N=O')):
        raise ValueError('ambiguous_nitroso_phenol_tautomer')
    rw=Chem.RWMol(neutral)
    for i in acid:
        a=rw.GetAtomWithIdx(i)
        a.SetFormalCharge(-1)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(True)
    m=rw.GetMol()
    Chem.SanitizeMol(m)
    if Chem.GetFormalCharge(m)!=int(row['ligand_charge']):
        raise ValueError('species_charge_mismatch')
    # Unassigned stereochemistry with more than an enantiomeric pair can imply
    # diastereomers. A single unassigned center in an achiral medium is harmless.
    unspecified=[i for i,lab in Chem.FindMolChiralCenters(m,includeUnassigned=True) if lab=='?']
    if len(unspecified)>1:
        raise ValueError('ambiguous_relative_stereochemistry')
    return Chem.AddHs(m)


def connected_assignment(m, ds):
    if len(ds)<=1:
        return True
    D=Chem.GetDistanceMatrix(m)
    reached={ds[0]}
    while True:
        nxt={b for a in reached for b in ds if D[a,b]<=5}
        if nxt<=reached:
            break
        reached|=nxt
    return len(reached)==len(ds)


def calculate(row):
    start=time.process_time()
    out={'lig_id':row['lig_id'],'family':row['family']}
    m=Chem.MolFromSmiles(row['smiles'])
    ds=potential_donors(m)
    out.update({k:len(m.GetSubstructMatches(p)) for k,p in PAT.items()})
    out['macrocycle']=int(any(len(r)>=9 for r in m.GetRingInfo().AtomRings()))
    out['ligand_charge']=int(row['ligand_charge'])
    out.update(topology(m,ds))
    trace={'lig_id':row['lig_id'],'smiles':row['smiles'],'potential_donors':ds,'candidates':[]}
    try:
        organic=organic_species(row,m)
        if not ds:
            raise ValueError('no_supported_donors')
        if not AllChem.MMFFHasAllMoleculeParams(organic):
            raise ValueError('MMFF_parameterization_failure')
        # Keep the same assignment pool at La and Dy, including over-dentate cases.
        assignments=list(itertools.combinations(ds,min(8,len(ds))))
        assignments=[a for a in assignments if connected_assignment(m,a)]
        if not assignments:
            raise ValueError('disconnected_donor_assignment')
        if len(assignments)>8:
            assignments=[assignments[i] for i in np.linspace(0,len(assignments)-1,8,dtype=int)]
        trace['assignments']=assignments
        ps=AllChem.ETKDGv3()
        ps.randomSeed=20261030
        ps.numThreads=1
        ps.pruneRmsThresh=-1
        confids=list(AllChem.EmbedMultipleConfs(organic,numConfs=4,params=ps))
        if not confids:
            raise ValueError('conformer_embedding_failure')
        prop=AllChem.MMFFGetMoleculeProperties(organic,mmffVariant='MMFF94')
        candidates={}
        for cid in confids:
            ff=AllChem.MMFFGetMoleculeForceField(organic,prop,confId=cid)
            x0=np.asarray(organic.GetConformer(cid).GetPositions()).reshape(-1)
            for aid,assignment in enumerate(assignments):
                donor_idx=np.array(assignment)
                # Offset avoids the non-differentiable d=0 start for one donor.
                center=x0.reshape(-1,3)[donor_idx].mean(axis=0)+np.array([0.,0.,1.])
                initial=np.r_[x0,center]
                for k in (20.,10.,40.):
                    for metal,radius in (('La',1.216),('Dy',1.083)):
                        targets=np.array([radius+(1.55 if organic.GetAtomWithIdx(i).GetSymbol()=='N' else 1.40) for i in assignment])
                        def objective(v):
                            coords=v[:-3].reshape(-1,3)
                            dx=coords[donor_idx]-v[-3:]
                            distances=np.linalg.norm(dx,axis=1)
                            delta=distances-targets
                            grad_r=k*delta[:,None]*dx/np.maximum(distances[:,None],1e-10)
                            restr=.5*k*float(delta@delta)
                            energy=ff.CalcEnergy(v[:-3].tolist())
                            grad=np.array(ff.CalcGrad(v[:-3].tolist())).reshape(-1,3)
                            np.add.at(grad,donor_idx,grad_r)
                            return energy+restr,np.r_[grad.ravel(),-grad_r.sum(axis=0)]
                        opt=minimize(objective,initial,jac=True,method='L-BFGS-B',
                                     options={'maxiter':1000,'gtol':1e-4,'ftol':1e-10})
                        e=ff.CalcEnergy(opt.x[:-3].tolist())
                        delta=np.linalg.norm(opt.x[:-3].reshape(-1,3)[donor_idx]-opt.x[-3:],axis=1)-targets
                        datum={'conformer':cid,'assignment':aid,'k':k,'metal':metal,
                               'converged':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),
                               'organic_energy_kcal':float(e),'restraint_energy_kcal':float(.5*k*(delta@delta)),
                               'objective_kcal':float(opt.fun),'max_gradient':float(abs(opt.jac).max()),
                               'distance_rms_error_A':float(np.sqrt(np.mean(delta**2))),
                               'coordinates':opt.x[:-3].reshape(-1,3).tolist(),'center':opt.x[-3:].tolist()}
                        trace['candidates'].append(datum)
                        candidates[cid,aid,k,metal]=datum
        settings={'primary':(20.,4),'soft':(10.,4),'stiff':(40.,4),'two_conformers':(20.,2)}
        selected={}
        for setting,(k,nconf) in settings.items():
            # A pair enters only if both radii converged from this same start.
            pool=[(cid,aid) for cid in confids[:nconf] for aid in range(len(assignments))
                  if all(candidates[cid,aid,k,met]['converged'] for met in ('La','Dy'))]
            if not pool:
                out['s_'+setting]=''
                selected[setting]={'status':'no_matched_converged_candidate'}
                continue
            best={met:min((candidates[cid,aid,k,met] for cid,aid in pool),key=lambda x:x['objective_kcal'])
                  for met in ('La','Dy')}
            out['s_'+setting]=best['Dy']['organic_energy_kcal']-best['La']['organic_energy_kcal']
            selected[setting]={met:{kk:vv for kk,vv in b.items() if kk not in ('coordinates','center')}
                               for met,b in best.items()}
        trace['selected']=selected
        out['descriptor_status']='ok' if out['s_primary']!='' else 'minimization_failure'
        out['n_conformers']=len(confids)
        out['n_assignments']=len(assignments)
        out['n_minimizations']=len(candidates)
    except (ValueError,RuntimeError) as e:
        out['descriptor_status']=str(e)
        trace['failure']=str(e)
    out['process_cpu_seconds']=time.process_time()-start
    dest=W/'geometry'/f"{row['lig_id']}.json"
    dest.parent.mkdir(exist_ok=True)
    dest.write_text(json.dumps(trace,indent=2)+'\n')
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'Numerical descriptor work must use Slurm'
    rows=read(W/'primary.tsv')
    start=time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        out=list(pool.map(calculate,rows))
    write(W/'features.tsv',out)
    costs={'wall_seconds':time.time()-start,'worker_cpu_seconds':sum(r['process_cpu_seconds'] for r in out),
           'job_id':os.environ['SLURM_JOB_ID'],'rows':len(out),'coverage':sum(r['descriptor_status']=='ok' for r in out)}
    (W/'structure_cost.json').write_text(json.dumps(costs,indent=2)+'\n')
    print(json.dumps(costs,indent=2),flush=True)


if __name__=='__main__':
    main()
