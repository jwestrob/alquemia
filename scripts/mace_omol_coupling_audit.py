"""Saved-state direct/cross electrostatic accounting, without a new score."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from affordable_environment import COULOMB_KCAL_A
from mace_hybrid import rotation

BOUNDS=(0.,6.,12.,18.,24.,30.,36.,float('inf'))


def cross_terms(physical,qm,env,coords=None):
    xyz=np.array([a['xyz_A'] for a in physical]) if coords is None else np.array(coords)
    qi=np.flatnonzero(qm);ei=np.flatnonzero(env)
    if set(qi)&set(ei):raise InvalidArtifact('QM and forcefield charge overlap')
    distances=np.linalg.norm(xyz[qi,None,:]-xyz[None,ei,:],axis=2)
    if np.min(distances)<1e-8:raise InvalidArtifact('charge/nucleus overlap')
    terms=np.zeros(len(physical));terms[ei]=COULOMB_KCAL_A*env[ei]*np.sum(qm[qi,None]/distances,axis=0)
    return terms


def audit(collection,agreement,output):
    from mace_omol_solvent import collect
    start=time.monotonic();cpu=time.process_time();r=read_json(collection);mp=verify(r['manifest']);actual=collect(mp)
    for value in (r,actual):verify(value['collection_implementation'])
    if (r['collection_implementation']['sha256']!=actual['collection_implementation']['sha256']
        or {k:v for k,v in r.items() if k!='collection_implementation'}!={k:v for k,v in actual.items() if k!='collection_implementation'}
        or r['status']!='complete' or not r['numerical_gate_pass']):
        raise InvalidArtifact('actual completed solver result differs')
    m=read_json(mp);h=read_json(verify(m['vacuum_hybrid']));cases={};checks=[]
    def check(name,error):checks.append({'name':name,'error_kcal_mol':error,'pass':abs(error)<=1e-7})
    for name,pin in m['states'].items():
        s=read_json(verify(pin));physical=s['physical_atoms'];coords=np.array([a['xyz_A'] for a in physical]);env=np.array(s['environment_charges_e'])
        center=coords[next(i for i,a in enumerate(physical) if a['id']=='metal')];distance=np.linalg.norm(coords-center,axis=1)
        terms={};endpoints={}
        for metal in ('Ca','La'):
            qm=np.array(s['endpoints'][metal]['QM_charges_e']);terms[metal]=cross_terms(physical,qm,env)
            endpoints[metal]={'direct_cross_kcal_mol':float(terms[metal].sum()),'per_atom_kcal_mol':terms[metal].tolist()}
            moved=(coords-center)@rotation().T+center+np.array([10.,-7.,3.])
            check(name+'_'+metal+'_rigid',float(cross_terms(physical,qm,env,moved).sum()-terms[metal].sum()))
        delta=terms['Ca']-terms['La'];paired=float(delta.sum());bins=[];residues={}
        for low,high in zip(BOUNDS[:-1],BOUNDS[1:]):
            mask=(distance>=low)&(distance<high)
            bins.append({'from_A':low,'to_A':high if np.isfinite(high) else None,'direct_R_kcal_mol':float(delta[mask].sum()),'environment_atom_count':int(np.count_nonzero(env[mask]))})
        for atom,value in zip(physical,delta):
            key=atom['id'].rsplit('/',1)[0];residues[key]=residues.get(key,0.)+float(value)
        check(name+'_bin_closure',sum(b['direct_R_kcal_mol'] for b in bins)-paired)
        check(name+'_residue_closure',sum(residues.values())-paired)
        check(name+'_paired_closure',endpoints['Ca']['direct_cross_kcal_mol']-endpoints['La']['direct_cross_kcal_mol']-paired)
        cross=r['cases'][name]['GB_cross_R_kcal_mol']
        cases[name]={'state':pin,'endpoints':endpoints,'direct_cross_R_kcal_mol':paired,'GB_cross_R_kcal_mol':cross,
            'direct_plus_GB_cross_R_kcal_mol':paired+cross,'MACE_context_R_model_kcal':h['cases'][name]['context_R_model_kcal'],
            'radial_bins':bins,'per_residue_R_kcal_mol':residues,'per_atom_R_kcal_mol':delta.tolist()}
    result={'protocol_id':'saved_QMFF_direct_cross_component_audit_v1','collection':record(collection),'agreement':record(agreement),
        'implementation':record(__file__),'Coulomb_constant':COULOMB_KCAL_A,'constant_source':record(Path(__file__).with_name('affordable_environment.py')),
        'cases':cases,'checks':checks,'checks_pass':all(c['pass'] for c in checks),'new_score':None,'calibrated_class':None,
        'claim':'direct and solvent cross-term accounting; not a unique decomposition of learned energy',
        'new_DFT_calls':0,'new_MACE_calls':0,'new_solver_calls':0,'baseline_changed':False,
        'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name)
    result['preserved_implementation']=record(out/Path(__file__).name)
    write_new(out/'result.json',result)
    print(json.dumps({'checks_pass':result['checks_pass'],'cases':{k:{f:v[f] for f in ('direct_cross_R_kcal_mol','GB_cross_R_kcal_mol','direct_plus_GB_cross_R_kcal_mol','MACE_context_R_model_kcal','radial_bins')} for k,v in cases.items()}},indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('collection','agreement','output'):p.add_argument('--'+k,required=True)
    a=p.parse_args();audit(a.collection,a.agreement,a.output)
