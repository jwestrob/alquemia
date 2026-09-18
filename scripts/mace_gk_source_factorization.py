"""Symmetric algebraic attribution of a validated GK source-self difference."""
from __future__ import annotations
import argparse
import itertools
import math
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new

FACTORS=('projected_charge','source_distances','effective_Born_radii')
PROTOCOL='saved_GK_source_symmetric_factorization_v1'


def audit(config,output):
    start=time.monotonic();cpu=time.process_time();cfg=read_json(config);verify(cfg['plan'])
    if cfg['protocol']!=PROTOCOL:raise InvalidArtifact('unsupported factorization scope')
    g=read_json(verify(cfg['accounting']));loc=read_json(verify(cfg['localization']));native=read_json(verify(cfg['native_validation']))
    if not all(x['checks_pass'] for x in (g,loc,native)) or not native['complete']:raise InvalidArtifact('actual component validation required')
    nm=read_json(verify(native['manifest']))
    if nm['analysis']!=cfg['accounting']:raise InvalidArtifact('native verification pertains to another analytical result')
    left,right=cfg['left'],cfg['right'];cases=(right,left);out={};checks=[];inventory=None
    for metal in ('Ca','La'):
        qs=[];ds=[];bs=[];electric=None
        for case in cases:
            e=g['cases'][case]['endpoints'][metal];ids=sorted(e['source_ids'])
            if inventory is None:inventory=ids
            if ids!=inventory:raise InvalidArtifact('source identities differ; factorization unsupported')
            order=[e['source_ids'].index(i) for i in ids]
            p=read_json(verify(e['static_output']));full=loc['cases'][case]['physical_ids'];index={pid:i for i,pid in enumerate(full)}
            xyz=np.array([p['parameters']['atoms'][index[i]]['xyz_A'] for i in ids]);q=np.array(e['charge_e'])[order];b=np.array(e['effective_Born_A'])[order]
            qs.append(q);ds.append(np.sum((xyz[:,None]-xyz[None,:])**2,axis=2));bs.append(b)
            if electric is not None and electric!=e['electric']:raise InvalidArtifact('native constants differ')
            electric=e['electric']
        if abs(math.fsum(qs[1])-math.fsum(qs[0]))>2e-5:raise InvalidArtifact('different source formal charge')
        def value(bits):
            q=qs[bits[0]];d=ds[bits[1]];b=bs[bits[2]];aa=b[:,None]*b[None,:]
            return .5*electric*(1.-78.3)/78.3*math.fsum((q[:,None]*q[None,:]/np.sqrt(d+aa*np.exp(-d/(2.455*aa)))).ravel())
        values={bits:value(bits) for bits in itertools.product((0,1),repeat=3)}
        contributions={k:0. for k in FACTORS};paths=[]
        for perm in itertools.permutations(range(3)):
            bits=[0,0,0];steps={}
            for i in perm:
                before=values[tuple(bits)];bits[i]=1;change=values[tuple(bits)]-before
                contributions[FACTORS[i]]+=change/6;steps[FACTORS[i]]=change
            paths.append(dict(order=[FACTORS[i] for i in perm],marginal_kcal=steps))
        errors={case:values[(j,j,j)]-g['cases'][case]['endpoints'][metal]['source_self_kcal'] for j,case in enumerate(cases)}
        errors['factor_closure']=sum(contributions.values())-(values[(1,1,1)]-values[(0,0,0)])
        checks.append(dict(metal=metal,errors_kcal=errors,pass_=max(abs(x) for x in errors.values())<=1e-7))
        out[metal]=dict(algebraic_counterfactual_energies_kcal={''.join(map(str,k)):v for k,v in values.items()},bit_order=list(FACTORS),bit_states={'0':right,'1':left},paths=paths,contributions_kcal=contributions,
            total_difference_kcal=values[(1,1,1)]-values[(0,0,0)],charge_sums_e={case:math.fsum(q) for case,q in zip(cases,qs)},
            per_source_changes={pid:dict(charge_change_e=float(qs[1][i]-qs[0][i]),Born_change_A=float(bs[1][i]-bs[0][i]),right_Born_A=float(bs[0][i]),left_Born_A=float(bs[1][i])) for i,pid in enumerate(ids)})
    paired={k:out['Ca']['contributions_kcal'][k]-out['La']['contributions_kcal'][k] for k in FACTORS}
    expected=g['contrasts'][left+'__minus__'+right]['source_self_kcal'];error=sum(paired.values())-expected
    checks.append(dict(name='paired_native_closure',error_kcal=error,pass_=abs(error)<=1e-7))
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,root/Path(__file__).name)
    result=dict(protocol=PROTOCOL,config=record(config),implementation=record(root/Path(__file__).name),left=left,right=right,endpoints=out,paired_contributions_kcal=paired,native_total_difference_kcal=expected,checks=checks,checks_pass=all(c['pass_'] for c in checks),
        new_scientific_calls=0,new_score=None,interpretation='symmetric attribution within this mathematical kernel; mixed arrays are algebraic counterfactuals, not computed quantum states',
        wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write_new(root/'result.json',result)
    if not result['checks_pass']:raise InvalidArtifact('factorization does not reproduce native endpoints')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=audit(a.config,a.output);print({k:r[k] for k in ('checks_pass','paired_contributions_kcal','wall_seconds')})
