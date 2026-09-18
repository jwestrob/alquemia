"""Analytical source-source component of the pinned native GK expression."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_density_gk_hybrid import accepted

PROTOCOL='native_GK_source_monopole_component_accounting_v1'


def component(collection,case,variant):
    r=read_json(collection);m=read_json(verify(r['manifest']));root=Path(r['manifest']['path']).parent
    if not r['complete'] or not r['numerical_pass']:raise InvalidArtifact('qualified actual output required')
    values={};born=None
    for metal in ('Ca','La'):
        t=next(t for t in m['static_tasks'] if t['case_id']==case and t['state']==metal and t['variant']==variant)
        item=accepted(t,root/'tasks'/t['task_id'])
        if item is None:raise InvalidArtifact('actual static GK missing')
        s=item[0];b=read_json(verify(t['boundary']));ix=np.array(t['frozen_indices'])-1
        moments=np.array([a['global_'] for a in s['moments']])[ix]
        if np.any(moments[:,1:]!=0):raise InvalidArtifact('source contains nonmonopole moments')
        if s['parameters']['global'][:5]!=[1.,78.3,2.455,.3,.09] or s['energies_kcal_mol']['polarization']!=0:
            raise InvalidArtifact('unsupported native model settings')
        radii=np.array([s['fields'][str(i+1)][0] for i in ix]);q=moments[:,0]
        if born is not None and not np.array_equal(radii,born):raise InvalidArtifact('endpoint effective radii differ')
        born=radii
        xyz=np.loadtxt(verify(t['xyz']),skiprows=1,usecols=(2,3,4))[ix]
        r2=np.sum((xyz[:,None]-xyz[None,:])**2,axis=2);aa=radii[:,None]*radii[None,:]
        f=np.sqrt(r2+aa*np.exp(-r2/(2.455*aa)))
        pair=.5*s['electric']*(1.-78.3)/78.3*q[:,None]*q[None,:]/f
        self_energy=math.fsum(pair.ravel());g=r['variants'][variant]['cases'][case]['endpoints'][metal]['components']['GK_permanent_transfer_kcal']
        values[metal]=dict(source_self_kcal=self_energy,source_environment_cross_kcal=g-self_energy,full_transfer_kcal=g,
            source_ids=[b['physical_ids'][i] for i in ix],charge_e=q.tolist(),effective_Born_A=radii.tolist(),
            metal_Born_A=float(radii[next(j for j,i in enumerate(ix) if b['physical_ids'][i]=='metal')]),
            static_output=record(item[1]),electric=s['electric'])
    contrasts={k:values['Ca'][k]-values['La'][k] for k in ('source_self_kcal','source_environment_cross_kcal','full_transfer_kcal')}
    error=contrasts['source_self_kcal']+contrasts['source_environment_cross_kcal']-contrasts['full_transfer_kcal']
    if abs(error)>1e-7:raise InvalidArtifact('GK component closure failed')
    return dict(endpoints=values,contrasts=contrasts,closure_error_kcal=error)


def audit(config,output):
    start=time.monotonic();cpu=time.process_time();c=read_json(config);verify(c['plan'])
    if c['protocol']!=PROTOCOL:raise InvalidArtifact('unsupported GK accounting')
    source=verify(c['native_energy_source']).read_text();bornsource=verify(c['native_born_source']).read_text()
    for fragment in ('fc = electric * 1.0d0 * (1.0d0-dwater)/(0.0d0+1.0d0*dwater)',
                     'expterm = exp(-r2/(gkc*rb2))','esym = ci*ck*gc(1)'):
        if fragment not in source:raise InvalidArtifact('native expression source mismatch')
    if 'maxborn = 30.0d0' not in bornsource:raise InvalidArtifact('native Born model differs')
    cases={};checks=[]
    for name,pin in c['collections'].items():
        p=verify(pin);a=component(p,name,'primary');b=component(p,name,'rigid')
        err=max(abs(a['endpoints'][metal][key]-b['endpoints'][metal][key]) for metal in ('Ca','La') for key in ('source_self_kcal','source_environment_cross_kcal','full_transfer_kcal'))
        checks.append(dict(case_id=name,rigid_error_kcal=err,pass_=err<=1e-7));cases[name]=a
    contrasts={}
    for left,right in c['contrasts']:
        contrasts[left+'__minus__'+right]={k:cases[left]['contrasts'][k]-cases[right]['contrasts'][k] for k in cases[left]['contrasts']}
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,root/Path(__file__).name)
    result=dict(protocol=PROTOCOL,config=record(config),implementation=record(root/Path(__file__).name),cases=cases,contrasts=contrasts,
        checks=checks,checks_pass=all(x['pass_'] for x in checks),new_scientific_calls=0,new_score=None,
        native_source_only_execution_status='not_executed; analytical decomposition of pinned native expression',
        wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write_new(root/'result.json',result)
    if not result['checks_pass']:raise InvalidArtifact('component rigid check failed')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=audit(a.config,a.output);print({k:r[k] for k in ('checks_pass','contrasts','wall_seconds')})
