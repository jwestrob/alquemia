"""Same-density charge sampling and validated GK source-self sensitivity."""
from __future__ import annotations
import argparse
import itertools
import math
from pathlib import Path
import shutil
import sys
import time
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from density_embedding import parse_potential
from mace_omol_charges import potential, quality, TOL
import mace_chelpg_sampling_panel as panel


def report(population, identity, output):
    start=time.monotonic();pop=read_json(population);ident=read_json(identity)
    mp=verify(pop['manifest']);panel.validate(mp);m=read_json(mp);c=read_json(verify(m['config']))
    im=panel.validate_identity(verify(ident['manifest']));original,oi=panel.qualified(c)
    if im['population']!=record(population):raise InvalidArtifact('identity collection belongs to another population')
    sources,_,_=panel.source_inventory(c)
    g=read_json(verify(c['gk_accounting']));native=read_json(verify(c['native_validation']))
    if not(g['checks_pass'] and native['checks_pass'] and native['complete']):raise InvalidArtifact('native GK validation unavailable')
    nm=read_json(verify(native['manifest']))
    if nm['analysis']!=c['gk_accounting']:raise InvalidArtifact('native validation incompatible')
    gcfg=read_json(verify(g['config']));kernels={};checks=[]
    for case in panel.CASES:
        collection=read_json(verify(gcfg['collections'][case])); cm=read_json(verify(collection['manifest']))
        kernels[case]={}
        for metal in ('Ca','La'):
            e=g['cases'][case]['endpoints'][metal];s=read_json(verify(e['static_output']))
            ct=next(t for t in cm['static_tasks'] if t['case_id']==case and t['state']==metal and t['variant']=='primary')
            b=read_json(verify(ct['boundary']));order=[b['physical_ids'].index(pid) for pid in e['source_ids']]
            coords=np.array([s['parameters']['atoms'][i]['xyz_A'] for i in order]);a=np.array(e['effective_Born_A'])
            d=np.sum((coords[:,None]-coords[None,:])**2,axis=2);aa=a[:,None]*a[None,:]
            kernel=.5*e['electric']*(1.-78.3)/78.3/np.sqrt(d+aa*np.exp(-d/(2.455*aa)))
            q=np.array(e['charge_e']);error=math.fsum((q[:,None]*q[None,:]*kernel).ravel())-e['source_self_kcal']
            if abs(error)>1e-7:raise InvalidArtifact('native source-self kernel mismatch')
            kernels[case][metal]=(kernel,e['source_ids'],coords)
            checks.append(dict(name='native_kernel_identity',case=case,metal=metal,error_kcal=error,pass_=True))
    rows={};arrays={}
    for case,metal,setting in itertools.product(panel.CASES,('Ca','La'),panel.SETTINGS):
        tid=case+'_'+metal+'_'+setting;source=sources[case+'_'+metal]
        row=dict(case_id=case,metal=metal,setting=setting,status='invalid_or_unavailable')
        try:
            if tid in original['rows']:
                r=original['rows'][tid];it=next(t for t in read_json(verify(oi['manifest']))['tasks'] if t['task_id']==tid);ir=oi['rows'][tid]
                default_error=r['default_charge_max_error_e']
                row['reused_qualification']=True
            else:
                t=next(t for t in m['tasks'] if t['task_id']==tid);r=panel.population_row(t,mp)
                if r!=pop['rows'][tid]:raise InvalidArtifact('population collection changed')
                it=next(t for t in im['tasks'] if t['task_id']==tid);ir=ident['rows'][tid];default_error=r['charge_change_from_original_max_e']
                row['reused_qualification']=False
            if not r['preliminary_identity_pass']:raise InvalidArtifact('energy or default-charge identity gate failed')
            if ir['status']!='computed' or ir['receipt']['returncode']:raise InvalidArtifact('new density query unavailable')
            for pin in (ir['receipt']['log'],ir['receipt']['resource_usage'],r['execution_receipt'],r['output']):verify(pin)
            points=np.loadtxt(verify(it['points']),skiprows=1);exact=parse_potential(verify(ir['potential']),points)
            error=float(np.max(abs(exact-np.array(read_json(verify(it['expected']))['potential_au']))))
            if error>1e-8 or error!=ir['maximum_identity_error_au']:raise InvalidArtifact('saved density changed')
            q=np.array(r['charge_e']);atoms=xyz(verify(source['quantum_task']['xyz']));coords=np.array([a[1:] for a in atoms])
            proj=read_json(verify(source['observation_task']['projection']));pq=np.array(proj['weights'])@q;pc=np.array(proj['coordinates_A'])
            cq=abs(math.fsum(pq)-math.fsum(q));dq=float(np.max(abs(pc.T@pq-coords.T@q)))
            fitted=potential(q,coords,points);projected=potential(pq,pc,points)
            fit=quality(fitted,exact);mapped=quality(projected,exact)
            conservation=cq<=TOL['projection_charge_e'] and dq<=TOL['projection_dipole_eA']
            kernel,ids,ncoords=kernels[case][metal];order=[proj['physical_ids'].index(pid) for pid in ids]
            if set(ids)!=set(proj['physical_ids']) or np.max(abs(pc[order]-ncoords))>1e-8:raise InvalidArtifact('native physical source geometry differs')
            qq=pq[order];energy=math.fsum((qq[:,None]*qq[None,:]*kernel).ravel())
            arrays[tid]=dict(points=points,exact=exact,fitted=fitted,projected=projected)
            row.update(status='computed',charge_e=q.tolist(),projected_charge_e=pq.tolist(),projection=source['observation_task']['projection'],
                       charge_sum_e=math.fsum(q),charge_change_from_original_max_e=default_error,source_self_kcal=energy,
                       original_source_self_kcal=g['cases'][case]['endpoints'][metal]['source_self_kcal'],
                       fit_quality=fit,projected_quality=mapped,projection_charge_error_e=cq,projection_dipole_error_eA=dq,
                       physical_charge_quality_pass=fit['pass'] and mapped['pass'] and conservation,
                       density_identity_error_au=error,energy_identity_error_hartree=r['energy_identity_error_hartree'],
                       population_receipt=r['execution_receipt'],density_query_receipt=ir['receipt'],density_potential=ir['potential'])
        except (ValueError,OSError,KeyError) as exc:row['reason']=str(exc)
        rows[tid]=row
    complete=all(r['status']=='computed' for r in rows.values());paired={};differences={};screens=[];pair_quality={}
    for case in panel.CASES:
        paired[case]={}
        for setting in panel.SETTINGS:
            ca,la=(rows[case+'_'+metal+'_'+setting] for metal in ('Ca','La'))
            paired[case][setting]=ca['source_self_kcal']-la['source_self_kcal'] if ca['status']==la['status']=='computed' else None
            if all(case+'_'+metal+'_'+setting in arrays for metal in ('Ca','La')):
                ac,al=(arrays[case+'_'+metal+'_'+setting] for metal in ('Ca','La'))
                if not np.array_equal(ac['points'],al['points']):raise InvalidArtifact('paired exterior probe coordinates differ')
                fq=quality(ac['fitted']-al['fitted'],ac['exact']-al['exact']);pq=quality(ac['projected']-al['projected'],ac['exact']-al['exact'])
                pair_quality[case+'_'+setting]=dict(fit_quality=fq,projected_quality=pq,pass_=fq['pass'] and pq['pass'])
    for left,right in itertools.combinations(panel.CASES[:3],2):
        differences[left+'__minus__'+right]={s:(paired[left][s]-paired[right][s] if paired[left][s] is not None and paired[right][s] is not None else None) for s in panel.SETTINGS}
    for group,values in {**paired,**differences}.items():
        for first,second in [('finer','finest'),('original','extent'),('original','finer')]:
            delta=values[second]-values[first] if values[first] is not None and values[second] is not None else None
            primary=(first,second)!=('original','finer')
            screens.append(dict(group=group,comparison=second+'__minus__'+first,change_kcal=delta,primary=primary,tolerance_kcal=.5 if primary else None,pass_=(abs(delta)<=.5 if primary and delta is not None else None)))
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,root/Path(__file__).name)
    dependencies={name:record(sys.modules[name].__file__) for name in ('affordable_common','density_embedding','mace_omol_charges','mace_chelpg_sampling_panel','mace_chelpg_sampling')}
    result=dict(protocol=panel.PROTOCOL+'_report',population=record(population),identity=record(identity),config=m['config'],implementation=record(root/Path(__file__).name),dependencies=dependencies,
                rows=rows,paired_source_self_R_kcal=paired,GGR_structural_differences_kcal=differences,checks=checks,sensitivity_screens=screens,
                complete=complete,all_density_identities_pass=complete,pair_quality=pair_quality,
                all_charge_quality_pass=complete and all(r['physical_charge_quality_pass'] for r in rows.values()) and all(r['pass_'] for r in pair_quality.values()),
                sensitivity_pass=complete and all(s['pass_'] for s in screens if s['primary']),new_score=None,biological_classification=None,baseline_changed=False,
                wall_seconds=time.monotonic()-start,scientific_scope='sampling dependence of source-self GK only; full hybrid not recomputed; no prediction or model promotion')
    write_new(root/'result.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--population',required=True);p.add_argument('--identity',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=report(a.population,a.identity,a.output)
    print({k:r[k] for k in ('complete','all_density_identities_pass','all_charge_quality_pass','sensitivity_pass','paired_source_self_R_kcal','GGR_structural_differences_kcal','wall_seconds')})
