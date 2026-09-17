"""Fit/project embedded endpoint densities and test actual coupling probes."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from density_embedding import parse_chelpg,parse_potential
from mace_omol_charges import MODEL as OLD_MODEL,TOL,projection,potential,quality,accepted
from mace_file_checks import cached_file_checks

PROTOCOL='embedded_CHELPG_source_projection_diagnostic_v1'
MODEL={**OLD_MODEL,'probes':'union_existing_exterior_and_actual_environment_sites','density_source':'native_embedded_permanent_protein_field'}
NULLS={'aqueous_affinity_score':None,'calibrated_class':None,'reference':None,'combined_gradient':None,'baseline_changed':False}


@cached_file_checks
def sources(quantum_result):
    from mace_responsive_field import collect
    q=read_json(quantum_result);mp=verify(q['manifest']);actual=collect(mp)
    for r in (q,actual):
        verify(r['collection_implementation']);verify(r['parser_implementation'])
    if ({k:v for k,v in q.items() if k not in ('collection_implementation','parser_implementation')}!={k:v for k,v in actual.items() if k not in ('collection_implementation','parser_implementation')}
            or any(q[k]['sha256']!=actual[k]['sha256'] for k in ('collection_implementation','parser_implementation'))
            or q['status']!='complete' or not q['variational_gate_pass']):
        raise InvalidArtifact('actual embedded quantum result/variational checks required')
    qm=read_json(mp);f=read_json(verify(qm['source']));fm=read_json(verify(f['manifest']));cm=read_json(verify(fm['charges']))
    return q,qm,f,fm,cm


def probe_inventory(old_charge,old_short):
    exterior=np.loadtxt(verify(old_charge['points']),skiprows=1);env=np.loadtxt(verify(old_short['points']),skiprows=1)
    allpoints=[];lookup={};groups=[]
    for points in (exterior,env):
        ids=[]
        for row in points:
            key=tuple(row)
            if key not in lookup:lookup[key]=len(allpoints);allpoints.append(row)
            ids.append(lookup[key])
        groups.append(ids)
    return np.array(allpoints),{'exterior_indices':groups[0],'environment_indices':groups[1],
        'environment_weights':old_short['weights'],'old_exterior_points':old_charge['points'],'old_environment_points':old_short['points']}


def same_projection(actual,saved):
    if ({k:v for k,v in actual.items() if k not in ('weights','caps')}!={k:v for k,v in saved.items() if k not in ('weights','caps')}
            or not np.allclose(actual['weights'],saved['weights'],atol=1e-12,rtol=0)
            or len(actual['caps'])!=len(saved['caps'])):return False
    return all({k:v for k,v in a.items() if k!='lambda'}=={k:v for k,v in b.items() if k!='lambda'}
               and abs(a['lambda']-b['lambda'])<=1e-12 for a,b in zip(actual['caps'],saved['caps']))


def prepare(quantum,output):
    q,qm,f,fm,cm=sources(quantum);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    paths={k:verify(v) for k,v in qm['implementation'].items()}
    paths['mace_responsive_field.py']=verify(q['collection_implementation'])
    paths['mace_omol_vacuum.py']=verify(q['parser_implementation'])
    for name in ('mace_responsive_charges.py','mace_omol_charges.py'):paths[name]=Path(__file__).with_name(name)
    pins={}
    for name,path in paths.items():shutil.copyfile(path,impl/name);pins[name]=record(impl/name)
    tasks=[]
    for t in qm['tasks']:
        tid=t['task_id'];oc=next(v for v in cm['tasks'] if v['task_id']==tid);os=next(v for v in fm['tasks'] if v['task_id']==tid)
        d=out/'inputs'/tid;d.mkdir(parents=True);points,groups=probe_inventory(oc,os)
        pp=d/'points_bohr.xyz';pp.write_text(str(len(points))+'\n'+''.join(' '.join(format(v,'.12f') for v in row)+'\n' for row in points))
        gp=d/'probe_groups.json';write_new(gp,groups);files={};sources_wf={}
        for wf_key,suffix in [('gbw','gbw'),('density','densities'),('density_info','densitiesinfo')]:
            src=Path(t['output_path']).parent/('endpoint.runtime.'+suffix);dst=d/src.name;sources_wf[wf_key]=record(src);shutil.copyfile(src,dst);files[wf_key]=record(dst)
        task={'task_id':tid,'case_id':t['case_id'],'metal':t['metal'],'charge':t['charge'],'xyz':t['xyz'],
              'source_receipt':q['rows'][tid]['receipt'],'source_output':q['rows'][tid]['output'],'source_mapping':oc['source_mapping'],
              'projection':oc['projection'],'points':record(pp),'probe_groups':record(gp),'files':files,'source_wavefunctions':sources_wf,
              'state':os['state']}
        tasks.append(task)
    m={'protocol_id':PROTOCOL,'quantum':record(quantum),'agreement':qm['agreement'],'tasks':tasks,'model':MODEL,'tolerances':TOL,
       'utilities':cm['utilities'],'implementation':pins,'python':cm['python'],'coupling_tolerance_kcal_mol':1.,
       'new_DFT_calls':0,'charge_utility_calls':8,'potential_utility_calls':8,'compute_budget':None,'wall_time_limit':None,**NULLS}
    for t in tasks:t['cache_key']=key(t,m)
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def key(t,m):return cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},**{k:m[k] for k in ('quantum','model','tolerances','utilities','implementation','coupling_tolerance_kcal_mol')}})


@cached_file_checks
def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['model']!=MODEL or m['tolerances']!=TOL or m['coupling_tolerance_kcal_mol']!=1.:
        raise InvalidArtifact('embedded fit/projection settings changed')
    q,qm,f,fm,cm=sources(verify(m['quantum']))
    if m['utilities']!=cm['utilities']:raise InvalidArtifact('charge utilities changed')
    for pin in [m['agreement'],m['python'],*m['utilities'].values(),*m['implementation'].values()]:verify(pin)
    if len(m['tasks'])!=8 or {t['task_id'] for t in m['tasks']}!={t['task_id'] for t in qm['tasks']}:
        raise InvalidArtifact('embedded fit inventory changed')
    for t in m['tasks']:
        old=next(v for v in qm['tasks'] if v['task_id']==t['task_id']);oc=next(v for v in cm['tasks'] if v['task_id']==t['task_id']);os=next(v for v in fm['tasks'] if v['task_id']==t['task_id'])
        if any(t[k]!=old[k] for k in ('case_id','metal','charge','xyz')) or t['source_receipt']!=q['rows'][t['task_id']]['receipt'] or t['source_output']!=q['rows'][t['task_id']]['output']:
            raise InvalidArtifact('embedded charge electronic source differs')
        if t['projection']!=oc['projection'] or t['source_mapping']!=oc['source_mapping'] or t['state']!=os['state']:
            raise InvalidArtifact('physical projection/boundary changed')
        actual=projection(read_json(verify(t['source_mapping'])),xyz(verify(t['xyz'])))
        saved=read_json(verify(t['projection']))
        if not same_projection(actual,saved):raise InvalidArtifact('source projection does not replay')
        points,groups=probe_inventory(oc,os)
        if read_json(verify(t['probe_groups']))!=groups or not np.allclose(np.loadtxt(verify(t['points']),skiprows=1),points,atol=1e-11,rtol=0):
            raise InvalidArtifact('combined probes differ')
        for k,pin in t['files'].items():
            verify(pin);verify(t['source_wavefunctions'][k])
            source=Path(old['output_path']).parent/Path(pin['path']).name
            if record(source)!=t['source_wavefunctions'][k] or pin['sha256']!=t['source_wavefunctions'][k]['sha256']:
                raise InvalidArtifact('embedded wavefunction/density changed')
        if key(t,m)!=t['cache_key']:raise InvalidArtifact('charge cache differs')
    return {'status':'pass','tasks':8,'utility_calls':16,'manifest':record(manifest)}


def report(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};arrays={};attempts=[]
    for t in m['tasks']:
        tid=t['task_id'];good=[]
        for a in sorted((mp.parent/'execution'/tid).glob('attempt_*')):
            r=accepted(a,t,mp);attempts.append({'task_id':tid,'accepted':r is not None,'directory':str(a),'receipt':record(a/'execution.json') if (a/'execution.json').exists() else None})
            if r:good.append((a,r))
        row={'status':'unavailable'}
        try:
            if not good:raise InvalidArtifact('no accepted embedded utility receipt')
            a,r=good[-1];atoms=xyz(verify(t['xyz']));points=np.loadtxt(verify(t['points']),skiprows=1);groups=read_json(verify(t['probe_groups']))
            ex,env=groups['exterior_indices'],groups['environment_indices'];weights=read_json(verify(groups['environment_weights']))
            charges=parse_chelpg(verify(r['chelpg']['log']).read_text(),[v[0] for v in atoms],t['charge'])
            exact=parse_potential(verify(r['potential']),points);fit=potential(charges,[v[1:] for v in atoms],points)
            proj=read_json(verify(t['projection']));pq=np.array(proj['weights'])@charges;mapped=potential(pq,proj['coordinates_A'],points)
            cq=abs(float(pq.sum()-charges.sum()));dq=float(np.max(np.abs(np.array(proj['coordinates_A']).T@pq-np.array([a[1:] for a in atoms]).T@charges)))
            direct={label:float(np.dot(values[env],weights['weights_e'])*HA_TO_KCAL) for label,values in [('exact',exact),('fitted',fit),('projected',mapped)]}
            row={'status':'computed','charge_e':charges.tolist(),'charge_sum_e':float(charges.sum()),'projected_charge_e':pq.tolist(),
                'fit_quality':quality(fit[ex],exact[ex]),'projected_quality':quality(mapped[ex],exact[ex]),
                'projection_charge_error_e':cq,'projection_dipole_error_eA':dq,'projection_conservation_pass':cq<=TOL['projection_charge_e'] and dq<=TOL['projection_dipole_eA'],
                'direct_coupling_kcal_mol':direct,'execution_receipt':record(a/'execution.json'),'potential':r['potential'],
                'projection':t['projection'],'points':t['points'],'probe_groups':t['probe_groups'],
                'utility_wall_seconds':r['chelpg']['wall_seconds']+r['vpot']['wall_seconds']}
            arrays[tid]={k:v.tolist() for k,v in [('exact',exact),('fitted',fit),('projected',mapped)]}
        except (OSError,ValueError,KeyError) as exc:row['reason']=str(exc)
        rows[tid]=row
    complete=all(r['status']=='computed' for r in rows.values());pairs={};checks=[]
    if complete:
        for name in sorted({t['case_id'] for t in m['tasks']}):
            t=next(t for t in m['tasks'] if t['case_id']==name);ex=read_json(verify(t['probe_groups']))['exterior_indices']
            a,b=(arrays[name+'_'+metal] for metal in ('Ca','La'));delta={k:np.array(a[k])-b[k] for k in a}
            direct={k:rows[name+'_Ca']['direct_coupling_kcal_mol'][k]-rows[name+'_La']['direct_coupling_kcal_mol'][k] for k in ('exact','fitted','projected')}
            error=direct['projected']-direct['exact'];checks.append({'name':name,'paired_projected_minus_exact_kcal_mol':error,'pass':abs(error)<=1.})
            pairs[name]={'fit_quality':quality(delta['fitted'][ex],delta['exact'][ex]),'projected_quality':quality(delta['projected'][ex],delta['exact'][ex]),'direct_R_kcal_mol':direct}
        error=(pairs['GGR_connected']['direct_R_kcal_mol']['projected']-pairs['GGR_connected']['direct_R_kcal_mol']['exact'])-(pairs['GGR_extended']['direct_R_kcal_mol']['projected']-pairs['GGR_extended']['direct_R_kcal_mol']['exact'])
        checks.append({'name':'GGR_partition','paired_projected_minus_exact_kcal_mol':error,'pass':abs(error)<=1.})
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'status':'complete' if complete else 'incomplete','rows':rows,'paired':pairs,'coupling_checks':checks,'attempts':attempts,
        'fit_quality_gate_pass':complete and all(r['fit_quality']['pass'] for r in [*rows.values(),*pairs.values()]),
        'projection_gate_pass':complete and all(r['projected_quality']['pass'] for r in [*rows.values(),*pairs.values()]) and all(r['projection_conservation_pass'] for r in rows.values()),
        'coupling_gate_pass':complete and all(c['pass'] for c in checks),'implementation':record(__file__),**NULLS}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'potentials.json',arrays);result['potential_arrays']=record(out/'potentials.json');write_new(out/'result.json',result)
    return {k:v for k,v in result.items() if k not in ('rows','attempts')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare');a.add_argument('--quantum',required=True);a.add_argument('--output',required=True)
    for op in ('dry-run','report'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output',required=op=='report')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.quantum,a.output)
    elif a.op=='dry-run':r=validate(a.manifest)
    else:r=report(a.manifest,a.output)
    print(json.dumps(r,indent=2))
