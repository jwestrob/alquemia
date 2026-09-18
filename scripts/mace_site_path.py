"""Finite coupled donor/metal path, reusing existing MACE and quantum runners."""
from __future__ import annotations
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms,write_xyz
from mace_site_kinematics import Kinematics
from mace_omol_gradient_run import CONFIG,check_parent,model
from mace_omol_ablation import ADAPTER as MASK,COMPONENT
from mace_omol_ablation_run import SEMANTICS

POLICY='matched_vacuum_hybrid_coupled_donor_path_v1'
STAGE='coupled_donor_path'
FRACTIONS={'p025':.25,'p050':.5,'p075':.75,'p100':1.}


def membership_at(site,positions):
    from affordable_response import source_key
    from ggr_sensitivity import membership
    from mace_mechanics import physical_id
    graph=site.graph;coords={}
    selected=site.repair['selected_site'];mk=graph.key(graph.locate(selected).key,selected['atom'])
    for key,atom in graph.atoms.items():
        meta=graph.meta[key];i=site.lookup.get(physical_id(meta));coords[source_key(meta)]=np.array(tuple(atom.pos)) if i is None else positions[i]
    return membership(graph,coords,positions[site.metal],mk)


@cached_file_checks
def prepare(projections,agreement,output):
    from mace_site_coordinates import SiteCoordinates
    start=time.monotonic();cpu=time.process_time();p=read_json(projections)
    if p['status']!='complete_gradient_projection' or len(p['rows'])!=16:raise InvalidArtifact('sixteen actual projected states required')
    cfg=read_json(verify(p['config']));a=read_json(verify(cfg['transfer_assessment']));parent=read_json(verify(a['preparation']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);states={};excluded={}
    for key,ref in p['rows'].items():
        row=read_json(verify(ref));s=row['source'];site=SiteCoordinates(read_json(verify(s['repair'])),read_json(verify(s['whole_preparation'])),verify(s['core_state']['xyz']))
        kin=Kinematics.from_site(site);dire=kin.direction(row['hybrid_projected_gradient']);q=np.array(dire['q']);d=out/key;d.mkdir()
        write_new(d/'kinematics.json',kin.data);points={};base=membership_at(site,site.positions);bad=[]
        for point,f in FRACTIONS.items():
            v=f*q;full,core,_,_=kin.evaluate(v);check=kin.check(v);member=membership_at(site,full)
            if not check['pass']:raise InvalidArtifact('physical coordinate implementation check failed: '+key)
            if member!=base:bad.append(point)
            ends={}
            for kind,positions in [('full',full),('core',core)]:
                e=s[kind+'_state'];old=xyz(verify(e['xyz']));xp=d/(point+'_'+kind+'.xyz');write_xyz(xp,[(r[0],*v) for r,v in zip(old,positions)])
                ends[kind]={'xyz':record(xp),'charge':e['charge'],'spin_multiplicity':1,'metal_index':e['metal_index'],'state':e['state']}
            points[point]={'q':v.tolist(),'fraction':f,'endpoints':ends,'geometry_check':check,'membership':member}
        if bad:excluded[key]={'status':'donor_membership_changed','points':bad,'original_membership':base}
        elif dire['status']!='prepared':excluded[key]={'status':dire['status']}
        states[key]={'source_projection':ref,'case_id':s['case_id'],'metal':s['metal'],'global_id':s['global_id'],'kinematics':record(d/'kinematics.json'),
                     'direction':dire,'points':points,'original_membership':base,'status':'unsupported' if key in excluded else 'prepared'}
    result={'protocol_id':POLICY,'projections':record(projections),'agreement':record(agreement),'states':states,'excluded':excluded,
            **{k:parent[k] for k in ('software','inventory','core_report','quantum_parent')},
            'new_initial_MACE_calls':4*(16-len(excluded)),'maximum_native_DFT_calls':16-len(excluded),'maximum_native_MACE_calls':2*(16-len(excluded)),
            'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'baseline_changed':False}
    write_new(out/'preparation.json',result);return {'status':'prepared','preparation':record(out/'preparation.json'),'excluded':excluded,'initial_MACE_calls':result['new_initial_MACE_calls']}


@cached_file_checks
def prepared(path):
    p=read_json(path);parent=read_json(verify(p['projections']))
    if p['protocol_id']!=POLICY or set(p['states'])!=set(parent['rows']) or len(p['states'])!=16:raise InvalidArtifact('physical path inventory changed')
    for key,s in p['states'].items():
        if s['source_projection']!=parent['rows'][key]:raise InvalidArtifact('center gradient source changed')
        row=read_json(verify(s['source_projection']));kin=Kinematics(read_json(verify(s['kinematics'])))
        direction=kin.direction(row['hybrid_projected_gradient'])
        if direction['status']!=s['direction']['status'] or not np.allclose(direction['q'],s['direction']['q'],rtol=0,atol=1e-12):raise InvalidArtifact('fixed physical direction changed')
        if (s['status']=='unsupported')!=(key in p['excluded']):raise InvalidArtifact('unsupported state omitted from denominator')
        for point,f in FRACTIONS.items():
            v=s['points'][point];q=f*np.array(s['direction']['q']);values=kin.evaluate(q)
            if not np.array_equal(v['q'],q) or not v['geometry_check']['pass']:raise InvalidArtifact('fixed path point/check differs')
            for kind,wanted in [('full',values[0]),('core',values[1])]:
                e=v['endpoints'][kind];src=row['source'][kind+'_state'];coords=xyz(verify(e['xyz']));original=xyz(verify(src['xyz']))
                if any(e[k]!=src[k] for k in ('charge','metal_index','state')) or e['spin_multiplicity']!=1:raise InvalidArtifact('path electronic state differs')
                if [a[0] for a in coords]!=[a[0] for a in original] or not np.allclose([a[1:] for a in coords],wanted,rtol=0,atol=5e-10):raise InvalidArtifact('path source geometry differs')
                if check_atoms(coords,e['charge'])!=e['state']:raise InvalidArtifact('path charge/parity differs')
    return p


def task(key,state,point,kind,energy_only):
    e=state['points'][point]['endpoints'][kind]
    return {**copy.deepcopy(e),'task_id':kind+'__'+key+'__'+point,'case_id':state['case_id'],'metal':state['metal'],'kind':kind,
            'variant':'primary','point':point,'source_projection':state['source_projection'],'q':state['points'][point]['q'],
            'energy_only':energy_only,'descriptor_gradient_experiment':POLICY,'derivative_backend':'checkpointed',
            'energy_component':COMPONENT,'charge_feature_adapter':MASK,'output_semantics':SEMANTICS,'capture_native_readout':True}


@cached_file_checks
def initial(preparation,output):
    from mace_omol import common,seal
    p=prepared(preparation);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);jobs=[]
    for group in sorted({s['global_id'] for s in p['states'].values() if s['status']=='prepared'}):
        _,d,m=common(verify(p['inventory']),verify(p['software']),verify(p['agreement']),out/group,STAGE)
        tasks=[task(k,s,point,'full',True) for k,s in p['states'].items() if s['global_id']==group and s['status']=='prepared' for point in FRACTIONS]
        m.update(protocol_id=POLICY,phase='initial',group=group,preparation=record(preparation),core_report=p['core_report'],gradient_settings=CONFIG,model=model(verify(p['software'])),tasks=tasks)
        check_parent(m);seal(d,m);jobs.append(record(d/'manifest.json'))
    r={'preparation':record(preparation),'jobs':jobs,'new_MACE_calls':p['new_initial_MACE_calls']};write_new(out/'initial.json',r);return r


def expected_tasks(m,p):
    if m['phase']=='initial':
        return [task(k,s,point,'full',True) for k,s in p['states'].items() if s['global_id']==m['group'] and s['status']=='prepared' for point in FRACTIONS]
    if m['phase']=='native':
        a=read_json(verify(m['selection']));verify(a['implementation'])
        if a['preparation']!=m['preparation'] or a['status']!='complete_selection':raise InvalidArtifact('native selection source differs')
        return [task(k,p['states'][k],s['point'],kind,False) for k,s in a['selected'].items() if s['point']!='center' for kind in ('full','core')]
    raise InvalidArtifact('unsupported coupled-path phase')


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);p=prepared(verify(m['preparation']));wanted=expected_tasks(m,p)
    if m['protocol_id']!=POLICY or m['stage']!=STAGE or m['software']!=p['software'] or m['model']!=model(verify(p['software'])) or m['gradient_settings']!=CONFIG or m['core_report']!=p['core_report']:raise InvalidArtifact('coupled-path method differs')
    check_parent(m)
    for ref in [m['agreement'],*m['implementation'].values()]:verify(ref)
    if len(m['tasks'])!=len(wanted):raise InvalidArtifact('finite path task inventory differs')
    for t,w in zip(m['tasks'],wanted):
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if payload!=w or t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('path task/state/cache differs')
    return {'status':'pass','tasks':len(wanted),'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        accepted=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp);attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
            if r is not None:accepted.append(r)
        rows[t['task_id']]=accepted[-1] if accepted else {'status':'unavailable','energy_eV':None}
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete','manifest':record(mp),'rows':rows,'attempts':attempts}


@cached_file_checks
def select(initial,output):
    job=read_json(initial);p=prepared(verify(job['preparation']));cols=[collect(verify(ref)) for ref in job['jobs']]
    if any(c['status']!='complete' for c in cols):raise InvalidArtifact('path energies incomplete')
    actual={k:r for c in cols for k,r in c['rows'].items()};selected={}
    for key,s in p['states'].items():
        if s['status']!='prepared':continue
        row=read_json(verify(s['source_projection']));src=row['source'];kin=Kinematics(read_json(verify(s['kinematics'])))
        dg=np.array(src['DFT']['gradient_kcal_mol_per_A'])-np.load(verify(src['MACE_core']['gradient']))*EV_TO_KCAL
        values={'center':{'prediction_kcal':0.,'learned_full_change_kcal':0.,'core_tangent_kcal':0.}}
        for point in FRACTIONS:
            full,core,_,_=kin.evaluate(s['points'][point]['q']);r=actual['full__'+key+'__'+point]
            delta=(r['energy_eV']-src['MACE_full']['energy_eV'])*EV_TO_KCAL;tangent=float(np.sum(dg*(core-kin.core)))
            values[point]={'prediction_kcal':delta+tangent,'learned_full_change_kcal':delta,'core_tangent_kcal':tangent,'MACE':r}
        point=min(values,key=lambda k:values[k]['prediction_kcal']);selected[key]={'point':point,'prediction_kcal':values[point]['prediction_kcal'],'all_points':values}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);sp=out/'selection_source.py';sp.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete_selection','protocol_id':POLICY,'preparation':job['preparation'],'initial':record(initial),'collections':cols,'selected':selected,'excluded':p['excluded'],'implementation':record(sp)}
    write_new(out/'result.json',result);return {'status':result['status'],'selected':{k:{x:v[x] for x in ('point','prediction_kcal')} for k,v in selected.items()},'excluded':p['excluded']}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for op,keys in [('prepare',('projections','agreement','output')),('initial',('preparation','output')),('validate',('manifest',)),('collect',('manifest',)),('select',('initial','output'))]:
        a=sub.add_parser(op)
        for key in keys:a.add_argument('--'+key,required=True)
    a=vars(parser.parse_args());op=a.pop('command');print(json.dumps(globals()[op](**a),indent=2))
