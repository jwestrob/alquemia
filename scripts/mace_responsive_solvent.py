"""Full fixed-boundary GB from responsive quantum charges and paired comparison."""
from __future__ import annotations
import argparse,copy,json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import accepted_attempt,write_xyz
from mace_omol_solvent import MODEL as OLD_MODEL,TOL,NULLS,CASES,specs,task_payload,geometry

SCHEMA='alquemia.mace_responsive_GB.v1'
PROTOCOL='embedded_QM_ff19SB_full_OBC2_POLAR_short_hybrid_v1'
MODEL={**OLD_MODEL,'vacuum_parent':None,'quantum_parent':'native_r2SCAN3c_embedded_fixed_permanent_protein_field',
       'charge_representation':'embedded_native_CHELPG_barycentric_QM_plus_identical_ff19SB_exterior',
       'energy_definition':'embedded_core_including_direct_field_plus_full_GB_plus_short_context'}


@cached_file_checks
def sources(charges):
    from mace_responsive_charges import validate,sources as charge_sources
    c=read_json(charges);mp=verify(c['manifest']);validate(mp);cm=read_json(mp)
    if c['status']!='complete' or not all(c[k] for k in ('fit_quality_gate_pass','projection_gate_pass','coupling_gate_pass')):
        raise InvalidArtifact('responsive charge representation failed/incomplete')
    q,qm,f,fm,oldcm=charge_sources(verify(cm['quantum']))
    g=read_json(verify(fm['GB']));gm=read_json(verify(g['manifest']));states={}
    for name in CASES:
        s=read_json(verify(gm['states'][name]));lookup={a['id']:i for i,a in enumerate(s['physical_atoms'])};env=np.array(s['environment_charges_e'])
        for metal in ('Ca','La'):
            tid=name+'_'+metal;r=c['rows'][tid];t=next(t for t in cm['tasks'] if t['task_id']==tid)
            receipt=read_json(verify(r['execution_receipt']))
            from mace_omol_charges import accepted
            if accepted(Path(r['execution_receipt']['path']).parent,t,mp)!=receipt:raise InvalidArtifact('embedded charge receipt invalid')
            proj=read_json(verify(t['projection']));pq=np.array(r['projected_charge_e'])
            from density_embedding import parse_chelpg
            fitted=parse_chelpg(verify(receipt['chelpg']['log']).read_text(),[a[0] for a in xyz(verify(t['xyz']))],t['charge'])
            if not np.array_equal(fitted,r['charge_e']) or not np.allclose(np.array(proj['weights'])@fitted,pq,atol=1e-12,rtol=0):
                raise InvalidArtifact('embedded charge arrays differ from actual native fit')
            if len(pq)!=len(proj['physical_ids']):raise InvalidArtifact('projected charge inventory differs')
            qmcharges=np.zeros(len(env))
            for pid,val in zip(proj['physical_ids'],pq):qmcharges[lookup[pid]]=val
            if set(proj['physical_ids'])!=set(s['projection_support_ids']) or np.any(env[np.flatnonzero(qmcharges)]!=0):
                raise InvalidArtifact('forcefield/QM charge overlap or changed support')
            e=s['endpoints'][metal];full=qmcharges+env
            if abs(float(qmcharges.sum())-t['charge'])>TOL['charge_e'] or abs(float(full.sum())-e['charge'])>TOL['charge_e']:
                raise InvalidArtifact('responsive charge closure failed')
            e.update(QM_charges_e=qmcharges.tolist(),full_charges_e=full.tolist(),actual_charge_receipt=r['execution_receipt'])
        states[name]=s
    return states,c,cm,q,qm,f,fm,g,gm


def experiment_specs():return [s for s in specs() if s[2]!='environment']


def prepare(charges,software,output):
    states,c,cm,q,qm,f,fm,g,gm=sources(charges);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    if record(software)!=gm['software']:raise InvalidArtifact('same validated GB software required')
    impl=out/'implementation';impl.mkdir();paths={k:verify(v) for k,v in cm['implementation'].items()}
    for name in ('mace_hybrid.py','mace_omol_solvent.py','mace_responsive_solvent.py','mace_gb.py','mace_global_prepare.py'):paths[name]=Path(__file__).with_name(name)
    pins={}
    for name,path in paths.items():shutil.copyfile(path,impl/name);pins[name]=record(impl/name)
    stored={};qfiles={}
    for name,s in states.items():
        d=out/'states'/name;d.mkdir(parents=True);write_new(d/'state.json',s);stored[name]=record(d/'state.json')
        for metal in ('Ca','La'):
            for component in ('full','QM'):
                p=d/(metal+'_'+component+'.npy');np.save(p,np.array(s['endpoints'][metal][component+'_charges_e']));qfiles[(name,metal,component)]=record(p)
    coords=out/'coordinates';coords.mkdir();tasks=[]
    for spec in experiment_specs():
        name,metal,component,variant,_=spec;p=coords/('_'.join(spec[:4])+'.xyz');write_xyz(p,geometry(states[name],metal,variant))
        tasks.append(task_payload(spec,states[name],stored[name],record(p),qfiles[(name,metal,component)]))
    m={'schema_version':SCHEMA,'protocol_id':PROTOCOL,'model':MODEL,'tolerances':TOL,'software':record(software),'agreement':qm['agreement'],
       'charges':record(charges),'quantum':cm['quantum'],'frozen_comparison':qm['source'],'states':stored,'tasks':tasks,'implementation':pins,
       'compute_budget':None,'wall_time_limit':None,'new_GB_calls':44,'reused_environment_GB_calls':4,**NULLS}
    for t in tasks:t['cache_key']=key(t,m)
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def key(t,m):return cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},**{k:m[k] for k in ('model','software','implementation','charges','quantum','frozen_comparison')}})


@cached_file_checks
def validate(manifest):
    m=read_json(manifest)
    if m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['model']!=MODEL or m['tolerances']!=TOL or len(m['tasks'])!=44:
        raise InvalidArtifact('responsive GB settings/inventory changed')
    for pin in [m['agreement'],m['software'],*m['implementation'].values()]:verify(pin)
    sm=read_json(verify(m['software']))
    for pin in [sm['python'],sm['requirements'],sm['backend_source_inventory'],*read_json(verify(sm['backend_source_inventory']))['files']]:verify(pin)
    states,c,cm,q,qm,f,fm,g,gm=sources(verify(m['charges']))
    if m['quantum']!=cm['quantum'] or m['frozen_comparison']!=qm['source'] or m['software']!=gm['software']:
        raise InvalidArtifact('responsive GB source/software mismatch')
    for name,pin in m['states'].items():
        if read_json(verify(pin))!=states[name]:raise InvalidArtifact('responsive physical charge state differs')
    for t,spec in zip(m['tasks'],experiment_specs()):
        name,metal,component,variant,_=spec;s=states[name];atoms=xyz(verify(t['xyz']));expected=geometry(s,metal,variant)
        if ([a[0] for a in atoms]!=[a[0] for a in expected] or not np.allclose([a[1:] for a in atoms],[a[1:] for a in expected],atol=1e-12,rtol=0)
                or not np.array_equal(np.load(verify(t['charges'])),s['endpoints'][metal][component+'_charges_e'])):
            raise InvalidArtifact('responsive solver coordinates/charges changed')
        wanted=task_payload(spec,s,m['states'][name],t['xyz'],t['charges'])
        if {k:v for k,v in t.items() if k!='cache_key'}!=wanted or key(t,m)!=t['cache_key']:raise InvalidArtifact('solver cache/state differs')
    return {'status':'pass','tasks':44,'manifest':record(manifest)}


def collect(manifest):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        good=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp);attempts.append({'task_id':t['task_id'],'directory':str(a),'accepted':r is not None,'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
            if r:good.append(r)
        rows[t['task_id']]=good[-1] if good else {'status':'unavailable'}
    complete=all(r['status']=='computed' for r in rows.values());checks=[];cases={};contrasts=[];partition=None
    def energy(name,metal,component='full',variant='primary'):return rows['_'.join((name,metal,component,variant))]['GB_reaction_kcal_mol']
    def check(name,error,tol):checks.append({'name':name,'error_kcal_scale':error,'tolerance':tol,'pass':abs(error)<=tol})
    if complete:
        states,c,cm,q,qm,f,fm,g,gm=sources(verify(m['charges']))
        for name in CASES:
            ends={}
            for metal in ('Ca','La'):
                base=energy(name,metal);tid=name+'_'+metal;old=f['cases'][name]['endpoints'][metal]
                check(tid+'_identity',energy(name,metal,variant='identity'),TOL['identity_kcal_mol'])
                for variant in ('rotate','translate')+(('repeat','reference') if name=='GGR_connected' else ()):
                    check(tid+'_'+variant,energy(name,metal,variant=variant)-base,TOL['numerical_kcal_mol'])
                env=g['cases'][name]['endpoints'][metal]['GB_env_kcal_mol'];selfpart=energy(name,metal,'QM')
                embedded=q['rows'][tid]['energy_hartree']*HA_TO_KCAL;direct=c['rows'][tid]['direct_coupling_kcal_mol']['exact']
                short=old['short_context_model_kcal'];hybrid=embedded+base+short
                ends[metal]={'embedded_DFT_kcal_mol':embedded,'GB_total_kcal_mol':base,'GB_QM_kcal_mol':selfpart,'GB_env_kcal_mol':env,
                    'GB_cross_kcal_mol':base-selfpart-env,'short_context_model_kcal':short,'hybrid_kcal_scale':hybrid,
                    'electronic_response_kcal_mol':q['rows'][tid]['electronic_response_kcal_mol'],'direct_embedded_density_kcal_mol':direct,
                    'intrinsic_core_kcal_mol':embedded-direct,'GB_change_from_frozen_kcal_mol':base-old['GB_kcal_mol'],
                    'change_from_frozen_hybrid_kcal_scale':hybrid-old['hybrid_kcal_scale']}
            components={k:ends['Ca'][k]-ends['La'][k] for k in ('embedded_DFT_kcal_mol','GB_total_kcal_mol','short_context_model_kcal')}
            r=sum(components.values());check(name+'_algebra',r-(ends['Ca']['hybrid_kcal_scale']-ends['La']['hybrid_kcal_scale']),TOL['algebra_kcal_scale'])
            check(name+'_GB_components',components['GB_total_kcal_mol']-sum(ends['Ca'][k]-ends['La'][k] for k in ('GB_QM_kcal_mol','GB_env_kcal_mol','GB_cross_kcal_mol')),TOL['algebra_kcal_scale'])
            for variant in ('rotate','translate')+(('repeat','reference') if name=='GGR_connected' else ()):
                check(name+'_'+variant+'_contrast',energy(name,'Ca',variant=variant)-energy(name,'La',variant=variant)-components['GB_total_kcal_mol'],TOL['numerical_kcal_mol'])
            response=ends['Ca']['electronic_response_kcal_mol']-ends['La']['electronic_response_kcal_mol'];gbchange=ends['Ca']['GB_change_from_frozen_kcal_mol']-ends['La']['GB_change_from_frozen_kcal_mol']
            check(name+'_response_update',r-f['cases'][name]['hybrid_R_kcal_scale']-response-gbchange,TOL['algebra_kcal_scale'])
            cases[name]={'endpoints':ends,'components_R':components,'hybrid_R_kcal_scale':r,'response_R_kcal_mol':response,
                'GB_update_R_kcal_mol':gbchange,'frozen_hybrid_R_kcal_scale':f['cases'][name]['hybrid_R_kcal_scale'],'evidence':f['cases'][name]['evidence']}
        partition=cases['GGR_connected']['hybrid_R_kcal_scale']-cases['GGR_extended']['hybrid_R_kcal_scale']
        for a in ('ALPHA_1F6S','ALPHA_6IP9'):
            for b in ('GGR_extended','GGR_connected'):
                delta=cases[a]['hybrid_R_kcal_scale']-cases[b]['hybrid_R_kcal_scale'];contrasts.append({'alpha':a,'GGR':b,'difference_kcal_scale':delta,'pass':delta>TOL['ordering_kcal_scale']})
    return {'protocol_id':PROTOCOL,'manifest':record(mp),'status':'complete' if complete else 'incomplete','rows':rows,'attempts':attempts,
        'checks':checks,'numerical_gate_pass':complete and all(c['pass'] for c in checks),'cases':cases,'contrasts':contrasts,
        'partition_shift_kcal_scale':partition,'partition_gate_pass':partition is not None and abs(partition)<=TOL['partition_kcal_scale'],
        'ordering_gate_pass':complete and all(c['pass'] for c in contrasts),'collection_implementation':record(__file__),**NULLS}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('charges','software','output'):a.add_argument('--'+k,required=True)
    for op in ('dry-run','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.charges,a.software,a.output)
    elif a.op=='dry-run':r=validate(a.manifest)
    else:r=collect(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','attempts','cases')},indent=2))
