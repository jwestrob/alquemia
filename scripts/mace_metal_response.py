"""DFT-anchored complete metal-translation response using existing MACE workers."""
from __future__ import annotations
import argparse
import copy
import itertools
import json
from pathlib import Path
import numpy as np
from affordable_common import (InvalidArtifact, HA_TO_KCAL, BOHR_TO_A, cache_key,
                               read_json, record, verify, write_new, xyz)
from mace_hybrid import check_atoms, write_xyz, accepted_attempt, EV_TO_KCAL
from mace_global_benchmark import snapshot, numerical_parent_gate
from mace_gb import MODEL as GB_MODEL

POLICY='DFT_anchored_MACE_GB_full_metal_translation_v1'
SCHEMAS={k:f'alquemia.mace_metal_response_{k}.v1' for k in ('core','short','gb')}
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')
TOL={'total_charge_e':1e-5,'energy_kcal_mol':.01,'force_max_eV_A':.001}
LIMIT=.20
UNAVAILABLE={'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,
             'entropy_correction_kcal_mol':None,'calibrated_class':None}


def grid():
    points={'center':[0.,0.,0.]}
    for label,h in [('coarse',.02),('fine',.01)]:
        for axes in [(0,),(1,),(2,),(0,1),(0,2),(1,2)]:
            for signs in itertools.product((-1,1),repeat=len(axes)):
                v=np.zeros(3)
                for i,s in zip(axes,signs):v[i]=s*h
                name=label+'_'+''.join('xyz'[i]+('p' if s>0 else 'm') for i,s in zip(axes,signs))
                points[name]=v.tolist()
    return points


def shift(rows,metal,u):
    indices=[i for i,a in enumerate(rows) if a[0]==metal]
    if len(indices)!=1:raise InvalidArtifact('exactly one selected metal required')
    i=indices[0];result=list(rows);result[i]=(metal,*(np.array(rows[i][1:])+u))
    return result,i


def rows_with_reuse(collection):
    c=read_json(collection);mp=verify(c['manifest']);m=read_json(mp)
    if c['status']!='complete':raise InvalidArtifact('incomplete parent collection')
    result=dict(c['rows'])
    for key,item in c.get('reused_rows',{}).items():result[key]=item['result']
    for r in result.values():
        if r['status']!='computed':raise InvalidArtifact('failed parent endpoint')
    return result,m


def prepare(assessment,agreement,output):
    from mace_mechanics import physical_model
    from ggr_sensitivity import membership
    from affordable_response import read_engrad
    a=read_json(assessment);p=read_json(verify(a['sources']['prepared']))
    cores,cm=rows_with_reuse(verify(a['sources']['core']))
    gbs,bm=rows_with_reuse(verify(a['sources']['gb']))
    shorts,sm=rows_with_reuse(verify(a['sources']['short']))
    gate=read_json(verify(sm['short_gate']))
    if not gate['numerical_checks_pass']:raise InvalidArtifact('short force prerequisite failed')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    source_copy=out/'preparation_source.py';source_copy.write_bytes(Path(__file__).read_bytes())
    mappings={};physical={};center={};jobs=[]
    for name in CASES:
        c=read_json(verify(p['cases'][name]));repair=read_json(verify(c['source_preparation']))
        graph,coords,pos,mk,_,_=physical_model(repair);inventory=membership(graph,coords,pos,mk)
        directory=out/name;directory.mkdir();states={};centers={}
        for metal in ('La','Ca'):
            key=f'{name}_center_{metal}';src=c['grids']['center']['endpoints'][metal]
            pin=a['gradient_sources'][key];gradient=read_engrad(verify(pin['engrad']))
            for ref in pin.values():verify(ref)
            rows=xyz(verify(src['xyz']))
            if rows!=xyz(verify(pin['xyz'])) or rows[0][0]!=metal:
                raise InvalidArtifact('DFT metal ordering/coordinates do not match')
            fullkey=f'full__{c["global_id"]}_center_{metal}'
            core_short=(shorts['core__'+key] if name.startswith('ALPHA') else
                        gate['rows']['core__'+key.removeprefix('GGR_')])
            fullsrc=p['physical_cases'][c['global_id']]['grids']['center']['endpoints'][metal]
            fullrows=xyz(verify(fullsrc['xyz']));_,fi=shift(fullrows,metal,[0,0,0])
            gj=(-np.load(verify(shorts[fullkey]['forces']))[fi]+
                 np.load(verify(core_short['forces']))[0])*EV_TO_KCAL
            gd=np.array(gradient['gradient_Ha_per_bohr'][0])*HA_TO_KCAL/BOHR_TO_A
            if len(rows)!=gradient['atom_count']:raise InvalidArtifact('gradient atom count mismatch')
            centers[metal]={'core':cores[key],'gb':gbs[key],'full_short':shorts[fullkey],
                            'core_short':core_short,'DFT_artifacts':pin,
                            'DFT_energy_Ha':gradient['energy_Ha'],'DFT_gradient_kcal_mol_A':gd.tolist(),
                            'J_gradient_kcal_mol_A':gj.tolist(),'core_state':src,'full_state':fullsrc,
                            'full_metal_index':fi,'old_axis_projection_fraction':abs(a['rows'][name]['gradient_components'][metal]['DFT_core'][0])/float(np.linalg.norm(gd))}
            for point,u in grid().items():
                if membership(graph,coords,pos+np.array(u),mk)!=inventory:
                    raise InvalidArtifact('typed donor membership changed at '+name+'/'+point)
                moved,idx=shift(rows,metal,u)
                if idx!=0:raise InvalidArtifact('QM metal must remain first')
                xp=directory/f'{point}_{metal}.xyz'
                if point=='center':xp.write_bytes(verify(src['xyz']).read_bytes())
                else:write_xyz(xp,moved)
                states.setdefault(point,{})[metal]={'xyz':record(xp),'charge':src['charge'],
                    'spin_multiplicity':1,'state':check_atoms(moved,src['charge'])}
        mappings[name]={'source_mapping':p['cases'][name],'global_id':c['global_id'],
                         'source_preparation':c['source_preparation'],'states':states}
        center[name]=centers
        gid=c['global_id']
        if gid not in physical:
            fd=out/gid if gid not in CASES else out/(gid+'_physical');fd.mkdir(exist_ok=True)
            physical[gid]={}
            for point,u in grid().items():
                for metal in ('La','Ca'):
                    src=centers[metal]['full_state'];rows=xyz(verify(src['xyz']));moved,_=shift(rows,metal,u)
                    xp=fd/f'full_{point}_{metal}.xyz'
                    if point=='center':xp.write_bytes(verify(src['xyz']).read_bytes())
                    else:write_xyz(xp,moved)
                    physical[gid].setdefault(point,{})[metal]={'xyz':record(xp),'charge':src['charge'],
                         'spin_multiplicity':1,'state':check_atoms(moved,src['charge'])}
    top={'policy_id':POLICY,'assessment':record(assessment),'agreement':record(agreement),
         'implementation':record(source_copy),'cases':mappings,'physical':physical,'centers':center,'grid':grid(),
         'new_calls':{'core':288,'gb':288,'short':216,'conditional_DFT_max':8,'conditional_short_max':16},**UNAVAILABLE}
    write_new(out/'preparation.json',top)
    for kind,parent,groups in [('core',cm,{k:v['states'] for k,v in mappings.items()}),('short',sm,physical)]:
        for name,states in groups.items():
            tasks=[]
            for point,ends in states.items():
                if point=='center':continue
                for metal,e in ends.items():
                    t={'task_id':f'{name}_{point}_{metal}','case_id':name,'point':point,'metal':metal,
                       'kind':'core' if kind=='core' else 'full','variant':'metal_translation',
                       'displacement_A':grid()[point],**e}
                    if kind=='short':t['energy_component']='interaction_energy'
                    tasks.append(t)
            job=out/(kind+'_'+name);job.mkdir()
            pins=snapshot(job,(*parent['implementation'],'mace_metal_response.py'))
            model=copy.deepcopy(parent['model']);model['preparation_policy']=POLICY
            m={'schema_version':SCHEMAS[kind],'protocol_id':POLICY+'_'+kind,'kind':kind,
               'preparation':record(out/'preparation.json'),'agreement':record(agreement),
               'parent_collection':a['sources']['core' if kind=='core' else 'short'],
               'short_gate':sm['short_gate'],'numerical_reference':parent['numerical_reference'],
               'model':model,'software':parent['software'],'implementation':pins,'tasks':tasks,'tolerances':TOL,**UNAVAILABLE}
            seal(job,m);validate(job/'manifest.json');jobs.append(record(job/'manifest.json'))
    write_new(out/'jobs.json',{'preparation':record(out/'preparation.json'),'jobs':jobs})
    return {'status':'prepared','jobs':jobs,'new_calls':top['new_calls']}


def seal(out,m):
    for t in m['tasks']:
        t.pop('cache_key',None)
        t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
    write_new(out/'manifest.json',m)


def validate(manifest):
    m=read_json(manifest);kind=m['kind'];p=read_json(verify(m['preparation']))
    if (m['schema_version']!=SCHEMAS[kind] or p['policy_id']!=POLICY or m['tolerances']!=TOL or
            p['grid']!=grid() or numerical_parent_gate(m)['status']!='pass'):
        raise InvalidArtifact('response protocol/grid/prerequisite changed')
    for ref in [m['agreement'],*m['implementation'].values(),p['assessment'],p['implementation']]:verify(ref)
    software=read_json(verify(m['software']))
    for ref in [software['python'],software['requirements'],*read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    _,parent=rows_with_reuse(verify(m['parent_collection']))
    expected=copy.deepcopy(parent['model']);expected['preparation_policy']=POLICY
    if kind=='gb':
        source=read_json(verify(m['source_mace_collection']));source_m=read_json(verify(source['manifest']))
        gate=read_json(verify(m['solver_validation']));gm=read_json(verify(gate['manifest']))
        if not gate['numerical_checks_pass'] or m['model']!=GB_MODEL or m['software']!=gm['software'] or source_m['preparation']!=m['preparation']:
            raise InvalidArtifact('GB model/source mismatch')
    elif m['model']!=expected or m['software']!=parent['software']:
        raise InvalidArtifact('scientific method differs from parent')
    names={t['case_id'] for t in m['tasks']}
    if len(names)!=1 or len(m['tasks'])!=72 or len({t['task_id'] for t in m['tasks']})!=72:
        raise InvalidArtifact('finite task inventory changed')
    name=next(iter(names));states=p['physical'][name] if kind=='short' else p['cases'][name]['states']
    expected_ids={f'{name}_{point}_{metal}' for point in grid() if point!='center' for metal in ('La','Ca')}
    if {t['task_id'] for t in m['tasks']}!=expected_ids:raise InvalidArtifact('task identity changed')
    for t in m['tasks']:
        e=states[t['point']][t['metal']]
        if any(t.get(k)!=v for k,v in e.items()) or t['displacement_A']!=grid()[t['point']]:
            raise InvalidArtifact('task state changed')
        rows=xyz(verify(t['xyz']));base=xyz(verify(states['center'][t['metal']]['xyz']))
        moved,_=shift(base,t['metal'],t['displacement_A'])
        if [r[0] for r in rows]!=[r[0] for r in moved] or not np.allclose([r[1:] for r in rows],[r[1:] for r in moved],atol=5e-10,rtol=0):
            raise InvalidArtifact('nonmetal coordinates/order or displacement changed')
        check_atoms(rows,t['charge'])
        if kind=='gb':
            r=source['rows'][t['task_id']]
            if t['source_density']!=r['density_coefficients'] or t['source_vacuum_energy_eV']!=r['energy_eV'] or t['solver']!='native' or t['platform']!='CUDA' or t['solvent_dielectric']!=78.5:
                raise InvalidArtifact('GB density or solvent changed')
            verify(t['source_density'])
        base={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':base,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('cache key changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def collect(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:
                if m['kind']=='gb':verify(r['serialized_system'])
                valid.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=valid[-1] if valid else {'status':'unavailable','energy_eV':None}
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete',
            'manifest':record(mp),'protocol_id':m['protocol_id'],'rows':rows,'attempts':attempts,**UNAVAILABLE}


def prepare_gb(collection,solver_validation,output):
    c=read_json(collection);m=read_json(verify(c['manifest']));validate(verify(c['manifest']))
    actual=collect(verify(c['manifest']))
    if c['status']!='complete' or actual['rows']!=c['rows'] or m['kind']!='core':
        raise InvalidArtifact('GB requires actual complete core results')
    gate=read_json(solver_validation);gm=read_json(verify(gate['manifest']))
    if not gate['numerical_checks_pass']:raise InvalidArtifact('GB prerequisite failed')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,m['implementation']);result=copy.deepcopy(m)
    for t in result['tasks']:
        t.pop('cache_key');r=c['rows'][t['task_id']]
        t.update(solver='native',platform='CUDA',solvent_dielectric=78.5,
                 source_density=r['density_coefficients'],source_vacuum_energy_eV=r['energy_eV'])
    result.update(kind='gb',schema_version=SCHEMAS['gb'],protocol_id=POLICY+'_gb',model=GB_MODEL,
                  software=gm['software'],implementation=pins,source_mace_collection=record(collection),
                  solver_validation=record(solver_validation))
    seal(out,result);return validate(out/'manifest.json')


def matrix(values,scale):
    h={'coarse':.02,'fine':.01}[scale];k=np.zeros((3,3));e0=values['center']
    def e(axes,signs):return values[scale+'_'+''.join('xyz'[i]+('p' if s>0 else 'm') for i,s in zip(axes,signs))]-e0
    for i in range(3):k[i,i]=(e((i,),(1,))+e((i,),(-1,)))/h**2
    for i,j in itertools.combinations(range(3),2):
        k[i,j]=k[j,i]=sum(s*t*e((i,j),(s,t)) for s,t in itertools.product((-1,1),repeat=2))/(4*h*h)
    if not np.isfinite(k).all():raise InvalidArtifact('nonfinite curvature')
    return k


def prediction(g,components):
    matrices={key:{s:matrix(v,s) for s in ('coarse','fine')} for key,v in components.items()}
    fine=sum(v['fine'] for v in matrices.values());coarse=sum(v['coarse'] for v in matrices.values())
    eig={s:np.linalg.eigvalsh(k).tolist() for s,k in [('fine',fine),('coarse',coarse)]}
    result={'matrices_kcal_mol_A2':{key:{s:k.tolist() for s,k in value.items()} for key,value in matrices.items()},
            'eigenvalues_kcal_mol_A2':eig,'gradient_kcal_mol_A':list(g),'status':'unstable_or_unresolved_curvature',
            'predicted_response_kcal_mol':None,'proposed_displacement_A':None}
    if min(eig['fine'])<=0 or min(eig['coarse'])<=0:return result
    u=-np.linalg.solve(fine,g);uc=-np.linalg.solve(coarse,g)
    errors={key:float(.5*u@(v['fine']-v['coarse'])@u) for key,v in matrices.items()}
    errors['combined']=float(.5*u@(fine-coarse)@u)
    residual=float(np.linalg.norm(fine@u+g)/max(np.linalg.norm(g),np.finfo(float).tiny))
    result.update(proposed_displacement_A=u.tolist(),displacement_norm_A=float(np.linalg.norm(u)),
                  condition_number=float(np.linalg.cond(fine)),solve_relative_residual=residual,
                  refinement_energy_change_kcal_mol=errors,refinement_displacement_change_A=float(np.linalg.norm(u-uc)),
                  exploratory_quadratic_change_kcal_mol=float(.5*np.array(g)@u))
    if residual>1e-10:result['status']='unreliable_linear_solve'
    elif max(map(abs,errors.values()))>.02 or np.linalg.norm(u-uc)>.01:result['status']='curvature_not_converged'
    elif np.linalg.norm(u)>LIMIT:result['status']='trust_region_exceeded'
    else:result.update(status='eligible_for_DFT_validation',predicted_response_kcal_mol=float(.5*np.array(g)@u))
    return result


def assess(prepared,collections,output):
    p=read_json(prepared);index=read_json(collections);data={k:{} for k in SCHEMAS};sources=[]
    for path in index['collections']:
        c=read_json(path);mp=verify(c['manifest']);m=read_json(mp);validate(mp)
        if c['status']!='complete' or collect(mp)['rows']!=c['rows'] or m['preparation']!=record(prepared):
            raise InvalidArtifact('incomplete/mismatched actual grid collection')
        data[m['kind']].update(c['rows']);sources.append(record(path))
    if {k:len(v) for k,v in data.items()}!={'core':288,'gb':288,'short':216}:
        raise InvalidArtifact('complete declared grid inventory required')
    rows={};scores={}
    for name,c in p['cases'].items():
        rows[name]={}
        for metal in ('La','Ca'):
            center=p['centers'][name][metal];values={'C':{'center':0.},'J':{'center':0.}}
            for point in grid():
                if point=='center':continue
                key=f'{name}_{point}_{metal}';fk=f'{c["global_id"]}_{point}_{metal}'
                core=data['core'][key];gb=data['gb'][key];full=data['short'][fk]
                values['C'][point]=(core['energy_eV']-center['core']['energy_eV'])*EV_TO_KCAL+(gb['GB_reaction_kcal_mol']-center['gb']['GB_reaction_kcal_mol'])
                values['J'][point]=((full['energy_eV']-center['full_short']['energy_eV'])-(core['energy_components_eV']['interaction_energy']-center['core']['energy_components_eV']['interaction_energy']))*EV_TO_KCAL
            gd=np.array(center['DFT_gradient_kcal_mol_A']);gj=np.array(center['J_gradient_kcal_mol_A'])
            rows[name][metal]={'prediction':prediction(gd+gj,values),'relative_grid_energies_kcal_mol':values,
                                'DFT_gradient_kcal_mol_A':gd.tolist(),'J_gradient_kcal_mol_A':gj.tolist()}
        r=(p['centers'][name]['Ca']['DFT_energy_Ha']-p['centers'][name]['La']['DFT_energy_Ha'])*HA_TO_KCAL
        deltas={m:rows[name][m]['prediction']['predicted_response_kcal_mol'] for m in ('La','Ca')}
        delta=None if None in deltas.values() else deltas['Ca']-deltas['La']
        scores[name]={'baseline_R_kcal_mol':r,'unvalidated_predicted_delta_R_kcal_mol':delta,
                     'unvalidated_predicted_R_kcal_mol':None if delta is None else r+delta,**UNAVAILABLE}
    result={'status':'complete','policy_id':POLICY,'preparation':record(prepared),'collections':sources,
            'implementation':record(__file__),'rows':rows,'scores':scores,'predictions_are_not_validated_scores':True}
    out=Path(output);out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    return {'status':'complete','endpoint_status':{k:{m:v['prediction']['status'] for m,v in row.items()} for k,row in rows.items()},'scores':scores}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for command,keys in [('prepare',('assessment','agreement','output')),('prepare-gb',('collection','solver-validation','output')),
                         ('assess',('prepared','collections','output')),('validate',('manifest',)),('collect',('manifest','output'))]:
        q=sub.add_parser(command)
        for key in keys:q.add_argument('--'+key,required=True)
    args=vars(parser.parse_args());command=args.pop('command')
    if command=='collect':
        output=args.pop('output');result=collect(**args);write_new(output,result)
        result={'status':result['status'],'completed':sum(r['status']=='computed' for r in result['rows'].values())}
    else:result={'prepare':prepare,'prepare-gb':prepare_gb,'assess':assess,'validate':validate}[command](**args)
    print(json.dumps(result,indent=2))
