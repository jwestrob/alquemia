"""Eight consumed water cells; symmetric energy/rigid response, no engines."""
import math

GEOMETRIES=('water_minus','water_central','water_plus','aquo8_central')
CELLS=tuple(q+'_'+m for q in GEOMETRIES for m in ('La','Dy'))
WEIGHTS={'energy':.25,'La':.25,'Dy':.25,'paired':.25}


def require(ok,message):
    if not ok:raise ValueError(message)


def evaluate(data,model,body):
    rows=data['rows'];require(set(rows)==set(model)==set(CELLS),'Exactly eight water cells required')
    ref={k:r['reference'] for k,r in rows.items()}
    for q in GEOMETRIES:
        a,b=rows[q+'_La'],rows[q+'_Dy']
        require(a['geometry']['positions_A']==b['geometry']['positions_A'] and a['bodies']==b['bodies'] and
                a['input']['bonds_zero_based']==b['input']['bonds_zero_based'],'Paired coordinates/bodies differ')
        for metal in ('La','Dy'):
            k=q+'_'+metal;r=rows[k];g=r['geometry']
            require(g['metal']==metal and g['charge']==3 and g['atoms']==(25 if q=='aquo8_central' else 4),'Water state/atoms differ')
            require(all(x['category'] in ('ion','water') for x in r['bodies']),'Nonwater body')
            for source in (ref,model):
                require(math.isfinite(source[k]['energy']) and len(source[k]['forces'])==g['atoms'] and
                    all(len(v)==3 and all(math.isfinite(x) for x in v) for v in source[k]['forces']),'Invalid E/F')
    def D(values,q):return values[q+'_Dy']['energy']-values[q+'_La']['energy']
    exchange={q:(D(model,q)-D(model,'aquo8_central'))-(D(ref,q)-D(ref,'aquo8_central')) for q in GEOMETRIES[:3]}
    work={m+'_'+s:(model['water_'+s+'_'+m]['energy']-model['water_central_'+m]['energy'])-
        (ref['water_'+s+'_'+m]['energy']-ref['water_central_'+m]['energy']) for m in ('La','Dy') for s in ('minus','plus')}
    paired_work={s:work['Dy_'+s]-work['La_'+s] for s in ('minus','plus')}
    vectors={'energy':[v*math.sqrt(.5/3) for v in exchange.values()]+[v*math.sqrt(.5/4) for v in work.values()]}
    controls={};per={};atomic={}
    for kind in ('La','Dy','paired'):
        vec=[];raw=[];per[kind]={}
        for q in GEOMETRIES:
            row=rows[q+'_La']
            if kind=='paired':
                mf=body.subtract(model[q+'_Dy']['forces'],model[q+'_La']['forces']);rf=body.subtract(ref[q+'_Dy']['forces'],ref[q+'_La']['forces'])
            else:mf=model[q+'_'+kind]['forces'];rf=ref[q+'_'+kind]['forces']
            weighted,components=body.rigid_residual(row,mf,rf);vec.extend(v/math.sqrt(4) for v in weighted);raw.extend(components)
            rms=2*math.sqrt(sum(v*v for v in weighted));mx=max(abs(v) for v in components)
            per[kind][q]={'rigid_RMS_kcal_mol_A':rms,'maximum_body_component_error_kcal_mol_A':mx,'passed':rms<=2 and mx<=5,
                'raw_body_components_kcal_mol_A':components,'model_body_response':body.response(row,mf),'reference_body_response':body.response(row,rf)}
            residual=body.subtract(mf,rf);atomic[q+'_'+kind]={'atomic_force_residual_kcal_mol_A':residual,
                'atom_vector_RMS_kcal_mol_A':math.sqrt(sum(x*x for v in residual for x in v)/len(residual)),'used_in_loss':False}
        vectors[kind]=vec
        controls[kind]={'category_balanced_rigid_RMS_kcal_mol_A':2*math.sqrt(sum(v*v for v in vec)),
            'max_body_component_error_kcal_mol_A':max(abs(v) for v in raw),'passed':all(r['passed'] for r in per[kind].values()),
            'every_fixture_required':True}
    weighted=[math.sqrt(WEIGHTS[k])*v for k in WEIGHTS for v in vectors[k]]
    energies={'balanced_exchange':exchange,'individual_radial_work':work,'paired_radial_work':paired_work}
    energy_pass=all(abs(x)<=1 for values in energies.values() for x in values.values())
    return {'objective':sum(v*v for v in weighted),'weighted_residuals':weighted,'blocks':{k:sum(v*v for v in vectors[k]) for k in WEIGHTS},
        'block_weights':WEIGHTS,'energy_subblock_weights':{'exchange':.5,'individual_radial_work':.5},'energy_errors_kcal_mol':energies,
        'energy_checks_passed':energy_pass,'controls':controls,'per_fixture':per,'atomic_diagnostics':atomic,
        'development_checks_passed':energy_pass and all(v['passed'] for v in controls.values()),'physical_accuracy_qualified':False,'production_released':False}
