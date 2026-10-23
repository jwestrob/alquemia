"""Ten consumed water cells; fixed engineering loss, no engine or file access."""
import math

GEOMETRIES = ('water_minus', 'water_central', 'water_plus', 'aquo8_central', 'aquo9_central')
EXCHANGE_GEOMETRIES = ('water_minus', 'water_central', 'water_plus', 'aquo9_central')
CELLS = tuple(q+'_'+metal for q in GEOMETRIES for metal in ('La','Dy'))
WEIGHTS = {'energy':.25, 'La':.25, 'Dy':.25, 'paired':.25}
RESIDUAL_COUNT = 413

def require(ok, message):
    if not ok:
        raise ValueError(message)

def finite_number(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)

def evaluate(data, model, body):
    rows = data['rows']
    require(isinstance(rows,dict) and isinstance(model,dict) and set(rows)==set(model)==set(CELLS), 'Exactly ten water cells required')
    ref = {}
    for q in GEOMETRIES:
        la,dy = rows[q+'_La'],rows[q+'_Dy']
        require(la['geometry']['positions_A']==dy['geometry']['positions_A'] and la['bodies']==dy['bodies'] and la['input']['bonds_zero_based']==dy['input']['bonds_zero_based'], 'Paired geometry/bodies/topology differ')
        atoms = 25 if q=='aquo8_central' else 28 if q=='aquo9_central' else 4
        waters = (atoms-1)//3
        for metal in ('La','Dy'):
            gid=q+'_'+metal; row=rows[gid]; g=row['geometry'];ref[gid]=row['reference']
            require(g['metal']==metal and g['charge']==3 and g['atoms']==atoms and len(g['positions_A'])==atoms, 'Water geometry/state/size differs')
            require(row['bodies']==body.bodies(row) and len(row['bodies'])==waters+1 and sum(b['category']=='water' for b in row['bodies'])==waters, 'Water body partition differs')
            for value in (ref[gid],model[gid]):
                require(isinstance(value,dict) and finite_number(value.get('energy')) and isinstance(value.get('forces'),(list,tuple)) and len(value['forces'])==atoms and all(isinstance(v,(list,tuple)) and len(v)==3 and all(finite_number(x) for x in v) for v in value['forces']), 'Missing or invalid energy/force value')
    def difference(values,q):
        return values[q+'_Dy']['energy']-values[q+'_La']['energy']
    exchange = {q:(difference(model,q)-difference(model,'aquo8_central'))-(difference(ref,q)-difference(ref,'aquo8_central')) for q in EXCHANGE_GEOMETRIES}
    work = {metal+'_'+sign:(model['water_'+sign+'_'+metal]['energy']-model['water_central_'+metal]['energy'])-(ref['water_'+sign+'_'+metal]['energy']-ref['water_central_'+metal]['energy']) for metal in ('La','Dy') for sign in ('minus','plus')}
    paired_work = {sign:work['Dy_'+sign]-work['La_'+sign] for sign in ('minus','plus')}
    vectors = {'energy':[v*math.sqrt(.5/4) for v in exchange.values()]+[v*math.sqrt(.5/4) for v in work.values()]}
    per,controls,atomic = {},{},{}
    for kind in ('La','Dy','paired'):
        vector=[];raw=[];per[kind]={}
        for q in GEOMETRIES:
            row=rows[q+'_La']
            if kind=='paired':
                mf=body.subtract(model[q+'_Dy']['forces'],model[q+'_La']['forces'])
                rf=body.subtract(ref[q+'_Dy']['forces'],ref[q+'_La']['forces'])
            else:
                mf=model[q+'_'+kind]['forces'];rf=ref[q+'_'+kind]['forces']
            projected,components=body.rigid_residual(row,mf,rf)
            vector.extend(v/math.sqrt(5) for v in projected);raw.extend(components)
            rms=2*math.sqrt(sum(v*v for v in projected));maximum=max(abs(v) for v in components)
            per[kind][q]={'rigid_RMS_kcal_mol_A':rms,'maximum_body_component_error_kcal_mol_A':maximum,'passed':rms<=2 and maximum<=5,'raw_body_components_kcal_mol_A':components,'model_body_response':body.response(row,mf),'reference_body_response':body.response(row,rf)}
            residual=body.subtract(mf,rf)
            atomic[q+'_'+kind]={'atomic_force_residual_kcal_mol_A':residual,'atom_vector_RMS_kcal_mol_A':math.sqrt(sum(x*x for v in residual for x in v)/len(residual)),'used_in_loss':False}
        vectors[kind]=vector
        controls[kind]={'category_balanced_rigid_RMS_kcal_mol_A':2*math.sqrt(sum(v*v for v in vector)),'max_body_component_error_kcal_mol_A':max(abs(v) for v in raw),'passed':all(v['passed'] for v in per[kind].values()),'every_fixture_required':True,'fixture_count':5}
    weighted=[math.sqrt(WEIGHTS[k])*v for k in WEIGHTS for v in vectors[k]]
    require(len(weighted)==RESIDUAL_COUNT and all(math.isfinite(v) for v in weighted), 'Residual count or finiteness differs')
    energies={'balanced_exchange':exchange,'individual_radial_work':work,'paired_radial_work':paired_work}
    energy_controls={group:{key:{'error_kcal_mol':value,'passed':abs(value)<=1} for key,value in vals.items()} for group,vals in energies.items()}
    energy_pass=all(v['passed'] for vals in energy_controls.values() for v in vals.values())
    old_energy_pass=all(abs(v)<=1 for q,v in exchange.items() if q!='aquo9_central') and all(abs(v)<=1 for vals in (work,paired_work) for v in vals.values())
    old_force_pass=all(per[kind][q]['passed'] for kind in per for q in GEOMETRIES[:4])
    w9_force_pass=all(per[kind]['aquo9_central']['passed'] for kind in per)
    return {'objective':sum(v*v for v in weighted),'weighted_residuals':weighted,'blocks':{k:sum(v*v for v in vectors[k]) for k in WEIGHTS},'block_weights':dict(WEIGHTS),'energy_subblock_weights':{'exchange':.5,'individual_radial_work':.5},'energy_errors_kcal_mol':energies,'energy_controls':energy_controls,'energy_checks_passed':energy_pass,'controls':controls,'per_fixture':per,'atomic_diagnostics':atomic,'development_checks_passed':energy_pass and all(v['passed'] for v in controls.values()),'declared_energy_controls':10,'available_energy_controls':10,'declared_force_controls':15,'available_force_controls':15,'residual_layout':{'energy':[0,8],'La':[8,143],'Dy':[143,278],'paired':[278,413],'geometry_order':list(GEOMETRIES),'exchange_order':list(EXCHANGE_GEOMETRIES),'individual_radial_order':list(work),'body_order':'Unchanged helper: category ion then water; body order; net xyz then centroid torque/1A xyz for waters.'},'development_partitions':{'old_eight':{'energy_controls':9,'force_controls':12,'energy_passed':old_energy_pass,'force_passed':old_force_pass,'passed':old_energy_pass and old_force_pass},'added_W9':{'energy_controls':1,'force_controls':3,'energy_passed':abs(exchange['aquo9_central'])<=1,'force_passed':w9_force_pass,'passed':abs(exchange['aquo9_central'])<=1 and w9_force_pass}},'residual_interpretation':'Correlated engineering residuals, not independent observations or an uncertainty likelihood.','physical_accuracy_qualified':False,'production_released':False}
