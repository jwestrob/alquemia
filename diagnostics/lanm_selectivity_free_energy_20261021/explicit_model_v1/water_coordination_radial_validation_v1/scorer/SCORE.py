"""Four fresh W9 radial cells; twelve fixed checks, no fitting or engines."""
import math

CELLS=tuple('aquo9_'+s+'_'+m for s in ('minus','plus') for m in ('La','Dy'))

def require(ok,message):
    if not ok:raise ValueError(message)

def evaluate(data,model,body,anchors):
    rows=data['rows'];require(set(rows)==set(CELLS) and set(model)==set(CELLS),'Exactly four declared validation slots required')
    refs={g:rows[g].get('reference') for g in CELLS}
    def available(v,n):
        if v is None or v.get('status') in ('unavailable','missing','failed'):return False
        require(isinstance(v.get('energy'),(float,int)) and not isinstance(v['energy'],bool) and math.isfinite(v['energy']),'Invalid available energy')
        require(len(v.get('forces',[]))==n and all(len(a)==3 and all(isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) for x in a) for a in v['forces']),'Invalid available forces')
        return True
    for gid in CELLS:
        r=rows[gid];g=r['geometry'];m=gid.rsplit('_',1)[1]
        require(g['atoms']==28 and g['charge']==3 and g['metal']==m and g['water_count']==9,'W9 validation chemistry differs')
        require(all(b['category'] in ('ion','water') for b in r['bodies']),'Nonwater body present')
    for sign in ('minus','plus'):
        a,b=(rows['aquo9_'+sign+'_'+m] for m in ('La','Dy'))
        require(a['geometry']['positions_A']==b['geometry']['positions_A'] and a['bodies']==b['bodies'] and a['input']['bonds_zero_based']==b['input']['bonds_zero_based'],'Paired geometry/body mismatch')
    valid={g:available(model[g],28) and available(refs[g],28) for g in CELLS}
    am=anchors.get('model',{});ar=anchors.get('reference',{})
    anchor_ok={m:available(am.get(m),28) and available(ar.get(m),28) for m in ('La','Dy')}
    unavailable=lambda missing:{'status':'unavailable','passed':None,'missing':missing}
    individual={};paired={};rigid={};count=0
    for sign in ('minus','plus'):
        for m in ('La','Dy'):
            gid='aquo9_'+sign+'_'+m
            if not valid[gid] or not anchor_ok[m]:individual[gid]=unavailable(([gid] if not valid[gid] else [])+(['anchor_'+m] if not anchor_ok[m] else []));continue
            mv=model[gid]['energy']-am[m]['energy'];rv=refs[gid]['energy']-ar[m]['energy'];err=mv-rv
            individual[gid]={'status':'complete','model_kcal_mol':mv,'reference_kcal_mol':rv,'error_kcal_mol':err,'passed':abs(err)<=1};count+=1
        a,b=('aquo9_'+sign+'_'+m for m in ('La','Dy'))
        missing=[g for g in (a,b) if individual[g]['status']!='complete']
        if missing:paired[sign]=unavailable(missing)
        else:
            err=individual[b]['error_kcal_mol']-individual[a]['error_kcal_mol']
            paired[sign]={'status':'complete','model_kcal_mol':individual[b]['model_kcal_mol']-individual[a]['model_kcal_mol'],
                'reference_kcal_mol':individual[b]['reference_kcal_mol']-individual[a]['reference_kcal_mol'],'error_kcal_mol':err,'passed':abs(err)<=1};count+=1
        rigid[sign]={}
        for kind in ('La','Dy','paired'):
            needed=(a,b) if kind=='paired' else ('aquo9_'+sign+'_'+kind,);missing=[g for g in needed if not valid[g]]
            if missing:rigid[sign][kind]=unavailable(missing);continue
            if kind=='paired':mf=body.subtract(model[b]['forces'],model[a]['forces']);rf=body.subtract(refs[b]['forces'],refs[a]['forces'])
            else:mf=model[needed[0]]['forces'];rf=refs[needed[0]]['forces']
            weighted,components=body.rigid_residual(rows[a],mf,rf);rms=2*math.sqrt(sum(v*v for v in weighted));mx=max(abs(v) for v in components)
            rigid[sign][kind]={'status':'complete','RMS_kcal_mol_A':rms,'maximum_component_kcal_mol_A':mx,'passed':rms<=2 and mx<=5,
                'model_body_response':body.response(rows[a],mf),'reference_body_response':body.response(rows[a],rf)};count+=1
    checks=list(individual.values())+list(paired.values())+[v for bykind in rigid.values() for v in bykind.values()]
    require(len(checks)==12,'Declared denominator changed')
    return {'individual_radial_work':individual,'paired_radial_work':paired,'rigid':rigid,'declared_checks':12,'available_checks':count,
        'passed_checks':sum(x['passed'] is True for x in checks),'all_validation_checks_passed':count==12 and all(x['passed'] is True for x in checks),
        'model_cells_available':sum(available(model[g],28) for g in CELLS),'reference_cells_available':sum(available(refs[g],28) for g in CELLS),
        'energy_tolerance_kcal_mol':1,'rigid_RMS_tolerance_kcal_mol_A':2,'maximum_component_tolerance_kcal_mol_A':5,
        'physical_accuracy_qualified':False,'automatic_fitting_release':False,'production_released':False}
