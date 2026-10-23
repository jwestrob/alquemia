"""Frozen per-cell consumed transfer comparisons; no engines or fitting."""
import math

def evaluate(data,model,body):
    rows=data['rows'];assert set(model)==set(rows) and len(rows)==32
    original={k:r['baseline']for k,r in rows.items()};ref={k:r['reference']for k,r in rows.items()}
    values={'baseline':original,'candidate':model,'reference':ref}
    for source in values.values():
        for k,v in source.items():
            assert math.isfinite(v['energy']) and len(v['forces'])==rows[k]['geometry']['atoms']
            assert all(len(a)==3 and all(math.isfinite(x)for x in a)for a in v['forces'])
    checks=[];works=[];exchange=[];addition=[];rigid=[];paths=[]
    def compare(v,group,name):
        out={'values_kcal_mol':v,'errors':{},'group':group,'id':name}
        for method in ('baseline','candidate'):
            error=v[method]-v['reference'];out['errors'][method]={'signed_kcal_mol':error,'passed':abs(error)<=1}
        out['absolute_error_change_kcal_mol']=abs(out['errors']['candidate']['signed_kcal_mol'])-abs(out['errors']['baseline']['signed_kcal_mol'])
        out['deteriorated']=out['absolute_error_change_kcal_mol']>0
        checks.append({'group':group,'id':name,'kind':'energy','baseline_passed':out['errors']['baseline']['passed'],'candidate_passed':out['errors']['candidate']['passed']})
        return out
    def delta(src,q):return src[q+'_Dy']['energy']-src[q+'_La']['energy']
    families=('water','hydrated_acetate_mono','hydrated_acetate_bi','hydrated_acetamide')
    for f in families:
        group='water_radial'if f=='water'else'hydrated_donor'
        for side in ('minus','plus'):
            per={}
            for metal in ('La','Dy'):
                v={m:s[f+'_'+side+'_'+metal]['energy']-s[f+'_central_'+metal]['energy']for m,s in values.items()};per[metal]=v
                works.append(compare(v,group,f+'_'+side+'_'+metal))
            works.append(compare({m:per['Dy'][m]-per['La'][m]for m in values},group,f+'_'+side+'_Dy_minus_La'))
        for metal in ('La','Dy'):
            path=data['paths'][f+'_'+metal];u=path['direction'];ids=path['moving_atoms'];record={'family':f,'metal':metal,'moving_atoms':ids,'direction':u,'total_endpoint_separation_A':.30,'small_step_qualification':False,'models':{}}
            for method,src in values.items():
                slopes=[-sum(sum(src[f+'_'+side+'_'+metal]['forces'][i][k]*u[k]for k in range(3))for i in ids)for side in ('minus','central','plus')]
                secant=(src[f+'_plus_'+metal]['energy']-src[f+'_minus_'+metal]['energy'])/.30
                record['models'][method]={'slopes_minus_central_plus_kcal_mol_A':slopes,'secant_kcal_mol_A':secant,'central_minus_secant_kcal_mol_A':slopes[1]-secant}
            paths.append(record)
        pair={'family':f,'metal':'Dy_minus_La','small_step_qualification':False,'models':{}}
        a,b=paths[-2:]
        for m in values:
            am,bm=a['models'][m],b['models'][m];pair['models'][m]={'slopes_minus_central_plus_kcal_mol_A':[y-x for x,y in zip(am['slopes_minus_central_plus_kcal_mol_A'],bm['slopes_minus_central_plus_kcal_mol_A'])],'secant_kcal_mol_A':bm['secant_kcal_mol_A']-am['secant_kcal_mol_A'],'central_minus_secant_kcal_mol_A':bm['central_minus_secant_kcal_mol_A']-am['central_minus_secant_kcal_mol_A']}
        paths.append(pair)
    for q in sorted({k.rsplit('_',1)[0]for k in rows if not k.startswith('matched_dry_') and not k.startswith('aquo8_')}):
        group='water_radial'if q.startswith('water_')else'hydrated_donor'
        c=compare({m:delta(s,q)-delta(s,'aquo8_central')for m,s in values.items()},group,q+'_balanced_exchange')
        c['per_metal_intermediate_energy_differences_kcal_mol']={m:{metal:s[q+'_'+metal]['energy']-s['aquo8_central_'+metal]['energy']for metal in ('La','Dy')}for m,s in values.items()}
        c['isolated_intermediates_are_not_affinities']=True;exchange.append(c)
    for f in ('acetate_mono','acetate_bi','acetamide'):
        addition.append(compare({m:delta(s,'hydrated_'+f+'_central')-delta(s,'matched_dry_'+f+'_central')for m,s in values.items()},'hydrated_donor',f+'_matched_water_addition'))
    for q in sorted({k.rsplit('_',1)[0]for k in rows}):
        group='dry_baseline'if q.startswith('matched_dry_')else'hydrated_donor'if q.startswith('hydrated_')else'fitted_centers'if q in ('water_central','aquo8_central')else'water_radial'
        row=rows[q+'_La']
        for kind in ('La','Dy','paired'):
            record={'geometry':q,'kind':kind,'group':group,'models':{}}
            def force(src):return body.subtract(src[q+'_Dy']['forces'],src[q+'_La']['forces'])if kind=='paired'else src[q+'_'+kind]['forces']
            target=force(ref)
            for method,src in (('baseline',original),('candidate',model)):
                predicted=force(src);weighted,raw=body.rigid_residual(row,predicted,target);rms=2*math.sqrt(sum(x*x for x in weighted));maximum=max(abs(x)for x in raw);atomic=body.subtract(predicted,target)
                record['models'][method]={'rigid_RMS_kcal_mol_A':rms,'maximum_body_component_error_kcal_mol_A':maximum,'passed':rms<=2 and maximum<=5,'raw_body_component_errors_kcal_mol_A':raw,'model_response':body.response(row,predicted),'reference_response':body.response(row,target),'atomic_force_residual_kcal_mol_A':atomic,'atom_vector_RMS_kcal_mol_A':math.sqrt(sum(x*x for v in atomic for x in v)/len(atomic))}
            a,b=record['models']['baseline'],record['models']['candidate'];record['RMS_error_change']=b['rigid_RMS_kcal_mol_A']-a['rigid_RMS_kcal_mol_A'];record['maximum_error_change']=b['maximum_body_component_error_kcal_mol_A']-a['maximum_body_component_error_kcal_mol_A'];record['deteriorated']=record['RMS_error_change']>0 or record['maximum_error_change']>0
            rigid.append(record);checks.append({'group':group,'id':q+'_'+kind,'kind':'rigid','baseline_passed':a['passed'],'candidate_passed':b['passed']})
    grouped={g:{'checks':len([v for v in checks if v['group']==g]),'baseline_passed':all(v['baseline_passed']for v in checks if v['group']==g),'candidate_passed':all(v['candidate_passed']for v in checks if v['group']==g),'failed_candidate':[v['id']for v in checks if v['group']==g and not v['candidate_passed']]}for g in ('fitted_centers','water_radial','hydrated_donor','dry_baseline')}
    return {'radial_work':works,'balanced_exchange':exchange,'matched_water_addition':addition,'rigid_response':rigid,'coarse_radial_diagnostics':paths,'development_checks':checks,'development_groups':grouped,'all_development_checks_passed':all(v['candidate_passed']for v in checks),'dry_baseline_failures_do_not_imply_water_parameter_deterioration':True,'independent_validation':False,'physical_accuracy_qualified':False,'protein_or_periodic_transfer_qualified':False}
