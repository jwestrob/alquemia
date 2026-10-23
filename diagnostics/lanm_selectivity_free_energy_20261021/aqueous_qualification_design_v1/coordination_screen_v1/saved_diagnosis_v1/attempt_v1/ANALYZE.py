#!/usr/bin/env python3
"""Saved W8/W9 native replay and balanced component/force accounting only."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('screen',HERE.parent/'RUN.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main():
    start,cpu=time.monotonic(),time.process_time();checks=[];pins={}
    def load(name,path):
        pins[name]=m.pin(path);return m.read(pins[name])
    def check(ok,name):
        checks.append({'name':name,'passed':bool(ok)})
        m.require(ok,name)
    actual=load('actual17670',m.W/'aqueous_qualification_design_v1/coordination_screen_v1/run_v1/execution/RESULT.json')
    bound=m.read(actual['bound']);pins['bound']=actual['bound'];ccfg=m.read(bound['candidate_config'])
    transfer=m.read(ccfg['transfer_result']);pins['transfer17612']=ccfg['transfer_result'];tb=m.read(transfer['bound']);data=m.read(tb['data']);pins['consumed32_metadata']=tb['data']
    old=load('original_panel',m.W/'explicit_model_v1/panel_v1/RESULT.json');oldrows={r['id']:r for r in old['rows']};oldcfg=m.read(old['config'])
    native=m.module(bound['native_helper'],'diag_native');body=m.module(bound['body_helper'],'diag_body');numerical=m.module(bound['numerical_helper'],'diag_num')
    check(actual['status']=='complete' and actual['job_id']=='17670' and actual['native_program_invocations']==4 and actual['comparison']['passed'] is False,'Actual failed17670 is complete and preserved')
    check(old['status']=='panel_complete' and old['implementation']==bound['native_helper'] and oldcfg['executables']==bound['executables'],'Original and candidate use exact same native parser/executables')
    rows={};models={'original':{},'frozen17605':{}};references={};admissions={}
    for n in (8,9):
        for metal in ('La','Dy'):
            gid=f'aquo{n}_central_{metal}';o=oldrows[gid]
            row=copy.deepcopy(data['rows'][gid] if n==8 else bound['rows'][gid]);rows[gid]=row
            ref=bound['anchor_reference'][metal] if n==8 else row['reference'];references[gid]=ref
            m.verify(ref['raw_result']);check(o['status']=='complete',gid+' old complete')
            for key in ('xyz','key','parameters'):m.verify(o['input'][key])
            check(o['input']['xyz']['sha256']==row['input']['xyz']['sha256'] and o['input']['bonds_zero_based']==row['input']['bonds_zero_based'] and o['input']['types']==row['input']['types'],gid+' exact typed coordinate/topology identity')
            og=copy.deepcopy(o['input']['source_geometry']);m.verify(og['xyz']);m.verify(row['geometry']['xyz']);check(og['xyz']['sha256']==row['geometry']['xyz']['sha256'],gid+' exact source geometry bytes');og['xyz']=row['geometry']['xyz'];check(og==row['geometry'],gid+' all geometry metadata agree after verified path alias')
            current=bound['anchor_model'][metal] if n==8 else actual['rows'][gid]['model']
            currentfiles=tb['inputs'][gid]['files'] if n==8 else row['files']
            for key in ('xyz','key','parameters'):m.verify(currentfiles[key]);check(currentfiles[key]['sha256']==(o['input'][key]['sha256'] if key!='key' else currentfiles[key]['sha256']),gid+' '+key+' admitted')
            check(m.verify(currentfiles['key']).read_text().startswith(m.verify(o['input']['key']).read_text()),gid+' original key prefix retained')
            versions={'original':{'energy':o['native']['total_energy_kcal_mol'],'forces':o['native']['force_kcal_mol_A'],'native':o['native']},'frozen17605':current}
            for version,v in versions.items():
                ni=v['native'];inp=copy.deepcopy(o['input']);inp['derivative_test']=ni['finite_difference'] is not None
                ap=m.verify(ni['analyze_output']);gp=m.verify(ni['gradient_output'])
                numerical.program_admission(ap,'analyze',native);numerical.program_admission(gp,'testgrad',native)
                replay=native.parse(ap,gp,inp)
                check(replay==ni and replay['total_energy_kcal_mol']==v['energy'] and replay['force_kcal_mol_A']==v['forces'],gid+' '+version+' exact raw native replay')
                tol=(len(ni['components'])+1)*.5e-8+1e-10
                check(abs(ni['printed_component_sum_remainder_kcal_mol'])<=tol,gid+' '+version+' component closure within printed precision')
                models[version][gid]=v;admissions[version+'__'+gid]={'analyze':ni['analyze_output'],'gradient':ni['gradient_output'],'component_closure_remainder':ni['printed_component_sum_remainder_kcal_mol'],'component_printing_bound':tol,'old_parameters':o['input']['parameters'],'candidate_files':currentfiles}
    ids=lambda n:(f'aquo{n}_central_La',f'aquo{n}_central_Dy')
    def contrast(v,component=None):
        def E(g):return v[g]['energy'] if component is None else v[g]['native']['components'][component]['energy_kcal_mol']
        return E(ids(9)[1])-E(ids(9)[0])-E(ids(8)[1])+E(ids(8)[0])
    refcontrast=contrast(references);compkeys=list(models['original'][ids(8)[0]]['native']['components']);energy={};force={};radial={}
    for version,values in models.items():
        energy[version]={'contrast_kcal_mol':contrast(values),'reference_contrast_kcal_mol':refcontrast,'error_kcal_mol':contrast(values)-refcontrast,'components':{k:contrast(values,k) for k in compkeys},'per_metal_W9_minus_W8_model':{metal:values[f'aquo9_central_{metal}']['energy']-values[f'aquo8_central_{metal}']['energy'] for metal in ('La','Dy')},'per_cell_energy_kcal_mol':{g:v['energy'] for g,v in values.items()}}
        check(abs(sum(energy[version]['components'].values())-energy[version]['contrast_kcal_mol'])<=4*(len(compkeys)+1)*.5e-8+1e-9,version+' balanced component closure')
        force[version]={};radial[version]={}
        for n in (8,9):
            la,dy=ids(n);row=rows[la]
            for kind in ('La','Dy','paired'):
                a=body.subtract(values[dy]['forces'],values[la]['forces']) if kind=='paired' else values[f'aquo{n}_central_{kind}']['forces']
                b=body.subtract(references[dy]['forces'],references[la]['forces']) if kind=='paired' else references[f'aquo{n}_central_{kind}']['forces']
                weighted,raw=body.rigid_residual(row,a,b);rms=2*math.sqrt(sum(v*v for v in weighted));maximum=max(abs(v) for v in raw)
                force[version][f'W{n}_{kind}']={'RMS_kcal_mol_A':rms,'maximum_component_kcal_mol_A':maximum,'passed':rms<=2 and maximum<=5}
            for metal in ('La','Dy'):
                gid=f'aquo{n}_central_{metal}';row=rows[gid];x=row['geometry']['positions_A'];mi=row['geometry']['metal_index'];records=[]
                ma=body.response(row,values[gid]['forces']);ra=body.response(row,references[gid]['forces'])
                for ordinal,(a,b) in enumerate(zip(ma,ra)):
                    if a['category']!='water':continue
                    oi=next(i for i in a['indices'] if row['geometry']['symbols'][i]=='O');vec=[x[oi][k]-x[mi][k] for k in range(3)];radius=math.sqrt(sum(z*z for z in vec));u=[z/radius for z in vec]
                    am=sum(z*t for z,t in zip(a['net_force'],u));ar=sum(z*t for z,t in zip(b['net_force'],u));records.append({'water_ordinal':ordinal,'oxygen_index_zero_based':oi,'radius_A':radius,'model_outward_force':am,'reference_outward_force':ar,'error_outward_force':am-ar})
                radial[version][gid]={'water_bodies':records,'summed_radial_force_error':sum(r['error_outward_force'] for r in records),'first_eight_sum_error':sum(r['error_outward_force'] for r in records[:8]),'ninth_error':records[8]['error_outward_force'] if n==9 else None,'metal_force_error':[a-b for a,b in zip(values[gid]['forces'][mi],references[gid]['forces'][mi])]}
    change={'contrast_change_kcal_mol':energy['frozen17605']['contrast_kcal_mol']-energy['original']['contrast_kcal_mol'],'component_changes_kcal_mol':{k:energy['frozen17605']['components'][k]-energy['original']['components'][k] for k in compkeys}}
    check(abs(energy['frozen17605']['error_kcal_mol']-actual['comparison']['energy']['error_kcal_mol'])<1e-8,'Independent energy arithmetic matches preserved17670')
    for kind in ('La','Dy','paired'):
        check(force['frozen17605']['W9_'+kind]=={k:actual['comparison']['rigid'][kind][k] for k in ('RMS_kcal_mol_A','maximum_component_kcal_mol_A','passed')},kind+' exact preserved force scorer replay')
    panel=m.read(m.read(m.read(bound['config'])['references']['aquo9_central_La'])['panel']);pure=[c for c in panel['cells'] if c['recipe']=='ecp_pbe0' and (c['geometry']['case_id'] in ('water','aquo8','aquo9'))]
    check(len(pure)==10,'Consumed pure-water reference inventory has exact10 cells')
    inventory=[]
    for c in pure:
        gid=c['geometry']['geometry_id'];inventory.append({'id':gid,'method':c['recipe'],'charge':c['charge'],'multiplicity_effective':c['multiplicity'],'multiplicity_physical':c['physical_multiplicity'],'ecp_core':c['ecp_core'],'explicit_electrons':c['explicit_electrons'],'geometry':c['geometry']['xyz'],'input':c['input'],'candidate_native_available':gid in models['frozen17605'] or gid in transfer['rows'],'reference_source':'four references admitted in coordination preparation' if gid in references else 'consumed radial training metadata; not new outcome','reference_raw_result':references[gid]['raw_result'] if gid in references else None})
    result={'status':'complete','implementation':m.pin(__file__),'specification':m.pin(HERE/'SPEC.md'),'input_pins':pins,'native_admissions':admissions,'checks':checks,'checks_passed':len(checks),'energy':energy,'fit_induced_change':change,'rigid_force':force,'radial_force':radial,'consumed_pure_water_inventory':inventory,'original_screen_passed':abs(energy['original']['error_kcal_mol'])<=1 and all(force['original']['W9_'+k]['passed'] for k in ('La','Dy','paired')),'reserved_donor_outcomes_opened':False,'native_calls':0,'model_refit':False,'physical_accuracy_qualified':False,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    out=m.W/'aqueous_qualification_design_v1/coordination_screen_v1/saved_diagnosis_v1';out.mkdir(parents=True,exist_ok=False);m.save(out/'RESULT.json',result)
    print(json.dumps({'result':m.pin(out/'RESULT.json'),'energy':energy,'change':change,'force':force,'checks':len(checks),'wall_seconds':result['wall_seconds']}))

if __name__=='__main__':main()
