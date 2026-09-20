"""One geometry-selected reference torsion using exact campaign center receipts."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,record,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL,write_xyz
from mace_site_kinematics import Kinematics
from accommodation_torsion_profiles import geometry_check
from compact_solvation import input_text,completed,diagnostics
CID='mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-1'
ROLE='extra_acidic_ligand_homolog'
PROTOCOL='geometry_selected_reference_extraAsp_grid_v1'


def prepare(preparation,mapping_result,solvent_manifest,comparison,agreement,output):
    source=read_json(preparation);case=next(x for x in source['cases'] if x['case_id']==CID)
    check=next(x for x in read_json(mapping_result)['rows'] if x['case_id']==CID)
    if not check['roles'][ROLE]['warning'] or check['mapping_status']!='supported':raise InvalidArtifact('geometry selection/mapping differs')
    sm=read_json(solvent_manifest);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    gp=read_json(verify(check['mapping_check']['mapping']));mode=check['mapping_check']['roles'][ROLE]['mode'];kin=Kinematics(gp['context'])
    index=next(i for i,m in enumerate(kin.modes) if m['id']==mode);points=[];tasks=[];centers={};model=None;software=None
    for metal in ('Ca','La'):
        quartet=[t for t in sm['all_tasks'] if t['case_id']==CID and t['representation']=='context' and t['metal']==metal]
        if len(quartet)!=2:raise InvalidArtifact('source center missing')
        e=quartet[0]['source_endpoint'];rows=xyz(verify(e['xyz']));mr=read_json(verify(e['native_MACE_receipt']));mm=read_json(verify(e['source_manifest']))
        mt=next(t for t in mm['tasks'] if t['task_id']==mr['task_id'])
        if mr['status']!='computed' or mr['manifest']!=e['source_manifest'] or mr['cache_key']!=mt['cache_key'] or mt['xyz']!=e['xyz'] or mt['charge']!=e['charge']:
            raise InvalidArtifact('native center state/receipt differs')
        if xyz(verify(case['representations']['context']['endpoints'][metal]['xyz']))!=rows:raise InvalidArtifact('center prepared geometry differs')
        if not np.allclose(kin.core,[a[1:] for a in rows],atol=1e-12,rtol=0):raise InvalidArtifact('physical center differs')
        if model is not None and model!=mm['model']:raise InvalidArtifact('metal-dependent model differs')
        model=mm['model'];software=mm['software'];low={}
        for t in quartet:
            pin=completed(solvent_manifest,t['task_id'])
            if pin is None:raise InvalidArtifact('q0 source GFN2 not complete')
            audit=diagnostics(pin,t)
            if audit['charge_sanity_status']!='pass':raise InvalidArtifact('q0 native charges unsupported')
            if verify(t['input']).read_text()!=input_text(t['charge'],1,t['medium'],'native'):raise InvalidArtifact('q0 recipe differs')
            low[t['medium']]={**pin,'audit':audit}
        centers[metal]={'xyz':e['xyz'],'charge':e['charge'],'MACE_receipt':e['native_MACE_receipt'],'MACE_manifest':e['source_manifest'],
            'MACE_energy_eV':mr['energy_eV'],'GFN2':low}
        for angle in (-.2,.2):
            q=np.zeros(len(kin.modes));q[index]=angle;coords=kin.evaluate(q)[1];tid=CID+'__'+str(angle).replace('-','m').replace('.','p')+'__'+metal
            pd=out/'points'/tid;pd.mkdir(parents=True);xp=pd/'core.xyz';write_xyz(xp,[(a[0],*p) for a,p in zip(rows,coords)])
            gc=geometry_check(kin,q,[a[0] for a in rows]);point={'task_id':tid,'case_id':CID,'metal':metal,'angle_radian':angle,'role':ROLE,
                'mode_id':mode,'q':q.tolist(),'xyz':record(xp),'charge':e['charge'],'multiplicity':1,'geometry_checks':gc,'mapping':check['mapping_check']['mapping']}
            if not gc['pass']:raise InvalidArtifact('new geometric overlap/map failure')
            points.append(point)
            for medium in ('vacuum','alpb'):
                ld=out/'low/tasks'/(tid+'__'+medium);ld.mkdir(parents=True);lx=ld/'core.xyz';shutil.copyfile(xp,lx)
                ip=ld/'endpoint.inp';ip.write_text(input_text(e['charge'],1,medium,'native'))
                tasks.append({'task_id':tid+'__'+medium,'point_id':tid,'case_id':CID,'case':CID,'metal':metal,'medium':medium,
                    'charge':e['charge'],'multiplicity':1,'xyz':record(lx),'input':record(ip),'output_path':str(ld/'endpoint.out'),'geometry_supported':True})
    bands=read_json(comparison)['frozen_bands']['context']
    design={'protocol_id':PROTOCOL,'case_id':CID,'agreement':record(agreement),'preparation':record(preparation),
        'mapping_result':record(mapping_result),'solvent_source_manifest':record(solvent_manifest),'source_conditioning_metal':case['source']['source_conditioning_metal'],
        'source_metal':case['source']['raw_source_metal'],'expected_class':case['expected_class'],'mode_id':mode,
        'centers':centers,'points':points,'model':model,'software':software,'frozen_band_source':record(comparison),'frozen_bands':bands,
        'bounded_grid_radians':[-.2,0.,.2],'new_calls':{'MACE':4,'GFN2':8,'DFT':0},'all_evidence_consumed':True,'baseline_changed':False}
    write_new(out/'design.json',design);impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        dp=impl/p.name;shutil.copyfile(p,dp);pins[p.name]=record(dp)
    lm={'protocol_id':PROTOCOL,'stage':'low','design':record(out/'design.json'),'agreement':record(agreement),'orca':sm['orca'],
        'implementation':pins,'tasks':tasks,'all_tasks':tasks,'execution_resources':{'mpi_ranks':8,'concurrent_tasks':8},
        'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']}}
    write_new(out/'low/manifest.json',lm)
    write_new(out/'manifest.json',{'protocol_id':PROTOCOL,'design':record(out/'design.json'),'implementation':pins,'low_manifest':record(out/'low/manifest.json')})
    from affordable_workflow import dry_run
    return dry_run(out/'low/manifest.json')


def decision(value,bands):
    if value is None:return None
    return 'Ca' if value<=bands['Ca_max'] else ('La' if value>=bands['La_min'] else 'indeterminate')


def analyze(manifest,mace_result,low_collection,output):
    m=read_json(manifest);d=read_json(verify(m['design']));mr=read_json(mace_result);lc=read_json(low_collection)
    if mr['manifest']!=record(manifest) or lc['manifest']!=m['low_manifest']:raise InvalidArtifact('new result receipts differ')
    mc={r['task_id']:r for r in mr['rows']};lr={r['task_id']:r for r in lc['rows']};points=[]
    for metal,c in d['centers'].items():
        transfer=(c['GFN2']['alpb']['energy_hartree']-c['GFN2']['vacuum']['energy_hartree'])*HA_TO_KCAL
        points.append({'metal':metal,'angle_radian':0.,'native_kcal_mol':c['MACE_energy_eV']*EV_TO_KCAL,
            'composite_kcal_mol':c['MACE_energy_eV']*EV_TO_KCAL+transfer,'solvent_transfer_kcal_mol':transfer,'reused':True,'status':'complete'})
    for p in d['points']:
        a=mc[p['task_id']];v,b=[lr[p['task_id']+'__'+medium] for medium in ('vacuum','alpb')]
        row={'metal':p['metal'],'angle_radian':p['angle_radian'],'native_kcal_mol':None,'composite_kcal_mol':None,'solvent_transfer_kcal_mol':None,'reused':False,'status':'unavailable'}
        if a['status']=='complete':row['native_kcal_mol']=a['energy_eV']*EV_TO_KCAL
        if all(x['status']=='complete' for x in (a,v,b)):
            tr=(b['energy_hartree']-v['energy_hartree'])*HA_TO_KCAL;row.update(status='complete',composite_kcal_mol=row['native_kcal_mol']+tr,solvent_transfer_kcal_mol=tr)
        points.append(row)
    models={}
    for model in ('native','composite'):
        key=model+'_kcal_mol';result={'status':'unavailable','R_grid_kcal_mol':None,'free_energy':None,'old_band_transfer_only':True}
        if all(p[key] is not None for p in points):
            by={metal:[p for p in points if p['metal']==metal] for metal in ('Ca','La')}
            center={metal:next(p[key] for p in rows if p['angle_radian']==0) for metal,rows in by.items()}
            minima={metal:min(rows,key=lambda p:p[key]) for metal,rows in by.items()}
            R0=center['Ca']-center['La'];Rg=minima['Ca'][key]-minima['La'][key]
            contrasts=[]
            for q in d['bounded_grid_radians']:
                energies={z:next(p[key] for p in by[z] if p['angle_radian']==q) for z in by}
                R=energies['Ca']-energies['La'];contrasts.append({'angle_radian':q,'R_kcal_mol':R,'delta_R_from_center_kcal_mol':R-R0,
                    'endpoint_work_kcal_mol':{z:energies[z]-center[z] for z in by},'old_band_transfer_call':decision(R,d['frozen_bands'][model])})
            result.update(status='available',R_grid_kcal_mol=Rg,R_center_kcal_mol=R0,delta_R_grid_kcal_mol=Rg-R0,
                selected_angles={z:minima[z]['angle_radian'] for z in by},endpoint_lowering_kcal_mol={z:minima[z][key]-center[z] for z in by},
                fixed_displacements=contrasts,old_band_transfer_center=decision(R0,d['frozen_bands'][model]),old_band_transfer_grid=decision(Rg,d['frozen_bands'][model]))
        models[model]=result
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'MACE_result':record(mace_result),'low_collection':record(low_collection),
        'analyzer':record(__file__),'points':points,'models':models,'expected_class':d['expected_class'],'source_conditioning_metal':d['source_conditioning_metal'],
        'frozen_bands':d['frozen_bands'],'all_evidence_consumed':True,'threshold_refitted':False,'new_DFT_calls':0,'baseline_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('prepare')
    for k in ('preparation','mapping_result','solvent_manifest','comparison','agreement','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    q=s.add_parser('analyze')
    for k in ('manifest','mace_result','low_collection','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());op=a.pop('op');r=globals()[op](**a);print(json.dumps({k:v for k,v in r.items() if k!='points'},indent=2))
if __name__=='__main__':main()
