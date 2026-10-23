#!/usr/bin/env python3
"""Saved-only La rigid-body discrepancy localization, no engine calls."""
import argparse,hashlib,importlib.util,json,math,os,time
from pathlib import Path
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):assert pin(p['path'])==p;return Path(p['path'])
def norm(a):return math.sqrt(sum(v*v for v in a))
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def run(config,out):
    start=time.monotonic();c=json.loads(Path(config).read_text());assert c['implementation']==pin(__file__)
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    verify(c['specification']);d=json.loads(verify(c['training']).read_text());assert d['status']=='training_complete' and d['complete']==28
    assert c['objective']==d['objective'];s=importlib.util.spec_from_file_location('la_objective',verify(c['objective']));o=importlib.util.module_from_spec(s);s.loader.exec_module(o)
    out.mkdir(parents=True,exist_ok=False);rows=[]
    for gid,r in d['rows'].items():
        if not gid.endswith('_La'):continue
        verify(r['baseline_result']);verify(r['reference_result'])
        geom=r['geometry'];x=geom['positions_A'];metal=geom['metal_index'];m=o.response(r,r['baseline']['forces']);q=o.response(r,r['reference']['forces'])
        bodies=[];cats=sorted(set(b['category'] for b in m))
        donor_ids=geom.get('donor_indices')
        if not donor_ids and gid=='water_central_La':donor_ids=[i for i,a in enumerate(geom['symbols']) if a=='O']
        ligand_axis=None
        if donor_ids:
            z=[sum(x[i][k] for i in donor_ids)/len(donor_ids)-x[metal][k] for k in range(3)];ligand_axis=[v/norm(z) for v in z]
        for a,b in zip(m,q):
            f=[u-v for u,v in zip(a['net_force'],b['net_force'])];t=[u-v for u,v in zip(a['torque_over_1A'],b['torque_over_1A'])]
            axis=ligand_axis
            if a['category']=='water':
                oi=next(i for i in a['indices'] if geom['symbols'][i]=='O');z=[x[oi][k]-x[metal][k] for k in range(3)];axis=[v/norm(z) for v in z]
            radial=sum(u*v for u,v in zip(f,axis)) if axis else None
            tang=[f[k]-radial*axis[k] for k in range(3)] if axis else None
            ncat=sum(z['category']==a['category'] for z in m)
            ss=sum(v*v for v in f)/3
            if a['category']!='ion':ss=(ss+sum(v*v for v in t)/3)/2
            bodies.append({'category':a['category'],'indices':a['indices'],'center_A':a['center_A'],
                'model_net_force_kcal_mol_A':a['net_force'],'reference_net_force_kcal_mol_A':b['net_force'],
                'force_error_vector_kcal_mol_A':f,'force_error_norm_kcal_mol_A':norm(f),
                'model_torque_kcal_mol':a['torque_over_1A'],'reference_torque_kcal_mol':b['torque_over_1A'],
                'torque_error_vector_kcal_mol':t,'torque_error_norm_kcal_mol':norm(t),
                'radial_outward_axis':axis,'radial_force_error_kcal_mol_A':radial,
                'tangential_force_error_vector_kcal_mol_A':tang,'tangential_force_error_norm_kcal_mol_A':norm(tang) if tang else None,
                'squared_rigid_error':ss,'weighted_squared_contribution':ss/(len(cats)*ncat)})
        group='hydrated_'+geom['distortion'] if gid.startswith('hydrated') else 'matched_dry' if gid.startswith('matched_dry') else 'water_controls'
        residual,_=o.rigid_residual(r,r['baseline']['forces'],r['reference']['forces']);rms=2*norm(residual)
        assert abs(rms*rms-sum(b['weighted_squared_contribution'] for b in bodies))<1e-10
        rows.append({'id':gid,'group':group,'category_balanced_rigid_RMS_kcal_mol_A':rms,
                     'geometry':geom,'bodies':bodies,'category_summaries':{cat:{
                         'mean_squared_rigid_error':sum(b['squared_rigid_error'] for b in bodies if b['category']==cat)/sum(b['category']==cat for b in bodies),
                         'outward_radial_error_values':[b['radial_force_error_kcal_mol_A'] for b in bodies if b['category']==cat],
                         'force_error_norm_values':[b['force_error_norm_kcal_mol_A'] for b in bodies if b['category']==cat]}
                         for cat in cats}})
    assert len(rows)==14
    for stage,ids in (('water',o.WATER),('donor',o.DONOR)):
        v=math.sqrt(sum(r['category_balanced_rigid_RMS_kcal_mol_A']**2 for r in rows if r['id'][:-3] in ids)/len(ids))
        assert abs(v-d['baseline_objectives'][stage]['controls']['La']['category_balanced_rigid_RMS_kcal_mol_A'])<1e-10
    rows.sort(key=lambda r:r['category_balanced_rigid_RMS_kcal_mol_A'],reverse=True)
    result={'status':'complete','implementation':pin(__file__),'config':pin(config),'training':c['training'],
        'complete_La_fixtures':14,'rows':rows,'new_molecular_calls':0,'local_threads':1,'validation_outcomes_opened':False,
        'radial_screen_unchanged':d['radial_derivative_screen'],'baseline_controls':d['baseline_objectives'],
        'bookkeeping_checks_passed':True,'fit_released':False,'physical_accuracy_qualified':False,'elapsed_seconds':time.monotonic()-start}
    save(out/'RESULT.json',result);print(json.dumps({'result':pin(out/'RESULT.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.config,a.output)
