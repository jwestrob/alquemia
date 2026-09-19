"""Native DFT checks of selected coupled water motions and interior minima."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from hydration_basin_coordinates import WaterCoordinates
from hydration_water_motion import CENTERS, METHOD, checked_gradient
from hydration_square import endpoint
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms, write_xyz

STAGE='hydration_coupled_validation'
PROTOCOL='native_r2scan3c_cpcm_cartesian_anchored_mace_coupled_check_v1'
TOL={'amplitude':.05,'rotation_limit_radian':.10,'parallel_cosine':.95,
     'energy_absolute':.02,'energy_relative_even':.25,'even_absolute':.005,'even_relative':.25,
     'gradient_absolute_work':.02,'gradient_relative_change':.25,
     'minimum_energy_absolute':.05,'minimum_energy_relative':.10,'minimum_gradient_max':.4}


def normalize_direction(frame,direction):
    v=np.asarray(direction,dtype=float).copy()
    if np.linalg.norm(v)<1e-10:raise InvalidArtifact('zero validation direction')
    maxvel=float(np.max(np.linalg.norm(np.einsum('aci,i->ac',frame.jacobian(),v),axis=1)))
    v/=maxvel
    rotation=.05*max(np.linalg.norm(v.reshape(-1,6)[:,3:],axis=1))
    if rotation>TOL['rotation_limit_radian']:v*=TOL['rotation_limit_radian']/rotation
    # Finite rotations can shift the max-motion atom; use both actual signs.
    for _ in range(8):
        largest=max(float(np.max(np.linalg.norm(frame.positions(sign*.05*v)-frame.initial,axis=1))) for sign in (-1,1))
        if largest<=.05+1e-13:break
        v*=.05/largest
    if largest>.050000001:raise InvalidArtifact('failed exact physical displacement normalization')
    i=int(np.argmax(abs(v)))
    if v[i]<0:v=-v
    return v


def choose_directions(frame,result):
    curve=result['initial_curvature'];m=np.array(curve['mass_matrix'])
    modes=np.array(curve['physical_eigenvectors']);a=modes[:,0];b=np.array(result['minimum']['q'])
    n=float(np.sqrt(b@m@b))
    if n<1e-6 or abs(float(a@m@b))/n>TOL['parallel_cosine']:b=modes[:,-1]
    return {'soft':normalize_direction(frame,a),'response':normalize_direction(frame,b)}


def prepare(collection,agreement,output,proposals_only=False):
    import mace_omol as omol
    from hydration_basin import collect
    c=read_json(collection);cm=read_json(verify(c['manifest']));prep=read_json(verify(cm['preparation']))
    if c!=collect(verify(c['manifest'])) or c['status']!='complete':raise InvalidArtifact('actual complete coupled collection required')
    rows={r['task_id']:r for r in c['rows']};out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    _,mo,mm=omol.common(verify(cm['inventory']),verify(cm['software']),agreement,out/'mace',STAGE)
    for name in ('hydration_basin_validation.py','hydration_basin.py','hydration_basin_coordinates.py','hydration_water_motion.py','hydration_square.py'):
        p=mo/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),p);mm['implementation'][name]=record(p)
    selections={};definitions=[]
    names=sorted(prep['centers']) if proposals_only else CENTERS
    if proposals_only and (prep.get('generation') not in (1,2) or not prep.get('previous_validation')):
        raise InvalidArtifact('proposal-only continuation requires actual native reanchor')
    for name in names:
        center=prep['centers'][name];r=read_json(verify(rows[name]['result']))
        frame=WaterCoordinates(xyz(verify(center['xyz'])),center['groups']);directions={} if proposals_only else choose_directions(frame,r)
        selections[name]={'center':center,'basin_result':rows[name]['result'],
                          'directions':{k:v.tolist() for k,v in directions.items()}}
        for kind,v in directions.items():
            for sign in (-1,1):definitions.append({'task_id':name+'__'+kind+('__m' if sign<0 else '__p'),
                'center_id':name,'kind':kind,'amplitude':sign*.05,'direction':v.tolist(),'q':(sign*.05*v).tolist()})
        kind='minimum' if r['minimum']['cheap_minimum_eligible'] else 'proposal'
        definitions.append({'task_id':name+'__'+kind,'center_id':name,'kind':kind,'amplitude':None,
                            'direction':None,'q':r['minimum']['q']})
    pp=out/'preparation.json';write_new(pp,{'protocol_id':PROTOCOL,'selection':selections,
        'tasks':definitions,'source_collection':record(collection),'agreement':record(agreement),'tolerances':TOL,
        'proposals_only':proposals_only,'previous_validation':prep.get('previous_validation')})
    d=out/'dft';d.mkdir();di=d/'implementation';shutil.copytree(mo/'implementation',di)
    dpins={p.name:record(p) for p in di.glob('*.py')};dt=[];mt=[]
    for t in definitions:
        center=selections[t['center_id']]['center'];f=WaterCoordinates(xyz(verify(center['xyz'])),center['groups'])
        td=d/'tasks'/t['task_id'];td.mkdir(parents=True);xp=td/'core.xyz';write_xyz(xp,f.rows(t['q']))
        ip=td/'endpoint.inp';ip.write_text(METHOD+f"\n* xyzfile {center['charge']} 1 core.xyz\n")
        a={**t,'case':center['case'],'metal':center['metal'],'charge':center['charge'],'multiplicity':1,
           'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),'engrad_path':str(td/'endpoint.engrad'),
           'preparation':record(pp)}
        a['cache_key']=cache_key(a|{'protocol_id':PROTOCOL});dt.append(a)
        mt.append({**t,'case_id':center['case'],'metal':center['metal'],'charge':center['charge'],'spin_multiplicity':1,
                   'metal_index':0,'variant':'primary','xyz':record(xp),'state':check_atoms(xyz(xp),center['charge']),
                   'energy_component':omol.COMPONENT,'energy_only':False,'preparation':record(pp)})
    mm.update(tasks=mt,preparation=record(pp),tolerances=TOL,protocol_id=PROTOCOL)
    mr=omol.seal(mo,mm)
    source=read_json(verify(read_json(verify(prep['source_collection']))['manifest']))
    dm={'protocol_id':PROTOCOL,'tasks':dt,'preparation':record(pp),'implementation':dpins,'agreement':record(agreement),
        'orca':source['orca'],'execution_policy':source['execution_policy'],
        'execution_resources':{'mpi_ranks':16,'concurrent_tasks':min(4,len(dt))},'baseline_changed':False}
    write_new(d/'manifest.json',dm)
    from affordable_workflow import dry_run
    return {'MACE':mr,'DFT':dry_run(d/'manifest.json')}


def validate(manifest):
    import mace_omol as omol
    m=read_json(manifest);p=read_json(verify(m['preparation']))
    if m['stage']!=STAGE or m['tolerances']!=TOL or p['tolerances']!=TOL or m['model']!=omol.model(verify(m['software'])):
        raise InvalidArtifact('changed validation physics')
    for pin in [p['agreement'],p['source_collection'],*m['implementation'].values()]:verify(pin)
    definitions={t['task_id']:t for t in p['tasks']}
    if len(m['tasks'])!=len(definitions) or {t['task_id'] for t in m['tasks']}!=set(definitions):raise InvalidArtifact('changed validation coverage')
    for t in m['tasks']:
        d=definitions[t['task_id']];c=p['selection'][t['center_id']]['center']
        f=WaterCoordinates(xyz(verify(c['xyz'])),c['groups'])
        actual=xyz(verify(t['xyz']));wanted=f.rows(t['q'])
        if (any(t[k]!=v for k,v in d.items()) or [a[0] for a in actual]!=[a[0] for a in wanted]
                or not np.allclose([a[1:] for a in actual],[a[1:] for a in wanted],atol=1e-12,rtol=0) or t['charge']!=c['charge']):
            raise InvalidArtifact('changed physical validation path')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('validation cache identity changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def collect_mace(manifest):
    validate(manifest);m=read_json(manifest);rows=[]
    for t in m['tasks']:
        valid=[(a,accepted_attempt(a,t,manifest)) for a in sorted((Path(manifest).parent/'execution'/t['task_id']).glob('attempt_*'))]
        valid=[(a,r) for a,r in valid if r is not None]
        if not valid:rows.append({'task_id':t['task_id'],'status':'unavailable'});continue
        a,r=valid[-1];rows.append({'task_id':t['task_id'],'status':'computed','result':record(a/'result.json'),
                                  'receipt':record(a/'receipt.json'),'energy_eV':r['energy_eV'],'forces':r['forces']})
    return {'manifest':record(manifest),'rows':rows,'status':'complete' if all(r['status']=='computed' for r in rows) else 'incomplete'}


def collect_dft(manifest,output):
    m=read_json(manifest);p=read_json(verify(m['preparation']));rows=[]
    for t in m['tasks']:
        try:
            c=p['selection'][t['center_id']]['center'];atoms=xyz(verify(t['xyz']))
            f=WaterCoordinates(xyz(verify(c['xyz'])),c['groups'])
            wanted=f.rows(t['q'])
            if ([a[0] for a in atoms]!=[a[0] for a in wanted] or
                    not np.allclose([a[1:] for a in atoms],[a[1:] for a in wanted],atol=1e-12,rtol=0)):
                raise InvalidArtifact('changed validation geometry')
            r=endpoint(record(t['output_path']),record(t['output_path']+'.execution.json'),t['xyz'],t['input'])
            gp=read_json(verify(r['receipt']))['artifacts']['engrad'];g=checked_gradient(gp,r,atoms)
            rows.append({'task_id':t['task_id'],'status':'complete','result':r,'gradient':gp,
                         'physical_gradient':f.gradient(t['q'],g).tolist()})
        except (OSError,ValueError,KeyError) as exc:rows.append({'task_id':t['task_id'],'status':'unavailable','failure':str(exc)})
    r={'manifest':record(manifest),'rows':rows,'status':'complete' if all(t['status']=='complete' for t in rows) else 'incomplete'}
    write_new(output,r);return {'status':r['status'],'completed':sum(t['status']=='complete' for t in rows)}


def compare(preparation,mace_collection,dft_collection,output):
    p=read_json(preparation);mc=read_json(mace_collection);dc=read_json(dft_collection)
    if mc!=collect_mace(verify(mc['manifest'])) or dc['status']!='complete' or mc['status']!='complete':
        raise InvalidArtifact('complete actual calculations required')
    dm=read_json(verify(dc['manifest']));mm=read_json(verify(mc['manifest']))
    if dm['preparation']!=record(preparation) or mm['preparation']!=record(preparation):raise InvalidArtifact('different comparison inputs')
    mr={r['task_id']:r for r in mc['rows']};dr={r['task_id']:r for r in dc['rows']};rows=[]
    for t in dm['tasks']:
        s=p['selection'][t['center_id']];c=s['center'];br=read_json(verify(s['basin_result']))
        f=WaterCoordinates(xyz(verify(c['xyz'])),c['groups']);q=np.array(t['q']);x=f.positions(q)
        d=dr[t['task_id']];r=d['result']
        if endpoint(r['output'],r['receipt'],t['xyz'],t['input'])!=r:raise InvalidArtifact('DFT output changed')
        gd=checked_gradient(d['gradient'],r,xyz(verify(t['xyz'])));gdq=f.gradient(q,gd)
        correction=np.array(br['cartesian_anchor_correction_eV_A'])
        force=np.load(verify(mr[t['task_id']]['forces']),allow_pickle=False)
        gpq=f.gradient(q,(-force+correction)*EV_TO_KCAL)
        prediction=((mr[t['task_id']]['energy_eV']-c['MACE_energy_eV'])+float(np.sum(correction*(x-f.initial))))*EV_TO_KCAL
        observed=(r['energy_hartree']-c['DFT_result']['energy_hartree'])*HA_TO_KCAL
        g0=f.gradient(np.zeros(f.dimension),checked_gradient(c['DFT_gradient'],c['DFT_result'],xyz(verify(c['xyz']))))
        row={'task_id':t['task_id'],'center_id':t['center_id'],'kind':t['kind'],'amplitude':t['amplitude'],
             'predicted_change_kcal_mol':prediction,'actual_change_kcal_mol':observed,'energy_error_kcal_mol':prediction-observed,
             'DFT_physical_gradient':gdq.tolist(),'predicted_physical_gradient':gpq.tolist(),
             'DFT_gradient_max':float(np.max(abs(gdq)))}
        if t['kind'] in ('soft','response'):
            v=np.array(t['direction']);row.update(DFT_projected_gradient=float(gdq@v),predicted_projected_gradient=float(gpq@v),
                                                 DFT_center_projected_gradient=float(g0@v))
        else:
            et=max(TOL['minimum_energy_absolute'],TOL['minimum_energy_relative']*abs(observed))
            row.update(energy_tolerance=et,energy_pass=abs(prediction-observed)<=et,
                       stationary_within_tolerance=bool(np.max(abs(gdq))<=TOL['minimum_gradient_max']))
        rows.append(row)
    checks=[]
    for name in ([] if p.get('proposals_only') else p['selection']):
        for kind in ('soft','response'):
            pair=[r for r in rows if r['center_id']==name and r['kind']==kind]
            if len(pair)!=2:raise InvalidArtifact('missing both signs')
            even=sum(r['actual_change_kcal_mol'] for r in pair)/2;pred=sum(r['predicted_change_kcal_mol'] for r in pair)/2
            et=max(TOL['energy_absolute'],TOL['energy_relative_even']*abs(even));ct=max(TOL['even_absolute'],TOL['even_relative']*abs(even))
            for r in pair:
                workerror=abs(r['predicted_projected_gradient']-r['DFT_projected_gradient'])*.05
                gt=max(TOL['gradient_absolute_work'],TOL['gradient_relative_change']*abs(r['DFT_projected_gradient']-r['DFT_center_projected_gradient'])*.05)
                r.update(energy_tolerance=et,energy_pass=abs(r['energy_error_kcal_mol'])<=et,
                         gradient_work_error=workerror,gradient_work_tolerance=gt,gradient_pass=workerror<=gt)
            checks.append({'center_id':name,'kind':kind,'DFT_even_kcal_mol':even,'predicted_even_kcal_mol':pred,
                           'even_error':pred-even,'even_tolerance':ct,'even_pass':abs(pred-even)<=ct,
                           'pass':abs(pred-even)<=ct and all(r['energy_pass'] and r['gradient_pass'] for r in pair)})
    result={'protocol_id':PROTOCOL,'preparation':record(preparation),'MACE_collection':record(mace_collection),
            'DFT_collection':record(dft_collection),'rows':rows,'direction_checks':checks,
            'all_direction_checks_pass':all(r['pass'] for r in checks) if checks else None,'status':'complete','baseline_changed':False,
            'occupancy_probabilities':None,'thermochemical_status':'basin_and_other_terms_not_yet_qualified',
            'implementation':record(__file__)}
    write_new(output,result)
    return {'status':'complete','all_direction_checks_pass':result['all_direction_checks_pass'],'direction_checks':checks}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,names in [('prepare',('collection','agreement','output')),('collect-dft',('manifest','output')),
                     ('compare',('preparation','mace-collection','dft-collection','output'))]:
        q=sub.add_parser(op)
        for name in names:q.add_argument('--'+name,required=True)
        if op=='prepare':q.add_argument('--proposals-only',action='store_true')
    a=vars(p.parse_args());op=a.pop('op').replace('-','_');print(json.dumps(globals()[op](**a),indent=2))
