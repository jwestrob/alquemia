"""Fixed-union context adapter for the existing scaled angular proposals/pool."""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import os
from pathlib import Path
import shutil
import time
import numpy as np
import scipy
import scipy.optimize._slsqp_py as scipy_slsqp
import scipy.optimize._slsqplib as scipy_kernel
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project,preview
from accommodation_nonlinear import WarmGPU
import adaptive_completion as completion
import adaptive_angular_proposals as angular
import nikasha_pool as pool

PROTOCOL='fixed_union_context_four_angular_native_OMOL_SLSQP_v1'
POOL_PROTOCOL='nikasha_union_adaptive_minimal_common_geometry_native_OMOL_GFN2_ALPB_v1'
SETTINGS=copy.deepcopy(completion.SETTINGS)
POOL_SETTINGS={**pool.SETTINGS,'candidate_order':['origin','adaptive_Ca','adaptive_La']}
CANONICAL_REUSE=('1H4I','4MAE','q9z4j7-pqq-la_model','q88jh5-pqq-la_model')


def declared_population(inputs,calibration):
    data=read_json(inputs);population=data.get('population','common8')
    if population=='common8':
        if len(data['cases'])!=8:raise InvalidArtifact('exact common8 required')
        return population,8
    if population!='canonical_remaining24':raise InvalidArtifact('undeclared union population')
    if data['calibration']!=record(calibration):raise InvalidArtifact('canonical union source collection differs')
    cal=read_json(calibration);canonical=[r for r in cal['rows'] if r['canonical_coordinate_match']]
    crystals=[r for r in cal['rows'] if not r['canonical_coordinate_match']]
    if len(canonical)!=25 or {r['case_id'] for r in crystals}!={'1H4I','4MAE','1KB0'}:raise InvalidArtifact('canonical25/crystal3 differs')
    expected=[{'case_id':r['root_case_id'] if r['canonical_coordinate_match'] else r['case_id'],
        'actual_union_case_id':r['case_id'],'known_class':r['expected_class'],'biological_group':r['biological_group'],
        'role':'calibration' if r['canonical_coordinate_match'] else 'consumed_crystal_transfer'} for r in cal['rows']]
    if data['all_cases']!=expected or data['cases']!=[r for r in expected if r['case_id'] not in CANONICAL_REUSE] or data['reuse_case_ids']!=list(CANONICAL_REUSE):
        raise InvalidArtifact('designated canonical/reused source membership changed')
    pilot=read_json(verify(data['pilot_collection']));pm=read_json(verify(pilot['manifest']))
    if pilot['protocol_id']!=POOL_PROTOCOL or pm['settings']!=POOL_SETTINGS:raise InvalidArtifact('reused pilot model policy differs')
    for cid in CANONICAL_REUSE:
        c=next(x for x in pilot['cases'] if x['case_id']==cid)
        row=next(x for x in expected if x['case_id']==cid)
        actual=next(x for x in cal['rows'] if x['case_id']==row['actual_union_case_id'])
        if c['source']['origin_row']!=actual or c['pool']['status']!='available':raise InvalidArtifact('reused canonical pool/origin differs')
    return population,24


def canonical_inputs(calibration,pilot,agreement,output):
    cal=read_json(calibration)
    rows=[{'case_id':r['root_case_id'] if r['canonical_coordinate_match'] else r['case_id'],
        'actual_union_case_id':r['case_id'],'known_class':r['expected_class'],'biological_group':r['biological_group'],
        'role':'calibration' if r['canonical_coordinate_match'] else 'consumed_crystal_transfer'} for r in cal['rows']]
    result={'population':'canonical_remaining24','calibration':record(calibration),'pilot_collection':record(pilot),
        'agreement':record(agreement),'all_cases':rows,'cases':[r for r in rows if r['case_id'] not in CANONICAL_REUSE],
        'reuse_case_ids':list(CANONICAL_REUSE),'canonical_denominator':25,'transfer_denominator':3,
        'new_sources':24,'new_searches':48,'maximum_cross_MACE_calls':48,'maximum_GFN2_calls':192,'new_q0_calls':0}
    write_new(output,result);declared_population(output,calibration);return result


def snapshot(out):
    out.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        q=out/p.name;shutil.copyfile(p,q);pins[p.name]=record(q)
    return pins


def lookup(preparation,crystals,calibration,transfer,inputs):
    prep=read_json(preparation);cr=read_json(crystals);cal=read_json(calibration);trans=read_json(transfer);common=read_json(inputs)
    _,count=declared_population(inputs,calibration)
    rows=[]
    for case in common['cases']:
        cid=case['case_id'];matches=[(pin,r) for pin,col in ((record(calibration),cal),(record(transfer),trans)) for r in col['rows']
            if r['case_id']==cid or (r.get('canonical_coordinate_match') and r.get('root_case_id')==cid)]
        if len(matches)!=1:raise InvalidArtifact('nonunique exact canonical/transfer source alias: '+cid)
        pin,row=matches[0];original=next(p for p in prep['cases']+cr['cases'] if p['case_id']==row['case_id'])
        if row['status']!='complete' or original['status']!='prepared':raise InvalidArtifact('required union source unavailable: '+cid)
        rows.append({'case_id':cid,'actual_union_case_id':row['case_id'],'base':case,'union':original,'collection':pin,'origin_row':row})
    if len(rows)!=count:raise InvalidArtifact('declared population differs')
    return prep,rows


def origin(ep,model):
    """Read actual legacy/native receipts; never fabricate a worker result."""
    r=read_json(verify(ep['native_MACE_receipt']));m=read_json(verify(r['manifest']))
    t=next(t for t in m['tasks'] if t['task_id']==r['task_id'])
    if m['model']!=model or r['status'] not in ('computed','complete') or r['energy_eV']!=ep['native_MACE_energy_eV']:
        raise InvalidArtifact('union origin native method/status/energy differs')
    mult=t.get('multiplicity',t.get('spin_multiplicity'))
    if t['charge']!=ep['charge'] or mult!=ep['multiplicity']:raise InvalidArtifact('union native origin state differs')
    a,b=xyz(verify(t['xyz'])),xyz(verify(ep['xyz']))
    if [x[0] for x in a]!=[x[0] for x in b] or not np.allclose([x[1:] for x in a],[x[1:] for x in b],atol=1e-12,rtol=0):
        raise InvalidArtifact('union native origin coordinates/order differ')
    force=np.load(verify(r['forces']),allow_pickle=False)
    if force.shape!=(len(a),3) or not np.isfinite(force).all():raise InvalidArtifact('archived forces unavailable')
    return r,force


def prepare(preparation,crystals,calibration,transfer,inputs,source,agreement,output):
    from second_shell_context import parent_state
    from coordination_preparation_context import geometry
    prep,rows=lookup(preparation,crystals,calibration,transfer,inputs);sm=read_json(source)
    population,count=declared_population(inputs,calibration)
    if sm['model']!=prep['config']['model'] or sm['settings']!=SETTINGS:raise InvalidArtifact('checkpoint or scaled optimizer differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];cases=[]
    for row in rows:
        cid=row['case_id'];u=row['union'];old=row['origin_row'];context=read_json(verify(old['source_preparation']))
        state=parent_state(u['original_core'],prep['config']['topology'],require_endpoint_receipts=False)
        points={};projections={};new=[]
        for z in ('Ca','La'):
            ep=old['native_endpoints'][z];native,forces=origin(ep,sm['model'])
            atoms=xyz(verify(ep['xyz']));core=xyz(verify(u['original_core']['endpoints'][z]['xyz']))
            geo=geometry(state,context,atoms,core);gp=out/'maps'/(cid+'__'+z+'.json');write_new(gp,geo)
            kin=Kinematics(geo['context']);zeros=np.zeros(len(kin.modes));coords,vectors,raw,normalized,lengths=project(kin.data,zeros,forces)
            low={s:old['solvent_endpoints'][z][s] for s in ('vacuum','alpb')}
            for s,endpoint in low.items():
                from compact_solvation import completed,diagnostics
                pin={k:endpoint[k] for k in ('energy_hartree','output','receipt','manifest','task_id')}
                if completed(verify(pin['manifest']),pin['task_id'])!=pin:raise InvalidArtifact('union q0 solvent receipt differs')
                lm=read_json(verify(pin['manifest']));lt=next(t for t in lm.get('all_tasks',lm['tasks']) if t['task_id']==pin['task_id'])
                if not pool.same_geometry(xyz(verify(lt['xyz'])),atoms) or lt['charge']!=ep['charge']:raise InvalidArtifact('union solvent q0 geometry/state differs')
                audit=diagnostics(pin,lt)
                if audit['charge_sanity_status']!='pass':raise InvalidArtifact('union q0 charge audit failed')
            point={'status':'complete','coordinate':ep['xyz'],'MACE':ep['native_MACE_receipt'],'MACE_eV':native['energy_eV'],
                'forces':native['forces'],'full_q':zeros.tolist(),'active_q_radian':[0.]*4,
                'source_receipt_format':'actual_legacy_native_OMOL','reused_scientific_origin':True}
            points[z]=point;projections[z]=(vectors,raw,normalized)
            new.append({'task_id':cid+'__'+z,'case_id':cid,'metal':z,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity'],
                'mapping':record(gp),'source_preparation':old['source_preparation'],'mode_count':len(kin.modes),'q0_status':'available',
                'q0':{'components':{'MACE_eV':native['energy_eV'],'GFN2_vacuum_hartree':low['vacuum']['energy_hartree'],
                                    'GFN2_ALPB_hartree':low['alpb']['energy_hartree']},'low':low,'native_MACE_receipt':ep['native_MACE_receipt']}})
        ca,la=new
        paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        if read_json(verify(ca['mapping']))!=read_json(verify(la['mapping'])):raise InvalidArtifact('paired physical maps differ')
        kin=Kinematics(read_json(verify(ca['mapping']))['context'])
        choice=preview(kin.modes,projections['Ca'][0],projections['Ca'][2],projections['La'][2]);ids=[r['id'] for r in choice['selected']]
        if len(ids)!=4:raise InvalidArtifact('four common independent modes unavailable')
        selected=[[m['id'] for m in kin.modes].index(x) for x in ids]
        for t in new:
            z=t['metal'];point={**points[z],'gradient_kcal_mol_rad':projections[z][1][selected].tolist()}
            t.update(active_indices=selected,active_mode_ids=ids,active_roles=['adaptive_physical_angular']*4,
                     selector=choice,origin_reuse={'proposal_receipt':points[z]['MACE'],'point':point,
                     'receipt_kind':'actual_archived_native_origin_not_a_proposal'})
            angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))]);tasks.append(t)
        cases.append({**row,'status':'prepared','selection':choice,'atom_count':len(xyz(verify(ca['xyz'])))})
    m={k:sm[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')}
    m.update(protocol_id=PROTOCOL,settings=SETTINGS,agreement=record(agreement),inputs=record(inputs),source=record(source),
        preparation=record(preparation),crystals=record(crystals),calibration=record(calibration),transfer=record(transfer),
        cases=cases,tasks=tasks,declared_case_ids=[r['case_id'] for r in rows],implementation=snapshot(out/'implementation'),
        optimizer_software={'version':scipy.__version__,'wrapper':record(scipy_slsqp.__file__),'kernel':record(scipy_kernel.__file__)},
        population=population,maximum_optimizer_starts=2*count,new_origin_calls=0,maximum_cross_MACE_calls=2*count,maximum_GFN2_calls=8*count,
        production_changed=False,new_DFT_calls=0,reference=None)
    mp=out/'manifest.json';write_new(mp,m);result=validate(mp);write_new(out/'PREFLIGHT.json',result);return result


def validate(manifest):
    m=read_json(manifest)
    population,count=declared_population(verify(m['inputs']),verify(m['calibration']))
    if m.get('population','common8')!=population or m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or len(m['tasks'])!=2*count:raise InvalidArtifact('fixed union/adaptive scope differs')
    for k in ('agreement','inputs','source','preparation','crystals','calibration','transfer','software','orca','cpu_executable','gpu_executable'):verify(m[k])
    for pin in m['implementation'].values():verify(pin)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer runtime differs')
    ids=[c['case_id'] for c in read_json(verify(m['inputs']))['cases']]
    if m['declared_case_ids']!=ids or {(t['case_id'],t['metal']) for t in m['tasks']}!={(c,z) for c in ids for z in ('Ca','La')}:
        raise InvalidArtifact('exact paired declared population required')
    for c in m['cases']:
        data=[]
        for z in ('Ca','La'):
            t=next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(c['case_id'],z))
            ep=c['origin_row']['native_endpoints'][z];native,forces=origin(ep,m['model'])
            low=c['origin_row']['solvent_endpoints'][z]
            expected_components={'MACE_eV':native['energy_eV'],'GFN2_vacuum_hartree':low['vacuum']['energy_hartree'],
                                 'GFN2_ALPB_hartree':low['alpb']['energy_hartree']}
            if t['q0']['components']!=expected_components or t['q0']['low']!=low:raise InvalidArtifact('actual origin components differ')
            if (t['xyz'],t['charge'],t['multiplicity'],t['source_preparation'])!=(ep['xyz'],ep['charge'],ep['multiplicity'],c['origin_row']['source_preparation']):raise InvalidArtifact('union state changed')
            k=Kinematics(read_json(verify(t['mapping']))['context']);point=t['origin_reuse']['point']
            coords,v,raw,normed,_=project(k.data,np.zeros(len(k.modes)),forces);data.append((k,v,normed))
            if (point['MACE']!=ep['native_MACE_receipt'] or point['forces']!=native['forces'] or point['MACE_eV']!=native['energy_eV'] or
                any(point['full_q']) or not np.array_equal(raw[t['active_indices']],point['gradient_kcal_mol_rad'])):raise InvalidArtifact('actual origin force reprojection differs')
            if not np.allclose(coords,[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('mapped q0 differs')
            if [k.modes[i]['id'] for i in t['active_indices']]!=t['active_mode_ids']:raise InvalidArtifact('selected map differs')
            angular.final_geometry(k,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
        choice=preview(data[0][0].modes,data[0][1],data[0][2],data[1][2])
        if choice!=c['selection'] or any(t['selector']!=choice for t in m['tasks'] if t['case_id']==c['case_id']):raise InvalidArtifact('fixed common selector changed')
    return {'status':'prepared','manifest':record(manifest),'cases':count,'optimizer_starts':2*count,'new_q0_calls':0,'molecular_calls_in_validation':0}


def execute(manifest):
    validate(manifest);m=read_json(manifest);root=Path(manifest).parent
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000):raise InvalidArtifact('declared H200 allocation required')
    lock=(root/'proposal.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);started=time.monotonic();gpu=None;error=None;results=[]
    try:
        gpu=WarmGPU(manifest)
        for task in m['tasks']:
            pin=completion.optimize(manifest,task,gpu);results.append(pin)
            print(json.dumps({'task_id':task['task_id'],'status':read_json(verify(pin))['status']}),flush=True)
    except Exception as exc:error=repr(exc)
    finally:
        if gpu is not None:
            try:gpu.close()
            except Exception as exc:error=(error+'; ' if error else '')+repr(exc)
        elapsed=time.monotonic()-started
        write_new(root/('EXECUTION_'+os.environ['SLURM_JOB_ID']+'.json'),{'manifest':record(manifest),'results':results,'error':error,
                  'wall_seconds':elapsed,'allocated_core_seconds':elapsed*32,'allocated_GPU_seconds':elapsed,'job_id':os.environ['SLURM_JOB_ID'],
                  'maximum_optimizer_starts':m['maximum_optimizer_starts'],'new_q0_calls':0,'GPU_command':gpu.command if gpu else None})
    if error:raise InvalidArtifact(error)
    return {'status':'completed','results':results}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);rows=[]
    for t in m['tasks']:
        p=Path(manifest).parent/'proposals'/t['task_id']/'result.json';r=read_json(p) if p.exists() else None
        if r and r['manifest']!=record(manifest):raise InvalidArtifact('proposal source differs')
        ok=bool(r and r['status']=='proposal_available')
        if ok:
            k=Kinematics(read_json(verify(t['mapping']))['context']);point=r['proposal'];q=np.array(point['full_q']);atoms=xyz(verify(point['coordinate']))
            if any(q[i]!=0 for i in set(range(len(q)))-set(t['active_indices'])):raise InvalidArtifact('unselected mode moved')
            angular.final_geometry(k,t,q[t['active_indices']],[a[0] for a in atoms])
            if not np.allclose(k.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('candidate map mismatch')
        rows.append({'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],'status':'candidate_available' if ok else 'unavailable',
                     'proposal_receipt':record(p) if r else None,'candidate':r['proposal'] if ok else None,
                     'reason':r.get('reason') if r else 'not_run','boundary_flag':r.get('boundary_flag') if r else None})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'endpoints':rows,'denominator':len(m['tasks']),'available':sum(r['status']=='candidate_available' for r in rows)}
    write_new(output,result);return result


def prepare_pool(proposals,agreement,output):
    """Only assemble cells; unchanged shared executors compute/collect them."""
    result=read_json(proposals);mp=verify(result['manifest']);validate(mp);m=read_json(mp)
    if result['protocol_id']!=PROTOCOL or len(result['endpoints'])!=len(m['tasks']):raise InvalidArtifact('proposal population differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);cases=[];tasks=[]
    for source in m['cases']:
        cid=source['case_id'];ts={z:next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(cid,z)) for z in ('Ca','La')}
        c={'case_id':cid,'status':'unavailable','reason':None,'source':source,'aliases':{},
           'candidates':[{'id':'origin','xyz':ts['Ca']['xyz']}],'matrix':{z:{'origin':{
               'status':'complete','components':ts[z]['q0']['components'],'xyz':ts[z]['xyz'],
               'MACE':ts[z]['q0']['native_MACE_receipt'],'low':ts[z]['q0']['low'],'reused':True}} for z in ('Ca','La')}}
        pending=[]
        try:
            points={};atoms_by_name={'origin':xyz(verify(ts['Ca']['xyz']))}
            for z in ('Ca','La'):
                r=next(r for r in result['endpoints'] if (r['case_id'],r['metal'])==(cid,z))
                if r['status']!='candidate_available':raise InvalidArtifact('required adaptive '+z+' unavailable: '+str(r.get('reason')))
                actual=read_json(verify(r['proposal_receipt']));point=r['candidate'];q=np.asarray(point['full_q']);k=Kinematics(read_json(verify(ts[z]['mapping']))['context'])
                if actual['manifest']!=record(mp) or actual['status']!='proposal_available' or actual['proposal']!=point:raise InvalidArtifact('candidate receipt differs')
                angular.final_geometry(k,ts[z],q[ts[z]['active_indices']],[a[0] for a in xyz(verify(ts[z]['xyz']))])
                atoms=xyz(verify(point['coordinate']))
                if not np.allclose(k.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('candidate coordinates differ')
                name='adaptive_'+z;match=next((n for n,a in atoms_by_name.items() if pool.same_geometry(a,atoms)),None)
                if match is None:
                    match=name;atoms_by_name[name]=atoms;c['candidates'].append({'id':name,'xyz':point['coordinate']})
                c['aliases'][name]={'representative':match,'source_xyz':point['coordinate'],'proposal_receipt':r['proposal_receipt']};points[z]=point
            for z in ('Ca','La'):
                t=ts[z]
                for candidate in c['candidates'][1:]:
                    name=candidate['id'];atoms=atoms_by_name[name];cross=[(z,*atoms[0][1:]),*atoms[1:]]
                    tid=cid+'__'+z+'__at_'+name;td=out/'cells'/tid;td.mkdir(parents=True);xp=td/'context.xyz'
                    pool.write_xyz(xp,cross)
                    if xyz(xp)!=cross:raise InvalidArtifact('cross coordinate serialization changed')
                    task={'task_id':tid,'case_id':cid,'metal':z,'candidate':name,'xyz':record(xp),
                        'charge':t['charge'],'multiplicity':t['multiplicity'],'source_mapping':t['mapping'],'source_preparation':t['source_preparation']}
                    if name==c['aliases']['adaptive_'+z]['representative'] and xyz(verify(points[z]['coordinate']))==cross:
                        task['native_reuse']=points[z]['MACE'];pool.native_reuse(task,m)
                    pending.append(task);c['matrix'][z][name]={'status':'pending','task_id':tid,'xyz':record(xp),'reused':False}
            c['status']='prepared';tasks.extend(pending)
        except (InvalidArtifact,KeyError,OSError,StopIteration) as exc:c['reason']=str(exc)
        cases.append(c)
    pm={k:m[k] for k in ('model','software','orca','cpu_python','gpu_python','resources')}
    pm.update(protocol_id=POOL_PROTOCOL,settings=POOL_SETTINGS,source=record(proposals),source_manifest=record(mp),agreement=record(agreement),
        declared_case_ids=m['declared_case_ids'],cases=cases,tasks=tasks,implementation=snapshot(out/'implementation'),
        population=m.get('population','common8')+'_union_adaptive_minimal',shard_count=1,new_MACE_cells=sum(not t.get('native_reuse') for t in tasks),
        maximum_new_GFN2_calls=2*len(tasks),new_DFT_calls=0,new_optimizations=0,GFN2_maxiter=500,
        numerical_policy_id='native_GFN2_MaxIter500_unchanged_convergence_v1',reference=None,production_changed=False)
    path=out/'manifest.json';write_new(path,pm);check=validate_pool(path);pool.low_prepare(path);write_new(out/'PREFLIGHT.json',check);return check


def validate_pool(manifest):
    m=read_json(manifest);source=read_json(verify(m['source']));sm=read_json(verify(m['source_manifest']));validate(verify(m['source_manifest']))
    count=len(sm['cases'])
    if (m['protocol_id']!=POOL_PROTOCOL or m['settings']!=POOL_SETTINGS or m['GFN2_maxiter']!=500 or
        source['manifest']!=m['source_manifest'] or m['declared_case_ids']!=sm['declared_case_ids'] or len(m['cases'])!=count or
        m['maximum_new_GFN2_calls']!=2*len(m['tasks']) or len(m['tasks'])>4*count or m['new_MACE_cells']>2*count):
        raise InvalidArtifact('fixed minimal union pool differs')
    for key in ('model','software','orca','cpu_python','gpu_python','resources'):
        if m[key]!=sm[key]:raise InvalidArtifact('union pool scientific runtime differs')
    verify(m['agreement'])
    for pin in m['implementation'].values():verify(pin)
    if len({t['task_id'] for t in m['tasks']})!=len(m['tasks']):raise InvalidArtifact('duplicate cell')
    for c,original in zip(m['cases'],sm['cases']):
        if c['case_id']!=original['case_id'] or c['source']!=original:raise InvalidArtifact('union pool source changed')
        cid=c['case_id'];names=[q['id'] for q in c['candidates']]
        if names[0]!='origin' or len(set(names))!=len(names) or not set(names)<=set(POOL_SETTINGS['candidate_order']):raise InvalidArtifact('candidate membership differs')
        for z in ('Ca','La'):
            t=next(t for t in sm['tasks'] if (t['case_id'],t['metal'])==(cid,z))
            origin=c['matrix'][z]['origin']
            if origin['components']!=t['q0']['components'] or origin['xyz']!=t['xyz'] or not origin['reused'] or origin['low']!=t['q0']['low']:raise InvalidArtifact('original source cell differs')
        if c['status']!='prepared':
            if any(t['case_id']==cid for t in m['tasks']):raise InvalidArtifact('unavailable pool gained cells')
            continue
        if set(c['aliases'])!={'adaptive_Ca','adaptive_La'}:raise InvalidArtifact('required endpoint proposal omitted')
        for z in ('Ca','La'):
            e=next(e for e in source['endpoints'] if (e['case_id'],e['metal'])==(cid,z))
            if e['status']!='candidate_available' or c['aliases']['adaptive_'+z]['proposal_receipt']!=e['proposal_receipt']:raise InvalidArtifact('successful proposal identity differs')
            receipt=read_json(verify(e['proposal_receipt']));point=e['candidate'];t=next(t for t in sm['tasks'] if (t['case_id'],t['metal'])==(cid,z))
            if receipt['proposal']!=point or receipt['manifest']!=m['source_manifest']:raise InvalidArtifact('physical proposal changed')
            q=np.asarray(point['full_q']);kin=Kinematics(read_json(verify(t['mapping']))['context'])
            if any(q[i]!=0 for i in set(range(len(q)))-set(t['active_indices'])):raise InvalidArtifact('unselected mode moved')
            angular.final_geometry(kin,t,q[t['active_indices']],[a[0] for a in xyz(verify(t['xyz']))])
            if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in xyz(verify(point['coordinate']))],atol=1e-12,rtol=0):raise InvalidArtifact('candidate map drift')
            alias=c['aliases']['adaptive_'+z]
            representative=next(v for v in c['candidates'] if v['id']==alias['representative'])
            if alias['source_xyz']!=point['coordinate'] or not pool.same_geometry(xyz(verify(representative['xyz'])),xyz(verify(point['coordinate']))):raise InvalidArtifact('candidate alias not numerical copy')
            if set(c['matrix'][z])!=set(names):raise InvalidArtifact('metals do not share same candidates')
            for name in names[1:]:
                cell=c['matrix'][z][name];task=next(t for t in m['tasks'] if t['task_id']==cell['task_id'])
                coords=xyz(verify(next(v['xyz'] for v in c['candidates'] if v['id']==name)))
                want=[(z,*coords[0][1:]),*coords[1:]]
                if (task['case_id'],task['metal'],task['candidate'],task['source_mapping'],task['source_preparation'])!=(cid,z,name,t['mapping'],t['source_preparation']):raise InvalidArtifact('cross source identity differs')
                if xyz(verify(task['xyz']))!=want or (task['charge'],task['multiplicity'])!=(t['charge'],t['multiplicity']):raise InvalidArtifact('cross geometry/state differs')
                if task.get('native_reuse'):pool.native_reuse(task,m)
    return {'status':'validated','manifest':record(manifest),'denominator':count,'prepared':sum(c['status']=='prepared' for c in m['cases']),
            'new_MACE_cells':m['new_MACE_cells'],'new_GFN2_calls':m['maximum_new_GFN2_calls'],'new_DFT_calls':0}


def compare(collection,union_reference,output):
    from collections import Counter
    from accommodation_nonlinear import contrast_components,relative_components
    from accommodation_folds_compare import decision
    from accommodation_fold_proposals import outcome
    data=read_json(collection);m=read_json(verify(data['manifest']));validate_pool(verify(data['manifest']))
    source=read_json(verify(m['source_manifest']));inputs=read_json(verify(source['inputs']))
    adaptive=read_json(verify(inputs['adaptive_reference']));union=read_json(union_reference)
    if union['collection']!=source['calibration']:raise InvalidArtifact('union reference/calibration source differs')
    bands={'released':adaptive['old_frozen_bands'],'adaptive':adaptive['variants']['operational']['bands'],'union':union['bands']}
    rows=[]
    for c in data['cases']:
        original=c['source'];base=original['base'];old=next(x for x in read_json(verify(base['pool_collection']))['cases'] if x['case_id']==c['case_id'])
        matrix=old['matrix'];names=['origin','adaptive_Ca','adaptive_La']
        old_minimal=pool.choose_rows(matrix,names)
        old_origin=contrast_components(matrix['Ca']['origin']['components'],matrix['La']['origin']['components'])
        union_origin=contrast_components(c['matrix']['Ca']['origin']['components'],c['matrix']['La']['origin']['components'])
        scores={'released_static':old_origin['composite_R_model_kcal_mol'],
            'original_adaptive_full':old['pool']['operational']['composite_R_model_kcal_mol'],
            'original_adaptive_minimal':old_minimal['operational']['composite_R_model_kcal_mol'] if old_minimal['status']=='available' else None,
            'union_static':union_origin['composite_R_model_kcal_mol'],
            'union_adaptive_minimal':c['pool']['operational']['composite_R_model_kcal_mol'] if c['pool']['status']=='available' else None}
        calls={ref:{name:{'decision':decision(value,limits),'outcome':outcome(decision(value,limits),base['known_class'])}
                    for name,value in scores.items()} for ref,limits in bands.items()}
        works=None
        if c['pool']['status']=='available':
            works={z:relative_components(c['matrix'][z][c['pool']['rows'][z]['operational_candidate']]['components'],c['matrix'][z]['origin']['components']) for z in ('Ca','La')}
            actual=works['Ca']['composite_kcal_mol']-works['La']['composite_kcal_mol']
            if abs(actual-(scores['union_adaptive_minimal']-scores['union_static']))>1e-6:raise InvalidArtifact('selected work sign/conversion differs')
        rows.append({'case_id':c['case_id'],'actual_union_case_id':original['actual_union_case_id'],
            'expected_class_for_reporting':base['known_class'],'biological_group':original['origin_row']['biological_group'],
            'raw_R_model_kcal_mol':scores,'old_band_transfer_only':calls,'union_selected_endpoint_work':works,
            'union_static_components':union_origin,'released_static_components':old_origin,'union_pool':c['pool'],
            'original_minimal_pool':old_minimal,'original_full_pool':old['pool'],'source_base_collection':base['pool_collection'],
            'pool_failure_reason':c['pool'].get('reason'),'new_reference':None})
    variants=list(rows[0]['raw_R_model_kcal_mol']);counts={};order={};spreads=[]
    for ref in bands:
        counts[ref]={name:dict(Counter(r['old_band_transfer_only'][ref][name]['outcome'] for r in rows)) for name in variants}
    for name in variants:
        ca=[r['raw_R_model_kcal_mol'][name] for r in rows if r['expected_class_for_reporting']=='Ca']
        la=[r['raw_R_model_kcal_mol'][name] for r in rows if r['expected_class_for_reporting']=='La']
        complete=all(v is not None for v in ca+la)
        order[name]={'complete':complete,'minimum_La_minus_maximum_Ca_model_kcal_mol':min(la)-max(ca) if complete else None,
            'correct_pairwise_directions':sum(a>b for a in la for b in ca) if complete else None,
            'pairwise_denominator':len(la)*len(ca),'biological_independence_claimed':False}
    for group in dict.fromkeys(r['biological_group'] for r in rows):
        members=[r for r in rows if r['biological_group']==group]
        if len(members)<2:continue
        spans={}
        for name in variants:
            vals=[r['raw_R_model_kcal_mol'][name] for r in members]
            spans[name]=max(vals)-min(vals) if all(v is not None for v in vals) else None
        spreads.append({'biological_group':group,'declared_sources':[r['case_id'] for r in members],'range_model_kcal_mol':spans})
    result={'protocol_id':POOL_PROTOCOL,'collection':record(collection),'analysis_implementation':record(__file__),
        'union_reference':record(union_reference),'adaptive_reference':inputs['adaptive_reference'],'frozen_bands':bands,
        'old_band_transfer_only':True,'new_calibration':None,'denominator':8,'rows':rows,'counts':counts,'raw_ordering':order,
        'prescribed_pair_spreads':spreads,'new_calls_in_comparison':0,'production_changed':False}
    write_new(output,result);return result


def canonical_reference(collection,inputs,agreement,output):
    """Freeze a distinct canonical25 reference; all three crystals stay outside."""
    from datetime import datetime,timezone
    from nikasha_pool_compare import extrema_reference
    from accommodation_folds_compare import decision
    from accommodation_fold_proposals import outcome
    selection=read_json(inputs);declared_population(inputs,verify(selection['calibration']))
    fresh=read_json(collection);fm=read_json(verify(fresh['manifest']));validate_pool(verify(fresh['manifest']))
    parent=read_json(verify(fm['source_manifest']))
    if parent['inputs']!=record(inputs) or parent['agreement']!=record(agreement) or fm['agreement']!=record(agreement):raise InvalidArtifact('canonical agreement/input differs')
    pilot=read_json(verify(selection['pilot_collection']));pm=read_json(verify(pilot['manifest']));validate_pool(verify(pilot['manifest']))
    if any(fm[k]!=pm[k] for k in ('protocol_id','settings','model','software','orca','GFN2_maxiter')):raise InvalidArtifact('reused pool method differs')
    if [c['case_id'] for c in fresh['cases']]!=[c['case_id'] for c in selection['cases']]:raise InvalidArtifact('fresh canonical denominator differs')
    rows=[]
    for src in selection['all_cases']:
        reused=src['case_id'] in CANONICAL_REUSE;col=pilot if reused else fresh
        c=next(c for c in col['cases'] if c['case_id']==src['case_id'])
        if c['source']['actual_union_case_id']!=src['actual_union_case_id'] or c['source']['origin_row']['expected_class']!=src['known_class']:raise InvalidArtifact('canonical source/state/label differs')
        # Recompute selections from actual component cells; an absent cell stays unavailable.
        if c['status']=='prepared':
            chosen=pool.choose_rows(c['matrix'],[q['id'] for q in c['candidates']])
            if chosen!=c['pool']:raise InvalidArtifact('collected pool differs from actual matrix')
        else:
            chosen=c['pool']
            if chosen['status']!='unavailable' or any(chosen[v] is not None for v in ('mathematical','operational')):
                raise InvalidArtifact('unsupported preparation was substituted with a score')
        scores={v:chosen[v]['composite_R_model_kcal_mol'] if chosen['status']=='available' else None for v in ('mathematical','operational')}
        rows.append({**src,'expected_class':src['known_class'],'status':chosen['status'],'scores':scores,
            'selected_candidates':{v:{z:chosen['rows'][z][v+'_candidate'] for z in ('Ca','La')} for v in scores} if chosen['status']=='available' else None,
            'pool':chosen,'collection':selection['pilot_collection'] if reused else record(collection),'reused_common8':reused})
    variants={}
    for variant in ('mathematical','operational'):
        canonical=[{**r,'R_model_kcal_mol':r['scores'][variant]} for r in rows if r['role']=='calibration']
        variants[variant]=extrema_reference(canonical,variant,'Nikasha_union_adaptive_minimal_canonical25_v1')
        values={z:[r['scores'][variant] for r in rows if r['role']=='calibration' and r['known_class']==z] for z in ('Ca','La')}
        if all(v is not None for vals in values.values() for v in vals):
            variants[variant]['class_spread']={z:max(v)-min(v) for z,v in values.items()}
        for r in rows:
            r.setdefault('own_reference',{})[variant]={'decision':decision(r['scores'][variant],variants[variant]['bands']) if variants[variant]['bands'] else 'unavailable',
                                                     'reference_status':variants[variant]['status']}
            r['own_reference'][variant]['outcome']=outcome(r['own_reference'][variant]['decision'],r['known_class'])
    result={'protocol_id':POOL_PROTOCOL,'reference_id':'Nikasha_union_adaptive_minimal_canonical25_v1','variants':variants,
        'inputs':record(inputs),'agreement':record(agreement),'collections':[selection['pilot_collection'],record(collection)],
        'denominator':28,'complete':sum(r['status']=='available' for r in rows),'rows':rows,
        'calibration_denominator':25,'crystal_transfer_denominator':3,'crystals_or_noncanonical_used_for_fit':False,
        'frozen_UTC':datetime.now(timezone.utc).isoformat(),'new_calls_in_analysis':0,'production_changed':False,
        'implementation':record(__file__),'model':fm['model'],'settings':fm['settings'],'optimizer_settings':parent['settings']}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('canonical_inputs')
    for name in ('calibration','pilot','agreement','output'):q.add_argument('--'+name,required=True)
    q=s.add_parser('canonical_reference')
    for name in ('collection','inputs','agreement','output'):q.add_argument('--'+name,required=True)
    q=s.add_parser('prepare')
    for name in ('preparation','crystals','calibration','transfer','inputs','source','agreement','output'):q.add_argument('--'+name,required=True)
    for name in ('validate','execute','collect'):
        q=s.add_parser(name);q.add_argument('--manifest',required=True)
        if name=='collect':q.add_argument('--output',required=True)
    q=s.add_parser('prepare_pool')
    for name in ('proposals','agreement','output'):q.add_argument('--'+name,required=True)
    q=s.add_parser('validate_pool');q.add_argument('--manifest',required=True)
    q=s.add_parser('compare')
    for name in ('collection','union_reference','output'):q.add_argument('--'+name,required=True)
    a=vars(p.parse_args());op=a.pop('op');result=globals()[op](**a)
    print(json.dumps({k:v for k,v in result.items() if k not in ('endpoints','results')},indent=2))
if __name__=='__main__':main()
