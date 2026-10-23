#!/usr/bin/env python3
"""Saved-only terminal review of the ten selected joint-water cells."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

GEOMETRIES = ('water_minus', 'water_central', 'water_plus', 'aquo8_central', 'aquo9_central')
CELLS = tuple(q+'_'+m for q in GEOMETRIES for m in ('La','Dy'))

def require(ok, message):
    if not ok:
        raise ValueError(message)

def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def verify(record):
    require(pin(record['path'])==record,'Changed pin: '+record['path'])
    return Path(record['path'])

def read(record):
    return json.loads(verify(record).read_text())

def module(record,name):
    spec=importlib.util.spec_from_file_location(name,verify(record))
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def close(a,b):
    return abs(a-b)<=1e-10*max(1.,abs(a),abs(b))

def native_authority(gid,value):
    """Join the producer's saved input, receipts and arrays to this request."""
    pv=value['provenance'];authority=read(pv['native_authority'])
    require(authority['status']=='complete','Native producer is unavailable: '+gid)
    if 'rows' in authority:
        require(pv['native_source_kind']=='17670_W9' and gid.startswith('aquo9_'),'Unexpected aggregate producer')
        row=authority['rows'][gid];native_value=row['model']
    else:
        require(authority['cell_id']==gid,'Native producer cell differs')
        row=authority;native_value=authority
    require(row['status']=='complete','Native producer row unavailable')
    for k in ('energy','forces','native'):
        require(native_value[k]==value[k],'Producer '+k+' differs: '+gid)
    inp=pv['effective_input'];receipts=row['receipts']
    require(len(receipts)==2,'Expected two native programs')
    for key in ('xyz','key','parameters'):
        require(verify(row['input'][key]) and row['input'][key]['sha256']==inp[key]['sha256'],'Producer input differs')
    for receipt,output_key,program in zip(receipts,('analyze_output','gradient_output'),('analyze','testgrad')):
        require(receipt['returncode']==0 and receipt['output']==value['native'][output_key],'Native receipt/output differs')
        verify(receipt['output'])
        require(Path(receipt['argv'][0]).name==program and receipt['argv'][1]=='system.xyz','Native invocation differs')
        require(all(receipt['threads'].get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'Native thread policy differs')
        for key,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
            require(pin(Path(receipt['cwd'])/name)==inp[key],'Runtime input differs from selected input')
    return {'authority':pv['native_authority'],'receipts':receipts}

def reference_replay(gid,row):
    """Recollect saved ORCA output/gradient using its frozen producer adapter."""
    rr=read(row['reference_result'])
    require(rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values()),'Reference unavailable')
    for k in ('output','engrad','input','geometry','panel','adapter'):verify(rr[k])
    adapter=module(rr['adapter'],'saved_reference_'+gid)
    panel=adapter.load_panel(verify(rr['panel']))
    cells=[c for c in panel['cells'] if c['cell_id']==rr['cell_id']]
    require(len(cells)==1,'Reference cell not unique')
    cell=cells[0]
    require(cell['geometry']['geometry_id']==gid and cell['geometry']['positions_A']==row['geometry']['positions_A'] and cell['geometry']['symbols']==row['geometry']['symbols'],'Reference physical geometry differs')
    raw=adapter.collect_one(panel,cell,verify(rr['output']).parent,rr['returncode'])
    require(all(rr[k]==v for k,v in raw.items()),'Raw reference collection differs: '+gid)
    require(raw['forces_kcal_mol_A']==row['reference']['forces'] and close(raw['energy_hartree']*627.509474,row['reference']['energy']),'Reference arrays differ')
    return {'result':row['reference_result'],'raw_pins':{k:rr[k] for k in ('output','engrad','input','geometry','panel','adapter')},'checks':raw['checks']}

def independent_score(data,model):
    """Direct Cartesian sums and explicit weights; no fitting/body helper calls."""
    ref={g:r['reference'] for g,r in data['rows'].items()}
    def delta(values,q):return values[q+'_Dy']['energy']-values[q+'_La']['energy']
    exchange={q:(delta(model,q)-delta(model,'aquo8_central'))-(delta(ref,q)-delta(ref,'aquo8_central'))
              for q in ('water_minus','water_central','water_plus','aquo9_central')}
    work={m+'_'+s:(model['water_'+s+'_'+m]['energy']-model['water_central_'+m]['energy'])-
          (ref['water_'+s+'_'+m]['energy']-ref['water_central_'+m]['energy']) for m in ('La','Dy') for s in ('minus','plus')}
    paired={s:work['Dy_'+s]-work['La_'+s] for s in ('minus','plus')}
    vector=[x*math.sqrt(.25*.5/4) for x in list(exchange.values())+list(work.values())]
    checks={}
    for kind in ('La','Dy','paired'):
        checks[kind]={}
        for q in GEOMETRIES:
            row=data['rows'][q+'_La'];xyz=row['geometry']['positions_A'];n=len(xyz);nw=(n-1)//3
            require(row['geometry']['metal_index']==0 and n in (4,25,28),'Unexpected water graph')
            expected=[{'indices':[0],'category':'ion'}]+[{'indices':list(range(i,i+3)),'category':'water'} for i in range(1,n,3)]
            require(row['bodies']==expected,'Unexpected body order')
            def errors(m):return [[a-b for a,b in zip(u,v)] for u,v in zip(model[q+'_'+m]['forces'],ref[q+'_'+m]['forces'])]
            if kind=='paired':
                err=[[d-l for d,l in zip(u,v)] for u,v in zip(errors('Dy'),errors('La'))]
            else:err=errors(kind)
            raw=err[0][:];weighted=[v/(2*math.sqrt(2*3)) for v in raw]
            for first in range(1,n,3):
                ids=range(first,first+3);center=[sum(xyz[i][k] for i in ids)/3 for k in range(3)]
                net=[sum(err[i][k] for i in ids) for k in range(3)];torque=[0.,0.,0.]
                for i in ids:
                    x,y,z=[xyz[i][k]-center[k] for k in range(3)];fx,fy,fz=err[i]
                    t=(y*fz-z*fy,z*fx-x*fz,x*fy-y*fx)
                    torque=[a+b for a,b in zip(torque,t)]
                raw.extend(net+torque);weighted.extend(v/(2*math.sqrt(2*nw*6)) for v in net+torque)
            rms=2*math.sqrt(sum(v*v for v in weighted));mx=max(abs(v) for v in raw)
            checks[kind][q]={'RMS':rms,'maximum':mx,'passed':rms<=2 and mx<=5}
            vector.extend(v*math.sqrt(.25/5) for v in weighted)
    require(len(vector)==413 and all(math.isfinite(x) for x in vector),'Invalid independent residual vector')
    energies={'balanced_exchange':exchange,'individual_radial_work':work,'paired_radial_work':paired}
    epass=all(abs(v)<=1 for group in energies.values() for v in group.values())
    return {'vector':vector,'objective':sum(v*v for v in vector),'energy_errors':energies,'force_checks':checks,
            'passed':epass and all(v['passed'] for group in checks.values() for v in group.values())}

def review(result_path,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    report={'status':'unavailable','implementation':pin(__file__),'result':pin(result_path),'new_native_calls':0,
            'complete':0,'error':None,'all_development_controls_passed':False,'all_selected_raw_replays_exact':False,
            'exact_selected_score_replayed':False,'physical_accuracy_qualified':False}
    try:
        require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'One numerical thread required')
        result=read(report['result']);cfg=read(result['config']);release=read(cfg['root_release'])
        require(result['job_id']=='17676' and result['stage']=='run' and result['status']=='bounded_joint_coordination_water_complete' and result['error'] is None,'No successful terminal fit application')
        require(result['candidate_frozen'] and result['original_results_unchanged'] and result['new_native_program_invocations']<=5120 and result['new_parameter_vectors']<=384 and result['completed_outer_steps']<=24,'Finite terminal scope differs')
        require(result['implementation']==cfg['implementation']==release['six_implementation'],'Driver/release binding differs')
        for rec in cfg.values():
            if isinstance(rec,dict) and set(rec)=={'path','sha256'}:verify(rec)
        ecfg=read(cfg['evaluator_config']);data=read(cfg['training_result']);model=result['selected']['model'];score=result['selected']['score']
        require(set(data['rows'])==set(model)==set(CELLS),'Selected ten-cell inventory differs')
        values=result['selected']['parameters'];require(len(values)==6,'Wrong parameter count')
        require(all(p['lower']<=v<=p['upper'] for p,v in zip(result['family']['parameters'],values)),'Selected parameters outside bounds')
        parser=module(ecfg['native_helper'],'terminal_native_parser');resolver=module(ecfg['resolver'],'terminal_pair_resolver')
        evidence={}
        for gid in CELLS:
            row=data['rows'][gid];value=model[gid];pv=value['provenance'];inp=pv['effective_input'];metal=row['input']['metal'];mt=500 if metal=='La' else 501;i=0 if metal=='La' else 2
            require(pv['requested_vector']==values,'Selected requested vector differs')
            producer=native_authority(gid,value)
            for rec in inp.values():verify(rec)
            for k in ('xyz','parameters'):require(inp[k]['sha256']==row['input'][k]['sha256'],'Physical source differs: '+gid+'/'+k)
            original=resolver.parse_parameters(verify(row['input']['parameters']).read_text(),verify(row['input']['key']).read_text())
            actual=resolver.parse_parameters(verify(inp['parameters']).read_text(),verify(inp['key']).read_text())
            pair=resolver.resolve_pair(actual,mt,36)
            require([pair['Rmin_pair_A'],pair['thole'],pair['epsilon_kcal_mol']]==values[i:i+2]+[values[4 if metal=='La' else 5]],'Selected effective parameters differ')
            for k in ('settings','atom','vdw','polarize','metal_multipoles'):require(original[k]==actual[k],'Nonpair setting changed')
            for k,allowed in (('vdwpr',(34,mt)),('polpair',(36,mt))):
                require({a:b for a,b in original[k].items() if a!=allowed}=={a:b for a,b in actual[k].items() if a!=allowed},'Nonwater pair changed')
            require(resolver.resolve_pair(original,mt,37)==resolver.resolve_pair(actual,mt,37),'Water hydrogen pair changed')
            native=value['native'];parsed=parser.parse(verify(native['analyze_output']),verify(native['gradient_output']),dict(row['input'],derivative_test=False))
            require(parsed==native and value['energy']==parsed['total_energy_kcal_mol'] and value['forces']==parsed['force_kcal_mol_A'],'Selected raw replay differs: '+gid)
            reference=reference_replay(gid,row)
            evidence[gid]={'input':inp,'native_authority':pv['native_authority'],'native_outputs':{k:native[k] for k in ('analyze_output','gradient_output')},'reference_result':row['reference_result'],'effective_pair':pair,'producer_join':producer,'raw_reference_replay':reference}
        obj=module(cfg['objective'],'terminal_objective');body=module(cfg['body_helper'],'terminal_body')
        require(obj.evaluate(data,model,body)==score,'Exact selected score replay differs')
        manual=independent_score(data,model)
        require(all(close(a,b) for a,b in zip(manual['vector'],score['weighted_residuals'])) and close(manual['objective'],score['objective']),'Independent413 arithmetic differs')
        require(manual['passed']==score['development_checks_passed']==result['development_checks_passed'],'Scientific acceptance differs')
        for kind in manual['force_checks']:
            for q,v in manual['force_checks'][kind].items():
                old=score['per_fixture'][kind][q]
                require(close(v['RMS'],old['rigid_RMS_kcal_mol_A']) and close(v['maximum'],old['maximum_body_component_error_kcal_mol_A']) and v['passed']==old['passed'],'Independent individual force gate differs')
        for group,errs in manual['energy_errors'].items():
            for q,v in errs.items():require(close(v,score['energy_errors_kcal_mol'][group][q]),'Independent energy control differs')
        report.update(status='selected_ten_raw_review_passed',selected_parameters=values,selected_cell_ids=list(CELLS),complete=10,
                      all_selected_raw_replays_exact=True,exact_selected_score_replayed=True,all_development_controls_passed=manual['passed'],
                      independent_arithmetic=manual,evidence=evidence,configuration=result['config'],release=cfg['root_release'],
                      stopping_reason=result['stopping_reason'],native_programs=result['new_native_program_invocations'],
                      last_Jacobian_origin_key=result['last_Jacobian_origin_key'],historical17670_review_loss_preserved=True)
    except Exception as exc:
        report.update(status='review_failed_or_unavailable',error=repr(exc))
    report.update(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    (out/'REVIEW.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('status','complete','all_development_controls_passed','error')}))
    return 0 if report['status']=='selected_ten_raw_review_passed' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();raise SystemExit(review(a.result,a.out))
