#!/usr/bin/env python3
"""Saved-only replay and component attribution of exactly steps92..172."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import time
import traceback

PROTOCOL='lbt_7ccn_trace_audit_v1'
PARTS=('Gelec','Ghb','Gsasa','Gshift')
COMPONENTS=('MACE_kcal_mol','fixed_ALPB_kcal_mol','composite_kcal_mol')
TOLERANCES={'coordinates_A_atol':1e-10,'derived_atol':1e-10,'derived_rtol':1e-10,'work_closure_kcal_mol_atol':1e-10}
FIXED={'run_result':'9159ea43cf76692fb9ad40acbf199c1a6b146acdf8dc21ff68a045350a764a10',
       'helper':'551d964badc41141cfc53a85d704135ec74c21f185e12b6cc50d578dc47910d2',
       'optimizer':'c5092f2480fbf393d6ca0d67733ebf2131580c922a62a67925dea8945eb5acb9',
       'initial_checkpoint':'317c2b24eff75a2af4a6fb7066f1754d1df1b9a284405a8aad2ea3b20d5c74cd'}


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve()
    with p.open('rb') as h: sha=hashlib.file_digest(h,'sha256').hexdigest()
    return {'path':str(p),'sha256':sha}


def verify(record):
    require(pin(record['path'])['sha256']==record['sha256'],'Changed pin: '+record['path']); return Path(record['path'])


def read(record): return json.loads(verify(record).read_text())


def write(path,value):
    with Path(path).open('x') as h: json.dump(value,h,indent=2,allow_nan=False); h.write('\n')
    return pin(path)


def module(record):
    spec=importlib.util.spec_from_file_location('saved_trace_completion',verify(record))
    obj=importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj


def inputs(cfg,out,result):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(cfg['protocol']==PROTOCOL and pin(__file__)['sha256']==cfg['implementation']['sha256']
            and cfg['tolerances']==TOLERANCES,'Frozen audit contract differs')
    for k,sha in FIXED.items(): require(cfg[k]['sha256']==sha,'Frozen input differs: '+k)
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    science=read(cfg['science_decision']); run=read(cfg['run_result']); rcfg=read(run['config']); search=read(cfg['search'])
    require(science['protocol']==PROTOCOL and science['status']=='declared_saved_only_mechanics_diagnosis'
            and all(cfg[k]==v for k,v in science['inputs'].items()),'Root fixed scope differs')
    require(run['status']=='complete' and run['stage']=='run' and run['slurm_job_id']=='15458'
            and run['implementation']==cfg['helper'] and rcfg['optimizer']==cfg['optimizer']
            and run['rows'][0]['search']==cfg['search'] and run['rows'][0]['checkpoint']==cfg['checkpoint']
            and run['rows'][0]['terminal_point']==cfg['terminal'] and run['rows'][0]['physical_state_admitted']
            and search['termination']=='accepted_step_limit' and not search['converged']
            and search['accepted_steps']==search['distinct_evaluations']==172 and len(search['trace'])==81
            and run['new_composite_evaluations']==run['new_MACE_calls']==run['new_ALPB_calls']==80
            and len(run['requests'])==80 and len(search['evaluated_points'])==80,'Exact completed80-transition trace required')
    h=module(cfg['helper']); inherited={}
    loaded,core,ev,initial,start,cp,evidence=h.prepare(rcfg,out,inherited)
    pre=read(run['preflight_result']); decision=read(run['run_decision'])
    require(run['preflight_result']==rcfg['preflight_result'] and run['run_decision']==rcfg['run_decision']
            and pre['status']=='preflight_passed' and run['preflight_evidence']==pre['preflight_evidence']==evidence
            and run['scientific_binding']==h.binding(rcfg)==decision['scientific_binding']
            and decision['status']=='approved_single_completion_attempt','Actual saved execution chain differs')
    require(cp==read(cfg['initial_checkpoint']) and cp['point']==search['trace'][0]['point']
            and search['start_point']==cp['initial_point'] and search['original_source']==evidence['original_source']
            and cfg['points']==[r['point'] for r in search['trace']]
            and len(cfg['checkpoints'])==81 and all(r['status']=='complete' for r in run['requests'])
            and search['evaluated_points']==[r['point'] for r in run['requests']]==cfg['points'][1:],
            'Point/request/checkpoint chain differs')
    require(ev.base.EV==23.06054783061903 and ev.base.HA==627.509474 and ev.base.BOHR==0.529177210903,
            'Original conversion constants differ')
    require(ev.engine.native is None and not ev.evaluation_requests and all(v==0 for v in ev.engine.cost.values()),'Saved helper made a call')
    def forbidden(*args,**kwargs): raise RuntimeError('Molecular evaluation forbidden in saved audit')
    ev.evaluate=forbidden; ev.engine.evaluate=forbidden
    result['input_evidence']={'run':cfg['run_result'],'initial_checkpoint':cfg['initial_checkpoint'],
        'final_checkpoint':cfg['checkpoint'],'initial_point':cp['initial_point'],'source':evidence['original_source'],
        'geometry_helper':evidence['geometry_helper'],'recipe':evidence['recipe'],'versions':evidence['versions'],
        'optimizer':cfg['optimizer'],'settings':core.SETTINGS,'tolerances':TOLERANCES}
    return loaded,core,ev,search,run,cp


def check_close(result,np,label,a,b,atol=None,rtol=None):
    a=np.asarray(a,dtype=float); b=np.asarray(b,dtype=float)
    atol=TOLERANCES['derived_atol'] if atol is None else atol
    rtol=TOLERANCES['derived_rtol'] if rtol is None else rtol
    require(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),'Invalid replay arrays: '+label)
    error=np.abs(a-b); maximum=float(np.max(error)) if error.size else 0.
    index=list(np.unravel_index(int(np.argmax(error)),error.shape)) if error.size else []
    entry={'label':label,'exact':bool(np.array_equal(a,b)),'max_abs_residual':maximum,
           'max_index':[int(i) for i in index],'atol':atol,'rtol':rtol,'pass':bool(np.all(error<=atol+rtol*np.abs(b)))}
    result['replay_checks'].append(entry); require(entry['pass'],'Replay mismatch: '+label)


def metric_replay(result,np,label,actual,saved):
    require(set(actual)==set(saved),'Metric keys differ: '+label)
    for k,v in actual.items():
        if isinstance(v,bool): require(v is saved[k],'Force convergence flag differs: '+label)
        else: check_close(result,np,label+'/'+k,v,saved[k])


def gradient_summary(np,core,parts):
    gm,gs,g=(parts[k] for k in COMPONENTS); nm=float(np.linalg.norm(gm)); ns=float(np.linalg.norm(gs)); ng=float(np.linalg.norm(g))
    dot=float(np.sum(gm*gs))
    return {'metrics':{k:core.force_metrics(v) for k,v in parts.items()},'MACE_dot_ALPB':dot,
        'dot_units':'(kcal/mol/Angstrom)^2','MACE_ALPB_cosine':dot/(nm*ns) if nm*ns else None,
        'total_L2_over_sum_component_L2':ng/(nm+ns) if nm+ns else None,
        'projections_on_negative_total_unit_gradient_kcal_mol_A':{k:float(np.sum(v*(-g/ng))) if ng else None for k,v in parts.items()}}


def energy_work(ev,a,b,pa,pb):
    value=ev.difference(b,a)
    alpb={k:(pb['native_energies']['ALPB_parts_Eh'][k]-pa['native_energies']['ALPB_parts_Eh'][k])*ev.base.HA for k in PARTS}
    return {'components_kcal_mol':value,'ALPB_parts_kcal_mol':alpb,
            'ALPB_parts_closure_residual_kcal_mol':math.fsum(alpb.values())-value['fixed_ALPB_kcal_mol']}


def attribution(np,core,atoms,parts):
    groups={}; rows=[]
    for i,a in enumerate(atoms):
        residue='metal' if a['id']=='metal' else '/'.join(a['id'].split('/')[:3])
        groups.setdefault(residue,[]).append(i)
        rows.append({'index':i,'physical_id':a['id'],'residue':residue,'resname':a.get('resname','La'),
            'force_kcal_mol_A':{k:(-v[i]).tolist() for k,v in parts.items()},
            'force_magnitudes_kcal_mol_A':{k:float(np.linalg.norm(v[i])) for k,v in parts.items()},
            'MACE_dot_ALPB':float(np.dot(parts[COMPONENTS[0]][i],parts[COMPONENTS[1]][i]))})
    residues=[{'residue':name,'indices':inds,'resname':atoms[inds[0]].get('resname','La'),
               'component_force_metrics':{k:core.force_metrics(v[inds]) for k,v in parts.items()}}
              for name,inds in groups.items()]
    return {'atoms':rows,'residues':residues,'residue_force_norms_are_not_energies':True,
        'atom_ranking_by_total_force':[r['physical_id'] for r in sorted(rows,key=lambda r:r['force_magnitudes_kcal_mol_A'][COMPONENTS[2]],reverse=True)]}


def audit(cfg,out,result):
    loaded,core,ev,search,run,initialcp=inputs(cfg,out,result)
    import numpy as np
    trace=search['trace']; points=[]; raw=[]; gradients=[]; checkpoints=[]; projected=[]
    for i,entry in enumerate(trace):
        step=92+i; record=cfg['points'][i]; p=read(record); point=ev.read_point(record); parts=ev.gradients(record)
        cp=read(cfg['checkpoints'][i]); checkpoints.append(cp)
        require(entry['accepted_steps']==entry['distinct_evaluations']==step and cp['point']==record
                and cp['accepted_steps']==cp['distinct_evaluations']==step and cp['limits']==search['limits']
                and cp['settings']==core.SETTINGS and cp['initial_point']==initialcp['initial_point']
                and cp['evaluation_requests']==initialcp['evaluation_requests']+run['requests'][:i]
                and cp['evaluated_points']==initialcp['evaluated_points']+search['evaluated_points'][:i]
                and cp['failure'] is None and entry['geometry']['pass'] and entry['geometry']['state_compatible'],
                'Recorded state/checkpoint chain differs at '+str(step))
        if i==0:
            require(cp['history']==initialcp['history'] and cp['trust_A']==initialcp['trust_A'],'Initial history/trust changed')
        require(p['recipe']==ev.read_point(cfg['points'][0])['physical_state'],'Point recipe differs')
        x=point['positions_A']; g=point['gradient_kcal_mol_A']; pg,rigid=core.remove_rigid(x,g)
        require(np.array_equal(parts[COMPONENTS[0]]+parts[COMPONENTS[1]],g),'Exact raw component/full gradient closure differs')
        metric_replay(result,np,str(step)+'/raw_force',core.force_metrics(g),entry['raw_full_force'])
        metric_replay(result,np,str(step)+'/projected_force',core.force_metrics(pg),entry['rigid_removed_force'])
        require(rigid['rigid_rank']==entry['rigid_projection']['rigid_rank'],'Rigid rank differs')
        for k in ('removed_L2','normalized_Gram_eigenvalues'):
            check_close(result,np,str(step)+'/rigid/'+k,rigid[k],entry['rigid_projection'][k])
        require(cp['trust_A']==entry['trust_A'],'Recorded trust differs')
        points.append(point); raw.append(p); gradients.append(parts); projected.append(pg)
        result['points'][i]={'step':step,'status':'available','point':record,'positions':p['positions'],
            'gradient':p['gradient_kcal_mol_A'],'MACE':p['MACE'],'ALPB':p['ALPB'],'checkpoint':cfg['checkpoints'][i],
            'exact_component_gradient_closure':True,'native_energies':p['native_energies'],
            'raw_force':core.force_metrics(g),'projected_force':core.force_metrics(pg),'rigid_projection':rigid,
            'trust_A':cp['trust_A'],'history_pairs':len(cp['history']),'component_gradients':gradient_summary(np,core,parts)}
    final=read(cfg['checkpoint'])
    require({k:v for k,v in final.items() if k!='trial_trace'}=={k:v for k,v in checkpoints[-1].items() if k!='trial_trace'}
            and final['trial_trace']==trace and cfg['points'][-1]==cfg['terminal'],'Final checkpoint numeric state differs')
    for i,entry in enumerate(trace[:-1]):
        step=92+i; cp=checkpoints[i]; nxt=checkpoints[i+1]; point=points[i]; other=points[i+1]
        x=point['positions_A']; g=point['gradient_kcal_mol_A']; pg=projected[i]
        history=[(np.asarray(h['s']),np.asarray(h['y']),h['rho']) for h in cp['history']]
        h0=float(np.dot(history[-1][0],history[-1][1])/np.dot(history[-1][1],history[-1][1])) if history else 1.
        d=loaded['world']['lbfgs'](pg.ravel(),history).reshape(x.shape); d,_=core.remove_rigid(x,d)
        before_fallback=float(np.sum(pg*d)); fallback=before_fallback>=-1e-14
        if fallback: history=[]; d=-pg
        maximum=float(np.linalg.norm(d,axis=1).max()); require(maximum>=1e-14,'Unexpected zero direction in accepted trace')
        scale=min(1.,cp['trust_A']/maximum); d*=scale; slope=float(np.sum(g*d)); require(slope<0,'Non-descent replay')
        require(len(entry['trials'])==1 and entry['trials'][0]['accepted'] and entry['trials'][0]['evaluation_requested'],
                'Actual single-trial denominator differs')
        trial=entry['trials'][0]; multiplier=.5**0
        require(trial['multiplier']==multiplier and trial['point']==cfg['points'][i+1],'Accepted transition identity differs')
        check_close(result,np,str(step)+'/predicted_positions',x+multiplier*d,other['positions_A'],TOLERANCES['coordinates_A_atol'],0.)
        check_close(result,np,str(step)+'/Armijo_limit',core.SETTINGS['armijo']*multiplier*slope,trial['Armijo_limit_kcal_mol'])
        work=energy_work(ev,point,other,raw[i],raw[i+1])
        require(work['components_kcal_mol']==trial['work'],'Exact native component work differs')
        require(work['components_kcal_mol']['composite_kcal_mol']<=trial['Armijo_limit_kcal_mol'],'Saved Armijo acceptance differs')
        check_close(result,np,str(step)+'/ALPB_work_closure',work['ALPB_parts_closure_residual_kcal_mol'],0.,TOLERANCES['work_closure_kcal_mol_atol'],0.)
        ng=projected[i+1]; s=(other['positions_A']-x).ravel(); y=(ng-pg).ravel(); sy=float(s@y)
        sn=float(np.linalg.norm(s)); yn=float(np.linalg.norm(y))
        accepted=sy>1e-10*np.linalg.norm(s)*np.linalg.norm(y) and sy>1e-16
        if accepted: history.append((s,y,1./sy)); history=history[-core.SETTINGS['history']:]
        require(len(history)==len(nxt['history']),'History length update differs')
        for j,((hs,hy,rho),saved) in enumerate(zip(history,nxt['history'])):
            for k,a,b in (('s',hs,saved['s']),('y',hy,saved['y']),('rho',rho,saved['rho'])):
                check_close(result,np,str(step)+'/history/'+str(j)+'/'+k,a,b)
        ratio=work['components_kcal_mol']['composite_kcal_mol']/(multiplier*slope); trust=cp['trust_A']
        if multiplier==1. and ratio>.75: trust=min(core.SETTINGS['maximum_trust_A'],2*trust)
        elif ratio<.25: trust=max(core.SETTINGS['minimum_trust_A'],.5*trust)
        require(trust==nxt['trust_A'],'Exact trust branch/update differs')
        dx=other['positions_A']-x; dl=float(np.linalg.norm(dx)); require(dl>0,'Zero accepted displacement')
        projections={end:{k:{'local_gradient_dot_displacement_kcal_mol':float(np.sum(v*dx)),
                              'gradient_on_unit_displacement_kcal_mol_A':float(np.sum(v*dx/dl))}
                          for k,v in parts.items()} for end,parts in (('before',gradients[i]),('after',gradients[i+1]))}
        result['transitions'][i]=dict(step_from=step,step_to=step+1,status='available',work=work,
            history_before=len(cp['history']),history_after=len(history),H0_scale=h0,descent_fallback=fallback,
            pre_fallback_projected_slope=before_fallback,pre_trust_max_atom_displacement_A=maximum,trust_scale=scale,
            trust_before_A=cp['trust_A'],trust_after_A=trust,predicted_full_step_work_kcal_mol=slope,
            actual_to_predicted_ratio=ratio,multiplier=multiplier,backtracks=0,
            displacement_L2_A=dl,displacement_max_atom_A=float(np.linalg.norm(dx,axis=1).max()),
            s_dot_y=sy,s_norm=sn,y_norm=yn,curvature_cosine=sy/(sn*yn) if sn*yn else None,
            secant_projected_curvature_kcal_mol_A2=sy/(sn*sn) if sn else None,history_pair_admitted=bool(accepted),
            component_projections=projections,projections_are_local_not_integrated_work=True)
    result['windows']=[]
    for left,right in ((0,20),(20,40),(40,60),(60,80),(0,80)):
        work=energy_work(ev,points[left],points[right],raw[left],raw[right]); trs=result['transitions'][left:right]
        for k,v in work['components_kcal_mol'].items():
            check_close(result,np,f'window{left}-{right}/{k}',math.fsum(t['work']['components_kcal_mol'][k] for t in trs),v,TOLERANCES['work_closure_kcal_mol_atol'],0.)
        check_close(result,np,f'window{left}-{right}/ALPB_parts',work['ALPB_parts_closure_residual_kcal_mol'],0.,TOLERANCES['work_closure_kcal_mol_atol'],0.)
        rows=result['points'][left:right+1]
        result['windows'].append(dict(step_from=92+left,step_to=92+right,transitions=right-left,work=work,
            initial_force=rows[0]['raw_force'],terminal_force=rows[-1]['raw_force'],
            minimum_RMS_force=min(r['raw_force']['RMS_atom_vector_kcal_mol_A'] for r in rows),
            maximum_RMS_force=max(r['raw_force']['RMS_atom_vector_kcal_mol_A'] for r in rows),
            total_path_L2_A=math.fsum(t['displacement_L2_A'] for t in trs),
            net_displacement_L2_A=float(np.linalg.norm(points[right]['positions_A']-points[left]['positions_A'])),
            fallback_count=sum(t['descent_fallback'] for t in trs),history_pair_admissions=sum(t['history_pair_admitted'] for t in trs),
            trust_scaled_steps=sum(t['trust_scale']<1. for t in trs),minimum_trust_scale=min(t['trust_scale'] for t in trs),
            minimum_H0_scale=min(t['H0_scale'] for t in trs),maximum_H0_scale=max(t['H0_scale'] for t in trs),
            minimum_secant_curvature=min(t['secant_projected_curvature_kcal_mol_A2'] for t in trs),
            maximum_secant_curvature=max(t['secant_projected_curvature_kcal_mol_A2'] for t in trs),
            minimum_actual_to_predicted_ratio=min(t['actual_to_predicted_ratio'] for t in trs),
            maximum_actual_to_predicted_ratio=max(t['actual_to_predicted_ratio'] for t in trs),
            initial_projected_force=rows[0]['projected_force'],terminal_projected_force=rows[-1]['projected_force']))
    require(result['windows'][-1]['work']['components_kcal_mol']==run['rows'][0]['work_from_checkpoint'],'Overall native work differs')
    result['endpoint_attribution']={}
    for label,i in (('initial',0),('terminal',80)):
        data=attribution(np,core,read(raw[i]['physical_source']),gradients[i]); data.update(step=92+i,point=cfg['points'][i])
        result['endpoint_attribution'][label]=write(out/(label.upper()+'_ATTRIBUTION.json'),data)
    require(ev.engine.native is None and not ev.evaluation_requests and all(v==0 for v in ev.engine.cost.values()),'Audit made a molecular call')
    result.update(status='complete',admitted_points=81,admitted_transitions=80,model_initialized=False,
        replay_passed=True,scientific_force_criteria_unchanged=True,
        interpretation={'gradient_units':'positive derivative kcal/mol/Angstrom; forces are negatives',
            'projection_warning':'Endpoint gradient dot displacement is a local linear prediction, not finite-difference qualification or integrated work.',
            'curvature_warning':'s,y and H0 are optimizer secant diagnostics, not a physical Hessian or proof of a stiffness defect.',
            'recommendation':'Use component work, cancellation and recorded scaling/history to distinguish continued model motion from a demonstrated implementation bottleneck. Do not extend automatically or label solvent wrong from its lowering alone.',
            'accuracy_claim':False})


def report(result):
    lines=['# Saved step92-to172 mechanics audit','',f"Status: **{result['status']}**. Zero molecular calls; no optimization.",'']
    if result['status']!='complete':
        return '\n'.join(lines+['Replay failed; preserve the mismatch before any continuation.',str(result.get('failure')),''])
    lines+=['All81 saved states and80 transitions admitted under frozen replay tolerances. Exact raw component-gradient closure retained.','',
            '| Steps | MACE work | ALPB work | Total work | RMS force before → after | Trust-scaled steps | History resets |',
            '|---|---:|---:|---:|---:|---:|---:|']
    for row in result['windows']:
        w=row['work']['components_kcal_mol'];lines.append(f"| {row['step_from']}–{row['step_to']} | {w['MACE_kcal_mol']:.6f} | {w['fixed_ALPB_kcal_mol']:.6f} | {w['composite_kcal_mol']:.6f} | {row['initial_force']['RMS_atom_vector_kcal_mol_A']:.6f} → {row['terminal_force']['RMS_atom_vector_kcal_mol_A']:.6f} | {row['trust_scaled_steps']} | {row['fallback_count']} |")
    lines+=['','Work: kcal/mol. Force: kcal/mol/Å. Windows are fixed before execution.','',
            '## Native solvent work components','', '| Steps | Gelec | Ghb | Gsasa | Gshift |','|---|---:|---:|---:|---:|']
    for row in result['windows']:
        p=row['work']['ALPB_parts_kcal_mol'];lines.append('| '+str(row['step_from'])+'–'+str(row['step_to'])+' | '+' | '.join(f'{p[k]:.6f}' for k in PARTS)+' |')
    lines+=['','Full per-step component norms, dot products, projections, H0 scales, curvature and replay residuals are in RESULT.json. Initial/terminal atom and residue attribution files include the metal separately.','',
            'Passing replay excludes the tested bookkeeping inconsistencies; it does not establish physical accuracy. Continued lowering, solvent dominance, trust scaling or nonconvergence alone do not identify a defective interaction. No further search or reference is released.','']
    work=result['windows'][-1]['work']['components_kcal_mol']
    if work['MACE_kcal_mol']>0 and work['fixed_ALPB_kcal_mol']<0:
        lines+=['The observed net lowering trades increased vacuum MACE energy for lower solvent energy. The next mechanics decision should address this component balance and the measured direction/scaling behavior, rather than assume another step allowance or a bookkeeping repair is needed. Solvent lowering alone is not evidence that solvent is wrong.','']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic()
    result=dict(protocol=PROTOCOL,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_points=81,declared_transitions=80,
        points=[{'step':s,'status':'unavailable'} for s in range(92,173)],
        transitions=[{'step_from':s,'step_to':s+1,'status':'unavailable'} for s in range(92,172)],replay_checks=[],
        new_molecular_calls=0,model_loads=0,optimizer_steps=0,geometry_trials=0,force_accuracy_qualified=False)
    try: audit(read(result['config']),out,result)
    except Exception as exc:result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc(),replay_passed=False)
    result['elapsed_seconds']=time.monotonic()-tick
    (out/'REPORT.md').write_text(report(result));result['report']=pin(out/'REPORT.md');write(out/'RESULT.json',result)
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
