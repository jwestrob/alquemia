#!/usr/bin/env python3
"""Same-model common5: original common3 plus two admitted continuation terminals.
CLI --stage pool|mm --config CONFIG --output FRESH; all numeric work requires Slurm.
Config: implementation,plan,pilot_result,search_pattern,frozen_scout_implementation,
remaining_implementation,scout_result,remaining_result. No search or parameter changes.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import time
import traceback

BASE={'path':'/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_coupled_response_20261016/MEX_CONTINUATION.py',
      'sha256':'c571f20ec76bc7bf7ecadb7b3f2c402caf0f1bdfcb57b6a84aaf4c3a2dd15712'}
REMAINING_SHA='dddc4d25fdd0934cb71d08a9a70d6b93dcf55d253dfe7f31576cc6aa7fc32446'
ALIASES=['origin','La_pilot','Dy_pilot','La_terminal','Dy_terminal']


def load_base():
    path=Path(BASE['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=BASE['sha256']:raise ValueError('Frozen evaluator differs')
    spec=importlib.util.spec_from_file_location('mex_pool_base',path);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base);return base


def parents(base,cfg,pilot,points,mech):
    require=base.require;np=mech.np
    require(cfg['frozen_scout_implementation']['sha256']==BASE['sha256'] and cfg['remaining_implementation']['sha256']==REMAINING_SHA,'Frozen search source differs')
    base.verify(cfg['frozen_scout_implementation']);remaining_module=base.module('mex_pool_remaining_gate',base.verify(cfg['remaining_implementation']))
    scout=base.read(cfg['scout_result']);dy=base.read(cfg['remaining_result'])
    for p,stage,sha in ((scout,'scout',BASE['sha256']),(dy,'remaining',REMAINING_SHA)):
        require(p['stage']==stage and p['implementation']['sha256']==sha,'Parent executed source/stage differs');base.verify(p['implementation'])
        pcfg=base.read(p['config']);require(pcfg['implementation']['sha256']==sha,'Parent config source differs')
        for key in ('plan','pilot_result','search_pattern'):
            base.verify(pcfg[key]);require(pcfg[key]==cfg[key],'Parent shared input differs: '+key)
        require(p['plan']==cfg['plan'] and p['pilot_result']==cfg['pilot_result'] and p['settings']==base.SETTINGS,'Parent declaration differs')
        pre=base.read(p['preflight']);require(pre['status']=='passed' and pre['new_molecular_calls']==0
            and pre['pilot_result']==cfg['pilot_result'] and pre['modes']==pilot['modes'] and len(pre['starts'])==4
            and {(r['proposer'],r['metal']) for r in pre['starts']}=={(a,b) for a in ('La','Dy') for b in ('La','Dy')}
            and all(r['passed'] and r['geometry']['passed'] for r in pre['starts']),'Parent preflight differs')
        require(p['status'] in ('complete','unavailable'),'Globally failed execution cannot supply new pool candidates')
    require(dy['scout_result']==cfg['scout_result'] and base.read(dy['config'])['scout_result']==cfg['scout_result'],'Dy scout gate lineage differs')
    rows={};slots=[]
    for metal,parent in (('La',scout),('Dy',dy)):
        require(len(parent['rows'])==1 and parent['rows'][0]['metal']==metal,'Parent search denominator differs')
        row=parent['rows'][0];search=base.read(row['search']);slot={'metal':metal,'alias':metal+'_terminal','search':row['search'],'status':'unavailable'}
        slots.append(slot)
        require(search['metal']==metal and search['settings']==base.SETTINGS and search['pilot_result']==cfg['pilot_result']
            and search['points']==row['points']==parent['points'] and search['requests']==parent['requests'],'Search state/receipt lineage differs')
        admitted=parent['status']=='complete' and row['candidate_admitted'] and search['candidate_admitted']
        if not admitted:slot['failure']='parent or search unavailable';continue
        require(row['status']==search['status']=='finite_proposal_available' and search['final_geometry']['passed']
            and search['accepted_steps']<=12 and search['new_own_composite_evaluations']<=24
            and search['initial_history_length']==0 and search['work_from_continuation_start']['composite_kcal_mol']<=1e-8,'Search admission differs')
        require(base.read(search['continuation_start'])==points[(metal,metal)],'Search did not start from exact own pilot')
        require(np.max(np.abs(np.asarray(search['coordinate_scaling_T'])-mech.T))<=1e-12,'Original whitening differs')
        point=base.read(search['candidate']);ph=base.read(point['physical']);y=np.asarray([a['xyz_A'] for a in ph]);q=np.asarray(search['full_q'])
        require(point['status']=='complete' and point['metal']==metal and [a['id'] for a in ph]==mech.ids
            and np.max(np.abs(mech.coordinates(q)-y))<=1e-12 and mech.check(q,y,q,y)['passed'],'Terminal coordinates/state/guard differs')
        require(search['candidate'] in parent['points'],'Terminal lacks parent point receipt')
        slot.update(status='available',candidate=search['candidate'],accepted_steps=search['accepted_steps'],termination=search['termination'],
            work_from_original_source=search['work_from_original_source'],work_from_continuation_start=search['work_from_continuation_start'],
            final_force_diagnostic=search['final_force_diagnostic'],final_geometry=search['final_geometry'])
        rows[metal]={'search':search,'parent':parent,'point':point,'q':q,'y':y}
    if len(rows)==2:remaining_module.scout_gate(base,cfg)
    return scout,dy,rows,slots


def summary(base,evaluator,menu,aliases):
    chosen=[]
    for alias in aliases:
        entry=next(r for r in menu if alias in r['aliases']);chosen.append((alias,entry))
    selected={}
    for metal in ('La','Dy'):
        alias,row=min(chosen,key=lambda p:p[1]['total_work_from_original_source_kcal_mol'][metal])
        selected[metal]={'geometry_alias':alias,'all_geometry_aliases':row['aliases'],'coordinate_sha256':row['coordinate_sha256'],
            'work_from_original_source_kcal_mol':row['total_work_from_original_source_kcal_mol'][metal],
            'point':row['cells'][metal]['point']}
    return {'status':'available','declared_menu':aliases,'selections':selected,
        'conditional_Dy_minus_La_accommodation_change_kcal_mol':selected['Dy']['work_from_original_source_kcal_mol']-selected['La']['work_from_original_source_kcal_mol'],
        'common_selected_work_kcal_mol':(selected['Dy']['work_from_original_source_kcal_mol']+selected['La']['work_from_original_source_kcal_mol'])/2}


def run(base,cfg,out,result):
    require=base.require;require(int(os.environ['SLURM_CPUS_PER_TASK'])>=2,'Two single-thread GFN workers require allocation')
    loaded=base.inputs(cfg,out);response,fcfg,world,mech,pilot,points=loaded;np=world[6]
    result.update(preflight=base.pin(out/'PREFLIGHT.json'),pilot_result=cfg['pilot_result'],modes=pilot['modes'],
        scout_result=cfg['scout_result'],remaining_result=cfg['remaining_result'],
        original_common3=pilot['comparison']['common3'],declared_menu=ALIASES,
        common5={'status':'unavailable','declared_menu':ALIASES},terminal_slots=[],new_terminal_cross_cells_attempted=0)
    scout,dy,terminals,slots=parents(base,cfg,pilot,points,mech);result['terminal_slots']=slots
    if len(terminals)!=2:
        result.update(status='unavailable',failure='Both declared terminal slots required; original common3 retained');return
    ev=base.Evaluator(cfg,out,result,loaded)
    for parent in (scout,dy):
        for record in parent['points']:
            p=base.read(record)
            if p['status']=='complete':ev.admit_saved(record)
    # Candidate membership is only the declared origins/pilots/admitted terminals.
    candidates=[('origin',np.zeros(7),mech.x)]
    for metal in ('La','Dy'):
        p=points[(metal,metal)];candidates.append((metal+'_pilot',np.asarray(p['q']),ev.coordinates(p)))
    for metal in ('La','Dy'):candidates.append((metal+'_terminal',terminals[metal]['q'],terminals[metal]['y']))
    groups={}
    for alias,q,y in candidates:
        key=ev.key('geometry',y).split(':',1)[1]
        if key not in groups:groups[key]={'aliases':[],'q':q,'y':y,'coordinate_sha256':key}
        groups[key]['aliases'].append(alias)
    menu=[];result['matrix']=menu;result['unique_physical_geometries']=len(groups)
    for group in groups.values():
        y,q=group['y'],group['q'];gc=mech.check(q,y,q,y);require(gc['passed'],'Pool candidate violates original-source guard')
        row={'aliases':group['aliases'],'coordinate_sha256':group['coordinate_sha256'],'full_q':q.tolist(),'geometry':gc,'cells':{},'status':'unavailable'};menu.append(row)
        components={};total={}
        for metal in ('La','Dy'):
            cached=ev.key(metal,y) in ev.cache
            if not cached:
                require(any(alias.endswith('_terminal') for alias in group['aliases']),'Old baseline cell unexpectedly missing')
                require(result['new_terminal_cross_cells_attempted']<2,'Terminal cross-cell limit exceeded')
                require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),'Allocated GPU required')
                result['new_terminal_cross_cells_attempted']+=1
            point,record=ev.evaluate(metal,q,y,gc);require(point['status']=='complete','Pool cell unavailable')
            parts=base.relative(response,world[0],point['native_energies'],response.origin_native(metal,world))
            components[metal]={k:v for k,v in parts.items() if k!='composite_kcal_mol'};total[metal]=parts['composite_kcal_mol']
            row['cells'][metal]={'point':record,'metal':metal,'reused':cached,'native_energies':point['native_energies'],
                'assembled_gradient':point['assembled_gradient'],'physical':point['physical']}
            require(np.array_equal(ev.coordinates(point),y),'Pool paired coordinate mismatch')
        row.update(status='complete',component_work_from_original_source_kcal_mol=components,
            total_work_from_original_source_kcal_mol=total,common_work_kcal_mol=(total['La']+total['Dy'])/2,
            Dy_minus_La_work_kcal_mol=total['Dy']-total['La'])
    common3=summary(base,ev,menu,ALIASES[:3]);common5=summary(base,ev,menu,ALIASES)
    # Reproduce the old numerical result before reporting any extension.
    old=pilot['comparison']['common3'];names={'origin':'origin','La_pilot':'La_proposal','Dy_pilot':'Dy_proposal'}
    for metal in ('La','Dy'):
        require(names[common3['selections'][metal]['geometry_alias']]==old['selections'][metal]['geometry']
            and abs(common3['selections'][metal]['work_from_original_source_kcal_mol']-old['selections'][metal]['work_from_origin_kcal_mol'])<=1e-12,'Original common3 selection/work replay differs')
    require(abs(common3['conditional_Dy_minus_La_accommodation_change_kcal_mol']-old['conditional_Dy_minus_La_accommodation_change_kcal_mol'])<=1e-12,'Original common3 contrast replay differs')
    result.update(status='complete',common3_replay=common3,common5=common5,
        common3_to_common5_Dy_minus_La_change_kcal_mol=common5['conditional_Dy_minus_La_accommodation_change_kcal_mol']-common3['conditional_Dy_minus_La_accommodation_change_kcal_mol'],
        points=ev.points,requests=ev.requests,
        parent_costs={m:{k:p[k] for k in ('new_composite_evaluations_attempted','new_MACE_calls_attempted','new_GFN2_calls_attempted',
            'new_MM_configurations_attempted','new_MM_group_queries_attempted','new_ALPB_calls_attempted','elapsed_seconds',
            'application_allocated_CPU_seconds','application_allocated_GPU_seconds')} for m,p in (('La',scout),('Dy',dy))},
        scope='Bounded seven-coordinate EF3 response with fixed three-Nd spectators; no stationary minimum or protein affinity claim')
    if ev.native is not None:result.update(native_model_calls=ev.native.calls,native_device=ev.native.device,model=ev.mcfg['model'],versions=ev.native.versions)
    require(result['new_composite_evaluations_attempted']<=2 and result['new_MACE_calls_attempted']<=2
        and result['new_GFN2_calls_attempted']<=4 and result['new_MM_configurations_attempted']<=2
        and result['new_MM_group_queries_attempted']<=8 and result['new_ALPB_calls_attempted']<=2,'Pool molecular ceiling exceeded')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',choices=('pool','mm'),default='pool')
    parser.add_argument('--config',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    base=load_base();base.require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    cfg=json.loads(args.config.read_text());base.require(base.pin(__file__)['sha256']==cfg['implementation']['sha256'],'Pool source differs');base.verify(cfg['plan'])
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'running','stage':args.stage,'implementation':base.pin(__file__),'config':base.pin(args.config),'plan':cfg['plan'],
        'base':BASE,'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'settings':base.SETTINGS,
        'new_composite_evaluations_attempted':0,'new_MACE_calls_attempted':0,'new_GFN2_calls_attempted':0,
        'new_MM_configurations_attempted':0,'new_MM_group_queries_attempted':0,'new_ALPB_calls_attempted':0,
        'new_native_model_loads':0,'observed_GFN2_child_CPU_seconds':0.,'new_DFT_calls':0,'new_wholeprotein_SCC_calls':0}
    try:
        if args.stage=='mm':base.mm_stage(cfg,out,result)
        else:run(base,cfg,out,result)
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-tick;result['process_CPU_seconds']=time.process_time()-cpu
    if args.stage=='pool':
        result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
        result['application_allocated_GPU_seconds']=result['elapsed_seconds']
    base.write(out/'RESULT.json',result);return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
