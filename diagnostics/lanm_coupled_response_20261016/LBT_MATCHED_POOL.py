#!/usr/bin/env python3
"""Rescore frozen18 geometries/36 states:4 exact origin reuses+32 native FF ALPB.
One CPU under Slurm. --config CONFIG --output FRESH. No new MACE/GFN/MM/search.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import struct
import time
import traceback


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def require(ok,message):
    if not ok:raise ValueError(message)


def read(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path']);return json.loads(Path(p['path']).read_text())


def write(p,value):
    with Path(p).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def coordinate_hash(x):
    return hashlib.sha256(b''.join(struct.pack('<d',0. if v==0. else v) for a in x for v in a)).hexdigest()


def inputs(cfg,out):
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Adapter differs');require(pin(cfg['plan']['path'])==cfg['plan'],'Plan differs')
    origin=read(cfg['matched_origin_result']);old=read(cfg['basis_pool_result']);inv=read(cfg['inventory'])
    require(origin['status']=='complete' and origin['slurm_job_id']=='12715' and origin['new_ALPB_calls_attempted']==4
        and origin['implementation']['sha256']==cfg['helper']['sha256'],'Matched origin/helper differs')
    require(old['status']=='complete' and old['slurm_job_id']=='12560' and old['declared_arm_states']==8
        and old['summary']['paired_arm_comparison_available'] and len(old['candidates'])==len(old['paired_matrix'])==18,'Complete frozen basis pool required')
    require(inv['basis_pool_result']==cfg['basis_pool_result'] and inv['matched_origin_result']==cfg['matched_origin_result']
        and inv['declared_geometries']==18 and inv['declared_metal_cells']==36 and len(inv['cells'])==36
        and inv['candidates']==old['candidates'],'Inventory denominator/source differs')
    require(pin(cfg['helper']['path'])['sha256']==cfg['helper']['sha256'],'Helper source differs')
    spec=importlib.util.spec_from_file_location('frozen_matched_lbt_v2',cfg['helper']['path']);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    h.verify(origin['implementation']);ocfg=read(origin['config']);d=out/'origin_admission';d.mkdir()
    bases,mex,scfg,build,native_source,pool_ok=h.inputs(ocfg,d)
    require(pool_ok and read(origin['preflight'])==json.loads((d/'PREFLIGHT.json').read_text()),'Origin chemistry/pool replay differs')
    require(h.verify(origin['build']['fortran']).read_text()==native_source,'Qualified wrapper source differs')
    binary=origin['build']['binary'];h.verify(binary)
    mc=read(read(ocfg['mace_result'])['config']);base={(c['source_id'],c['metal']):c for c in bases}
    candidates={c['candidate_id']:c for c in old['candidates']};pairs={r['candidate_id']:r['cells'] for r in old['paired_matrix']}
    indexed={(r['candidate_id'],r['metal']):r for r in inv['cells']}
    require(list(candidates)==['geometry_%02d'%i for i in range(18)] and len(indexed)==36
        and len({c['coordinate_sha256'] for c in candidates.values()})==18,'Frozen geometry identity/dedup differs')
    cells=[];origin_by={(r['source_id'],r['metal']):r for r in origin['rows']}
    for cid,candidate in candidates.items():
        sid=candidate['source_id'];require(sid in h.SOURCES,'Unexpected source')
        physical=read(candidate['physical_source']);reference=base[(sid,'La')]['physical']
        require(len(physical)==256 and [(a['id'],a['element'],a.get('resname')) for a in physical]
            ==[(a['id'],a['element'],a.get('resname')) for a in reference],'Candidate physical membership differs')
        px=h.vectors([a['xyz_A'] for a in physical],256)
        require(coordinate_hash(px)==candidate['coordinate_sha256'],'Saved physical coordinate hash differs')
        for metal in h.METALS:
            saved=pairs[cid][metal];i=indexed[(cid,metal)];b=base[(sid,metal)];mult=1 if metal=='La' else 6
            require(saved['status']=='accepted' and i['source_id']==sid and i['aliases']==candidate['aliases']
                and i['physical_source']==candidate['physical_source'] and i['coordinate_sha256']==candidate['coordinate_sha256'],'Saved matrix identity differs')
            obj=read(i['evidence_record']);raw=read(i['raw_MACE']);sy,x=h.xyz(i['xyz'])
            require(sy==b['symbols'] and h.coordinates_equal(x,px) and coordinate_hash(x)==candidate['coordinate_sha256'],'Saved XYZ/physical order differs')
            if i['kind']=='origin':
                ev=saved['evidence'][0];t=ev['origin_task'];require(i['evidence_record']==ev['MACE']==ocfg['mace_result'],'Origin MACE evidence differs')
                raw=next(r for r in raw['rows'] if r['task_id']==t['task_id'])
                require(t==b['task'] and i['xyz']==t['xyz'] and candidate['aliases']==[sid+':origin'],'Origin task differs')
            else:
                require(obj['MACE']==i['raw_MACE'] and obj['recipe']['model']==mc['model']
                    and obj['recipe']['MACE_versions']==mc['expected_versions'] and obj['recipe']['physical_charge']==-2
                    and obj['recipe']['physical_multiplicity']==mult and obj['recipe']['physical_ids']==[a['id'] for a in physical]
                    and obj['recipe']['symbols']==sy,'MACE model/state/identity differs')
                if i['kind']=='cross_cell':
                    require(saved['cell']==i['evidence_record'] and obj['status']=='accepted' and obj['task']==i['task_record'],'Cross receipt differs')
                    t=read(i['task_record']);require(t['xyz']==i['xyz'] and t['physical_source']==candidate['physical_source']
                        and t['source_id']==sid and t['metal']==metal and t['charge']==-2 and t['multiplicity']==mult,'Cross task differs')
                else:
                    require(i['kind']=='terminal_point' and saved['evidence'][0]['terminal_point']==i['evidence_record']
                        and obj['status']=='complete' and obj['source_id']==sid and obj['metal']==metal
                        and obj['xyz']==i['xyz'] and obj['physical_source']==candidate['physical_source'],'Own point differs')
                    h.verify(obj['positions'])
            h.state(raw['input_state_check'],256,-2,mult,255)
            require(math.isfinite(raw['energy_eV']) and raw['energy_eV']==saved['components']['MACE_eV'],'Raw MACE energy replay differs')
            c=dict(b,physical=physical,x=x,symbols=sy,MACE_eV=raw['energy_eV'],task=dict(b['task'],physical_source=candidate['physical_source'],xyz=i['xyz']))
            cell={'candidate_id':cid,'source_id':sid,'metal':metal,'MACE_evidence':i,'status':'unavailable','attempted':False,'reused':i['kind']=='origin'}
            if cell['reused']:
                r=origin_by[(sid,metal)];require(r['status']=='complete' and r['MACE_task_id']==b['task']['task_id'],'Origin reuse state differs')
                sol=read(r['ALPB']);h.pins_in(sol['artifacts']);start=read(sol['artifacts']['START.json'])
                require(h.xyz(start['xyz'])==(sy,x) and [float(v) for v in h.verify(start['charges']).read_text().split()]==c['q'],'Origin exact geometry/charges differ')
                require(sol['components_hartree']==r['native_energies']['ALPB_parts_Eh'] and sol['total_hartree']==r['native_energies']['ALPB_Eh']
                    and r['native_energies']['MACE_eV']==raw['energy_eV'],'Origin native components differ')
                cell.update(status='complete',ALPB=r['ALPB'],native_energies=r['native_energies'],component_energies_kcal_mol=r['component_energies_kcal_mol'],total_energy_kcal_mol=r['total_energy_kcal_mol'])
            cells.append((c,cell))
    require(sum(r['reused'] for _,r in cells)==4,'Four exact origins required')
    for source in h.SOURCES:
        require(sum(c['source_id']==source for c in candidates.values())==9,'Source-local denominator differs')
    write(out/'PREFLIGHT.json',{'status':'passed','molecular_calls':0,'geometries':18,'cells':36,'origin_reuses':4,'new_ALPB_ceiling':32,
        'basis_pool_result':cfg['basis_pool_result'],'matched_origin_result':cfg['matched_origin_result'],'inventory':cfg['inventory'],
        'original_origin_admission':h.pin(d/'PREFLIGHT.json'),'native_binary':binary,'MACE_model':mc['model'],'all_raw_MACE_replayed':True})
    return h,cells,mex,scfg,build,binary,old


def summarize(h,rows,mex,old):
    by={(r['candidate_id'],r['metal']):r for r in rows};candidates=old['candidates'];menus=[]
    for source in (*h.SOURCES,'global'):
        local=[c['candidate_id'] for c in candidates if source=='global' or c['source_id']==source]
        origin=[c['candidate_id'] for c in candidates if c['candidate_id'] in ('geometry_00','geometry_01') and c['candidate_id'] in local]
        oldmenus=old['summary']['global_menus'] if source=='global' else old['summary']['source_local_menus'][source]
        definitions={'origin':origin,**{name:oldmenus[name]['candidate_ids'] for name in ('prior','control','candidate')},'union':local}
        baseline=None
        for name,ids in definitions.items():
            r={'source_scope':source,'menu':name,'candidate_ids':ids,'status':'unavailable','declared_states':2*len(ids)};menus.append(r)
            if not all(by[(cid,m)]['status']=='complete' for cid in ids for m in h.METALS):continue
            selected={m:min((by[(cid,m)] for cid in ids),key=lambda x:x['total_energy_kcal_mol']) for m in h.METALS}
            v={m:selected[m]['native_energies'] for m in h.METALS};diff=h.difference(v['Dy'],v['La']);md=h.difference(mex['Dy'],mex['La'])
            rc={k:math.fsum([md[k],-diff[k]]) for k in h.NAMES}
            r.update(status='available',selected={m:selected[m]['candidate_id'] for m in h.METALS},
                exact_ties={m:[cid for cid in ids if by[(cid,m)]['total_energy_kcal_mol']==selected[m]['total_energy_kcal_mol']] for m in h.METALS},
                selected_native_energies=v,LBT_Dy_minus_La_components_kcal_mol=diff,R_components_kcal_mol=rc,R_kcal_mol=math.fsum(rc.values()))
            if name=='origin':baseline=r
            r['own_metal_work_from_origin_components_kcal_mol']={m:h.difference(v[m],baseline['selected_native_energies'][m]) for m in h.METALS}
            r['R_change_from_origin_kcal_mol']=math.fsum(rc[k]-baseline['R_components_kcal_mol'][k] for k in h.NAMES)
    return {'menus':menus,'old_summary_unchanged':old['summary_artifact'],'old_menu_selections':{
        s:{k:v['selected'] for k,v in (old['summary']['global_menus'] if s=='global' else old['summary']['source_local_menus'][s]).items()} for s in (*h.SOURCES,'global')},
        'balanced_exchange':'Mex-La + LBT-Dy -> Mex-Dy + LBT-La; positive R favors reactants',
        'limits':'Fixed old-objective candidate menu rescored; not a matched-FF search, population, affinity or classifier'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'One Slurm CPU required')
    cfg=json.loads(a.config.read_text());out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    result={'status':'running','protocol_id':'lbt_frozen18_matchedFF_rescore_v1','implementation':pin(__file__),'config':pin(a.config),'plan':cfg['plan'],
        'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'declared_geometries':18,'declared_cells':36,
        'new_ALPB_calls_attempted':0,'reused_origin_ALPB_cells':0,'reused_MACE_cells':0,'new_MACE_calls':0,'new_GFN2_calls':0,'new_MM_calls':0,'new_DFT_calls':0,'optimizer_steps':0,
        'rows':[{'candidate_id':'geometry_%02d'%i,'metal':m,'status':'unavailable','attempted':False} for i in range(18) for m in ('La','Dy')]}
    try:
        h,cells,mex,scfg,build,binary,old=inputs(cfg,out);result.update(preflight=pin(out/'PREFLIGHT.json'),rows=[r for _,r in cells],reused_origin_ALPB_cells=4,reused_MACE_cells=36)
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=scfg['library_directory'])
        for key in ('LD_PRELOAD','PYTHONHOME','PYTHONPATH'):env.pop(key,None)
        for c,r in cells:
            if r['reused']:continue
            try:
                for lib in build['linked_libraries']:h.verify(lib)
                result['new_ALPB_calls_attempted']+=1;r['attempted']=True;d=out/(r['candidate_id']+'_'+r['metal'])
                sol=h.solvent(binary,c,d,env);v={'MACE_eV':c['MACE_eV'],'field_Eh':0.,'MM_kcal_mol':0.,'ALPB_Eh':sol['total_hartree'],'ALPB_parts_Eh':sol['components_hartree']}
                comp=h.native_components(v);r.update(status='complete',ALPB=pin(d/'RESULT.json'),native_energies=v,component_energies_kcal_mol=comp,total_energy_kcal_mol=math.fsum(comp.values()))
            except Exception as exc:r.update(failure=repr(exc),traceback=traceback.format_exc())
        require(result['new_ALPB_calls_attempted']<=32,'Call ceiling exceeded');result['comparison']=summarize(h,result['rows'],mex,old)
        result['status']='complete' if all(r['status']=='complete' for r in result['rows']) else 'unavailable'
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-start;result['application_allocated_CPU_seconds']=result['elapsed_seconds'];write(out/'RESULT.json',result)
    lines=['# Frozen LBT pool under matched FF solvent','',f"Status: {result['status']}; new ALPB attempts {result['new_ALPB_calls_attempted']}/32.",'',
        '| Source scope | Menu | Status | La selected | Dy selected | R kcal/mol | Change from origin |','|---|---|---|---|---|---:|---:|']
    for r in result.get('comparison',{}).get('menus',[]):
        lines.append('| '+' | '.join(str(x) for x in (r['source_scope'],r['menu'],r['status'],r.get('selected',{}).get('La','—'),r.get('selected',{}).get('Dy','—'),r.get('R_kcal_mol','—'),r.get('R_change_from_origin_kcal_mol','—')))+' |')
    lines+=['','Positive R favors Mex-La+LBT-Dy. Full components and all36 cells remain in RESULT.json.',
        'Old-objective candidate geometries are rescored without new search; no population, affinity or classifier claim.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
