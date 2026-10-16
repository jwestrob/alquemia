#!/usr/bin/env python3
"""Frozen16 water slots: twelve new fixed-FF ALPB cells, four retained exclusions.
One Slurm CPU. --config CONFIG --output FRESH. No new MACE/GFN/MM/geometry.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import socket
import subprocess
import time
import traceback
import xml.etree.ElementTree as ET

SOURCES=('7CCO','7CCN');POSES=('inner_plane0','inner_plane90','outer_plane0','outer_plane90');METALS=('La','Dy')
WATER_IDS=['W/1//O','W/1//H1','W/1//H2']


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path']);return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def write(p,value):
    with Path(p).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def inputs(cfg,out):
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Adapter differs');verify(cfg['plan'])
    origin=read(cfg['matched_origin_result']);old=read(cfg['water_panel_result']);scout=read(cfg['water_scout_result']);prep=read(cfg['water_preparation_result'])
    require(origin['status']=='complete' and origin['slurm_job_id']=='12715' and origin['new_ALPB_calls_attempted']==4
        and origin['implementation']['sha256']==cfg['matched_helper']['sha256'],'Matched origin/helper differs')
    require(old['status']=='complete' and old['stage']=='remaining' and old['slurm_job_id']=='12670'
        and old['declared_cells']==16 and old['accepted_cells']==12 and old['geometry_rejected_cells']==4
        and old['all16_available'] is False and old['scout_result']==cfg['water_scout_result'],'Closed water panel differs')
    require(scout['status']=='complete' and scout['stage']=='scout' and scout['slurm_job_id']=='12664'
        and scout['scout_gate_passed'] and old['water_preparation_result']==scout['water_preparation_result']==cfg['water_preparation_result']
        and prep['status']=='complete' and prep['slurm_job_id']=='12565','Water source/scout differs')
    for parent in (origin,old,scout,prep):
        pc=read(parent['config']);verify(parent['implementation']);verify(pc['implementation']);verify(parent['plan']);verify(pc['plan'])
        require(parent['implementation']['sha256']==pc['implementation']['sha256'] and parent['plan']==pc['plan'],'Parent implementation/plan binding differs')
    spec=importlib.util.spec_from_file_location('matched_water_original_helper',verify(cfg['matched_helper']));h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    ocfg=read(origin['config']);d=out/'original_admission';d.mkdir();bases,mex,scfg,build,source,pool_ok=h.inputs(ocfg,d)
    require(read(origin['preflight'])==read(pin(d/'PREFLIGHT.json')) and verify(origin['build']['fortran']).read_text()==source,'Original FF/model/native replay differs')
    marker="if (n /= 256) error stop 'declared LBT physical scope requires 256 atoms'"
    require(source.count(marker)==1,'Native atom-count assertion differs')
    source=source.replace(marker,"if (n /= 259) error stop 'declared one-water LBT scope requires 259 atoms'")
    mc=read(read(ocfg['mace_result'])['config'])
    require(old['model']==scout['model']==mc['model'] and old['native_runtime']['versions']==scout['native_runtime']['versions']==mc['expected_versions']
        and prep['preparation_result']==ocfg['preparation'],'Native model/version/original source differs')
    export=read(cfg['mex_water_export']);ff=read(export['atoms']);require(cfg['tip3p_xml'] in export['forcefields'],'TIP3P provenance differs')
    tree=ET.parse(verify(cfg['tip3p_xml']));types={a.attrib['type']:float(a.attrib['charge']) for a in tree.findall('./NonbondedForce/Atom')}
    residues=tree.findall("./Residues/Residue[@name='HOH']/Atom")
    require([(a.attrib['name'],a.attrib['type']) for a in residues]==[('O','tip3p-O'),('H1','tip3p-H'),('H2','tip3p-H')],'TIP3P water membership differs')
    wq=[types[a.attrib['type']] for a in residues]
    require(wq==[-.834,.417,.417] and [a['charge_e'] for a in ff if a['id'] in ('A/301//O','A/301//H1','A/301//H2')]==wq,'TIP3P native charge replay differs')
    base={(c['source_id'],c['metal']):c for c in bases};expected=[s+'_'+p+'_'+m for s in SOURCES for p in POSES for m in METALS]
    require([t['task_id'] for t in prep['tasks']]==[r['task']['task_id'] for r in old['rows']]==expected
        and [r['source_id'] for r in prep['rows']]==list(SOURCES),'Frozen sixteen-slot order differs')
    sources={r['source_id']:r for r in prep['rows']};scout_rows={r['task']['task_id']:r for r in scout['rows']}
    cells=[];proof=[];paired={}
    for t,row in zip(prep['tasks'],old['rows']):
        sid,pose,metal=t['source_id'],t['pose_id'],t['metal'];sr=sources[sid];b=base[(sid,metal)]
        require(row['task']==t and t in sr['tasks'],'Water task binding differs')
        require(read(sr['result'])=={k:v for k,v in sr.items() if k!='result'},'Water source receipt differs')
        p=next(p for p in sr['poses'] if p['pose_id']==pose);ph=read(t['physical_source']);mapping=read(t['source_map']);graph=read(t['graph'])
        require(read(p['result'])=={k:v for k,v in p.items() if k!='result'} and p['admitted']==t['admitted']==p['contact_check']['pass']
            and p['original256_exactly_unchanged'] and t['physical_source']==p['physical_source']
            and t['source_map']==sr['source_map'] and t['graph']==sr['graph'],'Frozen pose/admission differs')
        require(len(ph)==len(mapping)==259 and ph[:256]==b['physical'] and [a['id'] for a in ph[256:]]==WATER_IDS
            and [a['element'] for a in ph[256:]]==['O','H','H'] and [a['source_atom'] for a in mapping[:256]]==b['physical']
            and [a['physical_id'] for a in mapping]==[a['id'] for a in ph],'Whole physical/map identity differs')
        require(graph['water_bonds']==[[256,257],[256,258]] and graph['metal_covalent_bonds']==[] and graph['water_formal_charge']==0
            and graph['water_template']==t['water_template']==prep['water_template']
            and t['original_source']==b['task']['physical_source'] and t['original_endpoint']==b['task'],'Water/source graph differs')
        h.pins_in(graph);verify(t['placement']);verify(t['water_template']);dp=read(sr['source_preparation'])
        require(graph['original_peptide_system']==dp['native_peptide_system'] and graph['original_peptide_bonds']==dp['covalent_bonds'],'Original FF graph lineage differs')
        mult=1 if metal=='La' else 6
        require(t['atom_count']==259 and t['charge']==-2 and t['metal_index']==255 and t['water_count']==1
            and t['water_indices']==[256,257,258] and t['multiplicity']==t['spin_multiplicity']==t['physical_multiplicity']==mult,'Water electronic state differs')
        sy,x=h.xyz(t['xyz']);require(sy==b['symbols']+['O','H','H'] and x==[a['xyz_A'] for a in ph],'Water XYZ/physical equality differs');h.vectors(x,259)
        key=(sid,pose)
        if key in paired:require(x==paired[key],'Paired water coordinates differ')
        else:paired[key]=x
        q=b['q']+wq;require(abs(math.fsum(q)+2)<=1e-6,'Full fixed charge closure differs')
        r={'task_id':t['task_id'],'source_id':sid,'pose_id':pose,'metal':metal,'task':t,'status':'unavailable','attempted':False,'old_cell_status':row['status']}
        rejected=sid=='7CCN' and pose.startswith('inner')
        require(t['admitted'] is not rejected,'Declared geometry exclusions differ')
        if rejected:
            require(row['status']=='geometry_rejected' and 'MACE' not in row,'Rejected cell was scored')
            r.update(status='geometry_rejected',reason=row['reason']);cells.append((None,r));continue
        require(row['status']=='complete' and read(row['cell'])=={k:v for k,v in row.items() if k!='cell'},'Saved energy CELL differs')
        raw=read(row['MACE']);tp=pin(Path(row['MACE']['path']).parent/'TASK.json');require(read(tp)==t,'Raw MACE TASK differs')
        h.state(raw['input_state_check'],259,-2,mult,255)
        require(raw['forces_eV_A'] is None and math.isfinite(raw['energy_eV']) and raw['energy_eV']==row['components']['MACE_eV'],'Raw energy-only MACE replay differs')
        if pose=='inner_plane0' and sid=='7CCO':require(row==scout_rows[t['task_id']],'Exact two-scout reuse differs')
        c={'task':t,'physical':ph,'symbols':sy,'x':x,'q':q,'MACE_eV':raw['energy_eV']}
        r.update(MACE=row['MACE'],MACE_task=tp,old_cell=row['cell']);cells.append((c,r));proof.append({'task_id':t['task_id'],'MACE':row['MACE'],'MACE_task':tp,'physical':t['physical_source'],'xyz':t['xyz']})
    require(len(cells)==16 and len(proof)==12,'Declared rescore denominator differs')
    write(out/'PREFLIGHT.json',{'status':'passed','molecular_calls':0,'declared_cells':16,'new_ALPB_ceiling':12,'retained_geometry_rejections':4,
        'original_admission':pin(d/'PREFLIGHT.json'),'water_panel_result':cfg['water_panel_result'],'MACE_model':mc['model'],
        'raw_MACE_replay':proof,'tip3p_xml':cfg['tip3p_xml'],'water_charges_e':wq,'water_geometry_changed':False})
    return h,cells,scfg,build,source,old


def solvent(h,binary,c,d,env):
    d.mkdir();xp=d/'physical.xyz';qp=d/'charges'
    xp.write_text('259\nFrozen whole LBT plus one water; native FF/TIP3P plus formal3\n'+''.join('%s %.17g %.17g %.17g\n'%(s,*v) for s,v in zip(c['symbols'],c['x'])))
    qp.write_text(''.join('%.17g\n'%q for q in c['q']))
    require(h.xyz(pin(xp))==(c['symbols'],c['x']) and [float(v) for v in qp.read_text().split()]==c['q'],'Input serialization differs')
    command=[str(verify(binary)),str(xp),str(qp)];tick=time.monotonic()
    write(d/'START.json',{'command':command,'xyz':pin(xp),'charges':pin(qp),'physical_source':c['task']['physical_source']})
    with (d/'stdout.txt').open('x') as so,(d/'stderr.txt').open('x') as se:proc=subprocess.run(command,cwd=d,env=env,stdout=so,stderr=se)
    write(d/'EXECUTION.json',{'command':command,'returncode':proc.returncode,'elapsed_seconds':time.monotonic()-tick})
    require(proc.returncode==0,'Native ALPB execution failed');text=(d/'stdout.txt').read_text()
    def values(label,n):
        rows=re.findall(r'^'+label+r'\s+(.+)$',text,re.M);require(len(rows)==1,'Native output field '+label)
        v=[float(s) for s in rows[0].split()];require(len(v)==n and all(math.isfinite(s) for s in v),'Invalid native field '+label);return v
    parts=dict(zip(h.PARTS,values('COMPONENTS',4)));charge=values('CHARGE_SUM',1)[0]
    require(abs(charge+2)<=1e-5 and values('ROUNDING_BOUNDS',2)==[0.,0.],'Native charge/rounding differs')
    h.vectors([[float(v) for v in line.split()] for line in (d/'gradient.dat').read_text().splitlines()],259)
    r={'components_hartree':parts,'total_hartree':math.fsum(parts.values()),'charge_e':charge,'rounding_bounds_hartree':[0.,0.],
       'cavity_global':values('CAVITY_GLOBAL',4),'gradient_units':'Hartree/Bohr; archived only','artifacts':{p.name:pin(p) for p in d.iterdir() if p.is_file()}}
    write(d/'RESULT.json',r);return r


def difference(h,a,b):
    v={'MACE_kcal_mol':(a['MACE_eV']-b['MACE_eV'])*h.EV,'fixed_FF_ALPB_kcal_mol':(a['ALPB_Eh']-b['ALPB_Eh'])*h.HA}
    return dict(v,composite_kcal_mol=math.fsum(v.values()))


def summarize(h,rows,old):
    by={(r['source_id'],r['pose_id'],r['metal']):r for r in rows};transfers=[];menus={}
    for sid in SOURCES:
        for plane in (0,90):
            pair={m:[by[(sid,s+'_plane'+str(plane),m)] for s in ('inner','outer')] for m in METALS}
            missing=[r['task_id'] for rs in pair.values() for r in rs if r['status']!='complete'];r={'source_id':sid,'plane_angle_degree':plane,'status':'unavailable','missing_cells':missing};transfers.append(r)
            if missing:continue
            own={m:difference(h,*[r['native_energies'] for r in pair[m]]) for m in METALS}
            v={m:[r['native_energies'] for r in pair[m]] for m in METALS}
            hc={'MACE_kcal_mol':math.fsum([v['Dy'][0]['MACE_eV']-v['Dy'][1]['MACE_eV'],-(v['La'][0]['MACE_eV']-v['La'][1]['MACE_eV'])])*h.EV,
                'fixed_FF_ALPB_kcal_mol':math.fsum([v['Dy'][0]['ALPB_Eh']-v['Dy'][1]['ALPB_Eh'],-(v['La'][0]['ALPB_Eh']-v['La'][1]['ALPB_Eh'])])*h.HA}
            hc['composite_kcal_mol']=math.fsum(hc.values())
            previous=next(x for x in old['summary']['transfers'] if x['source_id']==sid and x['plane_angle_degree']==plane)
            r.update(status='available',own_metal_inner_minus_outer=own,H_Dy_minus_La_transfer=hc,
                interpretation_at_1_kcal_mol='near_zero' if abs(hc['composite_kcal_mol'])<=1 else ('positive' if hc['composite_kcal_mol']>0 else 'negative'),
                old_selfconsistent_solvent_transfer=previous)
        missing=[by[(sid,p,m)]['task_id'] for p in POSES for m in METALS if by[(sid,p,m)]['status']!='complete']
        menu={'status':'unavailable','candidate_order':list(POSES),'missing_cells':missing};menus[sid]=menu
        if missing:continue
        relative={m:{p:difference(h,by[(sid,p,m)]['native_energies'],by[(sid,POSES[0],m)]['native_energies']) for p in POSES} for m in METALS}
        selected={m:min(POSES,key=lambda p:relative[m][p]['composite_kcal_mol']) for m in METALS}
        menu.update(status='available',selected=selected,same_metal_relative_energies=relative,
            exact_ties={m:[p for p in POSES if relative[m][p]['composite_kcal_mol']==relative[m][selected[m]]['composite_kcal_mol']] for m in METALS},
            selected_Dy_minus_La=difference(h,by[(sid,selected['Dy'],'Dy')]['native_energies'],by[(sid,selected['La'],'La')]['native_energies']))
    return {'transfers':transfers,'source_local_common_four_pose_menus':menus,'old_summary_unchanged':old['summary'],
        'interpretive_resolution_kcal_mol':1.,'all_four_transfer_contrasts_available':all(r['status']=='available' for r in transfers),
        'scope':'One fixed water on original sources only; no dry/wet minima, imputation, global pool, occupancy or affinity'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'One Slurm CPU required')
    cfg=json.loads(a.config.read_text());out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic()
    result={'status':'running','protocol_id':'lbt_closed_one_water_matched_FF_v1','implementation':pin(__file__),'config':pin(a.config),'plan':cfg['plan'],
        'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'declared_cells':16,'reused_MACE_cells':0,'new_ALPB_calls_attempted':0,
        'new_MACE_calls':0,'new_GFN2_calls':0,'new_MM_calls':0,'new_DFT_calls':0,'new_geometries':0,'optimizer_steps':0,
        'rows':[{'task_id':s+'_'+p+'_'+m,'source_id':s,'pose_id':p,'metal':m,'status':'not_admitted','attempted':False} for s in SOURCES for p in POSES for m in METALS]}
    try:
        h,cells,scfg,build,source,old=inputs(cfg,out);result.update(preflight=pin(out/'PREFLIGHT.json'),rows=[r for _,r in cells],reused_MACE_cells=12,water_panel_result=cfg['water_panel_result'])
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=scfg['library_directory'])
        for key in ('LD_PRELOAD','PYTHONHOME','PYTHONPATH'):env.pop(key,None)
        result['native_environment']={k:env.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','LD_LIBRARY_PATH','LD_PRELOAD','PYTHONHOME','PYTHONPATH')}
        fp=out/'LBT_WATER_FIXED_SOLVENT.f90';fp.write_text(source);binary=out/'lbt_water_fixed_solvent'
        command=scfg['compiler_command']+[str(fp),'-o',str(binary)]+scfg['link_arguments']
        with (out/'compile.out').open('x') as so,(out/'compile.err').open('x') as se:proc=subprocess.run(command,cwd=out,env=env,stdout=so,stderr=se)
        result['build']={'command':command,'returncode':proc.returncode,'fortran':pin(fp),'build_inputs':scfg['build_inputs'],'runtime_libraries':build['linked_libraries']}
        require(proc.returncode==0,'Native259 wrapper build failed');bp=pin(binary);result['build']['binary']=bp
        for c,r in cells:
            if c is None:continue
            try:
                for lib in build['linked_libraries']:verify(lib)
                result['new_ALPB_calls_attempted']+=1;r['attempted']=True;d=out/r['task_id'];s=solvent(h,bp,c,d,env)
                v={'MACE_eV':c['MACE_eV'],'ALPB_Eh':s['total_hartree'],'ALPB_parts_Eh':s['components_hartree']}
                r.update(status='complete',ALPB=pin(d/'RESULT.json'),native_energies=v,
                    component_energies_kcal_mol={'MACE_vacuum':v['MACE_eV']*h.EV,'fixed_FF_charge_ALPB':v['ALPB_Eh']*h.HA},partition_defined_zeros=['GFN2_field_difference','MM_ledger'])
            except Exception as exc:r.update(failure=repr(exc),traceback=traceback.format_exc())
        require(result['new_ALPB_calls_attempted']<=12,'Twelve-call ceiling exceeded')
        result.update(accepted_cells=sum(r['status']=='complete' for r in result['rows']),geometry_rejected_cells=4,all16_available=False,summary=summarize(h,result['rows'],old))
        result['status']='complete' if result['accepted_cells']==12 else 'unavailable'
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-tick;result['application_allocated_CPU_seconds']=result['elapsed_seconds'];write(out/'RESULT.json',result)
    lines=['# Matched-FF rescore of the closed one-water panel','',f"Status: {result['status']}; native ALPB attempts {result['new_ALPB_calls_attempted']}/12; four rejected slots retained.",'',
        '| Source | Plane | Status | H kcal/mol |','|---|---:|---|---:|']
    for r in result.get('summary',{}).get('transfers',[]):lines.append(f"| {r['source_id']} | {r['plane_angle_degree']} | {r['status']} | {r.get('H_Dy_minus_La_transfer',{}).get('composite_kcal_mol','unavailable')} |")
    lines+=['','All16 rows, raw259 MACE/native solvent components, source menus and old outcomes remain in RESULT.json.',
        'Negative H favors Dy relatively for inner transfer. Outer is not bulk water; no occupancy, dry/wet raw minimum or affinity claim.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
