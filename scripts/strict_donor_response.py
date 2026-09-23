"""Fixed three-context donor works with qualified strict native scalar stopping."""
from __future__ import annotations
import argparse,importlib.util,json,re,shutil
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,energy,HA_TO_KCAL
from mace_hybrid import EV_TO_KCAL
from strict_native_pool import recipe,collect_fresh
from compact_solvation import input_text,completed
from run_orca_task_manifest import load_manifest_tasks
from structure_informed_starts import scf_details
from union_adaptive import snapshot

PROTOCOL='native_GFN2_TolE1e10_archived_three_donor_paths_v1'
CASES=('4MAE','PQQSEQ_83440678cbbd658047c9','PQQSEQ_07ab500e3df76b30d71c')
ANGLES=(0.0,-.2,.2)
ROLE='extra_acidic_ligand_homolog'
FINAL_SHA='c6d19ebe2c4f5be178ee12551877de214a5ef8a6e6acbe69148091425ba787b4'


def data(pin):return read_json(verify(pin))

def same_xyz(a,b):
    xa=xyz(verify(a));xb=xyz(verify(b))
    if [r[0] for r in xa]!=[r[0] for r in xb] or np.max(np.abs(np.asarray([r[1:] for r in xa])-np.asarray([r[1:] for r in xb])))>1e-12:raise InvalidArtifact('exact atom/coordinate identity differs')


def inventory(components,strict_reuse):
    c=read_json(components)
    if c['analysis']['sha256']!=FINAL_SHA:raise InvalidArtifact('fixed original comparison differs')
    a=data(c['analysis']);m=data(a['manifest']);d=data(m['design']);low=data(a['low_collection']);lm=data(m['low_manifest'])
    df=data(a['DFT_collection']);native={r['point_id']:r for r in df['rows']};native.update(d['reused_DFT_centers'])
    mace={r['task_id']:r for r in data(a['MACE_result'])['rows']};old={r['task_id']:r for r in low['rows']}
    sr=read_json(strict_reuse);sm=data(sr['manifest']);reuse=[];rows=[];points=[]
    if sr['branch']!='fresh':raise InvalidArtifact('only fresh strict reuse allowed')
    for cid in CASES:
        for angle in ANGLES:
            pair=[]
            for z in ('Ca','La'):
                p=next(p for p in d['points'] if p['case_id']==cid and p['metal']==z and (p['point']=='origin' if angle==0 else p['role']==ROLE and p['angle_radian']==angle))
                actual=next(v for v in a['points'] if v['task_id']==p['task_id']);nm=mace[p['task_id']];nd=native[p['task_id']]
                if actual['status']!='complete' or nm['status']!='complete' or nd.get('status','complete')!='complete' or nd['energy_hartree']!=actual['DFT_hartree'] or nm['energy_eV']!=actual['MACE_eV']:raise InvalidArtifact('actual archived endpoint unavailable or differs')
                same_xyz(p['xyz'],nm['xyz']);dr=data(nd['receipt']);same_xyz(p['xyz'],dr['artifacts']['xyz'])
                inp=verify(dr['artifacts']['template_input']).read_text()
                if not re.search(r'^\* xyzfile '+str(p['charge'])+r' 1 ',inp,re.M) or energy(verify(nd['output']))!=nd['energy_hartree']:raise InvalidArtifact('native DFT state/energy mismatch')
                point={**p,'old':actual,'MACE':nm,'DFT':nd};points.append(point);pair.append(p)
                for medium in ('vacuum','alpb'):
                    tid=p['task_id']+'__'+medium;t=next(t for t in lm['all_tasks'] if t['task_id']==tid);v=old[tid];receipt=data(v['receipt']);same_xyz(p['xyz'],t['xyz'])
                    if v['status']!='complete' or v['charge_sanity_status']!='pass' or verify(t['input']).read_text()!=input_text(p['charge'],1,medium,'native') or energy(verify(v['output']))!=v['energy_hartree'] or not receipt['normal_termination'] or not receipt['scf_converged']:raise InvalidArtifact('old native input/receipt mismatch')
                    detail=scf_details(verify(v['output']).read_text())
                    src={'xyz':t['xyz'],'input':t['input'],'output':v['output'],'receipt':v['receipt'],'source_cell':v,'source_energy_hartree':v['energy_hartree'],
                         'charge':t['charge'],'multiplicity':1,'electron_count':detail['electrons'],'parameter_export':v['parameter_export'],'orca':receipt['orca_executable']}
                    row={'task_id':tid,'case_id':cid,'candidate':p['point'],'point_id':p['task_id'],'angle_radian':angle,'metal':z,'medium':medium,'charge':p['charge'],'multiplicity':1,'source':src}
                    matches=[v for v in sr['rows'] if v['case_id']==cid and v['candidate']==p['point'] and v['metal']==z and v['medium']==medium]
                    if matches:
                        r=matches[0];rt=next(t for t in sm['tasks'] if t['task_id']==r['task_id']);same_xyz(t['xyz'],rt['xyz']);actualpin=completed(verify(sr['manifest']),r['task_id'])
                        if r['status']!='complete' or r['observed_TolE_hartree']!=1e-10 or verify(rt['input']).read_text()!=recipe(p['charge'],medium,'fresh') or r['initial_guess']!=['SAD'] or actualpin!=r['actual'] or data(r['actual']['receipt'])['parallelism']['nprocs']!=1 or r['audit']['parameter_export']['sha256']!=src['parameter_export']['sha256']:raise InvalidArtifact('strict reuse incompatible')
                        reuse.append({**row,'actual_row':r,'collection':record(strict_reuse)})
                    else:rows.append(row)
            ca,la=pair;xc=xyz(verify(ca['xyz']));xl=xyz(verify(la['xyz']))
            if xc[0][0]!='Ca' or xl[0][0]!='La' or xc[1:]!=xl[1:] or xc[0][1:]!=xl[0][1:] or la['charge']-ca['charge']!=1:raise InvalidArtifact('paired source coordinates/state differ')
    if len(rows)!=32 or len(reuse)!=4 or len(points)!=18:raise InvalidArtifact('expected32fresh/4exactreuse/18points')
    return {'analysis':c['analysis'],'original_manifest':a['manifest'],'points':points,'fresh':rows,'reused':reuse}


def prepare(components,strict_reuse,agreement,output):
    inv=inventory(components,strict_reuse);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[]
    for r in inv['fresh']:
        d=out/'tasks'/r['task_id'];d.mkdir(parents=True);shutil.copyfile(verify(r['source']['xyz']),d/'core.xyz');(d/'endpoint.inp').write_text(recipe(r['charge'],r['medium'],'fresh'))
        tasks.append({**r,'seed_source':None,'gradient_requested':False,'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),
                      'output_path':str(d/'endpoint.out'),'active_seed_path':str(d/'endpoint.runtime.xtbw'),'active_gbw_path':str(d/'endpoint.runtime.gbw')})
    impl=snapshot(out/'implementation');m={'protocol_id':PROTOCOL,'branch':'fresh','tasks':tasks,'reused':inv['reused'],'points':inv['points'],
        'original_analysis':inv['analysis'],'original_manifest':inv['original_manifest'],'components':record(components),'strict_reuse':record(strict_reuse),'agreement':record(agreement),
        'orca':tasks[0]['source']['orca'],'implementation':impl,'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},
        'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
        'logical_cell_denominator':36,'fresh_calls':32,'reused_calls':4,'new_MACE_DFT_searches':0,'production_changed':False}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp)
    if m['protocol_id']!=PROTOCOL or m['branch']!='fresh' or len(m['tasks'])!=32 or len(m['reused'])!=4 or len(m['points'])!=18 or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}:raise InvalidArtifact('fixed scope/profile differs')
    for k in ('components','strict_reuse','agreement','original_analysis','original_manifest','orca'):verify(m[k])
    for pin in m['implementation'].values():verify(pin)
    required={(c,a,z,s) for c in CASES for a in ANGLES for z in ('Ca','La') for s in ('vacuum','alpb')}
    if {(t['case_id'],t['angle_radian'],t['metal'],t['medium']) for t in m['tasks']+m['reused']}!=required:raise InvalidArtifact('all36 source cells required')
    load_manifest_tasks(mp)
    for t in m['tasks']:
        if t['seed_source'] is not None or t['gradient_requested'] or t['charge']!=t['source']['charge'] or t['multiplicity']!=1 or verify(t['xyz']).read_bytes()!=verify(t['source']['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(t['charge'],t['medium'],'fresh'):raise InvalidArtifact('state/coordinates/recipe differs')
        d=Path(t['output_path']).parent
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('basename differs')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp'}:raise InvalidArtifact('nonfresh directory')
    return {'status':'validated','manifest':record(mp),'fresh_calls':32,'exact_reuses':4,'total_cells':36,'molecular_calls_in_validation':0}


def execute(manifest):
    path=Path(__file__).with_name('strict_native_pool.py');spec=importlib.util.spec_from_file_location('_strict_donor_executor',path);engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine);engine.validate=validate
    return engine.execute(manifest)


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def idx(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']} if p.exists() else {}
    before=idx('SEEDS_BEFORE.json');after=idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=collect_fresh(mp,t,before.get(t['task_id']),after.get(t['task_id']));r.update(point_id=t['point_id'],angle_radian=t['angle_radian'],reused=False)
        if r['status']=='complete':
            text=verify(r['actual']['output']).read_text();tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text);r['observed_TolE_hartree']=float(tol[0]) if len(tol)==1 else None
            if r['observed_TolE_hartree']!=1e-10 or data(r['actual']['receipt'])['parallelism']['nprocs']!=1:r.update(status='audit_failed',reason='effective tolerance/rank differs',energy_hartree=None)
        rows.append(r)
    rows.extend({**t['actual_row'],'point_id':t['point_id'],'angle_radian':t['angle_radian'],'reused':True,'reuse_collection':t['collection']} for t in m['reused'])
    index={(r['point_id'],r['medium']):r for r in rows};points=[]
    for p in m['points']:
        low={s:index[p['task_id'],s] for s in ('vacuum','alpb')};ok=all(r['status']=='complete' for r in low.values())
        transfer=(low['alpb']['energy_hartree']-low['vacuum']['energy_hartree'])*HA_TO_KCAL if ok else None
        points.append({'task_id':p['task_id'],'case_id':p['case_id'],'angle_radian':p['angle_radian'],'metal':p['metal'],'status':'complete' if ok else 'unavailable',
                       'MACE_eV':p['old']['MACE_eV'],'DFT_hartree':p['old']['DFT_hartree'],'old':p['old'],'strict_solvent_transfer_kcal_mol':transfer,
                       'strict_composite_kcal_mol':p['old']['MACE_eV']*EV_TO_KCAL+transfer if ok else None,'low':{s:r.get('actual') for s,r in low.items()}})
    pi={(p['case_id'],p['angle_radian'],p['metal']):p for p in points};works=[];diffs=[]
    for old in read_json(verify(m['components']))['rows']:
        endpoints={}
        for z in ('Ca','La'):
            p=pi[old['case_id'],old['angle_radian'],z];o=pi[old['case_id'],0.,z];ok=p['status']==o['status']=='complete'
            transfer=p['strict_solvent_transfer_kcal_mol']-o['strict_solvent_transfer_kcal_mol'] if ok else None
            mace=(p['MACE_eV']-o['MACE_eV'])*EV_TO_KCAL;dft=(p['DFT_hartree']-o['DFT_hartree'])*HA_TO_KCAL
            if abs(mace-old['endpoints'][z]['MACE'])>1e-7 or abs(dft-old['endpoints'][z]['DFT'])>1e-7:raise InvalidArtifact('archived work reconstruction differs')
            r={'case_id':old['case_id'],'angle_radian':old['angle_radian'],'metal':z,'status':'complete' if ok else 'unavailable','DFT':dft,'MACE':mace,
               'old_solvent':old['endpoints'][z]['solvent_transfer'],'old_composite':old['endpoints'][z]['composite'],
               'strict_solvent':transfer,'strict_composite':mace+transfer if ok else None,'strict_composite_error':mace+transfer-dft if ok else None}
            works.append(r);endpoints[z]=r
        ok=all(r['status']=='complete' for r in endpoints.values());diffs.append({'case_id':old['case_id'],'angle_radian':old['angle_radian'],'status':'complete' if ok else 'unavailable','old':old,
            'DFT':old['DFT_delta_R'],'MACE':old['MACE_delta_R'],'strict_solvent':endpoints['Ca']['strict_solvent']-endpoints['La']['strict_solvent'] if ok else None,
            'strict_composite':endpoints['Ca']['strict_composite']-endpoints['La']['strict_composite'] if ok else None,
            'strict_composite_error':endpoints['Ca']['strict_composite']-endpoints['La']['strict_composite']-old['DFT_delta_R'] if ok else None})
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'rows':rows,'points':points,'endpoint_works':works,'differentials':diffs,'complete_cells':sum(r['status']=='complete' for r in rows),
            'cell_denominator':36,'fresh_cells':32,'reused_cells':4,'endpoint_work_denominator':12,'differential_denominator':6,'response_accuracy_pass_threshold':None,'classifier_claim':None,'new_calls_in_collection':0}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','points','endpoint_works','differentials')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('components','strict_reuse','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        q=s.add_parser(op)
        for f in fields:q.add_argument('--'+f.replace('_','-'),type=Path,required=True)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
