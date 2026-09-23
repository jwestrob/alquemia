"""Calibrate standalone GFN2 static/minimal PQQ descriptors on declared references."""
from __future__ import annotations
import argparse
import concurrent.futures
import fcntl
import json
import os
from pathlib import Path
import shutil
import time
from functools import lru_cache
import standalone_xtb as backend
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
from nikasha_pool import choose_rows
from nikasha_pool_compare import extrema_reference,method_record

PROTOCOL='standalone_xtb671_canonical25_static_and_minimal_reference_v1'
CANDIDATES=('origin','adaptive_Ca','adaptive_La')
CRYSTALS=('1H4I','4MAE','1KB0')
@lru_cache(maxsize=None)
def cached(path):return read_json(path)
def source_task(cell,medium):
    low=cell['low'][medium];m=cached(str(verify(low['manifest'])))
    return next(t for t in m.get('all_tasks',m['tasks']) if t['task_id']==low['task_id'])

def prepare(source,pilot,designation,agreement,output):
    archive=read_json(source);original=read_json(designation);pr=read_json(pilot);pm=read_json(verify(pr['manifest']))
    if archive['denominator']!=30 or pr['complete']!=224 or pm['settings']!=backend.SETTINGS:raise InvalidArtifact('source protocol/scope differs')
    canonical={r['case_id']:r['expected_class'] for r in original['rows'] if r['representation']=='context' and r['role']=='calibration'}
    if len(canonical)!=25 or set(canonical)&set(CRYSTALS):raise InvalidArtifact('designated calibration differs')
    cases=[c for c in archive['cases'] if c['case_id'] in set(canonical)|set(CRYSTALS)]
    if len(cases)!=28:raise InvalidArtifact('missing/duplicate declared source')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        d=impl/p.name;shutil.copyfile(p,d);pins[p.name]=record(d)
    params=out/'parameters';params.mkdir();pp=verify(pm['parameter']);shutil.copyfile(pp,params/pp.name)
    previous={t['task_id']:t for t in pm['tasks']};index={(r['case_id'],r['candidate'],r['metal'],r['medium']):r for r in pr['rows'] if r['accuracy']=='0.02' and r['kind']=='pool'}
    tasks=[];cells=[];rows=[]
    for c in cases:
        cid=c['case_id'];label=canonical.get(cid,c['old_result']['expected_class'])
        if c['old_result']['expected_class']!=label:raise InvalidArtifact('source label differs from designation')
        rows.append({'case_id':cid,'expected_class_for_report_only':label,'role':'calibration' if cid in canonical else 'consumed_crystal_transfer','native_source':c})
        for q in CANDIDATES:
            rep=c.get('aliases',{}).get(q,{}).get('representative',q)
            for z in ('Ca','La'):
                sourcecell=c['matrix'][z][rep]
                if sourcecell['status']!='complete':raise InvalidArtifact('source native MACE/physical state missing')
                for medium in ('vacuum','alpb'):
                    src=source_task(sourcecell,medium);tid='__'.join((cid,q,z,medium));xp=sourcecell['xyz']
                    if xyz(verify(src['xyz']))!=xyz(verify(xp)):raise InvalidArtifact('source receipt geometry differs')
                    t={'task_id':tid,'kind':'pool','cell_id':tid,'case_id':cid,'candidate':q,'metal':z,'medium':medium,
                        'charge':src['charge'],'multiplicity':src['multiplicity'],'gradient_requested':False,'accuracy':'0.02',
                        'offset_radian':0. if q=='origin' else None,'source_xyz':xp,'source_task':src,'MACE_eV':sourcecell['components']['MACE_eV']}
                    prior=index.get((cid,q,z,medium));reuse=None
                    if prior:
                        old=previous[prior['task_id']]
                        if prior['status']!='complete' or any(old[k]!=t[k] for k in ('charge','multiplicity','accuracy','medium')) or xyz(verify(old['xyz']))!=xyz(verify(xp)):raise InvalidArtifact('pilot reuse incompatible')
                        reuse={'collection':record(pilot),'row':prior,'source_task':old}
                    if reuse:
                        cells.append({**t,'xyz':xp,'reuse':reuse});continue
                    d=out/'tasks'/tid;d.mkdir(parents=True);shutil.copyfile(verify(xp),d/'core.xyz');(d/'xcontrol').write_text(backend.control())
                    t.update(xyz=record(d/'core.xyz'),input=record(d/'xcontrol'),directory=str(d))
                    tasks.append(t);cells.append({**t,'reuse':None})
    m={k:pm[k] for k in ('executable','installed_parameter','package')}
    m.update(protocol_id=PROTOCOL,backend_settings=backend.SETTINGS,source=record(source),pilot=record(pilot),designation=record(designation),agreement=record(agreement),
        implementation=pins,parameter=record(params/pp.name),cases=rows,cells=cells,tasks=tasks,candidate_ids=list(CANDIDATES),
        case_denominator=28,cell_denominator=336,reused_cells=sum(bool(t['reuse']) for t in cells),maximum_new_calls=300,new_MACE_DFT_optimization_calls=0,
        calibrated_reference=None,production_changed=False)
    if len(cells)!=336 or len(tasks)!=300 or m['reused_cells']!=36:raise InvalidArtifact('frozen logical/reuse/new count differs')
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v

def validate(manifest,fresh=False):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['backend_settings']!=backend.SETTINGS or len(m['cells'])!=336 or len(m['tasks'])!=300 or m['candidate_ids']!=list(CANDIDATES):raise InvalidArtifact('frozen scope/backend differs')
    for k in ('executable','parameter','installed_parameter','package','source','pilot','designation','agreement'):verify(m[k])
    for p in m['implementation'].values():verify(p)
    if verify(m['parameter']).read_bytes()!=verify(m['installed_parameter']).read_bytes():raise InvalidArtifact('parameter changed')
    if {p.name for p in verify(m['parameter']).parent.iterdir()}!={'param_gfn2-xtb.txt'}:raise InvalidArtifact('parameter/rc override')
    source=read_json(verify(m['source']));lookup={c['case_id']:c for c in source['cases']};original=read_json(verify(m['designation']))
    designated={r['case_id']:r['expected_class'] for r in original['rows'] if r['representation']=='context' and r['role']=='calibration'}
    if {c['case_id'] for c in m['cases']}!=set(designated)|set(CRYSTALS):raise InvalidArtifact('calibration/crystal membership differs')
    index={(t['case_id'],t['candidate'],t['metal'],t['medium']):t for t in m['cells']}
    if len(index)!=336:raise InvalidArtifact('duplicate cell')
    for c in m['cases']:
        if c['native_source']!=lookup[c['case_id']]:raise InvalidArtifact('native source/failed history changed')
        expected=designated.get(c['case_id'],lookup[c['case_id']]['old_result']['expected_class'])
        if c['expected_class_for_report_only']!=expected or c['role']!=('calibration' if c['case_id'] in designated else 'consumed_crystal_transfer'):raise InvalidArtifact('label/role differs')
        for q in CANDIDATES:
            ca,la=[index[c['case_id'],q,z,'vacuum'] for z in ('Ca','La')]
            paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
    new=[]
    for t in m['cells']:
        c=lookup[t['case_id']];rep=c.get('aliases',{}).get(t['candidate'],{}).get('representative',t['candidate']);sc=c['matrix'][t['metal']][rep];src=source_task(sc,t['medium'])
        if t['MACE_eV']!=sc['components']['MACE_eV'] or t['source_task']!=src or t['accuracy']!='0.02' or t['multiplicity']!=1 or t['charge']!=src['charge'] or xyz(verify(t['xyz']))!=xyz(verify(sc['xyz'])):raise InvalidArtifact('physical/electronic input changed')
        if t['reuse']:
            reuse=t['reuse'];verify(reuse['collection']);actual=next(r for r in read_json(verify(reuse['collection']))['rows'] if r['task_id']==reuse['row']['task_id'])
            if actual!=reuse['row'] or actual['status']!='complete' or actual['accuracy']!='0.02':raise InvalidArtifact('not actual compatible reuse')
            verify(actual['receipt']);verify(actual['output'])
        else:
            new.append({k:v for k,v in t.items() if k!='reuse'})
            if verify(t['input']).read_text()!=backend.control():raise InvalidArtifact('recipe differs')
            directory=Path(t['directory']);directory.relative_to(Path(manifest).resolve().parent/'tasks')
            if fresh and {p.name for p in directory.iterdir()}!={'core.xyz','xcontrol'}:raise InvalidArtifact('nonfresh task')
    if new!=m['tasks']:raise InvalidArtifact('executable task subset differs')
    return {'status':'validated','manifest':record(manifest),'sources':28,'logical_cells':336,'actual_reuses':36,'new_calls':300,'new_MACE_DFT':0}

def execute(manifest):
    validate(manifest,True);m=read_json(manifest);root=Path(manifest).parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE',0))!=64:raise InvalidArtifact('64CPU allocation required')
    with (root/'execute.lock').open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);begin=time.monotonic();receipts=[];errors=[]
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
            futures={ex.submit(backend.execute_cell,manifest,m,t):t for t in m['tasks']}
            for future in concurrent.futures.as_completed(futures):
                try:receipts.append(future.result())
                except Exception as e:errors.append({'task_id':futures[future]['task_id'],'error':str(e)})
        elapsed=time.monotonic()-begin;write_new(root/'EXECUTION.json',{'manifest':record(manifest),'receipts':receipts,'errors':errors,'wall_seconds':elapsed,'allocated_core_seconds':elapsed*64,'job_id':os.environ['SLURM_JOB_ID'],'GPU_seconds':0})
    return {'receipts':len(receipts),'errors':len(errors)}

def collect(manifest,output):
    validate(manifest);m=read_json(manifest);rows=[]
    for t in m['cells']:
        row=backend.parse_task(manifest,t) if not t['reuse'] else {**t['reuse']['row'],'task_id':t['task_id'],'cell_id':t['cell_id']}
        rows.append({**row,'reused':bool(t['reuse']),'reuse':t['reuse']})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'rows':rows,'denominator':336,'complete':sum(r['status']=='complete' for r in rows),'new_attempts':300,'reused':36,'implementation':record(__file__)}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}

def calibrate(collection,output):
    r=read_json(collection);m=read_json(verify(r['manifest']));validate(verify(r['manifest']));idx={(t['case_id'],t['candidate'],t['metal'],t['medium']):t for t in r['rows']};states={t['task_id']:t for t in m['cells']};rows=[]
    for c in m['cases']:
        cid=c['case_id'];matrix={z:{} for z in ('Ca','La')}
        for z in matrix:
            for q in CANDIDATES:
                pair={s:idx[cid,q,z,s] for s in ('vacuum','alpb')};good=all(t['status']=='complete' for t in pair.values())
                matrix[z][q]={'status':'complete' if good else 'unavailable','components':{'MACE_eV':states[pair['vacuum']['task_id']]['MACE_eV'],
                    'GFN2_vacuum_hartree':pair['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':pair['alpb']['energy_hartree']} if good else None}
        static=choose_rows(matrix,('origin',));pool=choose_rows(matrix,CANDIDATES)
        values={'static':static['operational']['composite_R_model_kcal_mol'] if static['operational'] else None,
                'minimal':pool['operational']['composite_R_model_kcal_mol'] if pool['operational'] else None,
                'minimal_mathematical':pool['mathematical']['composite_R_model_kcal_mol'] if pool['mathematical'] else None}
        rows.append({'case_id':cid,'expected_class':c['expected_class_for_report_only'],'role':c['role'],'values':values,'matrix':matrix,'pool':pool,
                     'native_static':c['native_source']['old_result']['R0'],'native_adaptive':c['native_source']['pool']})
    variants={}
    for v in ('static','minimal','minimal_mathematical'):
        cal=[{'case_id':x['case_id'],'expected_class':x['expected_class'],'R_model_kcal_mol':x['values'][v]} for x in rows if x['role']=='calibration']
        variants[v]=extrema_reference(cal,v,'standalone_xtb671_canonical25_v1')
    for row in rows:row['decisions']={v:method_record(row['values'][v],variants[v]['bands'],row['expected_class']) for v in variants}
    result={'protocol_id':PROTOCOL,'collection':record(collection),'manifest':r['manifest'],'variants':variants,'rows':rows,'backend_settings':backend.SETTINGS,
        'reported_accuracy':'0.02','calibration_denominator':25,'crystal_denominator':3,'noncanonical_folds_used_for_calibration':False,
        'new_calls_in_calibration':0,'production_changed':False,'implementation':record(__file__),'candidate_ids':list(CANDIDATES)}
    write_new(output,result);return {v:{k:r[k] for k in ('status','bands','gap_model_kcal_mol')} for v,r in variants.items()}

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,names in {'prepare':('source','pilot','designation','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'calibrate':('collection','output')}.items():
        a=sub.add_parser(op)
        for n in names:a.add_argument('--'+n,required=True,type=Path)
        if op=='validate':a.add_argument('--fresh',action='store_true')
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
if __name__=='__main__':main()
