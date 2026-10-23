#!/usr/bin/env python3
"""One serial export per remaining fixture; numerical admission belongs to scout."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import time

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
SCOUT='aquo9_central_La__dimer_01'
MAX_CALLS=53

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):require(pin(q['path'])==q,'Changed pin: '+q['path']);return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def save(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def module(q,name):
    spec=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def tokens(text):return [s for raw in text.splitlines() if (s:=raw.split('#',1)[0].split())]
def render(key):
    lines=tokens(key); names=[s[0].lower() for s in lines]
    require(not set(names)&{'external-field','exfld-freq','response-matrix-export','polar-pair-response','exchpol','exchange-polar'},'Nonbaseline/export/field keyword present')
    additions=[]
    for name in ('polar-print','save-uinduce'):
        require(names.count(name)<=1 and all(len(s)==1 for s in lines if s[0].lower()==name),'Ambiguous reporting keyword')
        if name not in names:additions.append(name)
    additions.append('RESPONSE-MATRIX-EXPORT')
    return key+'\n# Reporting only: frozen remaining baseline operator inventory\n'+'\n'.join(additions)+'\n',additions
def helpers(cfg):
    return {k:module(cfg[k],'operator_inventory_'+k) for k in ('moment_adapter','induction_helper','native_helper','resolver','scout_helper')}
def configured(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Runner source/specification differs')
    require(cfg['maximum_native_programs']==53 and cfg['workers']==1,'Finite serial scope differs')
    for k in ('inventory','export_source_manifest','moment_adapter','induction_helper','native_helper','resolver','scout_helper'):verify(cfg[k])
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread required')

def source_rows(cfg,h):
    inv=read(cfg['inventory']);files=inv['artifacts']
    require(inv['status']=='source_inventory_complete_unreleased' and inv['counts']['total_baseline_exports']==54 and inv['first_gate']==SCOUT,'Inventory scope differs')
    require(len(inv['rows'])==54 and inv['remaining_after_actual_parity']==[k for k in inv['rows'] if k!=SCOUT],'Remaining denominator differs')
    data=read(files['development_data']);mr=read(files['model_moments']);mm=read(files['moment_manifest']);tr=read(files['trimer_model']);tm=read(files['trimer_manifest']);fm=read(files['field_manifest'])
    require(data['provenance']['candidate_result']==inv['candidate_result']==mr['candidate_result']==mm['candidate_result']==tm['candidate_result'],'Baseline model lineage differs')
    require(data['baseline_parameters']==inv['selected_parameters']==mm['selected_parameters']==tm['selected_parameters'],'Selected parameters differ')
    rows={}
    for cid,item in inv['rows'].items():
        q={k:files[v] for k,v in item['input'].items()};rr=[s.split() for s in verify(q['xyz']).read_text().splitlines()[1:] if s.strip()]
        types=[int(s[5]) for s in rr];positions=[[float(x) for x in s[2:5]] for s in rr]
        require(len(rr)==item['atoms'] and [int(s[0]) for s in rr]==list(range(1,len(rr)+1)),'Atom ordering differs')
        require(types[0] in (500,501) and types[1:]==[36,37,37]*((len(rr)-1)//3),'Water-only type inventory differs')
        key=verify(q['key']).read_text();prm=verify(q['parameters']).read_text();params=h['resolver'].parse_parameters(prm,key)
        masses=[params['atom'][t]['mass'] for t in types];mass=math.fsum(masses)
        com=[math.fsum(m*p[j] for m,p in zip(masses,positions))/mass for j in range(3)]
        row=dict(id=cid,atoms=len(rr),types=types,native_names=[params['atom'][t]['symbol'][:3] for t in types],positions_A=positions,
                 alpha_A3=[params['polarize'][t]['alpha_A3'] for t in types],masses_amu=masses,COM_A=com,charge=3,
                 kind='ion' if len(rr)==1 else item['category'],energy_replay_kcal_mol=None,source_input=q,
                 original_inventory_row=item,parent_geometry_id=item['parent'],parent_atom_indices=item['parent_atom_indices'])
        require(all(math.isfinite(a) and a>=1e-8 for a in row['alpha_A3']),'Unsupported native site alpha')
        if cid in mr['rows']:
            baseline=copy.deepcopy(mr['rows'][cid]['moments']);authority={'result':files['model_moments'],'cell':cid,'field':'moments'}
        elif cid in tr['rows']:
            baseline=copy.deepcopy(tr['rows'][cid]['parsed_moments']);authority={'result':files['trimer_model'],'cell':cid,'field':'parsed_moments'}
        else:
            old=data['baseline_model'][cid];baseline=h['moment_adapter'].parse_energy(verify(old['native']['analyze_output']))
            require(baseline['energy_kcal_mol']==old['energy'] and baseline['components']==old['native']['components'],'Original E/components authority differs')
            authority={'result':files['development_data'],'cell':cid,'field':'baseline_model','native_authority':old['provenance']['native_authority']}
        if cid in data['baseline_model']:
            old=data['baseline_model'][cid];require(baseline['energy_kcal_mol']==old['energy'] and baseline['components']==old['native']['components'],'Joined original-development baseline differs')
        row.update(baseline=baseline,baseline_authority=authority,baseline_u_available='induced_Debye' in baseline)
        rows[cid]=row
    require(sum(r['baseline_u_available'] for r in rows.values())==48,'Saved induced denominator differs')
    return inv,rows,fm['conversion']

def prepare(cfg,out):
    start=time.monotonic();cpu=time.process_time();configured(cfg);h=helpers(cfg);inv,rows,units=source_rows(cfg,h)
    for cid,row in rows.items():
        folder=out/'inputs'/cid;folder.mkdir(parents=True)
        for k,name in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['source_input'][k]),folder/name)
        rendered,added=render(verify(row['source_input']['key']).read_text());(folder/'system.key').write_text(rendered)
        row['input']={k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))};row['reporting_flags_added']=added
    m=dict(status='prepared_unreleased',implementation=cfg['implementation'],specification=cfg['specification'],inventory=cfg['inventory'],
           export_source_manifest=cfg['export_source_manifest'],scout_helper=cfg['scout_helper'],rows=rows,units=units,
           helpers={k:cfg[k] for k in ('moment_adapter','induction_helper','native_helper','resolver')},
           scout_id=SCOUT,remaining_cell_ids=inv['remaining_after_actual_parity'],maximum_native_programs=53,workers=1,
           snapshots_per_fixture=2,total_fixture_count=54,total_snapshot_count=108,new_native_calls=0,
           cost={'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'local_threads':1})
    save(out/'MANIFEST.json',m);return m

def context(cfg):
    configured(cfg);m=read(cfg['manifest']);inv=read(cfg['inventory'])
    for k in ('implementation','specification','inventory','export_source_manifest','scout_helper'):require(m[k]==cfg[k],'Prepared source binding differs: '+k)
    require(m['helpers']=={k:cfg[k] for k in ('moment_adapter','induction_helper','native_helper','resolver')},'Prepared helper implementation differs')
    require(m['status']=='prepared_unreleased' and len(m['rows'])==54 and list(m['rows'])==list(inv['rows']) and m['remaining_cell_ids']==inv['remaining_after_actual_parity'],'Prepared full denominator differs')
    require(m['scout_id']==SCOUT and m['maximum_native_programs']==53 and m['total_snapshot_count']==108 and m['workers']==1,'Prepared caps differ')
    for cid,row in m['rows'].items():
        require(row['source_input']=={k:inv['artifacts'][v] for k,v in inv['rows'][cid]['input'].items()},'Original input differs')
        for k in ('xyz','parameters'):require(verify(row['input'][k]).read_bytes()==verify(row['source_input'][k]).read_bytes(),'Physical bytes differ')
        key,flags=render(verify(row['source_input']['key']).read_text())
        require(verify(row['input']['key']).read_text()==key and row['reporting_flags_added']==flags,'Reporting-only delta differs')
    return m,helpers(cfg)

def release_gate(cfg,m,out,h):
    # The public scout admission below supplies all raw-build/native/math checks.
    require(all(isinstance(cfg.get(k),dict) for k in ('root_release','scout_review','scout_config','scout_result','build_receipt','analyze')),'Actual reviewed build/scout/root release remain unbound')
    release=read(cfg['root_release']);review=read(cfg['scout_review'])
    require(release['status']=='approved_remaining_baseline_matrix_exports','Root release status differs')
    for k in ('implementation','manifest','inventory','export_source_manifest','scout_helper','build_receipt','analyze','scout_config','scout_result','scout_review'):
        require(release[k]==cfg[k],'Root release binding differs: '+k)
    require(release['remaining_cell_ids']==m['remaining_cell_ids'] and release['maximum_native_programs']==53 and release['workers']==1,'Root finite inventory differs')
    require(release['execution_root']==str(out.resolve()) and release['new_fields'] is False and release['fit'] is False,'Root execution scope differs')
    require(review['status']=='passed_baseline_response_matrix_scout_review' and review['result']==cfg['scout_result'] and review['implementation']==cfg['scout_helper'],'Actual scout review differs')
    sc=read(cfg['scout_config']);sr=read(cfg['scout_result']);sm=read(sc['manifest'])
    require(sm['row']['id']==SCOUT and sm['units']==m['units'],'Scout fixture/units differ')
    for k in ('xyz','key','parameters'):
        require(verify(sm['row']['input'][k]).read_bytes()==verify(m['rows'][SCOUT]['source_input'][k]).read_bytes(),'Scout baseline input differs')
    admitted=h['scout_helper'].admit_result(sc,sm,sr)
    require(sc['build_receipt']==cfg['build_receipt'] and sc['analyze']==cfg['analyze'] and sc['source_manifest']==cfg['export_source_manifest'],'Scout native producer differs')
    require(admitted['status']=='admitted_response_matrix_scout','Scout numerical admission failed')
    return dict(status='complete',inherited=True,producer_result=cfg['scout_result'],producer_configuration=cfg['scout_config'],
                producer_review=cfg['scout_review'],admission=admitted)

def execute_remaining(cfg,m,out,h,reused):
    start=time.monotonic();cpu=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    result=dict(status='baseline_matrix_inventory_running',implementation=pin(__file__),config=cfg,manifest=cfg['manifest'],root_release=cfg['root_release'],
                rows={cid:{'status':'not_attempted'} for cid in m['rows']},native_program_invocations=0,maximum_native_programs=53,
                required_fixtures=54,required_snapshots=108,inherited_scout_native_programs=2,global_stop=False,error=None)
    result['rows'][SCOUT]=reused
    def persist():save(out/'RESULT.json',result)
    persist()
    for cid in m['remaining_cell_ids']:
        row=m['rows'][cid];record={'status':'running','receipts':[]};result['rows'][cid]=record
        try:
            require(not result['global_stop'] and result['native_program_invocations']<53,'Stop/native ceiling reached')
            folder=out/'cells'/cid;folder.mkdir(parents=True)
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),folder/name)
            record['runtime_input']={k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            save(folder/'START.json',{'implementation':cfg['implementation'],'cell_id':cid,'runtime_input':record['runtime_input'],
                                     'analyze':cfg['analyze'],'root_release':cfg['root_release'],'job_id':os.environ.get('SLURM_JOB_ID'),
                                     'threads':{k:os.environ.get(k) for k in THREADS}})
            record['start']=pin(folder/'START.json')
            result['native_program_invocations']+=1;persist()
            h['native_helper'].invoke([cfg['analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',record['receipts'])
            moment=h['moment_adapter'];output=folder/'analyze.out'
            parsed=moment.parse_moments(output,moment.induced_file(output),row,h['induction_helper'].induction)
            parity=h['scout_helper'].parity(parsed,row['baseline'],exact=False)
            require(parity['passed'],'Cached baseline parity failed')
            exports=sorted(folder.glob('response_matrix_*.dat'))
            require([p.name for p in exports]==['response_matrix_001.dat','response_matrix_002.dat'],'Exactly two export files required')
            names=re.findall(r'Response Matrix Written To\s*:\s*(\S+)',output.read_text())
            require(names==[p.name for p in exports],'Exactly two native export announcements required')
            numerical=h['scout_helper'].qualify_exports(exports,row,m['units'],parsed)
            require(numerical['passed'],'Matrix numerical qualification failed')
            record.update(status='complete',moments=parsed,baseline_parity=parity,numerical=numerical,exports=[pin(p) for p in exports])
        except Exception as e:
            result['global_stop']=True;record.update(status='failed',error=f'{type(e).__name__}: {e}')
            record['exports_present']=[pin(p) for p in (out/'cells'/cid).glob('response_matrix_*.dat')]
            result['error']={'cell_id':cid,'error':record['error']};persist();break
        persist()
    result['complete_fixtures']=sum(r['status']=='complete' for r in result['rows'].values())
    result['admitted_snapshots']=2*result['complete_fixtures']
    result['total_native_program_invocations_including_scout']=result['native_program_invocations']+2
    result['status']='baseline_matrix_inventory_complete' if result['complete_fixtures']==54 and result['native_program_invocations']==53 and result['error'] is None else 'baseline_matrix_inventory_incomplete'
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result['cost']={'wall_seconds':time.monotonic()-start,'parent_cpu_seconds':time.process_time()-cpu,
                    'child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'local_threads':1,
                    'inherited_scout_cost':'retained in producer_result; excluded from new runtime totals'}
    persist();return result

def run(cfg,out):
    start=time.monotonic();cpu=time.process_time()
    m,h=context(cfg);reused=release_gate(cfg,m,out,h)
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_NTASKS','1')=='1' and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Exactly one scheduled CPU/task required')
    with (out/'LAUNCH_LOCK.json').open('x') as handle:json.dump({'job_id':os.environ['SLURM_JOB_ID'],'configuration':cfg,'new_native_limit':53},handle,indent=2)
    result=execute_remaining(cfg,m,out,h,reused)
    result['cost']['total_wall_including_scout_admission_seconds']=time.monotonic()-start
    result['cost']['total_parent_cpu_including_scout_admission_seconds']=time.process_time()-cpu
    save(out/'RESULT.json',result);return result

def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    cfg=read(pin(a.config));out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:
        result=prepare(cfg,out) if a.mode=='prepare' else run(cfg,out)
        print(json.dumps({'status':result['status'],'out':str(out)}))
        if a.mode=='run' and result['status']!='baseline_matrix_inventory_complete':raise SystemExit(1)
    except Exception as e:
        save(out/'ADMISSION_FAILURE.json',{'status':'unavailable','error':f'{type(e).__name__}: {e}','new_native_calls':0 if not (out/'RESULT.json').exists() else None})
        raise

if __name__=='__main__':main()
