"""Frozen standalone static scoring on the225 declared structural repeat ledger."""
from __future__ import annotations
import argparse
import concurrent.futures
from collections import Counter
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
from nikasha_pool_compare import extrema_reference,method_record,aggregate
from accommodation_fold_proposals import summarize_rows

PROTOCOL='standalone_xtb671_static_primary225_transfer_v1'
@lru_cache(maxsize=None)
def cached(path):return read_json(path)

def reference_check(path):
 r=read_json(path)
 if r['protocol_id']!='standalone_xtb671_canonical25_static_and_minimal_reference_v1' or r['backend_settings']!=backend.SETTINGS or r['reported_accuracy']!='0.02' or r['noncanonical_folds_used_for_calibration']:raise InvalidArtifact('reference backend/calibration differs')
 v=r['variants']['static']
 if v['status']!='available' or extrema_reference(v['rows'],'static','standalone_xtb671_canonical25_v1')!=v:raise InvalidArtifact('static reference unavailable/changed')
 for k in ('manifest','collection','implementation'):verify(r[k])
 return r

def verify_origin(origin,model):
 nr=read_json(verify(origin['MACE']));nm=cached(str(verify(nr['manifest'])));nt=next(t for t in nm['tasks'] if t['task_id']==nr['task_id'])
 if nr['status']!='computed' or nr['energy_eV']!=origin['MACE_eV'] or nm['model']!=model:raise InvalidArtifact('native MACE origin receipt/method differs')
 if (nt['charge'],nt['spin_multiplicity'])!=(origin['charge'],origin['multiplicity']) or xyz(verify(nt['xyz']))!=xyz(verify(origin['xyz'])):raise InvalidArtifact('native origin geometry/state differs')

def prepare(inventory,reference,pilot,ledger,dft,agreement,output):
 inv=read_json(inventory);ref=reference_check(reference);prior=read_json(ledger);df=read_json(dft);pc=read_json(pilot);pm=read_json(verify(pc['manifest']))
 source=read_json(verify(inv['parent_source']));canonical=read_json(verify(read_json(verify(ref['manifest']))['source']));canonical_manifest=read_json(verify(canonical['manifest']))
 if source['model']!=canonical_manifest['model'] or len(prior['rows'])!=225 or {c['case_id'] for c in prior['rows']}!={c['case_id'] for c in df['rows']}:raise InvalidArtifact('reference model/225ledger identity differs')
 origins={(t['case_id'],t['metal']):t for t in inv['native_MACE_origins']};oldtasks={t['task_id']:t for t in pm['tasks']};reuseindex={(r['case_id'],r['metal'],r['medium']):r for r in pc['rows'] if r['accuracy']=='0.02' and r['kind']=='pool' and r['candidate']=='origin'}
 out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
 for p in Path(__file__).parent.glob('*.py'):
  d=impl/p.name;shutil.copyfile(p,d);pins[p.name]=record(d)
 params=out/'parameters';params.mkdir();pp=verify(pm['parameter']);shutil.copyfile(pp,params/pp.name)
 rows=[];cells=[];tasks=[]
 for row in prior['rows']:
  cid=row['case_id'];pair={z:origins.get((cid,z)) for z in ('Ca','La')};status='prepared' if all(pair.values()) else 'unavailable'
  rows.append({'case_id':cid,'status':status,'reason':None if status=='prepared' else 'source_preparation_unavailable','origins':pair,'prior_row':row})
  if status!='prepared':continue
  for z,o in pair.items():
   verify_origin(o,source['model'])
   for medium in ('vacuum','alpb'):
    tid='__'.join((cid,'origin',z,medium));t={'task_id':tid,'kind':'pool','cell_id':tid,'case_id':cid,'candidate':'origin','metal':z,'medium':medium,
      'charge':o['charge'],'multiplicity':o['multiplicity'],'gradient_requested':False,'accuracy':'0.02','offset_radian':0.,'source_xyz':o['xyz'],'MACE_eV':o['MACE_eV']}
    old=reuseindex.get((cid,z,medium));reuse=None
    if old:
     ot=oldtasks[old['task_id']]
     if old['status']!='complete' or any(ot[k]!=t[k] for k in ('charge','multiplicity','accuracy','medium')) or xyz(verify(ot['xyz']))!=xyz(verify(o['xyz'])):raise InvalidArtifact('pilot reuse incompatible')
     reuse={'collection':record(pilot),'row':old,'source_task':ot}
    if reuse:cells.append({**t,'xyz':o['xyz'],'reuse':reuse});continue
    d=out/'tasks'/tid;d.mkdir(parents=True);shutil.copyfile(verify(o['xyz']),d/'core.xyz');(d/'xcontrol').write_text(backend.control())
    t.update(xyz=record(d/'core.xyz'),input=record(d/'xcontrol'),directory=str(d));tasks.append(t);cells.append({**t,'reuse':None})
 m={k:pm[k] for k in ('executable','installed_parameter','package')};m.update(protocol_id=PROTOCOL,backend_settings=backend.SETTINGS,inventory=record(inventory),
  reference=record(reference),pilot=record(pilot),ledger=record(ledger),DFT=record(dft),agreement=record(agreement),implementation=pins,model=source['model'],
  parameter=record(params/pp.name),cases=rows,cells=cells,tasks=tasks,case_denominator=225,prepared_cases=208,cell_denominator=832,reused_cells=4,
  maximum_new_calls=828,new_MACE_DFT_optimization_calls=0,production_changed=False)
 if len(cells)!=832 or len(tasks)!=828 or sum(c['status']=='prepared' for c in rows)!=208:raise InvalidArtifact('finite input/new/reuse count differs')
 mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v

def validate(manifest,fresh=False):
 m=read_json(manifest);ref=reference_check(verify(m['reference']));prior=read_json(verify(m['ledger']));inv=read_json(verify(m['inventory']))
 if m['protocol_id']!=PROTOCOL or m['backend_settings']!=backend.SETTINGS or len(m['cases'])!=225 or len(m['cells'])!=832 or len(m['tasks'])!=828:raise InvalidArtifact('scope/backend differs')
 for k in ('parameter','installed_parameter','executable','package','pilot','DFT','agreement'):verify(m[k])
 for p in m['implementation'].values():verify(p)
 if verify(m['parameter']).read_bytes()!=verify(m['installed_parameter']).read_bytes() or {p.name for p in verify(m['parameter']).parent.iterdir()}!={'param_gfn2-xtb.txt'}:raise InvalidArtifact('parameter/rc override')
 original={r['case_id']:r for r in prior['rows']};origins={(t['case_id'],t['metal']):t for t in inv['native_MACE_origins']};expected=[]
 for c in m['cases']:
  cid=c['case_id']
  if c['prior_row']!=original[cid] or c['origins']!={z:origins.get((cid,z)) for z in ('Ca','La')}:raise InvalidArtifact('source identity/history differs')
  if c['status']=='prepared':
   ca,la=(c['origins'][z] for z in ('Ca','La'));paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
   for z in ('Ca','La'):
    verify_origin(c['origins'][z],m['model'])
    expected.extend((cid,z,s) for s in ('vacuum','alpb'))
  elif any(c['origins'].values()):raise InvalidArtifact('supported source excluded')
 if {(t['case_id'],t['metal'],t['medium']) for t in m['cells']}!=set(expected) or len(expected)!=832:raise InvalidArtifact('cell membership differs')
 pilotrows={r['task_id']:r for r in read_json(verify(m['pilot']))['rows']};new=[]
 for t in m['cells']:
  o=origins[t['case_id'],t['metal']]
  if (t['charge'],t['multiplicity'],t['MACE_eV'])!=(o['charge'],o['multiplicity'],o['MACE_eV']) or t['accuracy']!='0.02' or t['candidate']!='origin' or xyz(verify(t['xyz']))!=xyz(verify(o['xyz'])):raise InvalidArtifact('origin state changed')
  if t['reuse']:
   old=t['reuse']['row'];ot=t['reuse']['source_task']
   if old!=pilotrows[old['task_id']] or old['status']!='complete' or old['accuracy']!='0.02' or xyz(verify(ot['xyz']))!=xyz(verify(t['xyz'])) or ot['charge']!=t['charge']:raise InvalidArtifact('incompatible actual reuse')
   verify(old['receipt']);verify(old['output'])
  else:
   new.append({k:v for k,v in t.items() if k!='reuse'})
   if verify(t['input']).read_text()!=backend.control():raise InvalidArtifact('recipe differs')
   directory=Path(t['directory']);directory.relative_to(Path(manifest).resolve().parent/'tasks')
   if fresh and {p.name for p in directory.iterdir()}!={'core.xyz','xcontrol'}:raise InvalidArtifact('nonfresh input')
 if new!=m['tasks']:raise InvalidArtifact('executable subset differs')
 return {'status':'validated','manifest':record(manifest),'cases':225,'physical_sources':208,'unavailable_preparations':17,'logical_cells':832,'actual_reuses':4,'new_calls':828,'new_MACE_DFT':0}

def execute(manifest):
 validate(manifest,True);m=read_json(manifest);root=Path(manifest).parent
 if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE',0))!=64:raise InvalidArtifact('64CPU allocation required')
 with (root/'execute.lock').open('a') as f:
  fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);begin=time.monotonic();receipts=[];errors=[]
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
   fs={ex.submit(backend.execute_cell,manifest,m,t):t for t in m['tasks']}
   for future in concurrent.futures.as_completed(fs):
    try:receipts.append(future.result())
    except Exception as e:errors.append({'task_id':fs[future]['task_id'],'error':str(e)})
  elapsed=time.monotonic()-begin;write_new(root/'EXECUTION.json',{'manifest':record(manifest),'receipts':receipts,'errors':errors,'wall_seconds':elapsed,'allocated_core_seconds':elapsed*64,'job_id':os.environ['SLURM_JOB_ID'],'GPU_seconds':0})
 return {'receipts':len(receipts),'errors':len(errors)}

def collect(manifest,output):
 validate(manifest);m=read_json(manifest);rows=[]
 for t in m['cells']:
  row=backend.parse_task(manifest,t) if not t['reuse'] else {**t['reuse']['row'],'task_id':t['task_id'],'cell_id':t['cell_id']}
  rows.append({**row,'reused':bool(t['reuse']),'reuse':t['reuse']})
 result={'protocol_id':PROTOCOL,'manifest':record(manifest),'rows':rows,'denominator':832,'complete':sum(r['status']=='complete' for r in rows),'new_attempts':828,'reused':4,'implementation':record(__file__)}
 write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}

def compare(collection,output):
 r=read_json(collection);m=read_json(verify(r['manifest']));validate(verify(r['manifest']));ref=reference_check(verify(m['reference']));base=read_json(verify(m['ledger']));df=read_json(verify(m['DFT']))
 idx={(t['case_id'],t['metal'],t['medium']):t for t in r['rows']};didx={x['case_id']:x for x in df['rows']}
 bands={'released':base['bands']['released'],'native_adaptive':base['bands']['minimal_recovered'],'DFT':df['bands']['DFT'],'standalone_static':ref['variants']['static']['bands']}
 rows=[]
 for c in m['cases']:
  old=c['prior_row'];cid=c['case_id'];value=None;matrix={z:{} for z in ('Ca','La')};reasons=[]
  if c['status']=='prepared':
   for z in matrix:
    pair={s:idx[cid,z,s] for s in ('vacuum','alpb')};good=all(x['status']=='complete' for x in pair.values())
    matrix[z]['origin']={'status':'complete' if good else 'unavailable','components':{'MACE_eV':c['origins'][z]['MACE_eV'],
      'GFN2_vacuum_hartree':pair['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':pair['alpb']['energy_hartree']} if good else None}
    reasons.extend({'metal':z,'medium':s,'reason':x['reason']} for s,x in pair.items() if x['status']!='complete')
   pool=choose_rows(matrix,('origin',));value=pool['operational']['composite_R_model_kcal_mol'] if pool['operational'] else None
  else:pool={'status':'unavailable','reason':c['reason'],'operational':None,'mathematical':None};reasons.append(c['reason'])
  methods={'released':old['methods']['released'],'native_adaptive':old['methods']['minimal_recovered'],'DFT':didx[cid]['methods']['DFT'],
    'standalone_static':method_record(value,bands['standalone_static'],old['expected_class'])}
  rows.append({**{k:old[k] for k in ('case_id','biological_group','root_case_id','source_conditioning_metal','expected_class')},'methods':methods,'pool':pool,'matrix':matrix,'reasons':reasons,'preparation_status':c['status'],'historical_native_row':old,'historical_DFT_row':didx[cid]})
 ri={x['case_id']:x for x in rows};pools=[];triples=[]
 for group in base['pools']:
  if group['pool']=='balanced':
   methods={}
   for name,b in bands.items():
    values=[next(x for x in pools if x['root_case_id']==group['root_case_id'] and x['pool']==arm)['methods'][name]['R'] for arm in ('La4','Ca5')]
    value=sum(values)/2 if all(x is not None for x in values) else None
    methods[name]={**method_record(value,b,group['expected_class']),'Ca5_minus_La4':values[1]-values[0] if value is not None else None}
  else:methods={name:aggregate(group['members'],ri,name,b,group['expected_class']) for name,b in bands.items()}
  pools.append({**group,'methods':methods})
 for group in base['triples']:triples.append({**group,'methods':{name:aggregate(group['members'],ri,name,b,group['expected_class']) for name,b in bands.items()}})
 subsets={'all225':rows,'La100':[x for x in rows if x['source_conditioning_metal']=='La'],'Ca125':[x for x in rows if x['source_conditioning_metal']=='Ca'],'all100_triples':triples,**{arm:[x for x in pools if x['pool']==arm] for arm in ('La4','Ca5','balanced')}}
 counts={name:summarize_rows(values,bands) for name,values in subsets.items()};matched={}
 for name,values in subsets.items():
  matched[name]={}
  for other in ('released','native_adaptive','DFT'):
   common=[x for x in values if all(x['methods'][k]['decision']!='unavailable' for k in (other,'standalone_static'))]
   matched[name][other]={'common':len(common),'counts':summarize_rows(common,(other,'standalone_static')),'transitions':dict(Counter(x['methods'][other]['outcome']+'->'+x['methods']['standalone_static']['outcome'] for x in common))}
 result={'protocol_id':PROTOCOL,'collection':record(collection),'reference':m['reference'],'native_ledger':m['ledger'],'DFT_ledger':m['DFT'],'bands':bands,
  'rows':rows,'pools':pools,'triples':triples,'counts':counts,'matched':matched,'calibration_changed':False,'production_changed':False,'new_calls_in_comparison':0,'implementation':record(__file__)}
 write_new(output,result);return {'all225':counts['all225'],'matched':matched['all225']}

def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
 for op,names in {'prepare':('inventory','reference','pilot','ledger','dft','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'compare':('collection','output')}.items():
  a=sub.add_parser(op)
  for n in names:a.add_argument('--'+n,required=True,type=Path)
  if op=='validate':a.add_argument('--fresh',action='store_true')
 args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
if __name__=='__main__':main()
