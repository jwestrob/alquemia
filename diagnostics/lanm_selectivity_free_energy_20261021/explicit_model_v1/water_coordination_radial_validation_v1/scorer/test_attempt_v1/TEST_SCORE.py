#!/usr/bin/env python3
"""Synthetic scorer controls using geometry only; no validation E/F reads."""
import copy,importlib.util,json,math,sys,time,hashlib
from pathlib import Path
H=Path(__file__).resolve().parent;ROOT=H.parents[4]
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def pin(p):return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def main(out):
 t=time.monotonic();cpu=time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
 score=load(H/'SCORE.py','fresh_score');bodypath=H.parents[1]/'matched_dry_v1/OBJECTIVE.py';body=load(bodypath,'fresh_body')
 mp=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/water_coordination_radial_validation_v1/preparation_v1/MANIFEST.json';manifest=json.loads(mp.read_text())
 data={'rows':copy.deepcopy(manifest['rows'])};model={};anchors={'model':{},'reference':{}};checks=[]
 def ok(v,label):
  if not v:raise AssertionError(label)
  checks.append(label)
 def bad(fn,label):
  try:fn()
  except (ValueError,KeyError):checks.append(label)
  else:raise AssertionError(label)
 for m in ('La','Dy'):
  x={'energy':10. if m=='La' else 20.,'forces':[[0.,0.,0.] for _ in range(28)]}
  anchors['model'][m]=copy.deepcopy(x);anchors['reference'][m]=copy.deepcopy(x)
 for gid,row in data['rows'].items():
  row['bodies']=body.bodies(row);m=gid.rsplit('_',1)[1];value=copy.deepcopy(anchors['model'][m]);value['energy']+=-1 if '_minus_' in gid else 1
  row['reference']=copy.deepcopy(value);model[gid]=copy.deepcopy(value)
 r=score.evaluate(data,model,body,anchors);ok(r['available_checks']==12 and r['passed_checks']==12 and r['all_validation_checks_passed'],'Exact synthetic matching passes twelve')
 for gid in score.CELLS:
  m=copy.deepcopy(model);m[gid]=None;r=score.evaluate(data,m,body,anchors);ok(r['available_checks']==8 and not r['all_validation_checks_passed'],'Single missing cell retains denominator '+gid)
  m=copy.deepcopy(model);m[gid]['energy']+=1;r=score.evaluate(data,m,body,anchors);ok(r['all_validation_checks_passed'],'Energy threshold inclusive '+gid)
  m[gid]['energy']+=1e-6;r=score.evaluate(data,m,body,anchors);ok(not r['individual_radial_work'][gid]['passed'] and not r['all_validation_checks_passed'],'Energy threshold fails '+gid)
  m=copy.deepcopy(model);m[gid]['energy']=float('nan');bad(lambda:score.evaluate(data,m,body,anchors),'Reject nonfinite available '+gid)
  m=copy.deepcopy(model);m[gid]['forces']=m[gid]['forces'][:-1];bad(lambda:score.evaluate(data,m,body,anchors),'Reject wrong force shape '+gid)
 for sign in ('minus','plus'):
  m=copy.deepcopy(model);m['aquo9_'+sign+'_La']['energy']+=.75;m['aquo9_'+sign+'_Dy']['energy']-=.75;r=score.evaluate(data,m,body,anchors)
  ok(all(r['individual_radial_work']['aquo9_'+sign+'_'+x]['passed'] for x in ('La','Dy')) and not r['paired_radial_work'][sign]['passed'],'Paired work cannot hide two individually small errors '+sign)
  m=copy.deepcopy(model)
  for x in ('La','Dy'):m['aquo9_'+sign+'_'+x]['forces'][0][0]=10.
  r=score.evaluate(data,m,body,anchors);ok(r['rigid'][sign]['paired']['passed'] and not r['rigid'][sign]['La']['passed'] and not r['rigid'][sign]['Dy']['passed'] and not r['all_validation_checks_passed'],'Paired cancellation cannot rescue individuals '+sign)
  m=copy.deepcopy(model);m['aquo9_'+sign+'_La']['forces'][0][0]=5.;r=score.evaluate(data,m,body,anchors);ok(r['rigid'][sign]['La']['passed'],'Maximum component boundary5 inclusive '+sign)
  m['aquo9_'+sign+'_La']['forces'][0][0]+=1e-6;r=score.evaluate(data,m,body,anchors);ok(not r['rigid'][sign]['La']['passed'],'Maximum component beyond5 fails '+sign)
 m={g:None for g in score.CELLS};r=score.evaluate(data,m,body,anchors);ok(r['available_checks']==0 and r['declared_checks']==12 and not r['all_validation_checks_passed'],'All missing is unavailable')
 a=copy.deepcopy(anchors);a['model']['La']=None;r=score.evaluate(data,model,body,a);ok(r['available_checks']==8 and not r['all_validation_checks_passed'],'Missing anchor preserves force checks only')
 d=copy.deepcopy(data)
 for row in d['rows'].values():row.pop('reference')
 r=score.evaluate(d,model,body,anchors);ok(r['available_checks']==0 and not r['all_validation_checks_passed'],'References not yet present remain unavailable')
 m=copy.deepcopy(model);m['extra']=None;bad(lambda:score.evaluate(data,m,body,anchors),'Reject added slot')
 result={'status':'passed','source':pin(H/'SCORE.py'),'test_source':pin(Path(__file__)),'geometry_manifest':pin(mp),'body_helper':pin(bodypath),'checks':checks,'checks_passed':len(checks),'synthetic_values_only':True,'validation_outcomes_opened':False,'fit_outcomes_opened':False,'native_calls':0,'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1}
 (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','checks':len(checks),'result':pin(out/'RESULT.json')}))
if __name__=='__main__':main(sys.argv[1])
