"""Manufactured trace corruption + insertion/lineage controls; no native/FFT."""
from BINDING import *
from decimal import Decimal
import ast,copy,tempfile

def main():
 tick=time.monotonic();cpu=time.process_time();passed=[];failed=[]
 def check(name,fn):
  try:require(fn() is not False,name);passed.append(name)
  except Exception as e:failed.append(dict(name=name,error=repr(e)))
 def rejects(fn):
  try:fn()
  except (ValueError,KeyError):return True
  return False
 authority=read(pin(HERE/'INPUTS.json'));observer=load(pin(HERE/'OBSERVER.py'),'term_controls_observer')
 old=verify(authority['original_gradient']).read_text();new=observer.render(old)
 check('exact marked inverse',lambda:observer.inverse(new)==old)
 check('28 terms and all fixed form insertion lines',lambda:len(observer.TERMS)==28 and all(len(l)<=72 for block in new.split(observer.BEGIN)[1:] for l in block.split(observer.END)[0].splitlines()))
 check('all original native calls retain order',lambda:__import__('re').findall(r'call (\w+)',observer.inverse(new))==__import__('re').findall(r'call (\w+)',old))
 check('unique anchor rejects source drift',lambda:rejects(lambda:observer.render(old.replace('      einter = 0.0d0','      einter=0d0'))))
 check('new runtime sole addition',lambda:(HERE/'RUN.py').read_text().replace("            observer=load(m['observer'],'term_saved_parser')\n            trace=observer.parse(native/'term_trace.dat')\n            old=parser.parse_gradient(verify(row['baseline_gradient']),n=11958)\n            result['observer_checks']=observer.checks(trace,parsed,old,parser)\n            result['term_output']=pin(native/'term_trace.dat')\n",'')==verify(authority['original_native_runner']).read_text())
 helper=load(read(read(authority['declaration'])['baseline_manifest'])['collector'],'term_inherited_helpers')
 labels=['INITIAL']+[t[0] for t in observer.TERMS]+['FINAL']
 def token(x):return format(Decimal(x),'.17E')
 rows=[]
 for i,label in enumerate(labels):
  v=0 if i==0 else min(i,28);e=406 if i==29 else i
  rows.append(' '.join([str(i),label,'1',token(e),*[token(v if j==4 else 0) for j in range(9)]]))
 lines=['TERM_TRACE_V1','UNITS energy=kcal/mol virial=kcal/mol','VDW F BUFFERED-14-7',*rows,'FINAL_TERMS '+' '.join(token(i) for i in range(1,29)),'END snapshots=30 gradient_calls=1']
 with tempfile.TemporaryDirectory(prefix='term_observer_controls_') as temp:
  path=Path(temp)/'trace';path.write_text('\n'.join(lines)+'\n');trace=observer.parse(path)
  g=dict(E_G_virial_token_sha256='synthetic_equal',energy_token=trace['rows'][-1]['energy_token'],virial_tokens=trace['rows'][-1]['virial_tokens'],polarization_token=trace['rows'][19]['energy_token'],electric_token='1',dielec_token='1',role='BYPASS_OLD_ABSENT',state_id='central')
  check('manufactured full trace and closure',lambda:observer.checks(trace,g,g,helper)['virial_closure'])
  for name,index,replacement in [('missing termination',34,''),('wrong order',4,rows[2]),('wrong header',0,'BAD'),('wrong final energy',33,'FINAL_TERMS '+' '.join(token(0) for i in range(28))),('nonfinite',4,rows[1].replace(token(1),'NaN',1))]:
   def trial(index=index,replacement=replacement):
    corrupt=list(lines);corrupt[index]=replacement;path.write_text('\n'.join(corrupt)+'\n');return rejects(lambda:observer.parse(path))
   check(name,trial)
  for name,change in [('total token drift',dict(E_G_virial_token_sha256='different')),('polar token drift',dict(polarization_token='-99')),('state drift',dict(state_id='other'))]:
   check(name,lambda change=change:rejects(lambda:observer.checks(trace,g,dict(g,**change),helper)))
  changed=copy.deepcopy(trace);changed['rows'][2]['energy_token']='1000'
  check('energy closure failure',lambda:rejects(lambda:observer.checks(changed,g,g,helper)))
  changed=copy.deepcopy(trace);changed['rows'][-2]['virial_tokens'][4]='27'
  check('unassigned final virial rejected',lambda:rejects(lambda:observer.checks(changed,g,g,helper)))
  changed=copy.deepcopy(trace);changed['vdw_long_range_correction']=True
  check('VDW setting drift rejected',lambda:rejects(lambda:observer.checks(changed,g,g,helper)))
 for name in ('OBSERVER.py','BINDING.py','BUILD.py','RUN.py','PREPARE.py'):
  check('syntax '+name,lambda name=name:ast.parse((HERE/name).read_text()) is not None)
 result=dict(status='passed_term_observer_source_controls' if not failed else 'failed_term_observer_source_controls',implementation=pin(__file__),sources={n:pin(HERE/n) for n in ('OBSERVER.py','BINDING.py','BUILD.py','RUN.py','PREPARE.py')},passed=passed,failures=failed,native_calls=0,builds=0,scientific_replays=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);save(out/'RESULT.json',result);print(json.dumps(result));return bool(failed)
if __name__=='__main__':raise SystemExit(main())
