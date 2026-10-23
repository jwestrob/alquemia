"""Manufactured term slope closure and five-parent allocation controls only."""
import sys,json,time,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import COLLECT as c
from decimal import Decimal

def main():
 tick,cpu=time.monotonic(),time.process_time();passed=[];fail=[]
 def check(name,fn):
  try:c.require(fn() is not False,name);passed.append(name)
  except Exception as exc:fail.append(dict(name=name,error=repr(exc)))
 def rejected(fn):
  try:fn()
  except (ValueError,KeyError):return True
  return False
 inp=c.read(c.pin(c.HERE/'INPUTS.json'));h=c.load(inp['numerical_helper'],'manufactured_term_numerics')
 d=dict(kind='strain',L0_A=1.,strain_matrix=[[0.,0.,0.],[0.,1.,0.],[0.,0.,0.]])
 parsed={};traces={}
 for sid in c.IDS:
  q=Decimal(0) if sid=='central' else Decimal(sid.split('__')[1])*(-1 if sid.endswith('minus') else 1)
  rows=[dict(label='INITIAL',enabled=True,energy_token='0.00000000000000000E+0',virial_tokens=['0.00000000000000000E+0']*9)]
  for i in range(1,29):
   # Only two nonzero terms; slopes 2 and -1, total slope 1.
   a=2 if i==1 else (-1 if i==2 else 0);cum=2 if i==1 else 1
   rows.append(dict(label='term'+str(i),enabled=i<=2,energy_token=format(Decimal(a)*q,'.17E'),virial_tokens=[format(Decimal(cum if j==4 else 0),'.17E') for j in range(9)]))
  rows.append(dict(label='FINAL',enabled=True,energy_token=format(q,'.17E'),virial_tokens=list(rows[-1]['virial_tokens'])))
  parsed[sid]=dict(energy_token=format(q,'.17E'),virial_tokens=rows[-1]['virial_tokens']);traces[sid]=dict(rows=rows)
 x=c.decomposition(h,d,parsed,traces)
 check('two signed terms cancel to total derivative',lambda:x['total_original_checks']['passed'] and x['terms']['term1']['analytic_kcal_mol_A']==2 and x['terms']['term2']['analytic_kcal_mol_A']==-1 and x['central_virial_projection_closure'])
 check('all28 terms retained including disabled zeros',lambda:len(x['terms'])==28 and sum(t['enabled'] for t in x['terms'].values())==2)
 check('both stencils and closure exact',lambda:len(x['FD_component_closure'])==2 and all(y['error']==0 for y in x['FD_component_closure']))
 corrupt=copy.deepcopy(traces);corrupt[c.IDS[-1]]['rows'][1]['energy_token']='99'
 check('omitted/misassigned energy fails closure',lambda:rejected(lambda:c.decomposition(h,d,parsed,corrupt)))
 corrupt=copy.deepcopy(traces);corrupt['central']['rows'][-2]['virial_tokens'][4]='99'
 check('misassigned virial fails closure',lambda:rejected(lambda:c.decomposition(h,d,parsed,corrupt)))
 check('missing endpoint fails',lambda:rejected(lambda:c.decomposition(h,d,{k:v for k,v in parsed.items() if k!=c.IDS[-1]},traces)))
 subs=[dict(job_id=str(700+i)) for i in range(5)];lines=[]
 for sub in subs:
  for suffix in ('','.batch'):lines.append(sub['job_id']+suffix+'|test|COMPLETED|0:0|10|1|00:00:09')
 text='\n'.join(lines);account=c.scheduler(text,subs)
 check('five parents and no doublecount of steps',lambda:account['allocation_count']==5 and account['allocated_CPU_seconds']==50)
 check('four allocations rejected',lambda:rejected(lambda:c.scheduler(text,subs[:4])))
 check('missing parent rejected',lambda:rejected(lambda:c.scheduler('\n'.join(lines[1:]),subs)))
 check('failed parent rejected',lambda:rejected(lambda:c.scheduler(text.replace('COMPLETED','FAILED',1),subs)))
 check('duplicate parent rejected',lambda:rejected(lambda:c.scheduler(text+'\n'+lines[0],subs)))
 check('wrong CPU count rejected',lambda:rejected(lambda:c.scheduler(text.replace('|10|1|','|10|2|',1),subs)))
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
 result=dict(status='passed_term_collector_controls' if not fail else 'failed_term_collector_controls',implementation=c.pin(__file__),collector=c.pin(c.HERE/'COLLECT.py'),inputs=c.pin(c.HERE/'INPUTS.json'),passed=passed,failures=fail,actual_collection=False,native_calls=0,FFT_calls=0,factor_calls=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 c.save(out/'RESULT.json',result);print(json.dumps(result));return bool(fail)
if __name__=='__main__':raise SystemExit(main())
