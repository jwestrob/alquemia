"""Bounded manufactured anchor/scalar/serializer controls; zero native calls."""
from BINDING import *
from NUMERICS import *
from decimal import Decimal
import ast,copy,tempfile

def main():
 tick,cpu=time.monotonic(),time.process_time();checks_done=[];fail=[]
 def check(name,fn):
  try:require(fn() is not False,name);checks_done.append(name)
  except Exception as exc:fail.append(dict(name=name,error=repr(exc)))
 authority=read(pin(HERE/'INPUTS.json'));decl=read(authority['declaration']);base=read(read(decl['parent_build'])['manifest']);h=load(base['collector'],'refine_control_numerics')
 center=dict(gradient_tokens=['1.00000000000000000E+0','1.00000000000000000E+0','1.00000000000000000E+0','1.00000000000000000E+0','2.00000000000000000E+0','3.00000000000000000E+0'],virial_tokens=['1.00000000000000000E+1','0.00000000000000000E+0','0.00000000000000000E+0','0.00000000000000000E+0','2.00000000000000000E+1','0.00000000000000000E+0','0.00000000000000000E+0','0.00000000000000000E+0','3.00000000000000000E+1'])
 for axis,expected in enumerate((8,14,18)):
  S=[[0.,0.,0.] for i in range(3)];S[axis][axis]=1.;d=dict(kind='strain',strain_matrix=S,anchor_A=[1.,2.,3.],L0_A=1.)
  p=corrected_projection(h,center,d)
  check('axis'+str(axis)+' signed chain rule',lambda p=p,expected=expected:p['corrected']==expected and p['adjustment']<0)
  check('axis'+str(axis)+' literal derivative print allowance',lambda p=p:p['corrected_print']==p['raw_print']+p['adjustment_print'] and p['adjustment_print']>0)
  energies={q:(format(-Decimal(expected)*q,'.17E'),format(Decimal(expected)*q,'.17E')) for q in STEPS}
  corrected=checks(h,'manufactured',p['corrected'],p['corrected_print'],energies);raw=checks(h,'manufactured',p['raw'],p['raw_print'],energies)
  check('axis'+str(axis)+' raw failed corrected passed',lambda corrected=corrected,raw=raw:corrected['passed'] and not raw['passed'])
  check('axis'+str(axis)+' inherited two-step allowance',lambda corrected=corrected:corrected['two_step']['bound']==float(h.GATE+sum(Decimal.from_float(s['printing_bound']) for s in corrected['steps'])) or abs(corrected['two_step']['bound']-(.001+sum(s['printing_bound'] for s in corrected['steps'])))<1e-18)
 zero=corrected_projection(h,center,dict(kind='strain',strain_matrix=[[1.,0.,0.],[0.,0.,0.],[0.,0.,0.]],anchor_A=[0.,0.,0.],L0_A=1.))
 check('origin anchor preserves raw projection',lambda:zero['raw']==zero['corrected'] and zero['raw_print']==zero['corrected_print'])
 check('old helper global steps not mutated',lambda:h.STEPS==(Decimal('.0001'),Decimal('.00005')))
 check('exact12 new slots and no center',lambda:len(IDS)==len(set(IDS))==12 and 'central' not in IDS)
 check('unchanged native runner bytes',lambda:(HERE/'RUN.py').read_bytes()==verify(authority['baseline_runner']).read_bytes())
 for n in ('BINDING.py','RUN.py','PREPARE.py','NUMERICS.py','COLLECT.py'):
  check('syntax '+n,lambda n=n:ast.parse((HERE/n).read_text()) is not None)
 import numpy as np
 geom=load(authority['geometry_preparer'],'refine_control_geometry')
 with tempfile.TemporaryDirectory(prefix='strain_refinement_') as temp:
  p=Path(temp)/'in.xyz';p.write_text('2 sample\n49.324 49.324 49.324 90 90 90\n1 La 2.5405 1.8185 -1.86 900\n2 Na 8.12345 4.98765 3.14159 405\n')
  lines,rows,x,cell,ids=geom.parse_xyz(p)
  for axis in range(3):
   v=np.zeros_like(x);v[:,axis]=x[:,axis]-x[0,axis]
   for label in STEP_LABELS:
    for sign in (-1,1):
     q=sign*float(label);y=x+q*v;newcell=cell.copy();newcell[axis]*=1+q;dst=Path(temp)/(str(axis)+'_'+label+'_'+str(sign)+'.xyz');geom.write_xyz(dst,lines,rows,x,y,newcell,cell);_,_,z,c,keys=geom.parse_xyz(dst)
     check('serializer '+dst.stem,lambda z=z,c=c,y=y,newcell=newcell,keys=keys:np.array_equal(z,y) and np.array_equal(c,newcell) and keys==ids and np.array_equal(z[0],x[0]))
 result=dict(status='passed_anchored_strain_refinement_controls' if not fail else 'failed_anchored_strain_refinement_controls',implementation=pin(__file__),sources={n:pin(HERE/n) for n in ('BINDING.py','RUN.py','PREPARE.py','NUMERICS.py','COLLECT.py')},passed=checks_done,failures=fail,native_calls=0,builds=0,FFT_calls=0,factor_calls=0,actual_endpoint_evaluations=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);save(out/'RESULT.json',result);print(json.dumps(result));return bool(fail)
if __name__=='__main__':raise SystemExit(main())
