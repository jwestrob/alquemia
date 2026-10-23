"""Source/geometry preparation only, exact existing deformation and serializer."""
from BINDING import *
import numpy as np

def main():
 tick,cpu=time.monotonic(),time.process_time();a=read(pin(HERE/'INPUTS.json'));d=read(a['declaration']);build=read(d['parent_build']);base=read(build['manifest']);cm=read(d['central_manifest'])
 root=Path(build['manifest']['path']).parents[2]/'anchored_strain_refinement_v1';out=root/'preparation_v1';out.mkdir(parents=True,exist_ok=False)
 result=dict(status='failed_anchored_strain_refinement_preparation',implementation=pin(__file__),checks=[],error=None,native_calls=0,builds=0,FFT_calls=0,factorizations=0)
 def check(name,ok):result['checks'].append(dict(name=name,passed=bool(ok)));require(ok,name)
 try:
  check('One thread',all(os.environ.get(k)=='1' for k in THREADS));geom=load(a['geometry_preparer'],'strain_original_geometry')
  central=base['rows']['central']['inputs'];lines,rows,xyz,cell,ids=geom.parse_xyz(verify(central['xyz']));directions=read(d['directions'])
  check('Exact center and literal11958 map',len(rows)==11958 and ids[0][2]=='900' and len(cell)==6)
  chosen=[x for x in directions['directions'] if x['id'] in DIRECTIONS];check('Three original ordered directions',[x['id'] for x in chosen]==list(DIRECTIONS))
  m={k:base[k] for k in ('helpers','collector','environment','native_runtime','loader_command','candidate','composition','native_sources','directions','displaced_manifest')}
  m.update(status='prepared_twelve_anchored_strain_refinement_inputs',producer=pin(__file__),source_inputs=pin(HERE/'INPUTS.json'),root_declaration=a['declaration'],baseline_manifest=build['manifest'],build=d['parent_build'],build_root_review=a['build_root_review'],central_gradient=d['central_gradient'],central_result=cm['rows']['central']['baseline_result'],geometry_preparer=a['geometry_preparer'],execution_sources={n:pin(HERE/n) for n in ('BINDING.py','RUN.py')},original_wrapper=a['wrapper'],maximum_native_programs=12,rows={},resources=base['resources'])
  for direction in chosen:
   axis=direction['axis'];v=np.load(verify(direction['coefficients_dr_dq']),allow_pickle=False)
   expected=np.zeros_like(xyz);expected[:,axis]=xyz[:,axis]-xyz[0,axis]
   check(direction['id']+' original tangent/anchor/L0',direction['kind']=='strain' and direction['L0_A']==1. and np.array_equal(np.array(direction['anchor_A']),xyz[0]) and np.array_equal(v,expected))
   for label in STEP_LABELS:
    for sign,word in ((-1,'minus'),(1,'plus')):
     step=float(label);q=sign*step;sid=direction['id']+'__'+label+'__'+word;folder=out/'inputs'/sid;folder.mkdir(parents=True)
     new=xyz+q*v;newcell=cell.copy();newcell[axis]*=1+q;fp=folder/'system.xyz';geom.write_xyz(fp,lines,rows,xyz,new,newcell,cell)
     _,_,got,gcell,gids=geom.parse_xyz(fp)
     check(sid+' exact roundtrip/identity',np.array_equal(got,new) and np.array_equal(gcell,newcell) and gids==ids)
     check(sid+' original anchor/domain',np.array_equal(got[0],xyz[0]) and min(gcell[:3])>=44 and np.all(gcell[3:]==90) and 12<min(gcell[:3])/2)
     inputs=dict(central,xyz=pin(fp));m['rows'][sid]=dict(role='BYPASS_OLD_ABSENT',direction=direction['id'],half_step_A=step,sign=sign,inputs=inputs,source_inputs=inputs,output=str(root/'run_v1/tasks'/sid/'execution'))
  check('12 exact ordered states',tuple(m['rows'])==IDS);validate_manifest(m);save(out/'MANIFEST.json',m)
  result.update(status='anchored_strain_refinement_preparation_passed',manifest=pin(out/'MANIFEST.json'),source_inputs=pin(HERE/'INPUTS.json'))
 except Exception as exc:result['error']=repr(exc)
 result['cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1);save(out/'RESULT.json',result);print(json.dumps(result));return result['error'] is not None
if __name__=='__main__':raise SystemExit(main())
