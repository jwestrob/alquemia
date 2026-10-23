"""One-thread synthetic/source review only; never evaluates the real fixtures."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=Path.cwd()
W=ROOT/str(HERE.relative_to(ROOT)).replace('diagnostics/','workspaces/',1)
EXPECTED_SOURCE='2f6c5474dcc82308a65b8a4ae96d2d401630fcd783eb64cd376d4893f1386929'
EXPECTED_CFG='6ebbfc91'

def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
 start,cpu=time.monotonic(),time.process_time()
 checks=[]
 def ck(name,ok,evidence=None):
  checks.append({'name':name,'passed':bool(ok),'evidence':evidence})
  if not ok:raise AssertionError(name)
 src=pin(BASE/'BOUND.py');cfgp=pin(BASE/'CONFIG_TEMPLATE.json')
 ck('frozen_source',src['sha256']==EXPECTED_SOURCE,src)
 ck('frozen_configuration',cfgp['sha256'].startswith(EXPECTED_CFG),cfgp)
 cfg=json.loads(Path(cfgp['path']).read_text())
 for name,p in cfg['authorities'].items():ck('authority_'+name,pin(p['path'])==p,p)
 ck('template_implementation',cfg['implementation']==src)
 ck('baseline_fixed_engineering_recipe',cfg['ewald_inverse_A']==.4 and cfg['cutoff_A']==9. and cfg['custom_H']=='zero_baseline_only')
 ready=json.loads(Path(cfg['authorities']['handoff_ready']['path']).read_text())
 for name in ('inventory','inputs','namespace_and_representatives'):
  ck('handoff_'+name,pin(ready[name]['path'])==ready[name],ready[name])
 inv=json.loads(Path(ready['inventory']['path']).read_text())
 expected_params={'402':{'alpha':.837,'thole':.39},'403':{'alpha':.496,'thole':.39},
  '405':{'alpha':.120,'thole':.39},'414':{'alpha':4.,'thole':.39},
  '900':{'alpha':1.6597,'thole':.25},'901':{'alpha':.728,'thole':.139}}
 ck('exact_declared_type_parameters',cfg['type_parameters']==expected_params)
 bio=Path(cfg['authorities']['BIO']['path']).read_text().splitlines()
 records={int(t[1]):t for line in bio if (t:=line.split()) and t[0].lower()=='polarize'}
 for typ in (402,403,405,414):
  ck('BIO_alpha_thole_'+str(typ),[float(x) for x in records[typ][2:4]]==[expected_params[str(typ)]['alpha'],expected_params[str(typ)]['thole']])
 scales={t[0]:float(t[1]) for line in bio if (t:=line.split()) and t[0] in ('mutual-11-scale','mutual-12-scale','mutual-13-scale','mutual-14-scale')}
 ck('all_mutual_scales_one',len(scales)==4 and set(scales.values())=={1.})
 for c,f in zip(cfg['fixtures'],inv['fixtures']):
  ck('fixture_'+c['name'],c['xyz']==f['geometry'] and c['metal_type']==(900 if f['metal']=='La' else 901))
  ck('fixture_bytes_'+c['name'],pin(c['xyz']['path'])==c['xyz'])
  key=Path(f['key']['path']);ck('mapped_key_pin_'+c['name'],pin(key)==f['key'])
  polar={int(t[1]):t for line in key.read_text().splitlines() if (t:=line.split()) and t[0].lower()=='polarize'}
  t=c['metal_type'];ck('intrinsic_metal_'+c['name'],[float(x) for x in polar[t][2:4]]==[expected_params[str(t)]['alpha'],expected_params[str(t)]['thole']])
 ck('explicit_17676_polpair_not_old_key',cfg['metal_oxygen_polpair']=={'900':.17662704152555941,'901':.15457541153280463})
 s=importlib.util.spec_from_file_location('reviewed_bound',BASE/'BOUND.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
 # Exact damping cutoff at50, zero-gamma undamped path, including normalization.
 r=np.array([[1.,0.,0.],[1.,0.,0.]])
 tensor=m.tensor(r,np.ones(2),np.ones(2),np.array([50.,0.]),0.)
 ck('damp_50_and_gamma_zero_are_undamped',np.array_equal(tensor[0],np.diag([2.,-1.,-1.])) and np.array_equal(tensor[0],tensor[1]))
 # Both group order reversal and rectangular cross blocks occur in these fixtures.
 xyz=np.array([[0.,0.,0.],[1.7,.2,.1],[.3,2.1,.2],[.4,.1,2.3]])
 pairs=np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
 alpha=np.array([.6,.4,.8,.2]);ew=.4
 ts=m.tensor(xyz[pairs[:,1]]-xyz[pairs[:,0]],alpha[pairs[:,0]],alpha[pairs[:,1]],np.full(6,.39),ew)
 full=np.diag(np.repeat(1-alpha*4*ew**3/(3*np.sqrt(np.pi)),3))
 for (i,j),a in zip(pairs,ts):
  full[3*i:3*i+3,3*j:3*j+3]=-a;full[3*j:3*j+3,3*i:3*i+3]=-a.T
 eig=np.linalg.eigvalsh(full)
 for j,groups in enumerate(([[0,2],[1,3]],[[0,2,3],[1]])):
  b=m.bounds(4,pairs,ts,alpha,ew,groups)
  ck('interleaved_partition_enclosure_'+str(j),b['lower']<=eig[0]+1e-14 and b['upper']>=eig[-1]-1e-14,{'bounds':[b['lower'],b['upper']],'eigen_extrema':[float(eig[0]),float(eig[-1])]})
 # Independent scalar radial eigenvalue formula for normalized pair magnitudes.
 disp=np.array([[1.2,-.7,.9]])
 aa,bb,gamma=.6,.8,.39;r=float(np.linalg.norm(disp[0]));damp=gamma*r**3/np.sqrt(aa*bb)
 f3=1-np.exp(-damp);f5=1-(1+damp)*np.exp(-damp)
 expected=np.sort(np.sqrt(aa*bb)*np.array([-f3/r**3,-f3/r**3,(-f3+3*f5)/r**3]))
 got=np.linalg.eigvalsh(m.tensor(disp,np.array([aa]),np.array([bb]),np.array([gamma]),0.)[0])
 ck('normalized_damped_vacuum_eigenvalues',np.allclose(got,expected,rtol=0,atol=1e-15))
 author=W/'author_controls/RESULT.json';ac=json.loads(author.read_text())
 ck('author_tiny_controls',ac['status']=='passed_analytic_and_small_matrix_controls' and ac['controls']==4 and ac['new_native_calls']==0,pin(author))
 findings=[
  'dampewald recurrence gives bn1=(bn0+2a exp(-a²r²)/sqrt(pi))/r² and bn2=(3bn1+4a³ exp(-a²r²)/sqrt(pi))/r².',
  'dampthole initializes1 and changes only for nonzero gamma and damp<50. Native sixth-root pdamp algebra gives gamma*r³/sqrt(alpha_i*alpha_j); bitwise pdamp/cancellation parity is not claimed.',
  'Native mutual self field adds4a³/(3sqrt(pi))*u; B=I-sqrt(alpha)T sqrt(alpha) therefore subtracts alpha*self on its diagonal.',
  'Full Cartesian signs are assembled before each cross-block SVD norm. Opposite block entries use transposes; fixed water/ion grouping and Gershgorin bounds are coherent.',
  'Read_fixture enforces pinned bytes/cell/type counts, reciprocal two-H water bonds and unbonded ions; bounds separately rejects overlapping/incomplete blocks.',
  'This review binds the exact CONFIG_TEMPLATE values. Runtime authority hashes alone do not prove arbitrary changed scalar configurations are the reviewed recipe; only insertion of the review pin is allowed for the declared actual run.',
  'All results remain ordinary-double exploratory real+self bounds. A negative bound is not operator indefiniteness; a positive bound is not complete PME, native parity, nonzero-H, physical or sampling admission.'
 ]
 report={'status':'passed_real_bound_source_review','implementation':src,'configuration_template':cfgp,
  'reviewer_implementation':pin(Path(__file__)),'checks':checks,'checks_passed':len(checks),
  'author_controls':pin(author),'findings':findings,'native_calls':0,'actual_box_evaluations':0,'fit_outputs_read':0,
  'threads':1,'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,
  'complete_PME_stability_admitted':False,'rigorous_rounding_certificate':False}
 W.mkdir(parents=True,exist_ok=True)
 with (W/'REVIEW.json').open('x') as f:json.dump(report,f,indent=2,sort_keys=True);f.write('\n')
 with (W/'REPORT.txt').open('x') as f:f.write('Source review passed for the exact declared two-box H=0 exploratory real+self diagnostic.\n\n'+'\n'.join(findings)+'\n\nNo actual box, native program, fit artifact or candidate was evaluated.\n')
 print(json.dumps({'review':pin(W/'REVIEW.json'),'status':report['status'],'checks_passed':len(checks),'cpu_seconds':report['cpu_seconds']}))

if __name__=='__main__':main()
