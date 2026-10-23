"""Saved joins and four exact distance checks; no bound/native recomputation."""
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
ROOT=Path.cwd()
D=HERE.parent
W=ROOT/str(D.relative_to(ROOT)).replace('diagnostics/','workspaces/',1)
ACTUAL=W/'actual_v1'

def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):return json.loads(Path(p).read_text())

def main():
 start,cpu=time.monotonic(),time.process_time();checks=[]
 def ck(name,passed,evidence=None):
  checks.append({'name':name,'passed':bool(passed),'evidence':evidence})
  if not passed:raise AssertionError(name)
 def verify(name,p):ck(name,pin(p['path'])==p,p)
 rp=pin(ACTUAL/'RESULT.json');result=read(rp['path'])
 ck('exact_result',rp['sha256']=='7033c722b8fcfa6e8b51024f8d1b649460ecab2e724b2afe8af2e3cac7f72b4c')
 cfg=read(ACTUAL/'CONFIG.json');template=read(D/'CONFIG_TEMPLATE.json')
 release=read(ACTUAL/'ROOT_RELEASE.json');terminal=read(ACTUAL/'PROCESS_TERMINAL.json')
 for name,p in cfg['authorities'].items():verify('authority_'+name,p)
 for name in ('configuration','implementation','template','source_review'):verify('release_'+name,release[name])
 verify('result_configuration',result['configuration']);verify('terminal_result',terminal['result'])
 ck('actual_config_exact_reviewed_template_plus_receipt',{**template,'source_review':cfg['source_review']}==cfg)
 ck('implementation_join',release['implementation']==cfg['implementation'])
 review=read(cfg['source_review']['path'])
 ck('source_review_join',review['status']=='passed_real_bound_source_review' and review['implementation']==cfg['implementation'] and review['configuration_template']==release['template'])
 ck('root_release',release['status']=='released_fixed_baseline_real_bound_diagnostic' and release['configuration']==result['configuration'] and release['source_review']==cfg['source_review'])
 ck('terminal_success',terminal['status']=='complete' and terminal['exit_code']==0 and terminal['error'] is None)
 ck('stdout_exact_result',read(ACTUAL/'stdout.txt')==result)
 ck('stderr_empty',not (ACTUAL/'stderr.txt').read_text())
 ck('correct_complete_status',result['status']=='complete_exploratory_real_bounds')
 ck('scope_and_no_native',result['baseline_only_H_zero'] is True and result['rigorous_certificate'] is False and result['complete_PME_stability_admitted'] is False and result['new_native_calls']==terminal['new_native_calls']==release['new_native_calls']==0)
 ck('one_thread_login_no_scheduler',result['threads']==release['threads']==1 and terminal['allocated_scheduler_cpu_seconds']==0 and release['execution_mode']=='authorized_login_no_scheduler_allocation')
 ck('cost_ordering',0<result['cpu_seconds']<=terminal['child_cpu_seconds'] and 0<result['wall_seconds']<=terminal['wall_seconds'])
 ck('exact_two_fixture_denominator',[f['fixture'] for f in result['fixtures']]==cfg['fixtures'] and len(result['fixtures'])==2)
 geometry=[]
 for f in result['fixtures']:
  metal=f['fixture']['name'];a=f['atom_blocks'];b=f['water_ion_blocks']
  ck('counts_'+metal,f['pairs']==1815015 and b['blocks']==3998 and b['cross_blocks']==262878)
  ck('site_row_arithmetic_'+metal,abs((a['lower_diagonal']-a['lower_row_norm_sum'])-a['lower'])<1e-15 and abs((a['lower_diagonal']+a['lower_row_norm_sum'])-a['upper'])<1e-15)
  ck('water_row_arithmetic_'+metal,abs((b['lower_diagonal_min']-b['lower_offblock_norm_sum'])-b['lower'])<1e-15)
  ck('negative_not_certificate_'+metal,a['lower']<0 and b['lower']<0 and a['positive_real_lower_bound'] is False and b['positive_real_lower_bound'] is False and b['rigorous_rounding_certificate'] is False and b['complete_PME_stability_admitted'] is False)
  verify('geometry_'+metal,f['fixture']['xyz'])
  lines=Path(f['fixture']['xyz']['path']).read_text().splitlines()
  rows={int(t[0]):t for line in lines[2:] if (t:=line.split())}
  mids=[i for i,t in rows.items() if int(t[5])==f['fixture']['metal_type']]
  ck('one_metal_'+metal,len(mids)==1)
  mid=mids[0];m=[Fraction(x) for x in rows[mid][2:5]];L=Fraction('49.324')
  ids=[a['lower_atom_one_based']]+b['lower_atoms_one_based'];out=[]
  ck('worst_role_indices_'+metal,ids==[11956,6578,6579,6580] and [int(rows[i][5]) for i in ids]==[414,402,403,403])
  for i in ids:
   t=rows[i];dr=[Fraction(x)-y for x,y in zip(t[2:5],m)]
   shift=[int((x/L+Fraction(1,2))//1) for x in dr]
   nearest=[x-n*L for x,n in zip(dr,shift)]
   squared=sum(x*x for x in nearest);literal=sum(x*x for x in dr)
   ck('outside_real_cutoff_'+metal+'_'+str(i),squared>Fraction(81))
   with localcontext() as ctx:
    ctx.prec=60
    d2=Decimal(squared.numerator)/Decimal(squared.denominator)
    l2=Decimal(literal.numerator)/Decimal(literal.denominator)
    distance=d2.sqrt()
    out.append({'atom_one_based':i,'type':int(t[5]),'coordinate_tokens_A':t[2:5],
      'nearest_lattice_shift':shift,'literal_distance_squared_A2':str(l2),'minimum_image_distance_squared_A2':str(d2),
      'minimum_image_distance_A':str(distance),'display_distance_enclosure_A':{'lower':str(distance.next_minus()),'upper':str(distance.next_plus())},
      'exact_squared_distance_greater_than_81':True})
  geometry.append({'metal':metal,'metal_atom_one_based':mid,'metal_coordinate_tokens_A':rows[mid][2:5],'geometry':f['fixture']['xyz'],'worst_sites':out})
 ck('same_saved_extrema_both_metals',result['fixtures'][0]['atom_blocks']==result['fixtures'][1]['atom_blocks'] and result['fixtures'][0]['water_ion_blocks']==result['fixtures'][1]['water_ion_blocks'])
 ck('same_worst_geometry_both_metals',geometry[0]['worst_sites']==geometry[1]['worst_sites'])
 conclusion={
  'site_lower_bound':result['fixtures'][0]['atom_blocks']['lower'],
  'water_block_lower_bound':result['fixtures'][0]['water_ion_blocks']['lower'],
  'insufficient_real_space_certificate':True,'operator_indefiniteness_established':False,
  'full_PME_stability_admitted':False,
  'fixed_nonmetal_rows_unchanged_by_metal_only_diagonal_H':True,
  'scope':'At these fixed geometries, split, intrinsic-alpha normalization and fixed partitions, metal-only H changes only the metal diagonal block. The displayed nonmetal rows are unchanged. All four sites also lie outside the9A real cutoff, so their real cross blocks to the metal are zero. The unchanged negative rows prevent these unweighted Gershgorin bounds from becoming positive solely through metal stiffness. This is a limitation of this certificate, not of the complete operator or model stability.'}
 report={'status':'passed_real_bound_saved_result_review','result':rp,'configuration':pin(ACTUAL/'CONFIG.json'),'implementation':cfg['implementation'],
  'reviewer_implementation':pin(Path(__file__)),'release':pin(ACTUAL/'ROOT_RELEASE.json'),'source_review':cfg['source_review'],
  'process_terminal':pin(ACTUAL/'PROCESS_TERMINAL.json'),'checks':checks,'checks_passed':len(checks),'geometry_check':geometry,'conclusion':conclusion,
  'cost':{'process_wall_seconds':terminal['wall_seconds'],'child_cpu_seconds':terminal['child_cpu_seconds'],'in_program_wall_seconds':result['wall_seconds'],'in_program_cpu_seconds':result['cpu_seconds'],'allocated_scheduler_cpu_seconds':0},
  'review_wall_seconds':time.monotonic()-start,'review_cpu_seconds':time.process_time()-cpu,'threads':1,
  'new_native_calls':0,'whole_bound_recomputations':0,'fit_outputs_read':0,'candidate_or_geometry_changes':0}
 out=W/'independent_actual_review_v1';out.mkdir(exist_ok=True)
 with (out/'REVIEW.json').open('x') as f:json.dump(report,f,indent=2,sort_keys=True);f.write('\n')
 txt='Saved real-bound audit passed; neither fixed partition certifies a positive lower bound.\n\n'
 txt+='Both metals: atom bound -1.735670190813629 at Cl11956; water/ion block bound -1.1305549478256092 at water6578-6580.\n'
 txt+='Metal distances: Cl35.79974442A; water O21.77727994A, H20.86013471A and21.76946703A. Exact squared distances exceed81A²; no metal cross block exists at9A.\n\n'
 txt+=conclusion['scope']+'\n\n'
 txt+='No instability or full-PME claim. Exact inputs, template-only review-pin insertion, terminal/stdout and declared one-thread cost join; no whole-bound replay or native call.\n'
 txt+=f"Actual process cost {terminal['wall_seconds']:.6f}s wall/{terminal['child_cpu_seconds']:.6f}s child CPU; no scheduler allocation.\n"
 with (out/'REPORT.txt').open('x') as f:f.write(txt)
 print(json.dumps({'status':report['status'],'review':pin(out/'REVIEW.json'),'checks_passed':len(checks),'review_cpu_seconds':report['review_cpu_seconds']}))

if __name__=='__main__':main()
