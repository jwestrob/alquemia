#!/usr/bin/env python3
"""Ten saved Dy files, exact stage comparison; no corrected target or native call."""
import argparse,hashlib,importlib.util,json,os,re,time
from decimal import Decimal,getcontext
from pathlib import Path
getcontext().prec=50
D=Decimal
HERE=Path(__file__).resolve().parent

def pin(p):
 p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def checked(p):
 assert pin(p['path'])==p,'Changed source '+p['path'];return Path(p['path']).read_text()
def load(p):return json.loads(checked(p))
def mod(p,name):
 checked(p);s=importlib.util.spec_from_file_location(name,p['path']);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main(args):
 tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(args.config.read_text());assert not args.out.exists()
 assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'))
 assert cfg['implementation']==pin(__file__)
 diag=mod(cfg['diagnosis_source'],'frozen_dy_stage_helpers');native=mod(cfg['new_adapter'],'frozen_matching_grid_adapter');bound=diag.bound
 old=load(cfg['old_collection']);new=load(cfg['new_collection']);prior=load(cfg['old_diagnosis']);review=load(cfg['new_review']);checked(cfg['declaration'])
 assert prior['collection']==cfg['old_collection'] and prior['implementation']==cfg['diagnosis_source']
 assert review['collection']==cfg['new_collection'] and review['implementation']==cfg['new_adapter'] and review['review_passed'] is True and review['numerical_checks_passed'] is False
 assert review['failed_numerical_checks']==['Dy__odd_slope_h','Dy__odd_slope_half']
 assert old['complete']==10 and new['complete']==5 and new['declared']==5 and new['response']['declared_checks']==6
 assert old['response']['field_direction_lab']==new['response']['field_direction_lab'];direction=[D(str(x)) for x in old['response']['field_direction_lab']]
 checked(native.FROZEN_INPUTS['grid_parser']);checks={};rows={};calculations={};inputs={};excerpts={}
 def check(name,v):checks[name]=bool(v);assert v,name
 for label,col in [('original',old),('matching_grid',new)]:
  rows[label]={};inputs[label]={};excerpts[label]={}
  for tag in diag.TAGS:
   cid='pbe0__Dy__dry_mono_central__Ft_'+tag;q=col['rows'][cid]
   check(label+'/'+tag+'/native_admitted',q['status']=='complete' and q['returncode']==0 and all(q['parsed']['checks'].values()))
   checked(q['result']);checked(q['start']);text=checked(q['output']);inp=checked(q['input']);prop=checked(q['properties']['property_file']);inputs[label][tag]=inp
   check(label+'/'+tag+'/raw_property_tokens',all(t in prop for group in q['properties']['raw_tokens']['dipoles'].values() for t in group) and q['properties']['raw_tokens']['electronic'] in prop)
   # Exact existing saved-diagnosis last-iteration/final-stage extraction.
   lines=text.splitlines();end=[i for i,l in enumerate(lines) if 'SCF CONVERGED AFTER' in l];check(label+'/'+tag+'/one_SCF_terminal',len(end)==1);end=end[0]
   itrows=re.findall(r'^\s*(\d+)\s+(-?\d+\.\d{10,})\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*$', '\n'.join(lines[:end]),re.M)
   check(label+'/'+tag+'/SCF_rows_present',bool(itrows));last=itrows[-1];epre=D(last[1]);posttok=q['properties']['raw_tokens']['electronic'];epost=D(posttok)
   native_e=re.findall(r'^\s*Total Energy\s*:\s*([-\d.]+)\s+Eh',text,re.M)[-1]
   check(label+'/'+tag+'/final_electronic_join',abs(epost-D(native_e))<=bound(posttok)+bound(native_e))
   corr=re.findall(r'Exchange energy change after final integration\s*:\s*([-\d.]+)\s+Eh',text);check(label+'/'+tag+'/one_final_integration',len(corr)==1)
   delta=epost-epre;check(label+'/'+tag+'/stage_join',abs(delta-D(corr[0]))<=bound(last[1])+bound(posttok)+bound(corr[0]))
   mu=sum(t*D(x) for t,x in zip(direction,q['properties']['raw_tokens']['dipoles']['total']))
   d4=D(q['properties']['raw_tokens']['D4']);total=D(q['properties']['raw_tokens']['final'])
   check(label+'/'+tag+'/D4_closure',abs(total-epost-d4)<=D(str(q['properties']['energy_stage_closure_bound_hartree'])))
   grids=native.grid.grid_fingerprints(text);xc=native.dry.xc_metadata(text)
   check(label+'/'+tag+'/XC_metadata_replay',xc==q['parsed']['dft_xc_grid_metadata'])
   if label=='matching_grid':check(label+'/'+tag+'/grid_metadata_replay',grids==q['parsed']['physical_grid_fingerprints'])
   check(label+'/'+tag+'/origin_and_charge',q['properties']['origin_Bohr']==[0.,0.,0.] and q['cell']['charge']==2 and q['cell']['expected_electrons']==40)
   row=dict(result=q['result'],input=q['input'],output=q['output'],property=q['properties']['property_file'],field_au=q['cell']['field_au'],last_iteration_row=last,SCF_iterations=int(last[0]),energy_pre_final_hartree=epre,energy_final_electronic_hartree=epost,final_integration_delta_hartree=delta,printed_delta_hartree=corr[0],final_total_hartree=total,D4_hartree=d4,projected_mu_au=mu,grid_fingerprints=grids,XC_metadata=xc,printed_scf_criteria=q['parsed']['printed_scf_criteria'],native_stop=q['parsed']['native_scf_stop_signals'],native_mode=q['parsed']['native_scf_mode'],native_elapsed_seconds=q['elapsed_seconds'],allocated_CPU_seconds_proxy=q['allocated_elapsed_CPU_seconds_proxy'])
   if label=='original':
    for k in ('energy_pre_final_hartree','energy_final_electronic_hartree','final_integration_delta_hartree','projected_mu_au','D4_hartree'):
     check(label+'/'+tag+'/prior_diagnosis_'+k,row[k]==D(prior['rows'][tag][k]['decimal']))
   rows[label][tag]=row;excerpts[label][tag]={'first_line':max(0,end-8)+1,'text':'\n'.join(lines[max(0,end-8):end+22])}
  check(label+'/field_independent_D4',len({v['D4_hartree'] for v in rows[label].values()})==1)
  check(label+'/common_grids',all(v['grid_fingerprints']==rows[label]['zero']['grid_fingerprints'] for v in rows[label].values()))
  check(label+'/inputs_only_field',len({re.sub(r'^\s*EField\s+[^\n]+\n','',v,flags=re.M|re.I) for v in inputs[label].values()})==1)
  mu0=rows[label]['zero']['projected_mu_au'];calculations[label]={}
  for step,h in [('h',D('.001')),('half',D('.0005'))]:
   plus=rows[label]['plus_'+step];minus=rows[label]['minus_'+step];zero=rows[label]['zero'];r={}
   for stage,key in [('final','energy_final_electronic_hartree'),('pre_final_diagnostic','energy_pre_final_hartree')]:
    r['odd_'+stage]=-(plus[key]-minus[key])/(2*h);r['chi_'+stage]=-(plus[key]+minus[key]-2*zero[key])/h**2
   r.update(mu0=mu0,odd_final_minus_mu0=r['odd_final']-mu0,odd_pre_final_minus_mu0_diagnostic=r['odd_pre_final_diagnostic']-mu0,odd_final_integration_contribution=r['odd_final']-r['odd_pre_final_diagnostic'],chi_final_integration_contribution=r['chi_final']-r['chi_pre_final_diagnostic'],chi_dipole=(plus['projected_mu_au']-minus['projected_mu_au'])/(2*h),even_mu_numerator=plus['projected_mu_au']+minus['projected_mu_au']-2*mu0,printed_correction_odd_bound=(bound(plus['printed_delta_hartree'])+bound(minus['printed_delta_hartree']))/(2*h))
   # Existing scalar estimates include floating arithmetic; use their original print/binary propagation bounds.
   saved=col['response']['per_metal']['Dy']['estimates'][step]
   for k,sk in [('odd_final','odd_energy_moment'),('chi_final','chi_energy'),('chi_dipole','chi_moment'),('mu0','zero_moment')]:
    check(label+'/'+step+'/'+k,abs(float(r[k])-saved[sk])<=saved['roundoff_bounds'][sk])
   r['original_odd_gate_passes']=abs(r['odd_final_minus_mu0'])<=D('.0001')
   check(label+'/'+step+'/failed_gate_preserved',r['original_odd_gate_passes']==col['response']['checks']['Dy__odd_slope_'+step] and r['original_odd_gate_passes'] is False)
   calculations[label][step]=r
 changes={}
 for tag in diag.TAGS:
  before,after=rows['original'][tag],rows['matching_grid'][tag]
  check(tag+'/exact_grid_only_input',inputs['matching_grid'][tag].count(native.METHOD_BLOCK)==1 and inputs['matching_grid'][tag].replace(native.METHOD_BLOCK,'',1)==inputs['original'][tag])
  check(tag+'/exact_physical_grid_intervention',after['grid_fingerprints']==dict(before['grid_fingerprints'],grid2=before['grid_fingerprints']['grid3']) and after['grid_fingerprints']['grid2']==after['grid_fingerprints']['grid3'])
  check(tag+'/unchanged_XC_PostSCFGGA',before['XC_metadata']==after['XC_metadata'])
  check(tag+'/same_stopping_mode',before['native_stop']==after['native_stop'] and before['native_mode']==after['native_mode'])
  changes[tag]={k:after[k]-before[k] for k in ('energy_pre_final_hartree','energy_final_electronic_hartree','final_integration_delta_hartree','projected_mu_au')}
  changes[tag]['SCF_iteration_difference']=after['SCF_iterations']-before['SCF_iterations']
 stage_changes={step:{k:calculations['matching_grid'][step][k]-calculations['original'][step][k] for k in calculations['original'][step] if isinstance(calculations['original'][step][k],D)} for step in ('h','half')}
 model=old['model_comparison'];load(model['model_result']);load(model['model_review']);mc=model['model_center_estimates']['Dy'];scale={}
 for step in ('h','half'):
  mm=D(mc[step]['mu_E_decimal']);gap=mm-rows['matching_grid']['zero']['projected_mu_au'];res=abs(calculations['matching_grid'][step]['odd_final_minus_mu0'])
  scale[step]=dict(model_energy_conjugate_mu_au=mm,model_minus_new_zero_property_mu_au=gap,absolute_gap_over_odd_residual=abs(gap)/res,odd_residual_fraction_of_model_gap=res/abs(gap),model_minus_old_zero_property_mu_au=mm-rows['original']['zero']['projected_mu_au'],meaning='Descriptive scale only; neither residual nor recipe difference is a certified energy/moment error bound or target admission.')
 criteria={}
 for label,rr in rows.items():
  names=[x['name'] for x in rr['zero']['printed_scf_criteria']];criteria[label]={}
  for name in names:
   vv=[next(x for x in r['printed_scf_criteria'] if x['name']==name) for r in rr.values()]
   criteria[label][name]=dict(min=min(x['value'] for x in vv),max=max(x['value'] for x in vv),tolerances=sorted({x['tolerance'] for x in vv}),passed=sum(x['passes'] for x in vv),total=5)
  criteria[label]['SCF_iteration_counts']=[x['SCF_iterations'] for x in rr.values()]
 check('sources_unchanged',pin(cfg['new_collection']['path'])==cfg['new_collection'] and pin(cfg['old_collection']['path'])==cfg['old_collection'])
 result=dict(status='saved_grid_intervention_comparison_complete_both_Dy_panels_unqualified',implementation=pin(__file__),configuration=pin(args.config),declaration=cfg['declaration'],authorities={k:v for k,v in cfg.items() if isinstance(v,dict) and set(v)=={'path','sha256'}},checks=checks,check_count=len(checks),rows=rows,calculations=calculations,new_minus_original_by_endpoint=changes,new_minus_original_derived=stage_changes,printed_criteria_summary=criteria,model_scale=scale,matching_grid_actual_numerical_checks=new['response']['checks'],original_Dy_checks={k:v for k,v in old['response']['checks'].items() if k.startswith('Dy__')},new_recipe_native_and_collector_cost_reuse=review['scheduler_costs'],new_recipe_total_scheduler_CPU_seconds=review['total_scheduler_CPU_seconds'],new_native_calls=0,raw_outputs_read=10,corrected_target_admitted=False,Dy_target_qualified=False,paired_target_qualified=False,matching_La_released=False,new_recovery_released=False,unique_cause_identified=False,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 args.out.mkdir(parents=True,exist_ok=False);(args.out/'RESULT.json').write_text(json.dumps(diag.plain(result),indent=2)+'\n');(args.out/'RAW_EXCERPTS.json').write_text(json.dumps(excerpts,indent=2)+'\n')
 print(json.dumps(dict(result=pin(args.out/'RESULT.json'),check_count=len(checks),calculations=diag.plain(calculations),scale=diag.plain(scale),changes=diag.plain(changes),criteria=criteria,cost=result['cost']),indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();preexisted=args.out.exists();tick,cpu=time.monotonic(),time.process_time()
 try:main(args)
 except Exception as exc:
  if not preexisted:
   args.out.mkdir(parents=True,exist_ok=True);(args.out/'FAILURE.json').write_text(json.dumps(dict(status='failed_saved_grid_intervention_comparison',implementation=pin(__file__),configuration=pin(args.config),error=repr(exc),new_native_calls=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1)),indent=2)+'\n')
  raise
