#!/usr/bin/env python3
"""Collect two terminal, pinned campaigns. No kernel, optimizer or native calls."""
import argparse,csv,hashlib,importlib.util,io,json,math,os,platform,resource,subprocess,time,traceback
from pathlib import Path

ROLES=('shared12_control','metal15')
CORE_HASH={'shared12_control':'2aab54fe47c6c1f283fbdcda6d75e6303a441e0b2a3e2d2483ce4a18cd3c5fd0',
           'metal15':'a2fd7953f5cfb724ae6c92a2b878e0f7fd154fa76cb623f7dd12fbeda9b46b83'}
EXTRA_HASH='3355c4ac575144b964978736efb35974f10b45e0525cc937666d9b9efcc52f85'
SUPPORT_HASH='813b5a0aa081f502f3c247d798674d811b0e0239aa049f8ce946ffab76a2838e'
VERIFIED={}

def pin(path):
 p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def checked(p):
 if not isinstance(p,dict) or set(p)!={'path','sha256'} or pin(p['path'])!=p:raise ValueError('changed or invalid pin '+str(p))
 VERIFIED[p['path']]=p['sha256'];return Path(p['path'])
def read(p):return json.loads(checked(p).read_text())
def require(ok,message):
 if not ok:raise ValueError(message)
def write(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def load_extra(p):
 require(p['sha256']==EXTRA_HASH,'exact report-only helper');path=checked(p)
 spec=importlib.util.spec_from_file_location('joint_saved_report_extra',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def load_support(p):
 require(p['sha256']==SUPPORT_HASH,'exact source-bound support diagnostic');path=checked(p)
 spec=importlib.util.spec_from_file_location('joint_saved_support',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def allocation():
 job=os.environ.get('SLURM_JOB_ID','');require(job.isdigit() and os.environ.get('SLURM_CPUS_PER_TASK')=='1','one-CPU Slurm collection required')
 for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
  require(os.environ.get(key)=='1','one-thread setting '+key)
 raw=subprocess.run(['/usr/local/bin/scontrol','show','job','-o',job],capture_output=True,text=True,check=True).stdout
 row=dict(x.split('=',1) for x in raw.split() if '=' in x)
 for key,value in {'Partition':'standard','QOS':'standard_expand','Account':'banfield','NumCPUs':'1','NumNodes':'1','NumTasks':'1'}.items():require(row.get(key)==value,'allocation '+key)
 return {'job_id':job,'host':platform.node(),'scontrol':raw,'allocated_CPUs':1}

def accounting(p,job):
 if p is None:return {'status':'unbound','allocated_CPU_seconds':None,'scope':'No scheduler cost inferred from summed process CPU.'}
 raw=checked(p).read_text();rows=list(csv.DictReader(io.StringIO(raw),delimiter='|'))
 parent=[r for r in rows if r.get('JobIDRaw',r.get('JobID'))==str(job)]
 require(len(parent)==1,'one exact scheduler parent row');r=parent[0]
 seconds=int(r['ElapsedRaw']);cpus=int(r['AllocCPUS']);require(cpus==24,'campaign parent has 24 allocated CPUs')
 cost=seconds*cpus
 if r.get('CPUTimeRAW'):require(int(r['CPUTimeRAW'])==cost,'parent allocated CPU cost closure')
 return {'status':'saved','source':p,'parent':r,'allocated_CPU_seconds':cost,
         'scope':'Only exact parent row; .batch/.extern excluded. Measured TotalCPU retained separately as printed.'}

def model(role,entry):
 r=read(entry['result']);require(r['status'] in ('completed','completed_no_stable_candidate') and r['error'] is None,'terminal campaign required '+role)
 cfg=read(r['configuration']);src=cfg['sources'];release=read(r['release'])
 require(cfg['role']==role and src['runner']==r['implementation'],'campaign role/source join')
 require(release['status']=='released_joint_demonstration_v2','declared campaign release')
 require(release['configuration']==r['configuration'] and release['implementation']==r['implementation'] and release['output']==cfg['output'],'root release joins')
 require(Path(entry['result']['path']).parent==Path(cfg['output']),'terminal result directory')
 for p in src.values():checked(p)
 require(src['core']['sha256']==CORE_HASH[role],'exact qualified model core')
 for name in ('new_native_calls','new_quantum_calls','new_geometry_searches'):require(r[name]==0,'saved-only campaign '+name)
 settings=read(src['settings']);prov=read(r['data_provenance']);starts=read(r['starts']);grad=read(r['gradient_check'])
 require(len(starts['starts'])==16 and grad['status']=='passed','fixed starts and completed derivative check')
 require(r['operator_count']==64 and prov['operator_count']==64 and len(prov['fixture_pins'])==64,'complete operator denominator')
 require(prov['provenance']==read(src['data_inputs']),'saved DATA provenance')
 for p in prov['fixture_pins'].values():checked(p)
 checked(r['symmetry']);n=12 if role==ROLES[0] else 15
 require(len(settings['coefficients'])==n and len(settings['coefficient_scales'])==n,'coefficient dimension')
 baseline=read(r['baseline']);require(baseline['coefficients']==[0.]*n,'literal zero baseline')
 rows=r['start_results'];require(sorted(x['index'] for x in rows)==list(range(16)),'all sixteen start outcomes required')
 outcomes=[];eligible=[];counts={key:0 for key in ('H_SLSQP','h_BVLS','V_BVLS','joint_SLSQP')};byindex={}
 for wrapper in sorted(rows,key=lambda x:x['index']):
  if 'result' not in wrapper:
   outcomes.append({'index':wrapper['index'],'wrapper':wrapper,'raw_solver_records':None,'gap':'Worker returned no pinned terminal record; failure preserved.'});continue
  q=read(wrapper['result']);byindex[q['index']]=(q,wrapper['result'])
  require(q['index']==wrapper['index'] and q['status']==wrapper['status'] and q['stage']=='terminal','terminal start identity')
  solvers={k:q.get(k) for k in ('H_solver','h_solver','V_solver','joint_solver')}
  for key,field in [('H_SLSQP','H_solver'),('h_BVLS','h_solver'),('V_BVLS','V_solver'),('joint_SLSQP','joint_solver')]:counts[key]+=int(field in q)
  candidates={name:{'status':c['status'],'score':c.get('assessment',{}).get('score'),'coefficients':c.get('coefficients'),
                    'at_exact_box_boundary':c.get('at_exact_box_boundary'),'SPD_margins':c.get('SPD_margins')}
              for name,c in q.get('candidates',{}).items()}
  if q['status']=='complete_stable_candidate':
   stage=q['selected_stage'];c=q['candidates'][stage]
   require(c['status']=='complete_stable_candidate' and q['assessment']==c['assessment'] and q['coefficients']==c['coefficients'],'selected stage joins')
   permitted=[(v['assessment']['score'],0 if k=='conditional' else 1,k) for k,v in q['candidates'].items() if v['status']=='complete_stable_candidate']
   require(min(permitted)[2]==stage,'recorded within-start selection rule')
   require(wrapper['score']==q['assessment']['score'],'recorded terminal scalar score')
   eligible.append((q['assessment']['score'],q['index']))
  outcomes.append({'index':q['index'],'result':wrapper['result'],'status':q['status'],'stage':q['stage'],
                   'selected_stage':q.get('selected_stage'),'score':q.get('assessment',{}).get('score'),
                   'solver_attempts':q.get('solver_attempts',[]),'raw_solver_records':solvers,'candidates':candidates,
                   'error':q.get('error'),'joint_error':q.get('joint_error'),'failed_linear_stage':q.get('failed_linear_stage'),
                   'conditional_objective_join':q.get('conditional_objective_join'),'cost':q.get('cost')})
 require(counts==r['returned_solver_records'],'returned solver accounting')
 selected=ablation=None
 if r['selected_candidate'] is not None:
  pick=r['selected_candidate'];require(eligible and min(eligible)==(pick['score'],pick['start_index']),'terminal recorded campaign selection')
  q,qp=byindex[pick['start_index']];require(qp==pick['result'] and q['selected_stage']==pick['stage'],'selected source record')
  selected=q['candidates'][pick['stage']];require(selected['coefficients']==pick['coefficients'],'selected exact coefficients')
  ablation=read(r['V_only_ablation']);require(ablation['refitted'] is False,'unfitted V-only ablation')
  v=[x if settings['coefficient_scales'][i]==100 else 0. for i,x in enumerate(pick['coefficients'])]
  require(ablation['coefficients']==v,'V-only preserves each selected V and zeros h/H')
 else:require(not eligible and r['status']=='completed_no_stable_candidate','no hidden eligible candidate')
 for label,c in [('baseline',baseline),('selected',selected),('V_only',ablation)]:
  if c is None:continue
  require(len(c['coefficients'])==n and all(isinstance(x,(int,float)) and math.isfinite(x) and abs(x)<=s for x,s in zip(c['coefficients'],settings['coefficient_scales'])),'literal parameter box '+label)
  require(set(c['predictions'])==set(prov['fixture_pins']),'all64 saved predictions '+label)
  for id,row in c['predictions'].items():
   require(row['fixture']==prov['fixture_pins'][id],'prediction fixture pin '+id)
   require(row['minimum_normalized_eigenvalue']>=settings['normalized_SPD_floor'],'saved normalized SPD floor '+id)
 return {'result_pin':entry['result'],'result':r,'configuration':cfg,'sources':src,'settings':settings,'provenance':prov,
         'starts':starts,'baseline':baseline,'selected':selected,'V_only':ablation,'all_start_outcomes':outcomes,
         'accounting':accounting(entry.get('scheduler_accounting'),r['allocation']['job_id'])}

def search_coverage(outcomes):
 groups={'returned_H_endpoints':{},'complete_conditional_points':{}}
 for row in outcomes:
  h=((row.get('raw_solver_records') or {}).get('H_solver') or {}).get('solution')
  c=row.get('candidates',{}).get('conditional',{})
  for key,point in [('returned_H_endpoints',h),('complete_conditional_points',c.get('coefficients') if c.get('status')=='complete_stable_candidate' else None)]:
   if point is None:continue
   identity=json.dumps(point,separators=(',',':'));groups[key].setdefault(identity,{'point':point,'start_indices':[]})['start_indices'].append(row['index'])
 return {key:{'distinct_literal_count':len(value),'groups':list(value.values())} for key,value in groups.items()}|{
  'interpretation':'Exact stored floating-point tuples, without rounding or basin clustering. Sixteen starts do not establish sixteen independent basins or global coverage.'}

def report_bundle(m):
 inputs=read(m['sources']['data_inputs'])
 for p in inputs.values():checked(p)
 ds=read(inputs['water_dataset'])
 for p in ds.values():
  if isinstance(p,dict) and set(p)=={'path','sha256'}:checked(p)
 cycles=read(ds['cycles']);base=read(ds['baseline']);static=read(ds['static'])
 bundle={'pair_cycles':cycles['pair'],'water_parent_components_and_cycles':cycles['parent'],
         'water_common_basis_components':cycles['common_basis_decomposition'],'units':base['units']}
 return bundle,static

def view(c,settings,bundle,static,E):
 if c is None:return None
 pred={id:{'energy':r['energy_kcal_mol'],'mu':r['mu_Debye']} for id,r in static.items()}
 pred.update({id:{'energy':r['energy_kcal_mol'],'mu':r['mu_conjugate_Debye']} for id,r in c['predictions'].items()})
 coefficients=[{'name':name,'value':value,'bound':bound,'exact_bound_hit':abs(value)==bound}
               for name,value,bound in zip(settings['coefficients'],c['coefficients'],settings['coefficient_scales'])]
 return {'assessment':c['assessment'],'coefficients':coefficients,'parent_components':E.report(pred,bundle),
         'response':{'minimum_normalized_eigenvalue':min(r['minimum_normalized_eigenvalue'] for r in c['predictions'].values()),
                     'maximum_saved_solve_residual':max(r['solve_residual'] for r in c['predictions'].values())},
         'donor_channels':{id:{k:r[k] for k in ('mu_conjugate_Debye','mu_native_physical_u_Debye','fixture')}
                           for id,r in c['predictions'].items() if id in ('matched_dry_acetate_mono_central_La','matched_dry_acetate_mono_central_Dy')}}

def fmt(x):
 if x is None:return 'NA'
 if isinstance(x,(list,tuple)):return '('+', '.join(fmt(v) for v in x)+')'
 if isinstance(x,(int,float)):return f'{x:.8g}'
 return str(x)
def text_report(result):
 views=result['views'];main=[k for k in ('baseline','shared12','metal15') if views.get(k) is not None];lines=[]
 def line(x=''):lines.append(x)
 def table(title,header,rows):
  line(title);line(' | '.join(header))
  for row in rows:line(' | '.join(fmt(x) for x in row))
  line()
 line('FIXED SHORT-RANGE MODEL COMPARISON — CONSUMED DEVELOPMENT')
 line('Saved predictions only. Errors are prediction minus target; signed moment vectors are retained.');line()
 table('Five-family RMS (unit scales are fitting weights, not accuracy gates)',
       ['model','score','chi a0^3','mu D','E2 kcal/mol','E3 kcal/mol','slope kcal/mol/A'],
       [[name,v['assessment']['score'],*[v['assessment']['family_RMS_in_declared_units'][k] for k in ('chi','mu','E2','E3','slope')]] for name,v in views.items() if v is not None])
 first=views['baseline']['assessment']
 table('Individual E3 energies: all eight cycles; W9 rows were fitted', ['cycle','target']+[k+' prediction/error' for k in main],
       [[id,row['energy_kcal_mol']['target']]+[[views[k]['assessment']['trimers'][id]['energy_kcal_mol'][f] for f in ('prediction','error')] for k in main] for id,row in first['trimers'].items()])
 table('Individual signed mu3 vectors (Debye)', ['cycle','target xyz']+[k+' xyz/error xyz' for k in main],
       [[id,row['mu_Debye']['target']]+[[views[k]['assessment']['trimers'][id]['mu_Debye'][f] for f in ('prediction','error')] for k in main] for id,row in first['trimers'].items()])
 table('Paired Dy−La trimer contrasts', ['cycle','quantity','target']+[k+' prediction/error' for k in main],
       [[id,quantity,row[quantity]['target']]+[[views[k]['assessment']['trimer_Dy_minus_La'][id][quantity][f] for f in ('prediction','error')] for k in main] for id,row in first['trimer_Dy_minus_La'].items() for quantity in ('energy_kcal_mol','mu_Debye')])
 table('Five conditional susceptibility modes (source a0^3)', ['mode','target']+[k+' prediction/error' for k in main],
       [[id,row['target']]+[[views[k]['assessment']['chi'][id][f] for f in ('prediction','error')] for k in main] for id,row in first['chi'].items()])
 table('La donor slopes at both saved steps (kcal/mol/A)', ['direction/step A','target']+[k+' prediction/error' for k in main],
       [[direction+'/'+tag,row['target']]+[[views[k]['assessment']['donor_slopes'][direction][tag][f] for f in ('prediction','error')] for k in main] for direction,steps in first['donor_slopes'].items() for tag,row in steps.items()])
 table('Full signed donor conjugate moment vectors (Debye; Dy report-only)', ['metal','target xyz']+[k+' xyz/error xyz' for k in main],
       [[metal,first['individual_signed_moments']['donor/'+metal+'_acetate']['target']]+[[views[k]['assessment']['individual_signed_moments']['donor/'+metal+'_acetate'][f] for f in ('prediction','error')] for k in main] for metal in ('La','Dy')])
 for name in main:
  extra=views[name]['parent_components']
  table(name+' parent C/P/B: electronic common-basis targets, D4 excluded', ['parent/component','E prediction/target/error kcal/mol','mu prediction/target/error xyz D'],
        [[parent+'/'+component,[r['energy_kcal_mol'][f] for f in ('prediction','target','error')],[r['mu_Debye'][f] for f in ('prediction','target','error')]] for parent,parts in extra['parents'].items() for component,r in parts.items()])
  table(name+' four-parent contrast Dy9−La9−Dy8+La8', ['component','E prediction/target/error','mu prediction/target/error xyz'],
        [[component,[r['energy_kcal_mol'][f] for f in ('prediction','target','error')],[r['mu_Debye'][f] for f in ('prediction','target','error')]] for component,r in extra['four_parent_contrast']['components'].items()])
  table(name+' coefficients (literal bounds)', ['name','value','absolute bound','at exact bound'],
        [[c['name'],c['value'],c['bound'],c['exact_bound_hit']] for c in views[name]['coefficients']])
 for role,m in result['models'].items():
  table(role+' all16 raw solver outcomes (unsuccessful feasible endpoints remain labelled)',
        ['start','state','H status/success','h status/rank','V status/rank','joint status/success','selected stage','score'],
        [[r['index'],r.get('status',r.get('wrapper',{}).get('status')),
          *[[((r.get('raw_solver_records') or {}).get(stage) or {}).get(a),((r.get('raw_solver_records') or {}).get(stage) or {}).get(b)] for stage,a,b in [('H_solver','status','success'),('h_solver','status','rank'),('V_solver','status','rank'),('joint_solver','status','success')]],
          r.get('selected_stage'),r.get('score')] for r in m['all_start_outcomes']])
  line(role+' allocated parent CPU-s: '+fmt(m['scheduler_accounting']['allocated_CPU_seconds']))
  line(role+' distinct literal returned H endpoints / complete conditional points: '+fmt([m['search_coverage'][k]['distinct_literal_count'] for k in ('returned_H_endpoints','complete_conditional_points')]))
 support=result['support_diagnostic']
 table('Unbounded stiffness support: optimistic chi floors (ordinary double, source a0^3)',
       ['mode','donor O count','baseline','target','current H-support floor','all donor-O Cartesian floor'],
       [[id,r['donor_count'],r['baseline_projected_chi_source_a0_cubed'],r['reference_chi_source_a0_cubed'],r['optimistic_unbounded_stiffness_floor_source_a0_cubed'],r['full_oxygen_support']['optimistic_unbounded_chi_floor_source_a0_cubed']] for id,r in support['rows'].items()])
 table('Signed moment residual perpendicular to accessible response directions (Debye)',
       ['mode/reference','current-support xyz','current-support norm','all donor-O norm'],
       [[id+'/'+r['reference_id'],r['perpendicular_component_Debye'],r['irreducible_norm_for_current_support_Debye'],r['irreducible_norm_full_oxygen_support_Debye']] for id,mode in support['rows'].items() for r in mode['moment_targets']])
 line('Support floors allow unbounded stiffness beyond fitted boxes. SVD cutoffs/singular values and exact donor/frame maps are retained in RESULT.json. They are not rigorous interval certificates or unique quantum-mechanism assignments.');line()
 line();line('LIMITS / UNRESOLVED GAPS')
 for gap in result['interpretation_limits']:line('- '+gap)
 return '\n'.join(lines)+'\n'

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);parser.add_argument('--out',required=True);args=parser.parse_args()
 out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
 result={'status':'started','error':None,'implementation':pin(__file__),'configuration':pin(args.config),'new_native_calls':0,'new_quantum_calls':0,'new_response_energy_evaluations':0,'new_optimizer_calls':0};code=1
 try:
  cfg=json.loads(Path(args.config).read_text());require(cfg['implementation']==pin(__file__) and cfg['output']==str(out),'collector identity/output')
  require(set(cfg['models'])==set(ROLES),'exact two-model denominator');result['allocation']=allocation();E=load_extra(cfg['report_extra'])
  models={role:model(role,cfg['models'][role]) for role in ROLES};a,b=[models[k] for k in ROLES]
  for key in ('data','data_inputs','endpoint_map','kernel','declaration'):require(a['sources'][key]==b['sources'][key],'same scientific authority '+key)
  require(a['provenance']==b['provenance'] and a['starts']['starts']==b['starts']['starts'],'same64 data and16 fixed starts')
  require(a['baseline']['predictions']==b['baseline']['predictions'] and a['baseline']['assessment']==b['baseline']['assessment'],'exact common zero-baseline predictions')
  bundle,static=report_bundle(a);views={'baseline':view(a['baseline'],a['settings'],bundle,static,E)}
  for role,label in zip(ROLES,('shared12','metal15')):
   m=models[role];views[label]=view(m['selected'],m['settings'],bundle,static,E);views[label+'_V_only']=view(m['V_only'],m['settings'],bundle,static,E)
  result['views']=views;result['baseline_identity']='Both saved baseline predictions and assessments are literally equal; coefficient zero-vector dimensions differ.'
  result['models']={role:{'result':m['result_pin'],'configuration':m['result']['configuration'],'sources':m['sources'],
                         'settings':m['settings'],'selected_candidate':m['result']['selected_candidate'],'all_start_outcomes':m['all_start_outcomes'],
                         'search_coverage':search_coverage(m['all_start_outcomes']),
                         'returned_solver_records':m['result']['returned_solver_records'],'scheduler_accounting':m['accounting'],
                         'producer_process_cost':m['result']['cost'],'completed_worker_CPU_seconds':m['result']['completed_worker_CPU_seconds']}
                    for role,m in models.items()}
  S=load_support(cfg['support_diagnostic']);result['support_diagnostic_source']=cfg['support_diagnostic']
  result['support_diagnostic']=S.analyze(a['sources'],a['baseline']['predictions'])
  result['interpretation_limits']=[
   'Both models use the same consumed64-geometry data. More flexible energy interpolation is not independent validation or a global capacity proof.',
   'All32 finite starts and raw solver statuses are retained. Shared endpoints are counted exactly; sixteen starts do not imply sixteen independent basins. Feasible unsuccessful endpoints do not establish convergence or global optimality.',
   'Parent component errors, Dy−La differences, and net C=P+B cancellation are distinct; these cluster cycles are not exchange free energies.',
   'Water energies use electronic common-parent-basis targets without D4; donor relative energies/slopes use their original PBE0-D4 protocol.',
   'Conjugate moment uses (u+v)/2; native physical-u moment is recorded separately. Charged-system origin and reference-state qualifications remain inherited.',
   'Only five susceptibility modes are targeted, not a qualified full response tensor. Failed Dy susceptibility and unavailable Dy canonical torque are excluded.',
   'Both donor slope steps are reported; only half-step La slopes are fitted. Dy donor moment is report-only, with no Dy mechanical fitting.',
   'Water full coordinate forces are unavailable from these frozen-operator predictions; no force preservation, periodic dynamics or production qualification follows.',
   'All coefficient bounds and response-stability constraints are unchanged. Fitting scales are not physical accuracy gates; no tolerance was relaxed.',
   'Tightened W9far references and canonical donor results remain separate sensitivities; original targets were not replaced.']
  result.update(status='completed_saved_joint_model_comparison',all_start_outcomes_count=sum(len(m['all_start_outcomes']) for m in models.values()),report_extra=cfg['report_extra'],verified_source_and_data_files=VERIFIED)
  write(out/'RESULT.json',result);(out/'REPORT.txt').write_text(text_report(result));result['report']=pin(out/'REPORT.txt');code=0
 except Exception as exc:result.update(status='failed_saved_comparison',error=repr(exc),traceback=traceback.format_exc())
 finally:
  result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'own_scheduler_cost':'root captures terminal parent sacct'}
  write(out/'RESULT.json',result)
 return code
if __name__=='__main__':raise SystemExit(main())
