#!/usr/bin/env python3
"""Single saved BF review; no evaluator construction, molecular call or solver."""
import copy,hashlib,importlib.util,json,math,os,re,time
from pathlib import Path
T0,C0=time.monotonic(),time.process_time()
ROOT=Path.cwd(); REL=Path('explicit_model_v1/radial_drive_shape_v1/boundary_fraction_witness_v1')
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL
OUT=W/'independent_actual_review_v1/actual_v1'
checks=[];pins={}
def check(ok,name):
    if not ok:raise AssertionError(name)
    checks.append(name)
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(q):
    p=Path(q['path']);v=pin(p);check(v==q,'pin '+str(p));pins[str(p)]=q;return p
def read(q):return json.loads(verify(q).read_text())
def load(q,name):
    s=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def close(a,b,name):check(abs(a-b)<=1e-12*max(1.,abs(a),abs(b)),name)
def save(name,x):
    p=OUT/name;check(not p.exists(),'new output '+name);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');return pin(p)
check(all(os.environ[k]=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')),'one review thread')
rpin=pin(W/'run_v1/execution/RESULT.json');r=read(rpin);cfg=read(r['configuration']);release=read(cfg['root_release']);source=read(cfg['source_review']);prep=read(cfg['preparation']);witness=read(r['witness']);v=read(r['vector'])
check(r['status']=='completed_boundary_fraction_native_test' and r['error'] is None and r['complete_native_admission'],'terminal complete native result')
check(r['implementation']==cfg['implementation']==source['implementation']==release['implementation'],'exact frozen implementation')
for k in ('implementation','adapter','core','resource_helper','scheduler_wrapper','parent_evaluator','native_source_manifest'):verify(cfg[k])
check(source['status']=='passed_boundary_fraction_witness_source_review' and source['configuration_fingerprint']==release['configuration_fingerprint']==r['configuration_fingerprint']==prep['configuration_fingerprint'],'source preparation release fingerprint')
check(source['witness']==release['witness']==prep['witness']==r['witness'] and source['preparation']==cfg['preparation'],'unchanged admitted witness')
for k in ('adapter','root_declaration','witness_inputs','fit_configuration','core','resource_helper','scheduler_wrapper','witness_specification','limits','output_directory','source_review','preparation'):check(release[k]==cfg[k],'release '+k)
fp=digest({k:x for k,x in cfg.items() if k not in ('root_release','source_review','preparation','fingerprint')});check(fp==r['configuration_fingerprint'],'configuration digest')
check(witness['fraction']==.0625 and witness['eligible'] and witness['solver_success'] is False,'one eligible1/16 AlmostSolved witness')
check(r['witness_parameters']==witness['parameters']==v['parameters'] and digest(v['parameters'])==r['witness_key']==witness['key'],'literal point unchanged')
check(r['solver_status_retained']=='AlmostSolved' and r['solver_success'] is False and r['original_fit_stop']=='no_actual_safeguarded_improvement','original solver and stop retained')
check(r['native_programs']==56 and all(r[k]==0 for k in ('solver_calls','new_derivative_probes','reference_calls','heldout_reads')),'56 native zero new solver/probes/reference')
check(r['selected_candidate'] is None and not r['automatic_promotion'] and not r['automatic_extension'],'no selection or extension')
parent=read(cfg['parent_configuration']);m=read(parent['manifest']);domain=read(parent['domain']);check(set(v['cells'])==set(m['rows']) and len(v['cells'])==46,'exact46 cells')
J=load(pin(D.parent.parent/'joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/JOURNAL.py'),'bf_saved_journal')
# Journal is also pinned by the admitted original adapter package.
old_adapter=load(cfg['original_adapter'],'bf_saved_original_adapter');check(pin(old_adapter.J.__file__)==pin(J.__file__),'qualified journal identity')
jr=r['journal_reconciliation'];jpath=verify(jr['journal']);ledger=read(jr['ledger']);replay=J.replay(jpath)
check(replay['ledger']==ledger and not replay['partial_suffix_bytes'] and not replay['unresolved_attempt_numbers'],'complete journal exact final ledger')
check(jr['final_authoritative_snapshot'] and jr['new_vectors']==1 and jr['new_native_programs']==jr['returned_programs']==56 and jr['imported_cells']==0,'final reconciled counts')
side=read(pin(Path(jr['ledger']['path']+'.journal.json')));check(side['ledger']==jr['ledger'] and side['sequence']==replay['sequence'] and side['last_event_sha256']==replay['last_event_sha256'],'final sidecar')
check(ledger['vector_keys']==[r['witness_key']] and not ledger['stopped'] and len(ledger['requests'])==1 and ledger['requests'][0]['status']=='complete','one successful request no stop')
check(ledger['native_program_invocations']==ledger['native_returncodes_observed']==len(ledger['native_attempts'])==56,'56 durable successful attempts')
check(ledger['limits']==cfg['limits']==r['limits'] and ledger['requests'][0]['required_new_programs']==56,'same hard budgets')
attempts={(a['cell_id'],a['program']):a for a in ledger['native_attempts']};check(len(attempts)==56 and all(a['status']=='returned_success' and a['returncode_observed'] for a in attempts.values()),'unique successful attempts')
N=load(parent['native_helper'],'bf_saved_native_parser');programs=0;gradients=0;raw_energy={};cache_pins=[]
for cid,c in v['cells'].items():
    row=m['rows'][cid];check(c['status']=='complete' and c['cell_id']==cid and not c['inherited_native'] and c['job_id']=='20734','new admitted cell '+cid)
    check(c['first_requested_vector']==v['parameters']==c['request_provenance']['vector'] and c['root_release']==cfg['root_release'],'cell candidate/release '+cid)
    recordpin=pin(Path(c['receipts'][0]['cwd']).parent/'RESULT.json');record=read(recordpin);cache_pins.append(recordpin)
    check({k:x for k,x in c.items() if k!='request_provenance'}==record,'immutable raw cell '+cid)
    check(c['cache_key']==digest(c['cache_identity']) and c['cache_identity']['key_sha256']==c['program_inputs']['analyze']['key']['sha256'],'cell cache key '+cid)
    vv=v['parameters'];i=0 if row['metal']=='La' else 4;j=8 if row['metal']=='La' else 10;aa,bb,u,w=vv[i:i+4]
    expected=[aa,bb,domain['S']*(u+w)/4,domain['S']*(u-w)/4,vv[j],vv[j+1]];e=vv[12 if row['metal']=='La' else 13]
    check(c['render']['physical_values_a_b_c_d_R_epsilon']==expected and c['render']['radial_drive']==e,'physical radial mapping '+cid)
    expected_programs=['analyze','testgrad'] if row['needs_gradient'] else ['analyze']
    check(list(c['program_inputs'])==expected_programs and len(c['receipts'])==len(expected_programs),'program denominator '+cid)
    for prog,receipt in zip(expected_programs,c['receipts']):
        programs+=1;inp=c['program_inputs'][prog];folder=Path(receipt['cwd']);args=['system.xyz','EM'] if prog=='analyze' else ['system.xyz','Y','N']
        check(receipt['argv']==[cfg['executables'][prog]['path']]+args and receipt['returncode']==0 and receipt['timeout_seconds']==300,'native argv return timeout '+cid+'/'+prog)
        check(all(x=='1' for x in receipt['threads'].values()) and receipt['executor']['job_id']=='20734' and receipt['executor']['native_workers']==12,'truthful native allocation '+cid+'/'+prog)
        check(attempts[cid,prog]['folder']==str(folder) and attempts[cid,prog]['cache_key']==c['cache_key'],'attempt receipt join '+cid+'/'+prog)
        for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):check(verify(inp[k])==folder/name,'literal native input '+cid+'/'+prog+'/'+k)
        check(inp['xyz']['sha256']==row['source_input']['xyz']['sha256'] and inp['parameters']['sha256']==row['source_input']['parameters']['sha256'],'unchanged geometry parameters '+cid+'/'+prog)
        key=Path(inp['key']['path']).read_text();mt=500 if row['metal']=='La' else 501
        check(('POLAR-LOCAL-RESPONSE '+str(mt)+' 36 '+' '.join(format(x,'.17g') for x in expected[:4])) in key and ('POLAR-RADIAL-DRIVE '+str(mt)+' 36 '+format(e,'.17g')) in key,'literal native response '+cid+'/'+prog)
        output=verify(receipt['output']);text=output.read_text();check(output==folder/(prog+'.out') and 'Tinker' in text and 'Version 26.2' in text and not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b',text),'raw native completion '+cid+'/'+prog)
        energies=re.findall(r'Total Potential Energy\s*:\s*('+N.NUMBER+')',text);check(len(energies)==1 and float(energies[0])==c['energy'],'raw energy '+cid+'/'+prog)
        raw_energy[cid]=float(energies[0])
    if row['needs_gradient']:
        parsed=N.parse(Path(c['receipts'][0]['output']['path']),Path(c['receipts'][1]['output']['path']),dict(atoms=row['atoms'],derivative_test=False))
        check(parsed==c['model']['native'] and parsed['force_kcal_mol_A']==c['model']['forces']==v['model'][cid]['forces'],'inherited analytic force parse '+cid);gradients+=1
    for q in c.get('exports',[]):verify(q)
check(programs==56 and gradients==10,'56raw programs ten analytic force pairs')
# Reconstruct the saved scalar gates and component cycles, without any matrix solve.
radcfg=read(cfg['radial_configuration']);sources=read(radcfg['inputs'])['sources'];core=load(sources['core'],'bf_saved_physical_core')
physical=core.physical(v['score']);check(physical['all25'],'original25 physical reconstruction')
actual=r['actual_metrics'];inc=r['incumbent_metrics']
for name,error in zip(physical['energy_names'],physical['energy']):check(float(error)==actual['original25']['energy/'+name]['error_kcal_mol'],'energy gate scalar '+name)
for name,fr in physical['force'].items():
    a=actual['original25']['force/'+name];check(fr['RMS']==a['rigid_RMS_kcal_mol_A'] and fr['maximum']==a['maximum_body_component_error_kcal_mol_A'] and a['passed']==(fr['RMS']<=2 and fr['maximum']<=5),'force gate scalar '+name)
dec=read(sources['decomposition']);inv=read(cfg['reuse_inventory']);fixed=dict(static={k:x['energy_kcal_mol'] for k,x in m['static'].items()},terms={k:x['model_terms'] for k,x in dec['parents'].items()},references={k:x['reference_total_kcal_mol'] for k,x in dec['paired_Dy_minus_La'].items()})
errors,parents=core.components(v,fixed);check(dict(zip(core.COMPONENT_ORDER,errors.tolist()))==actual['component_errors_kcal_mol'],'component signed errors')
for kind,ix in [('P',(0,2)),('B',(1,3))]:
    score=math.fsum(float(errors[i])**2 for i in ix)/2;close(score,actual['J_'+kind],'J '+kind);close(math.sqrt(score),actual['component_RMS'][kind],'RMS '+kind)
    close(score/inv['fixed_authorities']['q_denominators']['J_'+kind],actual['component_squared_error_ratios'][kind],'original denominator '+kind)
check(actual['q']==max(actual['component_squared_error_ratios'].values()) and r['actual_q_improvement']==inc['q']-actual['q'] and r['strict_q_improvement_tolerance']==1e-10*max(1,abs(inc['q'])),'literal q and improvement')
for k,entry in actual['fixed_electrical_safeguards'].items():
    check(entry['score']==v['score']['block_objectives'][core.SCORE_KEYS[k]] and entry['ceiling']==inv['fixed_safeguards']['actual_scores'][k] and entry['passed']==(entry['score']<=entry['ceiling']),'fixed electrical gate '+k)
    close(r['actual_margins']['electrical'][k],entry['ceiling']-entry['score'],'electrical margin '+k)
C=v['response']['C_dipoles_Debye'];signed=actual['signed_W9_z_safeguard'];check(signed['model_Dy_minus_La_z_Debye']==math.fsum([C['aquo9_central_Dy'][2],-C['aquo9_central_La'][2]]) and signed['model_minus_reference_Debye']==math.fsum([C['aquo9_central_Dy'][2],-C['aquo9_central_La'][2],-signed['reference_Dy_minus_La_z_Debye']]),'signed W9 raw response')
check(signed['absolute_ceiling_Debye']==inv['fixed_safeguards']['W9_absolute_ceiling_Debye'] and signed['reference_Dy_minus_La_z_Debye']==inv['fixed_safeguards']['W9_reference_Debye'] and signed['passed']==(abs(signed['model_minus_reference_Debye'])<=signed['absolute_ceiling_Debye']),'unchanged signed W9 gate')
for kind,key in [('pair_mu','pair_Debye'),('C_mu','C_Debye'),('chi','chi_a0_cubed')]:close(math.sqrt(actual['original_block_scores'][kind]),actual['electrical_residual_RMS'][key],'absolute electrical RMS '+key)
for k,g in actual['original25'].items():
    expected={'absolute_energy_margin_kcal_mol':1-abs(g['error_kcal_mol'])} if k.startswith('energy/') else {'RMS_margin_kcal_mol_A':2-g['rigid_RMS_kcal_mol_A'],'maximum_margin_kcal_mol_A':5-g['maximum_body_component_error_kcal_mol_A']}
    check(r['actual_margins']['physical'][k]==expected and all(x>=0 for x in expected.values()),'literal physical margin '+k)
check(r['witness_scientific_passed']==(actual['all_declared_safeguards_pass'] and actual['q']<inc['q']-r['strict_q_improvement_tolerance']),'one point actual decision')
check(actual['all25_original_physical_gates'] and actual['physical_constraints_pass'] and actual['all_declared_safeguards_pass'] and not actual['old_score_ceilings_enforced'] and not actual['electrical_physical_acceptance_established'],'correct scientific gate semantics')
accounting=pin(OUT/'SACCT.psv');rows=[x.split('|') for x in verify(accounting).read_text().splitlines()];parentrow=next(x for x in rows if x[0]=='20734');batchrow=next(x for x in rows if x[0]=='20734.batch')
check(parentrow[2:6]==['COMPLETED','0:0','13','13'] and batchrow[2]=='COMPLETED','actual scheduler terminal')
cost=dict(native_programs=56,new_vectors=1,new_cells=46,wall_seconds=r['cost']['wall_seconds'],controller_CPU_seconds=r['cost']['controller_CPU_seconds'],child_user_CPU_seconds=r['cost']['child_user_CPU_seconds'],child_system_CPU_seconds=r['cost']['child_system_CPU_seconds'],controller_plus_child_CPU_seconds=sum(r['cost'][k] for k in ('controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds')),scheduler_parent_elapsed_seconds=int(parentrow[4]),scheduler_parent_allocated_CPUs=int(parentrow[5]),scheduler_allocated_CPU_seconds=int(parentrow[4])*int(parentrow[5]),scheduler_parent_total_CPU=parentrow[6],batch_MaxRSS=batchrow[9],sacct=accounting,per_program_CPU_not_summed=True)
check(r['allocation']['job_id']=='20734' and r['allocation']['allocated_CPUs']==13 and r['allocation']['native_workers']==12 and r['allocation']['scheduler_fields']['QOS']=='standard_expand','actual original resource envelope')
summary=dict(incumbent_q=inc['q'],actual_q=actual['q'],absolute_q_reduction=r['actual_q_improvement'],relative_q_reduction=r['actual_q_improvement']/inc['q'],predicted_q=witness['prediction']['actual_component_q'],actual_component_RMS=actual['component_RMS'],component_errors_kcal_mol=actual['component_errors_kcal_mol'],electrical_RMS=actual['electrical_residual_RMS'],signed_W9=signed,La8=actual['original25']['force/La/aquo8_central'],La8_max_margin=r['actual_margins']['physical']['force/La/aquo8_central']['maximum_margin_kcal_mol_A'],balanced_W9=actual['original25']['energy/balanced_exchange/aquo9_central']['error_kcal_mol'],actual_minus_predicted=r['actual_minus_predicted'])
checks_pin=save('CHECKS.json',dict(checks=checks,pins=list(pins.values()),count=len(checks)))
review=dict(status='passed_boundary_fraction_witness_actual_review',result=rpin,configuration=r['configuration'],implementation=r['implementation'],adapter=cfg['adapter'],source_review=cfg['source_review'],root_release=cfg['root_release'],preparation=cfg['preparation'],witness=r['witness'],vector=r['vector'],selected_candidate=None,witness_key=r['witness_key'],fraction=.0625,complete_native_admission=True,all25_physical_passed=True,fixed_electrical_and_signed_safeguards_passed=True,strict_q_improvement=True,witness_scientific_passed=True,original_fit_stop_preserved=True,solver_status='AlmostSolved',solver_success=False,automatic_promotion=False,automatic_extension=False,counts=dict(new_vectors=1,cells=46,programs=56,gradients=10,solver_calls=0,probes=0,reference_calls=0,imported_cells=0,journal_events=replay['sequence']),journal=jr,checks=checks_pin,summary=summary,cost=cost,review_cost=dict(wall_seconds=time.monotonic()-T0,CPU_seconds=time.process_time()-C0,threads=1,native_calls=0,solver_calls=0),scope='New56 output/input/receipt hashes; ten raw analytic force parser replays; original physical/component/electrical scalar gates and journal. No original fit/raw46/stencil replay, no scalar certificate recomputation, no susceptibility solve or eigenspectrum replay.',interpretation='The original1/8 fraction floor missed one actual feasible improving point along its saved step8 direction. The La8 force maximum remains only1.06e-5 below its unchanged5 limit. This is finite consumed development, not global convergence, model accuracy, transfer qualification or a reason to change fixed periodic/reference candidates.')
save('REVIEW.json',review);print(json.dumps(dict(status=review['status'],checks=len(checks),review=pin(OUT/'REVIEW.json'),summary=summary,cost=cost),indent=2))
