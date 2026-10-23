#!/usr/bin/env python3
"""Frozen quantum10 source and literal preparation review; no native or test-suite replay."""
import ast,hashlib,json,math,re,time
from pathlib import Path
ROOT=Path.cwd();REL=Path('reference_chemistry_v1/dry_donor_susceptibility_v1');DB=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL;WB=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL;OUT=WB/'independent_review_v1';OUT.mkdir(exist_ok=False);checks=[];t=time.monotonic();cpu=time.process_time()
def pin(p):
 p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def ck(q,n):
 checks.append(dict(check=n,passed=bool(q)));assert q,n
def verify(p):ck(pin(p['path'])==p,'pin:'+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
ready=json.loads((DB/'READY.json').read_text())
for v in ready.values():
 if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
p=read(ready['panel']);prep=read(ready['preparation']);cfg=read(ready['configuration']);inv=read(p['source_inventory']);controls=read(ready['controls']);modelreview=read(p['model_actual_review']);model=read(inv['model40_result'])
ck(p['adapter']==ready['implementation']==cfg['implementation'],'source joins')
ck(p['configuration']==ready['configuration'] and cfg['source_inventory']==p['source_inventory']==ready['source_inventory'] and p['specification']==ready['specification'],'configuration joins')
for k,v in cfg.items():
 if isinstance(v,dict) and set(v)=={'path','sha256'}:
  verify(v)
  if k in p:ck(p[k]==v,'helper bind:'+k)
ck(prep['status']=='passed_static_preflight' and prep['panel']==ready['panel'] and prep['declared']==10 and prep['new_native_calls']==0 and prep['old_SP_zeros_reused']==0 and prep['old_central_raw_replays']==2,'actual zero-call preparation')
ck(controls['status']=='passed_saved_controls' and controls['implementation']==p['adapter'] and controls['panel']==ready['panel'] and controls['new_native_calls']==0 and controls['check_count']==len(controls['checks'])==37 and all(x is True or isinstance(x,dict) and x['passed'] for x in controls['checks'].values()),'37 saved author controls inspected only')
ck(p['centers']==inv['centers'] and p['field_direction_lab']==inv['field_direction_lab'],'original centers/fixed tangent')
ck(len(p['cells'])==p['maximum_native_calls']==10 and p['ranks']==4 and p['nodes']==1,'finite ten/4MPI/one node')
ck(p['scout_cell_ids']==inv['first_stage_ids']==ready['first_stage_ids'] and p['remaining_cell_ids']==inv['remaining_stage_ids']==ready['remaining_stage_ids'] and len(p['scout_cell_ids'])==len(p['remaining_cell_ids'])==5 and not set(p['scout_cell_ids'])&set(p['remaining_cell_ids']),'La5 then Dy5')
ck(p['centers']['La']['geometry']['positions_A']==p['centers']['Dy']['geometry']['positions_A'] and abs(math.fsum(x*x for x in p['field_direction_lab'])-1)<1e-14,'paired unchanged geometry/unit direction')
ck(modelreview['status']=='passed_dry_acetate_electrical_actual_review' and modelreview['result']==inv['model40_result'] and model['complete_cells']==40 and model['summary']['numerical_checks_passed'],'admitted model40 separate authority')
fields={'zero':0.,'minus_h':-.001,'minus_half':-.0005,'plus_half':.0005,'plus_h':.001};seen=set()
for c in p['cells']:
 cid=c['cell_id'];metal=c['metal'];center=p['centers'][metal];g=center['geometry'];f=fields[c['field_tag']];seen.add((metal,c['field_tag']))
 ck(c['positions_A']==g['positions_A'] and c['symbols']==g['symbols'] and c['source_geometry']==g['xyz'],cid+'/exact original frame')
 coords=verify(c['source_geometry']).read_text().splitlines()[2:];ck([[float(v) for v in l.split()[1:4]] for l in coords if l.strip()]==c['positions_A'],cid+'/literal XYZ')
 ck(c['charge']==2 and c['multiplicity']==1 and c['expected_electrons']==40 and c['expected_effective_nuclear_charges']==[11,6,6,8,8,1,1,1],cid+'/charged ECP state')
 ck(c['expected_basis_functions']==189 and c['expected_basis_functions_by_center']==[47,31,31,31,31,6,6,6] and c['expected_aux_J_basis_functions']==(580 if metal=='La' else 598),cid+'/AO AuxJ')
 ck(not any(c['ghost_mask']) and c['expected_ecp_core']==(46 if metal=='La' else 55) and c['physical_multiplicity_of_parent']==(1 if metal=='La' else 6),cid+'/real metal effective/physical states')
 ck(c['field_au']==f and c['field_vector_au']==[f*x for x in p['field_direction_lab']],cid+'/nominal signed field')
 ck(c['explicit_field_origin']==c['explicit_property_origin']==[0,0,0] and c['new_single_point'] and c['old_zero_energy_reuse'] is False,cid+'/new SP zero and lab origin')
 old=verify(center['zero_EnGrad_input']).read_text();ck(c['source_input']==center['zero_EnGrad_input'],cid+'/original input authority')
 original=re.sub(r'^%maxcore \d+\s*$', '%maxcore 3000',old,flags=re.M);original=re.sub(r'^%pal nprocs \d+ end\s*$', '%pal nprocs 4 end',original,flags=re.M);ck(original.count(' EnGrad')==1,cid+'/original EnGrad inventory');original=original.replace(' EnGrad','',1)
 actual=verify(c['template_input']).read_text();field=' EField '+','.join(format(x,'.17g') for x in c['field_vector_au'])+'\n EFieldOrigin 0.0,0.0,0.0\n';prop='%elprop\n Dipole true\n Origin 0.0,0.0,0.0\nend\n'
 ck(actual.count(field)==1 and actual.count(prop)==1,cid+'/exact field/property directives');ck(actual.replace(field,'',1).replace(prop,'',1)==original,cid+'/sole scientific delta')
 ck('PBE0' in actual and 'D4' in actual and 'RIJCOSX' in actual and 'DefGrid3' in actual and 'VeryTightSCF' in actual and 'NoAutostart' in actual and '* xyzfile 2 1 geometry.xyz' in actual,cid+'/unchanged Hamiltonian and stops')
 ck(prep['records'][cid]['only_declared_delta'] and prep['records'][cid]['new_SP_zero']==(f==0),cid+'/preflight record')
ck(seen=={(m,f) for m in ('La','Dy') for f in fields},'complete exact10 geometry-field product')
for metal,center in p['centers'].items():
 for key in ['zero_EnGrad_input','zero_EnGrad_result','zero_EnGrad_output','zero_EnGrad_property','zero_EnGrad_start','original_panel','original_adapter','original_orca']:verify(center[key])
 for v in center['assets'].values():verify(v)
 context=p['old_EnGrad_context'][metal];ck(context['context_only'] and context['not_SP_zero'] and context['result']==center['zero_EnGrad_result'],'old center contextual only:'+metal)
 ck(center['model_central_field_estimates']==model['summary']['states'][metal+'_central']['steps'],'saved model field values:'+metal)
source=verify(p['adapter']).read_text();fn={n.name:ast.get_source_segment(source,n) for n in ast.parse(source).body if isinstance(n,ast.FunctionDef)}
ck("q=dry.parse(text,c,returncode)" in fn['parse'] and "declared_field_echo" in fn['parse'] and "declared_origins_echo" in fn['parse'],'strict dry native parser plus exact echoes')
ck("origin_Bohr=[0.,0.,0.]" in fn['field_properties'] and "Effective-ECP nuclear origin closure differs" in fn['field_properties'] and "Property/native dipole token closure differs" in fn['field_properties'],'lab zero and all dipole component closures')
ck("st=array_tokens(prop.section(raw,'SCF_Energy'),'totalEnergy',1)[0]" in fn['field_properties'] and "energy_hartree=s" in fn['field_properties'] and 'closure=e-s-d' in fn['field_properties'],'electronic estimator/final electronic+D4 split')
ck("old.estimate(vv,vv['zero'])" in fn['score'] and 'old.consistency(est,.5,odd=True)' in fn['score'] and "{'Dy':1,'La':-1}" in fn['score'] and '),1.)' in fn['score'],'unchanged two-step estimators and6+6+4 gates')
ck("d4checks[m]=all" in fn['score'] and "passed=all(checks.values()) and all(d4checks.values())" in fn['score'],'fixed-geometry D4 check participates in admission')
ck("if r['stage']=='remaining'" in fn['release_contract'] and "s['complete']==5" in fn['release_contract'] and "passed_dry_donor_susceptibility_scout_review" in fn['release_contract'],'actual La5 review precedes Dy5')
ck('claim(root,out,cid' in fn['run'] and 'stop(root,cid,error)' in fn['run'] and "native_program_invocations=attempted" in fn['run'],'immutable finite native attempt/global failure ledger')
ck("b.verify(r['properties']['property_file'])" in fn['admit_result'] and "props==r['properties']" in fn['admit_result'],'saved property hash and full replay')
ck('complete==len(ids)' in fn['collect'] and "response['passed']" in fn['collect'] and "'_unqualified'" in fn['collect'],'missing/failed/consistency failure retained')
wrapper=verify(ready['run_wrapper']).read_text();ck('#SBATCH --nodes=1' in wrapper and '#SBATCH --ntasks=4' in wrapper and '#SBATCH --cpus-per-task=1' in wrapper and not re.search(r'^#SBATCH\s+--(?:mem|time)',wrapper,re.M),'one node/default memory and time/one thread per rank')
ck("math.floor(.75*mem/4)" in fn['run'] and "all(start['threads'][k]=='1'" in fn['admit_result'],'allocation-derived MaxCore and runtime threading')
save(OUT/'CHECKS.json',checks);receipt={'status':'passed_dry_donor_susceptibility_source_review','implementation':ready['implementation'],'panel':ready['panel'],'preparation':ready['preparation'],'ready':pin(DB/'READY.json'),'specification':ready['specification'],'source_inventory':ready['source_inventory'],'configuration':ready['configuration'],'controls':ready['controls'],'model_actual_review':p['model_actual_review'],'review_implementation':pin(__file__),'checks':pin(OUT/'CHECKS.json'),'check_count':len(checks),'new_native_calls':0,'author_test_suite_runs':0,'existing_author_checks_read':37,'declared_cells':10,'new_matching_SP_zeros':2,'stages':[5,5],'response_check_denominator':16,'blockers':[],'limitations':['Native SCF stopping criteria and final-grid/property consistency remain actual admission questions.','Decimal/binary bounds and numerical resolutions are not physical accuracy bounds.','Old EnGrad zeros are contextual only; two new matching SP zeros are required.','Comparison uses global uniform-field projection and cannot uniquely identify local-field mechanism.','Collector cost proxies include admitted rows; complete scheduler costs belong to root.'],'review_cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1}}
save(OUT/'REVIEW.json',receipt);save(DB/'independent_review_v1/REVIEW.json',receipt);print(json.dumps({'review':pin(OUT/'REVIEW.json'),'checks':len(checks)},indent=2))
