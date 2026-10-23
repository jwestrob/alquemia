#!/usr/bin/env python3
"""One metadata inventory and saved zero-step check; no native/solver/Jacobian build."""
from pathlib import Path
import json,hashlib,math,time,os,socket
T=time.monotonic();C=time.process_time();HERE=Path(__file__).resolve().parent
ROOT=next(x for x in HERE.parents if (x/'workspaces/lanm_selectivity_free_energy_20261021').exists())
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021';W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021';RD='explicit_model_v1/radial_drive_shape_v1';F='explicit_model_v1/joint_component_feasibility_v1/component_fit_v1';J='explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1'
def pin(p):return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):return json.loads(p.read_text())
def verify(p):
 x=Path(p['path']);assert pin(x)==p;return read(x)
assert not (HERE/'INVENTORY.json').exists()
rootdec=read(HERE/'ROOT_DESIGN_DECLARATION.json');admit=verify(rootdec['terminal_trimer_admission']);tr=verify(admit['result']);candidate=verify(admit['candidate'])
fitpath=W/RD/'iterative_fit_v1/run_v1/execution/RESULT.json';fit=read(fitpath);step=fit['iterations'][-1];proposal=verify(step['proposal']);p=verify(proposal['linearization']);layout=verify(proposal['row_layout'])
assert step['outer_step']==proposal['outer_step']==8 and step['origin_key']==proposal['origin_key']==admit['selected_key']==rootdec['fixed_candidate_key']
assert step['derivative_keys']==proposal['derivative_keys'] and len(set(step['derivative_keys']))==28
assert p['center']==candidate['parameters'] and p['raw_center_q']==fit['selected_metrics']['q']==0.13804357390849725
trinv=verify(tr['inventory']);cy={k:v for k,v in trinv['cycles'].items() if k.startswith('aquo9_') and k.endswith('__far')};assert len(cy)==2
needed=sorted({k for v in cy.values() for k in v['model_terms']});dynamic=sorted(set(needed)&set(candidate['cells']));static=sorted(set(needed)&set(trinv['static_terms']))
assert len(needed)==14 and len(dynamic)==6 and len(static)==8 and set(dynamic+static)==set(needed)
rows=[]
for i,k in enumerate(step['derivative_keys']):
 f=fitpath.parent/'vectors'/(k+'.json');v=read(f)
 assert v['status']=='complete' and len(v['parameters'])==14 and all(v['cells'][c]['status']=='complete' for c in dynamic)
 rows.append({'index':i,'coordinate':i//2,'side':'minus' if i%2==0 else 'plus','key':k,'artifact':pin(f),'parameters':v['parameters'],'required_dynamic_fields':{c:{'json_path':['cells',c],'energy_field':'energy','moment_field':'mu_Debye','status':v['cells'][c]['status']} for c in dynamic},'raw_native_replay_performed':False})
for j,span in enumerate(p['spans']):assert (rows[2*j+1]['parameters'][j]-rows[2*j]['parameters'][j])/p['width'][j]==span
new='frozen_radial14_terminal5d992';groups=[]
for label,x in [('La',tr['comparison']['cycles']['aquo9_central_La__far']),('Dy',tr['comparison']['cycles']['aquo9_central_Dy__far']),('Dy_minus_La',tr['comparison']['paired_Dy_minus_La']['aquo9__far'])]:
 m=x['models'][new]
 groups.append({'name':'E3/'+label,'dimension':1,'r':[m['energy_error_kcal_mol']],'stored_scale':abs(m['energy_error_kcal_mol']),'units':'kcal/mol'})
 groups.append({'name':'mu3/'+label,'dimension':3,'r':m['dipole_error_Debye'],'stored_scale':m['dipole_error_norm_Debye'],'units':'Debye'})
assert len(groups)==6 and all(g['stored_scale']>0 for g in groups)
# Only the zero witness from saved center fields; no probe scientific values used.
q0=p['raw_center_q'];zero={};ratios={}
for name,b in p['blocks'].items():
 if name=='old_EF':continue
 ratio=math.fsum(float(v)*float(v) for v in b['r'])/b['denominator'];ratios[name]=ratio
 zero['block/'+name]=(q0 if name in ('P_energy','B_energy') else 1)-ratio
for name,v in zip(p['energy']['names'],p['energy']['r']):zero['energy/'+name]=1-abs(v)
for name,b in p['forces'].items():
 zero['force/RMS/'+name]=4-math.fsum(float(w)*float(v)*float(v) for w,v in zip(b['weights'],b['r']))
 zero['force/max/'+name]=5-max(abs(v) for v in b['r'])
zero['signed_W9']=p['signed_W9_z']['absolute_ceiling_Debye']-abs(p['signed_W9_z']['r_Debye'])
for j,(z,lo,hi) in enumerate(zip(p['z'],p['bounds_lower'],p['bounds_upper'])):
 zero['box/lower/'+str(j)]=z;zero['box/upper/'+str(j)]=1-z
 zero['trust/lower/'+str(j)]=.025;zero['trust/upper/'+str(j)]=.025
for g in groups:
 g['reconstructed_norm']=math.sqrt(math.fsum(v*v for v in g['r']));g['stored_scale_minus_reconstructed_norm']=g['stored_scale']-g['reconstructed_norm']
 g['normalized_zero_margin']=1-math.sqrt(math.fsum((v/g['reconstructed_norm'])**2 for v in g['r']))
zero_report={'scope':'finite zero-step scalar check only, not constructed conic/Jacobian or solver feasibility admission','d':[0.]*14,'t':0.,'literal_q_ceiling':q0,'old_auxiliary_q_not_used':p['q_current'],'ratios':ratios,'q_reconstructed':max(ratios[k] for k in ('P_energy','B_energy')),'margins':zero,'negative_margins':{k:v for k,v in zero.items() if v<0},'new_groups':groups,'sqrt_q':math.sqrt(q0),'sqrt_q_squared_minus_literal_q':math.sqrt(q0)**2-q0,'padding_added':False,'future_actual_conic_zero_check_required':True}
config=read(D/RD/'iterative_fit_v1/CONFIG_PREPARED.json');localinputs=verify(config['inputs'])
sources={k:config[k] for k in ('core','implementation','math','solver','settings','environment','targets','reuse_inventory','parent_evaluator')}
sources['attribution_math']=localinputs['sources']['attribution_math'];sources['parent_core']=localinputs['sources']['core'];sources['trimer_cycle_producer']=tr['implementation']
sources['original_cycle_helper']=pin(D/F/'physical_gate_component_fit_v1/postselection_regression_v1/RUN.py')
sources['nonduplicate_inventory']=pin(D/'explicit_model_v1/joint_component_feasibility_v1/inventory_review/REPORT.txt')
for x in sources.values():assert pin(Path(x['path']))==x
obj={'status':'source_design_inventory_no_solve','declaration':pin(HERE/'ROOT_DESIGN_DECLARATION.json'),'design':pin(HERE/'DESIGN.txt'),'inventory_implementation':pin(Path(__file__)),'trimer_admission':rootdec['terminal_trimer_admission'],'trimer_result':admit['result'],'trimer_inventory':tr['inventory'],'candidate':admit['candidate'],'fit_result':pin(fitpath),'fit_root_admission':pin(W/RD/'iterative_fit_v1/final_comparison_v1/actual_v1/ROOT_ACTUAL_REVIEW.json'),'proposal':step['proposal'],'linearization':proposal['linearization'],'conic_data':proposal['conic_data'],'row_layout':proposal['row_layout'],'origin_key':step['origin_key'],'ordered_derivative_vectors':rows,'actual_spans':p['spans'],'coordinate_lower':p['lower'],'coordinate_upper':p['upper'],'cycle_dictionaries':cy,'required_dynamic_cell_ids':dynamic,'static_terms':{k:trinv['static_terms'][k] for k in static},'groups':groups,'literal_q_ceiling':q0,'sources':sources,'old_SOC_names':[x['name'] for x in layout if x['kind']=='SOC'],'new_six_groups_absent_from_old_objective_and_constraints':True,'zero_step':zero_report,'new_native_calls':0,'new_probes':0,'new_reference_calls':0,'solver_calls':0,'new_Jacobians_constructed':0,'raw46_replays':0,'execution_release':None,'inherited_costs_not_recharged':{'full_fit':fit['cost'],'step8_proposal':proposal['new_solver_cost'],'trimer':tr['cost']}}
(HERE/'INVENTORY.json').write_text(json.dumps(obj,indent=2)+'\n')
receipt={'status':'metadata_inventory_and_zero_scalar_check_complete','inventory':pin(HERE/'INVENTORY.json'),'implementation':pin(Path(__file__)),'cost':{'wall_seconds':time.monotonic()-T,'CPU_seconds':time.process_time()-C,'hostname':socket.gethostname(),'PID':os.getpid(),'job_id':None,'threads':1},'negative_inherited_zero_margins':zero_report['negative_margins'],'negative_new_group_normalized_zero_margins':{g['name']:g['normalized_zero_margin'] for g in groups if g['normalized_zero_margin']<0},'new_native_calls':0,'solver_calls':0,'probe_Jacobians_constructed':0,'raw46_replays':0}
(HERE/'INVENTORY_RESULT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
