#!/usr/bin/env python3
"""One saved-scalar extraction; no raw native parsing or scientific re-admission."""
import hashlib,json,os,time,socket
from pathlib import Path
start=time.monotonic(); cpu=time.process_time(); here=Path(__file__).resolve().parent
assert not (here/'SAVED_VALUES.json').exists()
def pin(p):return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
dec=json.loads((here/'DECLARATION.json').read_text()); p=Path(dec['inputs']['result']['path'])
assert pin(p)==dec['inputs']['result']
r=json.loads(p.read_text()); assert r['status']=='radial14_consumed_trimer_regression_complete' and r['error'] is None
c=r['comparison']; assert len(c['cycles'])==8 and len(c['paired_Dy_minus_La'])==4
base='original_pre_response_baseline'; old='feca12_physical_gate_terminal'; new='frozen_radial14_terminal5d992'
paired=c['paired_Dy_minus_La']; cycles=c['cycles']
summary={
 'individual_E3_absolute_error_worse_than_baseline':sum(abs(x['models'][new]['energy_error_kcal_mol'])>abs(x['models'][base]['energy_error_kcal_mol']) for x in cycles.values()),
 'individual_E3_absolute_error_worse_than_feca12':sum(abs(x['models'][new]['energy_error_kcal_mol'])>abs(x['models'][old]['energy_error_kcal_mol']) for x in cycles.values()),
 'individual_mu3_error_norm_worse_than_feca12':sum(x['models'][new]['dipole_error_norm_Debye']>x['models'][old]['dipole_error_norm_Debye'] for x in cycles.values()),
 'individual_mu3_error_norm_better_than_baseline':sum(x['models'][new]['dipole_error_norm_Debye']<x['models'][base]['dipole_error_norm_Debye'] for x in cycles.values()),
 'paired_E3_absolute_error_worse_than_baseline':sum(abs(x['models'][new]['energy_error_kcal_mol'])>abs(x['models'][base]['energy_error_kcal_mol']) for x in paired.values()),
 'paired_E3_absolute_error_worse_than_feca12':sum(abs(x['models'][new]['energy_error_kcal_mol'])>abs(x['models'][old]['energy_error_kcal_mol']) for x in paired.values()),
 'paired_mu3_error_norm_better_than_baseline':sum(x['models'][new]['dipole_error_norm_Debye']<x['models'][base]['dipole_error_norm_Debye'] for x in paired.values()),
 'paired_mu3_error_norm_better_than_feca12':sum(x['models'][new]['dipole_error_norm_Debye']<x['models'][old]['dipole_error_norm_Debye'] for x in paired.values()),
 'paired_E3_model_negative_reference_positive':sum(x['models'][new]['Dy_minus_La_E3_kcal_mol']<0<x['reference_electronic_Dy_minus_La_E3_kcal_mol'] for x in paired.values())}
obj={'status':'saved_scalar_interpretation_only_not_numerical_admission','declaration':pin(here/'DECLARATION.json'),'implementation':pin(Path(__file__).resolve()),'result':pin(p),'comparison':c,'summary':summary,'terminal_metrics':r['terminal_metrics'],'reference_limits':r['reference']['historical_reference']['group_metadata'],'native_result_counts':r['actual_counts'],'native_result_cost':r['cost'],'new_native_calls':0,'raw_outputs_reparsed':0,'new_solver_calls':0}
(here/'SAVED_VALUES.json').write_text(json.dumps(obj,indent=2)+'\n')
lines=['EXACT SAVED COMPARISON DISPLAY','', 'Stage labels: baseline = original_pre_response_baseline; feca12 = feca12_physical_gate_terminal; radial14 = frozen_radial14_terminal5d992.','Values below only format already reported results. Full precision remains in SAVED_VALUES.json and source RESULT.','', 'PAIRED E3: Dy-minus-La, kcal/mol','case            reference   baseline     feca12    radial14   residual(base/feca/radial)']
for k,x in paired.items():
 ms=x['models']; vals=[x['reference_electronic_Dy_minus_La_E3_kcal_mol']]+[ms[m]['Dy_minus_La_E3_kcal_mol'] for m in (base,old,new)]
 lines.append(f'{k:15s} '+ ' '.join(f'{v:10.6f}' for v in vals)+'   '+ '/'.join(f'{ms[m]["energy_error_kcal_mol"]:+.6f}' for m in (base,old,new)))
lines += ['', 'PAIRED SIGNED mu3: Dy-minus-La vector (x,y,z), Debye. Norms below are error norms, not substitutes for signed vectors.']
for k,x in paired.items():
 lines.append(k)
 lines.append(' reference '+str(x['reference_Dy_minus_La_mu3_Debye']))
 for short,m in [('baseline',base),('feca12',old),('radial14',new)]:
  v=x['models'][m];lines.append(f' {short:8s} '+str(v['Dy_minus_La_mu3_Debye'])+f'  error norm {v["dipole_error_norm_Debye"]:.9f}')
lines += ['', 'INDIVIDUAL CYCLE ERROR: signed E3 residual kcal/mol and signed-vector mu3 error norm Debye','case                              E3 baseline/feca12/radial14             mu3 norm baseline/feca12/radial14']
for k,x in cycles.items():
 ms=x['models']; lines.append(f'{k:33s} '+ '/'.join(f'{ms[m]["energy_error_kcal_mol"]:+.6f}' for m in (base,old,new))+'       '+ '/'.join(f'{ms[m]["dipole_error_norm_Debye"]:.6f}' for m in (base,old,new)))
lines += ['', 'COUNTS (descriptive comparisons; no new acceptance threshold)']+[f'{k}: {v}' for k,v in summary.items()]
(here/'TABLES.txt').write_text('\n'.join(lines)+'\n')
receipt={'status':'completed_saved_scalar_extraction','implementation':pin(Path(__file__).resolve()),'declaration':pin(here/'DECLARATION.json'),'result':pin(p),'saved_values':pin(here/'SAVED_VALUES.json'),'tables':pin(here/'TABLES.txt'),'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'hostname':socket.gethostname(),'PID':os.getpid(),'job_id':None,'numerical_threads':1},'native_calls':0,'solver_calls':0,'raw_native_reparses':0}
(here/'EXTRACTION.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'summary':summary,'cost':receipt['cost'],'receipt':pin(here/'EXTRACTION.json')}))
