"""Read only the five declared energy-settling failures; no new SCF calls."""
import json,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,HA_TO_KCAL
from compact_solvation import diagnostics
from structure_informed_starts import scf_details

W=ROOT/'workspaces/precision_pool_continuation_20260923/run_v1'
c=read_json(W/'COMPARISON.json');inv=read_json(verify(c['inventory']));source={tuple(r[k] for k in ('case_id','candidate','metal','medium')):r for r in inv['rows']}
stages=[read_json(verify(c['stage'+str(i)])) for i in (1,2)]
look=[{tuple(r[k] for k in ('case_id','candidate','metal','medium')):r for r in s['rows']} for s in stages]
rows=[]
for cell in c['cells']:
    if cell['energy_settling_pass']:continue
    key=tuple(cell[k] for k in ('case_id','candidate','metal','medium'));src=source[key];steps=[];charges=[]
    for stage in (0,1,2):
        actual=src['source_cell'] if stage==0 else look[stage-1][key]['actual']
        text=verify(actual['output']).read_text();detail=scf_details(text)
        audit=diagnostics(actual,src) if stage==0 else look[stage-1][key]['audit']
        block=text.split('Iteration    Energy (Eh)')[-1].split('**** Energy Check signals convergence ****')[0]
        parsed=[]
        for line in block.splitlines():
            fields=line.split()
            if len(fields)==7 and fields[0].isdigit():
                try:parsed.append({'cycle':int(fields[0]),'energy_hartree':float(fields[1]),'delta_E_printed':float(fields[2]),'RMSDP':float(fields[3]),'MaxDP':float(fields[4]),'Omega':float(fields[5]),'time_seconds':float(fields[6])})
                except ValueError:pass
        if not parsed or parsed[-1]['cycle']!=detail['cycles']:raise RuntimeError('SCF iteration trace incomplete')
        q=np.asarray(audit['mulliken_charges_e']);charges.append(q)
        steps.append({'stage':stage,'actual':actual,'energy_hartree':cell[('cold' if stage==0 else 'stage'+str(stage))+'_hartree'],
            'cycles':detail['cycles'],'energy_check_signals_convergence':'Energy Check signals convergence' in text,
            'reported_check_mode':next((line.strip() for line in text.splitlines() if 'ConvCheckMode' in line),None),
            'metal_charge_e':float(q[0]),'fractional_orbitals_printed':detail['printed_fractional_orbital_count'],
            'convergence':{k:detail[k] for k in ('energy','max_density','rms_density')},'trace':parsed,
            'last5_energy_range_kcal_mol':(max(r['energy_hartree'] for r in parsed[-5:])-min(r['energy_hartree'] for r in parsed[-5:]))*HA_TO_KCAL})
    changes=[]
    for i,j in ((0,1),(1,2),(0,2)):
        d=charges[j]-charges[i]
        changes.append({'from':i,'to':j,'energy_change_kcal_mol':(steps[j]['energy_hartree']-steps[i]['energy_hartree'])*HA_TO_KCAL,
            'metal_charge_change_e':float(d[0]),'maximum_atomic_charge_change_e':float(np.max(np.abs(d))),
            'charge_L2_difference_e':float(np.linalg.norm(d))})
    rows.append({k:cell[k] for k in ('case_id','candidate','metal','medium')}|{'steps':steps,'changes':changes})
result={'comparison':record(W/'COMPARISON.json'),'selection':'exact five cells failing the frozen pass1_to_pass2_0.1_kcal_gate','rows':rows,'denominator':5,'new_molecular_calls':0,'implementation':record(__file__)}
write_new(W/'FAILING_TRACE_AUDIT_v1.json',result)
for r in rows:
    print(r['case_id'],r['candidate'],r['medium'])
    print('cycles,metalq,MAXDP,last5range:',[(x['cycles'],x['metal_charge_e'],x['convergence']['max_density']['residual'],x['last5_energy_range_kcal_mol']) for x in r['steps']])
    print('changes:',r['changes'])
