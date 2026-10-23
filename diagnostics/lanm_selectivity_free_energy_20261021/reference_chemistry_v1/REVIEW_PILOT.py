#!/usr/bin/env python3
"""Read-only native-output supplement; no molecular calls or admission changes."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
TASKS=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/pilot_v1/tasks'
PATTERN=re.compile(r'Last\s+(.+?)\s+\.{3}\s+([-+\d.eE]+)\s+Tolerance\s*:\s*([-+\d.eE]+)')

def pin(path):
    path=Path(path)
    return {'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);rows=[];pending=[]
    for directory in sorted(TASKS.iterdir()):
        output=directory/'execution/endpoint.out';result=directory/'execution/RESULT.json'
        if not result.exists():
            pending.append({'cell_id':directory.name,'native_result_present':False});continue
        r=json.loads(result.read_text());text=output.read_text();numeric=[]
        for name,value,tolerance in PATTERN.findall(text):
            numeric.append({'name':name.strip(),'value':float(value),'tolerance':float(tolerance),'passes_printed_tolerance':abs(float(value))<=float(tolerance)})
        cycles=re.findall(r'SCF CONVERGED AFTER\s+(\d+)\s+CYCLES',text)
        signals=re.findall(r'\*+\s*(.*?)\s+(?:signals|signifies) convergence\s*\*+',text,re.I)
        row={'cell_id':directory.name,'native_result':pin(result),'output':pin(output),'native_status':r['status'],
             'native_checks':r['checks'],'all_native_checks_pass':all(r['checks'].values()),
             'numeric_scf_criteria':numeric,'six_criteria_extracted':len(numeric)==6,
             'all_printed_residuals_pass':len(numeric)==6 and all(x['passes_printed_tolerance'] for x in numeric),
             'scf_cycles':int(cycles[-1]) if cycles else None,'termination_signals':signals,
             'spin_squared_values':r['spin_squared_values'],'energy_hartree':r['energy_hartree'],
             'elapsed_seconds':r['elapsed_seconds'],'allocated_cpu_seconds':r['allocated_cpu_seconds'],
             'physical_state_qualified':r['physical_state_qualified'],'method_accuracy_qualified':r['method_accuracy_qualified']}
        rows.append(row)
    result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'saved_output_review',
            'rows':rows,'not_yet_terminal':pending,'molecular_calls':0,'local_threads':1,'source':pin(__file__),
            'command':sys.argv,'interpretation':'Native protocol admission and strict printed residual compliance are distinct. No native gate changed; neither establishes physical-state or method accuracy.'}
    (a.out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Native pilot: supplemental SCF review','',
        'Saved output only. The frozen native adapter and its original admission remain unchanged. Numeric values below supplement the adapter\'s empty YES/NO residual field.','',
        '| Cell | Native checks | Energy | Max density | RMS density | DIIS | Orbital gradient | Orbital rotation |',
        '|---|---|---|---|---|---|---|---|']
    for r in rows:
        values=[f"{x['value']:.4g} ({'pass' if x['passes_printed_tolerance'] else 'fail'})" for x in r['numeric_scf_criteria']]
        lines.append('| '+r['cell_id']+' | '+str(r['all_native_checks_pass'])+' | '+' | '.join(values)+' |')
    lines+=['','Tolerances in order: 1e-9, 1e-8, 1e-9, 1e-8, 2e-6, 2e-6 (native atomic units).',
        '',f"Completed {len(rows)}/8; no terminal result yet for "+', '.join(x['cell_id'] for x in pending)+'.',
        '', 'The four frozen-f cells admit the remaining frozen-f panel as diagnostic coverage under the existing protocol. They do not qualify strict SCF/force accuracy. Explicit-f Dy admission remains separate; do not extrapolate the La results to Dy.',
        '', 'No physical-state, spin-orbit, basis-limit, correlated-reference, derivative-consistency or experimental selectivity qualification follows from this review.']
    (a.out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'reviewed':len(rows),'pending':len(pending),'native_checks_pass':[r['cell_id'] for r in rows if r['all_native_checks_pass']]}))

if __name__=='__main__':main()
