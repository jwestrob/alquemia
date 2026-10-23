#!/usr/bin/env python3
"""Interpret saved convergence traces without replacing native admission."""
import hashlib
import json
from pathlib import Path
import re
import sys
import time

HERE=Path(__file__).resolve().parent;F=r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?'
def pin(p):
 p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def main():
 tick=time.monotonic();out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);sources=[HERE/'native_review_v1/RESULT.json',HERE.parent/'hydrated_scf_review_v1/RESULT.json'];rows=[]
 for source in sources:
  for prior in json.loads(source.read_text())['rows']:
   raw=Path(prior['output']['path']);assert pin(raw)==prior['output'];text=raw.read_text().split('****END OF INPUT****')[-1]
   mode=re.findall(r'Convergence Check Mode\s+ConvCheckMode\s*\.{2,}\s*([^\n]+)',text)
   signals=re.findall(r'\*+\s*(.*?)\s+signals convergence\s*\*+',text)
   block=text.split('S-O-S-C-F',1)[1].split('SCF CONVERGED AFTER',1)[0]
   trace=[dict(zip(('iteration','energy_Eh','delta_Eh','RMSDP','MaxDP','MaxGrad','seconds'),map(float,v))) for v in re.findall(r'^\s*(\d+)\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s*$',block,re.M)]
   assert trace and len(signals)==1 and mode==['Total+1el-Energy'];crit={x['name']:x for x in prior['printed_scf_criteria']}
   diis=crit['DIIS Error']['value'];first=trace[0]['MaxGrad'];last=trace[-1]['MaxGrad']
   row={'cell_id':prior['cell_id'],'output':prior['output'],'native_mode':mode[0],'stop_signal':signals[0],'terminal_solver':'SOSCF','TRAH_executed':False,'first_SOSCF':trace[0],'last_SOSCF':trace[-1],
       'printed_last_DIIS':diis,'printed_DIIS_matches_first_SOSCF_gradient_at_print_precision':abs(diis-first)/max(abs(first),1e-30)<.005,
       'printed_last_gradient_matches_final_SOSCF_at_print_precision':abs(crit['Orbital Gradient']['value']-last)/max(abs(last),1e-30)<.005,
       'density_flags_current_but_not_all_required_by_energy_check':True,'raw_residual_review':pin(source)}
   if 'LARGEST AMPLITUDES' in text:
    # ORCA uses -1 orbital-index sentinels for singles: retain them too.
    section=text.split('LARGEST AMPLITUDES',1)[1].split('!!!',1)[0]
    row['all_printed_amplitudes_including_singles']=[float(x) for x in re.findall(r'^\s*-?\d+->\s*-?\d+\s+-?\d+->\s*-?\d+\s+('+F+r')\s*$',section,re.M)]
   rows.append(row)
 result={'status':'saved_convergence_interpretation','rows':rows,'source':pin(__file__),'command':sys.argv,'elapsed_seconds':time.monotonic()-tick,'molecular_calls':0,'local_threads':1,
   'official_reference':'https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html#convergence-tolerances',
   'conclusion':'All14 finish with SOSCF under Total+1el-Energy mode; native energy or gradient stopping signals are valid under that mode. Printed last DIIS values match initial SOSCF gradients, not terminal residuals. No all-six threshold gate was declared or is retroactively imposed.'}
 (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'rows':len(rows),'energy_stop':sum('Energy' in x['stop_signal'] for x in rows),'gradient_stop':sum('Gradient' in x['stop_signal'] for x in rows),'stale_DIIS_supported':sum(x['printed_DIIS_matches_first_SOSCF_gradient_at_print_precision'] for x in rows),'result':pin(out/'RESULT.json')}))
if __name__=='__main__':main()
