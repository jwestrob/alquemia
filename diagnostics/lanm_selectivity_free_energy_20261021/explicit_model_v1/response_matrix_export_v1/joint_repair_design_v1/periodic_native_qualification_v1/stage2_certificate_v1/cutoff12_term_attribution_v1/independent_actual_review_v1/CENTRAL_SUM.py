from pathlib import Path
import json
from decimal import Decimal,localcontext
T=Path('workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/periodic_native_qualification_v1/stage2_certificate_v1/cutoff12_term_attribution_v1')
r=json.loads((T/'collector_v1/run_v1/execution/RESULT.json').read_text())
with localcontext() as c:
 c.prec=70;s=[Decimal(0)]*3;n=0
 for l in Path(r['rows']['central']['gradient_output']['path']).read_text().splitlines():
  x=l.split()
  if x and x[0]=='TOTAL_GRADIENT':
   for j in range(3):s[j]+=Decimal(x[j+2])
   n+=1
 print('saved_sum',n,list(map(str,s)))
 print('anchor_y_correction',-Decimal('1.8185')*s[1])
