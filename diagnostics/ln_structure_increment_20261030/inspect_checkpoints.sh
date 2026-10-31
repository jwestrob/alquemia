#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
workspaces/ln_structure_increment_20261030/chemprop_env/bin/python - <<'PY'
from pathlib import Path
import torch,json
w=Path('workspaces/ln_structure_increment_20261030/upstream/models')
out={}
for model in ['stability_constant_model_M1','stability_constan_protonation_constant_model_M4']:
 p=w/model/'best_model/fold_0/model_0/model.pt'
 s=torch.load(p,map_location='cpu')
 args=vars(s['args'])
 out[model]={'args':args,'data_scaler':s.get('data_scaler'),'features_scaler':s.get('features_scaler')}
Path('workspaces/ln_structure_increment_20261030/checkpoint_metadata.json').write_text(json.dumps(out,indent=2,default=str))
for k,v in out.items(): print(k, {a:v['args'].get(a) for a in ('hidden_size','depth','features_size','features_path','epochs','features_generator','ffn_hidden_size','ffn_num_layers','smiles_columns')})
PY
