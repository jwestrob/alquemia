"""Serialization checks on real molecular artifacts; no scientific calculations."""
from pathlib import Path
import copy,importlib.util,json
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
MODULE=ROOT/'diagnostics/metal_environment_response_20260926/coupled_scaffold_20260928/result_io.py'
spec=importlib.util.spec_from_file_location('scaffold_result_io',MODULE);io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
RESULT=ROOT/'workspaces/metal_environment_response_20260926/coupled_scaffold_classical_v1/halfstep_receipt_recovery_v2/RESULT.json'
def test_real_result_numpy_roundtrip(tmp_path):
 r=json.loads(RESULT.read_text());expected=copy.deepcopy(r)
 for row in r['results'].values():
  for check in row['checks']:check['halfstep_passed']=np.bool_(check['halfstep_passed']);check['analytic']=np.float64(check['analytic'])
 r['energy_force_queries']=np.int64(r['energy_force_queries']);p=tmp_path/'result.json';io.write_new(p,r);assert json.loads(p.read_text())==expected

def test_real_force_array_checkpoint(tmp_path):
 base=ROOT/'workspaces/metal_environment_response_20260926/scaffold_H_repair_v1/Hans8DQ2';r=json.loads((base/'RECOVERED.json').read_text());f=json.loads((base/'forces.json').read_text());p=tmp_path/'configuration.json';io.write_configuration(p,energies={'preparation_kcal_mol':np.float64(r['actual_preserved_evidence']['after_energy_kcal_mol'])},forces=np.array(f['after_kcal_mol_A']),source='actual_saved_Hans_repair');assert json.loads(p.read_text())['forces']==f['after_kcal_mol_A']

def test_malformed_actual_result_leaves_no_partial(tmp_path):
 r=json.loads(RESULT.read_text());r['results']['La']['checks'][0]['analytic']=np.float64(float('nan'));p=tmp_path/'bad.json'
 with pytest.raises(ValueError):io.write_new(p,r)
 assert not p.exists()

def test_native_format_and_immutable(tmp_path):
 r=json.loads(RESULT.read_text());p=tmp_path/'result.json';io.write_new(p,r);original=p.read_bytes();assert original==(json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
 with pytest.raises(FileExistsError):io.write_new(p,r)
 assert p.read_bytes()==original
