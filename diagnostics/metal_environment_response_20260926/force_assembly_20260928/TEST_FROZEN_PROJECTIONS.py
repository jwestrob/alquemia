"""Real saved assembly algebra, not molecular force qualification."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
w=ROOT/'workspaces/metal_environment_response_20260926/force_assembly_frozen_f_v1'
d=json.load(open(w/'PROJECTION_DESIGN.json'))
r=json.load(open(Path(__file__).with_name('FROZEN_F_PROJECTION_RESULT.json')))
def test_real_mode_normalization():
 for c in d['configurations'].values():
  for m in c['modes'].values():
   u=np.array(m['unit_tangent']);t=np.array(m['raw_tangent'])
   assert abs(np.linalg.norm(u)-1)<1e-12
   assert np.max(abs(t-u*m['tangent_norm']))<1e-12

def test_actual_component_assembly_and_translation():
 for p in (w/'assembled').glob('*.json'):
  a=json.load(open(p));e=np.array(a['electronic']['gradient_kcal_mol_A']);cl=np.array(a['classical_gradient_kcal_mol_A']);g=np.array(a['total_gradient_kcal_mol_A'])
  assert np.array_equal(g,e+cl)
  assert np.array_equal(-g,a['total_force_kcal_mol_A'])
  assert np.linalg.norm(a['electronic']['mapping_translation_conservation_residual_kcal_mol_A'])<1e-10

def test_actual_differential_component_accounting():
 for row in r['rows'].values():
  for mode in row['modes'].values():
   c=mode['components'];key='Dy_minus_La_gradient_load_kcal_mol_A'
   assert abs(c['total'][key]-c['electronic'][key]-c['classical'][key])<1e-10
   assert c['total']['Dy_minus_La_force_load_kcal_mol_A']==-c['total'][key]
