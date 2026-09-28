"""Prepared molecular fixtures and corrupt copies only; no molecular evaluations."""
import copy,json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from metal_environment_lady_embedded_prepare import state,validate,input_text,asset_text,prepare
from metal_environment_lanm_reference import check_config
PREP=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1'

@pytest.mark.parametrize('source,electrons,charge',[('Hans8DQ2',806,-1),('Hans8FNR',816,-1),('Mex8FNS',832,0)])
def test_real_source_states_and_maps(source,electrons,charge):
 c=check_config(read_json(PREP/source/'INPUTS.json'))
 assert c['configurations']['A']['pointcharges']['sha256']==c['configurations']['B']['pointcharges']['sha256']
 for conf in c['configurations'].values():
  for metal,e in conf['endpoints'].items():
   s=state(verify(e['xyz']),metal,e['charge']);assert s['charge']==charge;assert s['explicit_electrons']==electrons;assert s['effective_multiplicity']==1;assert s['physical_multiplicity']==e['multiplicity']

def test_frozen_scout_and_actual_input():
 p=ROOT/'workspaces/metal_environment_response_20260926/lady_frozen_embedded_prepared_v1/Hans8DQ2/manifest.json';m=read_json(p);assert validate(p)['status']=='dry_run_pass'
 for t in m['tasks']:
  text=verify(t['input']).read_text();assert 'DoEQ false' in text and '%pointcharges "environment.pc"' in text;assert 'def2-TZVP' in text and 'gCP' not in text;assert f'* xyzfile {t["charge"]} 1 core.xyz' in text

def test_corrupted_real_spin_and_charge_rejected():
 c=read_json(PREP/'Mex8FNS/INPUTS.json');bad=copy.deepcopy(c);bad['configurations']['A']['endpoints']['Dy']['multiplicity']=1
 with pytest.raises(InvalidArtifact):check_config(bad)
 e=c['configurations']['A']['endpoints']['Dy']
 with pytest.raises(InvalidArtifact):state(verify(e['xyz']),'Dy',-1)

def test_wrong_actual_element_basis_rejected():
 m=read_json(ROOT/'workspaces/metal_environment_response_20260926/lady_frozen_embedded_prepared_v1/Hans8DQ2/manifest.json')
 with pytest.raises(InvalidArtifact):asset_text('La',m['assets']['Dy']['basis'],m['assets']['La']['aux'])

def test_real_origin_subset_and_missing_displacement(tmp_path):
 m=read_json(ROOT/'workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2/manifest.json')
 out=tmp_path/'la_origin'
 result=prepare(verify(m['inputs']),verify(m['agreement']),out,m['assets'],8,1,origin_metals=['La'])
 assert result['status']=='dry_run_pass'
 sub=read_json(out/'manifest.json')
 assert sub['protocol_id'].endswith('_origin_subset_v1')
 assert [t['task_id'] for t in sub['tasks']]==['La_A']
 assert sub['tasks'][0]['xyz']['sha256']==next(t for t in m['tasks'] if t['task_id']=='La_A')['xyz']['sha256']
 # Corrupt a real prepared manifest: a claimed second origin must not pass.
 sub['origin_metals']=['La','Dy'];(out/'manifest.json').write_text(json.dumps(sub))
 with pytest.raises(InvalidArtifact):validate(out/'manifest.json')

def test_subset_rejects_duplicate_metals_before_preparation(tmp_path):
 with pytest.raises(InvalidArtifact):prepare(None,None,tmp_path/'unused',{},8,1,origin_metals=['La','La'])
 with pytest.raises(InvalidArtifact):prepare(None,None,tmp_path/'unused',{},8,2,origin_metals=['La'])

def test_real_two_origin_trah_recovery_preserves_states(tmp_path):
 m=read_json(ROOT/'workspaces/metal_environment_response_20260926/lady_repaired_origins_Hans8FNR_v1/manifest.json')
 out=tmp_path/'trah_origins'
 result=prepare(verify(m['inputs']),verify(m['agreement']),out,m['assets'],8,2,origin_metals=['La','Dy'],origin_solver='TRAH')
 assert result['status']=='dry_run_pass'
 recovery=read_json(out/'manifest.json')
 assert recovery['protocol_id'].endswith('_origin_subset_v1_trah_v1')
 assert recovery['solver']=='TRAH'
 for old,new in zip(m['tasks'],recovery['tasks']):
  assert old['electronic_state']==new['electronic_state']
  assert old['xyz']['sha256']==new['xyz']['sha256']
  assert old['pointcharges']['sha256']==new['pointcharges']['sha256']
  assert ' TRAH\n' in verify(new['input']).read_text()

def test_seeded_repaired_origin_keeps_target_coordinates(tmp_path):
 from metal_environment_orbital_seed import stage_seed
 base=ROOT/'workspaces/metal_environment_response_20260926'
 source=base/'lady_frozen_embedded_hans_v2/manifest.json'
 target=base/'lady_repaired_la_origin_prepared_v1/manifest.json'
 stage_seed(source,source.parent/'FINAL_COLLECTION.json','La_A',target,'La_A',tmp_path/'seed')
 m=read_json(target);out=tmp_path/'seeded'
 result=prepare(verify(m['inputs']),verify(m['agreement']),out,m['assets'],8,1,origin_metals=['La'],seed_records={'La_A':tmp_path/'seed/SEED.json'})
 assert result['status']=='dry_run_pass'
 new=read_json(out/'manifest.json');t=new['tasks'][0]
 assert new['protocol_id'].endswith('_moread_v1')
 assert t['xyz']['sha256']==m['tasks'][0]['xyz']['sha256']
 assert 'Guess MORead' in verify(t['input']).read_text()
 assert 'Guess PModel' not in verify(t['input']).read_text()
 assert t['cache_key']!=m['tasks'][0]['cache_key']
 # A corrupted real seed target cannot silently seed another geometry.
 s=read_json(tmp_path/'seed/SEED.json');s['target_task_id']='Dy_A'
 (tmp_path/'seed/SEED.json').write_text(json.dumps(s))
 with pytest.raises(InvalidArtifact):validate(out/'manifest.json')
