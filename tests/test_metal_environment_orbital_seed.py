"""Real converged Hans fixtures; malformed copies exercise rejection only."""
import copy
import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from metal_environment_orbital_seed import validate_seed, stage_seed
from affordable_common import InvalidArtifact, read_json, record
ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT/'workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2'
MAN = FIX/'manifest.json'
COL = FIX/'FINAL_COLLECTION.json'
pytestmark = pytest.mark.skipif(not COL.exists(), reason='real completed embedded Hans artifacts unavailable')


@pytest.mark.parametrize('metal', ['La', 'Dy'])
def test_real_same_metal_seed_and_copy(metal, tmp_path):
    r = stage_seed(MAN, COL, metal+'_A', MAN, metal+'_B', tmp_path/'seed')
    assert r['source_gbw']['sha256'] == r['staged_gbw']['sha256']
    assert 0.08 < r['coordinate_changes']['max_displacement_A'] < 0.09
    assert not r['field_changes']['changed']
    assert r['target_energy'] is None and r['target_convergence'] is None
    assert r['electronic_state']['physical_multiplicity'] == (1 if metal == 'La' else 6)
    with pytest.raises(FileExistsError):
        stage_seed(MAN, COL, metal+'_A', MAN, metal+'_B', tmp_path/'seed')


def test_cross_metal_rejected():
    with pytest.raises(InvalidArtifact, match='incompatible metal'):
        validate_seed(MAN, COL, 'La_A', MAN, 'Dy_B')


@pytest.mark.parametrize('change', ['state', 'method', 'basis', 'mapping', 'executable'])
def test_corrupted_target_rejected(change, tmp_path):
    m = copy.deepcopy(read_json(MAN)); t = next(t for t in m['tasks'] if t['task_id'] == 'Dy_B')
    if change == 'state': t['electronic_state']['physical_multiplicity'] = 1
    elif change == 'method': m['method'] += ' HF'
    elif change == 'basis': m['assets']['Dy']['basis'] = m['assets']['La']['basis']
    elif change == 'executable': m['orca']['sha256'] = '0'*64
    else:
        mapping = read_json(t['core_mapping']['path']); mapping[1]['source_id'] += '_corrupted'
        p = tmp_path/'bad_map.json'; p.write_text(json.dumps(mapping)); t['core_mapping'] = record(p)
    p = tmp_path/'bad_manifest.json'; p.write_text(json.dumps(m))
    with pytest.raises(InvalidArtifact): validate_seed(MAN, COL, 'Dy_A', p, 'Dy_B')


def test_unaccepted_source_rejected(tmp_path):
    c = read_json(COL); c['rows']['Dy_A']['status'] = 'invalid'
    p = tmp_path/'bad_collection.json'; p.write_text(json.dumps(c))
    with pytest.raises(InvalidArtifact, match='not accepted'):
        validate_seed(MAN, p, 'Dy_A', MAN, 'Dy_B')
