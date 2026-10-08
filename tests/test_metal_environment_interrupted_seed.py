"""Real shutdown artifacts; changed copies below are malformed-input tests only."""
import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import InvalidArtifact, read_json, record
from metal_environment_interrupted_seed import validate_interrupted, stage
from metal_environment_orbital_seed import validate_seed

ROOT = Path(__file__).resolve().parents[1]
W = ROOT/'workspaces/metal_environment_response_20260926/lady_repaired_origins_Hans8FNR_trah_v1'
D = ROOT/'diagnostics/metal_environment_response_20260926/recovery_20261007'
MAN, EV, INV = W/'manifest.json', D/'RESULT.json', D/'ORBITAL_INVENTORY.json'
pytestmark = pytest.mark.skipif(not INV.exists(), reason='real interrupted checkpoint inventory unavailable')


@pytest.mark.parametrize('metal', ['La', 'Dy'])
def test_real_interrupted_checkpoint(metal, tmp_path):
    r = stage(MAN, MAN, metal+'_A', EV, INV, tmp_path/'seed')
    assert r['staged_gbw']['sha256'] == r['source_gbw']['sha256']
    assert r['target_energy'] is None and r['target_convergence'] is None
    assert not r['binary_integrity_qualified']
    with pytest.raises(FileExistsError):
        stage(MAN, MAN, metal+'_A', EV, INV, tmp_path/'seed')


@pytest.mark.parametrize('change', ['state', 'field', 'basis', 'method', 'solver'])
def test_changed_target_rejected(change, tmp_path):
    m = read_json(MAN); t = m['tasks'][0]
    if change == 'state': t['electronic_state']['physical_multiplicity'] = 99
    elif change == 'field': t['pointcharges']['sha256'] = '0'*64
    elif change == 'basis': m['assets']['La']['basis'] = m['assets']['Dy']['basis']
    elif change == 'method': m['method'] += ' HF'
    else: m['solver'] = 'default'
    p = tmp_path/'malformed.json'; p.write_text(json.dumps(m))
    with pytest.raises(InvalidArtifact): validate_interrupted(MAN, p, 'La_A', EV, INV)


def test_existing_accepted_seed_gate_stays_strict():
    with pytest.raises(InvalidArtifact, match='not accepted'):
        validate_seed(MAN, W/'COLLECTION_1220342.json', 'La_A', MAN, 'La_A')


def test_corrupted_checkpoint_pin_rejected(tmp_path):
    i = read_json(INV); i['files']['La']['sha256'] = '0'*64
    p = tmp_path/'bad_inventory.json'; p.write_text(json.dumps(i))
    with pytest.raises(InvalidArtifact): validate_interrupted(MAN, MAN, 'La_A', EV, p)


PREP = ROOT/'workspaces/metal_environment_response_20260926/lady_FNR_interrupted_continuation_20261007_v1/manifest.json'


@pytest.mark.skipif(not PREP.exists(), reason='real interrupted continuation preparation unavailable')
def test_prepared_continuation_and_unlabelled_seed_rejected(tmp_path):
    from metal_environment_lady_embedded_prepare import validate
    assert validate(PREP)['status'] == 'dry_run_pass'
    m = read_json(PREP)
    assert m['initial_guess'] == 'MORead' and m['seed_kind'] == 'interrupted'
    assert m['protocol_id'].endswith('_interrupted_v1')
    del m['seed_kind']
    m['protocol_id'] = m['protocol_id'].removesuffix('_interrupted_v1')
    p = tmp_path/'unlabelled.json'; p.write_text(json.dumps(m))
    with pytest.raises(InvalidArtifact, match='seed kind declaration differs'):
        validate(p)
