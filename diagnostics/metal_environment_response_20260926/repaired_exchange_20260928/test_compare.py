"""Actual completed Hans fixtures and explicitly missing/corrupted declarations."""
import importlib.util
import json
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('repaired_exchange_compare', HERE/'COMPARE.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ROOT = HERE.parents[2]


def declaration(tmp_path, missing=False):
    w = ROOT/'workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2'
    c = dict(sources=['Hans8DQ2'], hans_sources=[], mex_source=None, limitations=['Real archived parser fixture'], cells=[])
    for metal in ('La', 'Dy'):
        c['cells'].append(dict(source='Hans8DQ2', source_id='Hans8DQ2', metal=metal,
                               task_id=metal+'_A', manifest=str(w/'manifest.json'),
                               collection=str(tmp_path/'absent.json' if missing else w/'FINAL_COLLECTION.json')))
    path=tmp_path/'config.json';path.write_text(json.dumps(c));return path


def test_actual_archived_sign_and_units(tmp_path):
    result=module.compare(declaration(tmp_path))
    assert result['complete_cells']==result['declared_cells']==2
    expected=(-5732.680720450197-(-5727.211985065892))*module.HA_TO_KCAL
    assert result['raw_Dy_minus_La_kcal_mol']['Hans8DQ2']==pytest.approx(expected,abs=1e-10)
    assert result['affinity'] is None and result['classification'] is None


def test_missing_is_unavailable_not_zero(tmp_path):
    result=module.compare(declaration(tmp_path,missing=True))
    assert result['complete_cells']==0 and result['declared_cells']==2
    assert result['raw_Dy_minus_La_kcal_mol']['Hans8DQ2'] is None


def test_duplicate_real_cell_rejected(tmp_path):
    path=declaration(tmp_path);c=json.loads(path.read_text());c['cells'].append(c['cells'][0]);path.write_text(json.dumps(c))
    with pytest.raises(module.InvalidArtifact):module.compare(path)
