"""Collector acceptance/rejection coverage for metal_environment_lady_embedded.collect().

FND-2 (review_01_v1/FINDINGS.json): collect()'s MORead-fallback rejection,
TRAH-iteration-evidence rejection, and point-charge evidence check had zero
fixture tests on either branch. These tests exercise collect() against real
archived accepted Hans8DQ2 MORead/La_A and TRAH-scout/Dy_A artifacts (accept
branch, read-only, in place) and against explicitly-marked corrupted copies
of those same real outputs (reject branch). Corrupted copies retarget only
the task's output_path/engrad_path and recompute cache_key with the real
cache_key() function from affordable_common, so validate()'s identity check
stays genuine; every other field (xyz, input, pointcharges, orbital_seed,
basis/aux assets, implementation pins) still points at the real, untouched
archived files. No Hans8FNR production artifact is read or modified.
"""
import json
import hashlib
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if not (ROOT / 'workspaces').is_dir():
    # workspaces/ (raw products) is gitignored and only materialises in the
    # canonical checkout, never in a fresh worktree; scripts/ there is also
    # the only copy with the untracked production executors (e.g.
    # run_orca_task_manifest.py) that ggr_sensitivity.executed() needs.
    ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/'
                'lanthanide_binding/on_density_scanner/alchemical_bvs')
sys.path.insert(0, str(ROOT / 'scripts'))

from affordable_common import read_json, cache_key  # noqa: E402
from metal_environment_lady_embedded_prepare import METHOD  # noqa: E402
from metal_environment_lady_embedded import collect  # noqa: E402

LA_REAL = ROOT / 'workspaces/metal_environment_response_20260926/lady_repaired_la_seeded_v1'
DY_REAL = ROOT / 'workspaces/metal_environment_response_20260926/lady_repaired_trah_scout_v1'

pytestmark = pytest.mark.skipif(
    not (LA_REAL / 'manifest.json').exists() or not (DY_REAL / 'manifest.json').exists(),
    reason='real archived accepted MORead/La and TRAH/Dy fixtures unavailable',
)


def _retargeted_manifest(real_dir, dest, *, corrupt_output=None):
    """Copy the real manifest with only task output/engrad redirected into dest.

    cache_key is recomputed with the production cache_key() function over the
    modified task dict, so validate()'s identity check is genuinely satisfied
    rather than bypassed. Every other field keeps pointing at the real,
    untouched archived artifact.
    """
    dest.mkdir(parents=True)
    m = read_json(real_dir / 'manifest.json')
    t = dict(m['tasks'][0])
    task_id = t['task_id']

    # run_orca_task_manifest._resolve_under requires input/xyz/output/engrad to
    # resolve under the manifest's own directory, so these four (only) must be
    # physically copied alongside the retargeted manifest. The collector also
    # reads endpoint.runtime.pcgrad beside the output, so copy that real file.
    # Everything else
    # (pointcharges, core/boundary mapping, orbital_seed, initial_gbw, basis
    # assets) is validated by plain verify() against its real recorded path
    # and stays untouched at its real archived location.
    real_output = Path(t['output_path'])
    new_output = dest / real_output.name
    shutil.copyfile(real_output, new_output)
    if corrupt_output is not None:
        new_output.write_text(corrupt_output(new_output.read_text()))
    t['output_path'] = str(new_output)

    if t.get('engrad_path'):
        real_engrad = Path(t['engrad_path'])
        new_engrad = dest / real_engrad.name
        shutil.copyfile(real_engrad, new_engrad)
        t['engrad_path'] = str(new_engrad)

    shutil.copyfile(real_output.with_name('endpoint.runtime.pcgrad'),
                    dest / 'endpoint.runtime.pcgrad')

    real_input = Path(t['input']['path'])
    new_input = dest / real_input.name
    shutil.copyfile(real_input, new_input)
    t['input'] = {'path': str(new_input), 'sha256': t['input']['sha256']}

    real_xyz = Path(t['xyz']['path'])
    new_xyz = dest / real_xyz.name
    shutil.copyfile(real_xyz, new_xyz)
    t['xyz'] = {'path': str(new_xyz), 'sha256': t['xyz']['sha256']}

    bare = {k: v for k, v in t.items() if k != 'cache_key'}
    t['cache_key'] = cache_key(dict(task=bare, protocol=m['protocol_id'], method=METHOD,
                                     inputs=m['inputs'], assets=m['assets'], orca=m['orca'],
                                     implementation=m['implementation']))
    m['tasks'] = [t]

    new_manifest = dest / 'manifest.json'
    new_manifest.write_text(json.dumps(m))

    real_receipt = Path(str(real_output) + '.execution.json')
    receipt = json.loads(real_receipt.read_text())
    receipt['manifest']['sha256'] = hashlib.sha256(new_manifest.read_bytes()).hexdigest()
    receipt['artifacts']['output']['path'] = str(new_output)
    receipt['artifacts']['output']['sha256'] = hashlib.sha256(new_output.read_bytes()).hexdigest()
    receipt['artifacts']['template_input']['path'] = str(new_input)
    receipt['artifacts']['xyz']['path'] = str(new_xyz)
    if t.get('engrad_path'):
        receipt['artifacts']['engrad']['path'] = t['engrad_path']
    (dest / (new_output.name + '.execution.json')).write_text(json.dumps(receipt))

    return new_manifest, task_id


def test_real_accepted_moread_la_endpoint_collects(tmp_path):
    """Real archived Hans8DQ2 La_A MORead endpoint, read in place, is accepted."""
    result = collect(LA_REAL / 'manifest.json')
    row = result['rows']['La_A']
    assert row['status'] == 'complete'
    assert row['orbital_initialization']['target_energy_reused'] is False


def test_real_accepted_trah_dy_endpoint_collects(tmp_path):
    """Real archived Hans8DQ2 Dy_A TRAH-scout endpoint, read in place, is accepted."""
    result = collect(DY_REAL / 'manifest.json')
    row = result['rows']['Dy_A']
    assert row['status'] == 'complete'
    assert row['energy_hartree'] is not None


def test_relocated_moread_la_endpoint_collects(tmp_path):
    """Uncorrupted relocated control must pass before corruptions are meaningful."""
    manifest, task_id = _retargeted_manifest(LA_REAL, tmp_path / 'la_control')
    row = collect(manifest)['rows'][task_id]
    assert row['status'] == 'complete'
    assert row['orbital_initialization']['target_energy_reused'] is False


def test_relocated_trah_dy_endpoint_collects(tmp_path):
    """Relocation retains the actual accepted Dy point-charge gradient."""
    manifest, task_id = _retargeted_manifest(DY_REAL, tmp_path / 'dy_control')
    row = collect(manifest)['rows'][task_id]
    assert row['status'] == 'complete'
    assert row['pointcharge_gradient']['path'] == str(tmp_path / 'dy_control' / 'endpoint.runtime.pcgrad')


def test_silent_pmodel_fallback_rejected(tmp_path):
    """A completed job whose log shows the PModel fallback, not an orbital read, is rejected."""
    manifest, task_id = _retargeted_manifest(
        LA_REAL, tmp_path / 'la_fallback',
        corrupt_output=lambda text: text.replace(
            'Guess MOs are being read from file: initial.gbw',
            'INITIAL GUESS: MODEL POTENTIAL',
        ),
    )
    result = collect(manifest)
    row = result['rows'][task_id]
    assert row['status'] == 'invalid'
    assert 'orbital-read' in row['reason'] or 'fallback' in row['reason']


def test_missing_moread_file_evidence_rejected(tmp_path):
    """A completed job whose log never names initial.gbw at all is rejected."""
    manifest, task_id = _retargeted_manifest(
        LA_REAL, tmp_path / 'la_no_evidence',
        corrupt_output=lambda text: text.replace('initial.gbw', 'initial_renamed.gbw'),
    )
    result = collect(manifest)
    row = result['rows'][task_id]
    assert row['status'] == 'invalid'
    assert 'orbital-read' in row['reason']


def test_missing_trah_evidence_rejected(tmp_path):
    """A completed TRAH-solver job whose log has no TRAH/NR macro-iteration lines is rejected."""
    manifest, task_id = _retargeted_manifest(
        DY_REAL, tmp_path / 'dy_no_trah',
        corrupt_output=lambda text: text.replace('(TRAH MAcro)', '(SOSCF  Iter)').replace(
            '(NR MAcro)', '(SOSCF  Iter)'
        ),
    )
    result = collect(manifest)
    row = result['rows'][task_id]
    assert row['status'] == 'invalid'
    assert 'TRAH' in row['reason']


def test_point_charge_count_mismatch_rejected(tmp_path):
    """A completed job whose two point-charge read confirmations disagree with each other
    and with environment.pc's own declared count is rejected, even though environment.pc
    itself is untouched and still hashes correctly."""
    manifest, task_id = _retargeted_manifest(
        DY_REAL, tmp_path / 'dy_pc_mismatch',
        corrupt_output=lambda text: text.replace(
            'ok (1696 point charges)', 'ok (1695 point charges)', 1,
        ),
    )
    result = collect(manifest)
    row = result['rows'][task_id]
    assert row['status'] == 'invalid'
    assert 'permanent field differs' in row['reason']
