#!/usr/bin/env python3
"""Future Slurm-only saved-artifact controls for the frozen continuation draft.

No configuration, continuation decision, or admitted molecular observation is
created here. A real bound configuration must first pass the unchanged complete
inputs() replay. No native runner, model evaluator, or submission API is called.
"""
import argparse
import copy
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import sys
import time
import traceback

PROTOCOL = 'lbt_whole256_continuation_draft_negative_controls_v1'
ADAPTER_SHA = '06e40d37a44c16446a062bbc479b6a49db70c72e6d16057a2fd7f7278916ee0f'
SPEC_SHA = 'e23187e044eacef6f0d02e355c0189dfad2638107be97b4438a0068a5f56c0c1'
CASES = (
    ('missing_scout', 'Required evidence pin missing'),
    ('changed_scout_identity', 'Admitted actual13051 required'),
    ('altered_decision_scope', 'Continuation decision scope/source differs'),
    ('altered_decision_cost', 'Actual scout force/cost evidence required'),
    ('forbidden_before_La', 'Only the three remaining cells are executable; before-La repeat forbidden'),
    ('supplement_wrong_geometry', 'Saved ALPB physical provenance differs'),
    ('supplement_wrong_metal', 'Saved ALPB coordinates/charge/metal differ'),
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(record):
    require(isinstance(record, dict) and set(record) == {'path', 'sha256'}, 'Invalid evidence pin')
    require(pin(record['path']) == record, 'Evidence path/hash differs')
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return pin(path)


def deny_external_processes(attempts):
    """Irreversible process tripwire; unchanged saved readers need no children."""
    forbidden = {'subprocess.Popen', 'os.system', 'os.posix_spawn', 'os.exec', 'os.spawn'}

    def audit(event, args):
        if event in forbidden:
            attempts.append({'event': event})
            raise RuntimeError('External process forbidden in saved-only validation: ' + event)

    sys.addaudithook(audit)


def rejection(row, action, expected, checks=lambda: True):
    """An arbitrary exception is never accepted as the intended rejection."""
    try:
        action()
    except Exception as exc:
        row['observed_exception'] = type(exc).__name__
        row['observed_message'] = str(exc)
        row['status'] = 'passed' if type(exc) is ValueError and str(exc) == expected and checks() else 'failed'
        if row['status'] == 'failed':
            row['traceback'] = traceback.format_exc()
    else:
        row.update(status='failed', observed_message='Invalid evidence was accepted')


def validate(cfg, output, result):
    require(cfg['implementation']['sha256'] == ADAPTER_SHA
            and cfg['specification']['sha256'] == SPEC_SHA, 'Frozen draft/specification differs')
    spec = importlib.util.spec_from_file_location('whole_continuation_validation_target', verify(cfg['implementation']))
    verify(cfg['specification'])
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)

    # This is the single full positive replay, using only genuine admitted
    # artifacts and the real root decision. It also verifies raw native output,
    # gradients, original13049 state/decks/solvent and supplemental13142 ALPB.
    world = adapter.inputs(cfg)
    result['positive_replay'] = dict(
        status='passed', entrypoint='CONTINUATION_DRAFT.inputs',
        scout_result=cfg['scout_result'], scout_collection_result=cfg['scout_collection_result'],
        continuation_decision=cfg['continuation_decision'],
        original_preflight=cfg['original_preflight'], origin_receipt=world['origin']['receipt'],
        original_solvent_inventory={key: row['status'] for key, row in world['original_solvent'].items()},
        supplemental_after_Dy=world['solvent']['after_Dy']['receipt'],
        scientific_binding=adapter.binding(cfg))
    rows = {row['case']: row for row in result['controls']}

    def fixture(case, original, changed, changes):
        """Alter actual copies only; originals and their admission remain intact."""
        path = output / case
        path.mkdir(exist_ok=True)
        metadata = dict(test_only=True, not_scientific_evidence=True, case=case,
                        original=original, changes=changes)
        changed['_validation_fixture'] = metadata
        record = write(path / 'INVALID_EVIDENCE.json', changed)
        rows[case]['fixture'] = record
        rows[case]['fixture_provenance'] = write(path / 'FIXTURE_PROVENANCE.json', metadata)
        return record

    def config_control(case, changed_cfg):
        path = output / case
        path.mkdir(exist_ok=True)
        record = write(path / 'INVALID_CONFIG.json', changed_cfg)
        rows[case]['test_config'] = record
        changed_cfg['_config_pin'] = record
        rejection(rows[case], lambda: adapter.inputs(changed_cfg), rows[case]['expected_message'])

    case = 'missing_scout'
    test = copy.deepcopy(cfg)
    test.pop('_config_pin', None)
    test['scout_result'] = None
    config_control(case, test)

    case = 'changed_scout_identity'
    scout = read(cfg['scout_result'])
    require(scout['selected_tasks'] == ['before_La'], 'Positive scout selection changed')
    scout['selected_tasks'] = ['after_La']
    record = fixture(case, cfg['scout_result'], scout,
                     {'selected_tasks': {'original': ['before_La'], 'invalid': ['after_La']}})
    test = copy.deepcopy(cfg)
    test.pop('_config_pin', None)
    test['scout_result'] = record
    config_control(case, test)

    case = 'altered_decision_scope'
    decision = read(cfg['continuation_decision'])
    require(decision['allow_retry'] is False, 'Positive decision scope changed')
    decision['allow_retry'] = True
    record = fixture(case, cfg['continuation_decision'], decision,
                     {'allow_retry': {'original': False, 'invalid': True}})
    test = copy.deepcopy(cfg)
    test.pop('_config_pin', None)
    test['continuation_decision'] = record
    config_control(case, test)

    case = 'altered_decision_cost'
    decision = read(cfg['continuation_decision'])
    original = decision['measured_cost']['native_receipt']
    require(original != cfg['scout_result'], 'Distinct real artifact required for cost control')
    decision['measured_cost']['native_receipt'] = cfg['scout_result']
    record = fixture(case, cfg['continuation_decision'], decision,
                     {'measured_cost.native_receipt': {'original': original, 'invalid': cfg['scout_result']}})
    test = copy.deepcopy(cfg)
    test.pop('_config_pin', None)
    test['continuation_decision'] = record
    config_control(case, test)

    case = 'forbidden_before_La'
    forbidden_dir = output / case
    rows[case]['entrypoint'] = 'CONTINUATION_DRAFT.prepare_task'
    rejection(rows[case], lambda: adapter.prepare_task(world, 'before_La', forbidden_dir),
              rows[case]['expected_message'], lambda: not forbidden_dir.exists())
    rows[case]['native_directory_created'] = forbidden_dir.exists()

    # The full inputs() path hard-pins the original13142 collection and afterDy
    # receipt. Changing those would stop at an earlier immutable-parent gate.
    # Exercise its unchanged semantic ALPB reader directly using the admitted
    # afterDy cell and real, unaltered wrong-geometry / wrong-metal receipts.
    for case, substitute in (('supplement_wrong_geometry', 'before_Dy'),
                             ('supplement_wrong_metal', 'after_La')):
        receipt = world['whole'][substitute]['ALPB']
        verify(receipt)
        require(receipt != world['solvent']['after_Dy']['receipt'], 'Supplement control must change receipt')
        cell = dict(world['cells']['after_Dy'], row=dict(world['cells']['after_Dy']['row'], ALPB=receipt))
        rows[case].update(entrypoint='immutable REFERENCE.admit_alpb (helper-level semantic control)',
                          requested_cell='after_Dy', substituted_real_cell=substitute,
                          substituted_receipt=receipt, original_receipt=world['solvent']['after_Dy']['receipt'])
        rejection(rows[case], lambda cell=cell: world['old'].admit_alpb(
            cell, world['ocfg'], world['np'], world['model']), rows[case]['expected_message'])

    result['original_solvent_inventory_unchanged'] = (
        world['cells']['after_Dy']['row']['ALPB'] is None
        and world['original_solvent']['after_Dy']['status'] == 'unavailable')
    require(result['original_solvent_inventory_unchanged'], 'Validation mutated original solvent inventory')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Validation requires Slurm')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    attempts = []
    result = dict(protocol=PROTOCOL, status='execution_failure', implementation=pin(__file__),
                  config=pin(args.config), slurm_job_id=os.environ['SLURM_JOB_ID'],
                  hostname=socket.gethostname(), started_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                  positive_replay={'status': 'not_run'},
                  controls=[dict(case=case, expected_message=expected, status='not_run') for case, expected in CASES],
                  external_process_attempts=attempts,
                  new_molecular_calls=dict(DFT=0, MACE=0, GFN=0, MM=0, ALPB=0),
                  limitation='Saved-artifact gate validation only; no scientific accuracy or new native execution qualification.')
    deny_external_processes(attempts)
    try:
        cfg = read(result['config'])
        cfg['_config_pin'] = result['config']
        validate(cfg, output, result)
        require(not attempts, 'An external process was attempted')
        require(all(row['status'] == 'passed' for row in result['controls']), 'At least one intended rejection failed')
        result['status'] = 'complete'
    except Exception as exc:
        result.update(status='execution_failure', error=repr(exc), traceback=traceback.format_exc())
    result['application_wall_seconds'] = time.monotonic() - started
    result['completed_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
    write(output / 'RESULT.json', result)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
