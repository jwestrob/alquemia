#!/usr/bin/env python3
"""Slurm-only saved-matrix regression; no calculators or invented energies."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    if isinstance(p, dict):
        assert pin(p['path']) == p, ('changed pin', p)
        p = p['path']
    return json.loads(Path(p).read_text())


def close(a, b):
    assert math.isfinite(a) and math.isfinite(b) and abs(a-b) <= 1e-8, (a, b)


def main(config, output):
    assert os.environ.get('SLURM_JOB_ID'), 'actual Slurm required before analysis'
    c = read(config)
    assert pin(__file__) == c['test'] and pin(c['matrix']['path']) == c['matrix']
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    began = time.monotonic()
    receipt = {'status': 'failed', 'config': pin(config), 'job_id': os.environ['SLURM_JOB_ID'],
               'new_molecular_calls': 0, 'tests': []}
    try:
        spec = importlib.util.spec_from_file_location('exaa_test_matrix', c['matrix']['path'])
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        saved = read(c['targets']); final = read(c['reference'])
        reference = {r['record_id']: r for r in final['rows']}
        targets = [r for r in saved['rows'] if r['target']]
        assert len(targets) == 9 and len({r['record_id'] for r in targets}) == 9
        records = list(targets)
        # An actual marginal-energy example, if present, tests origin retention.
        marginal = next((r for r in saved['rows'] if r['status'] == 'available' and
            any(reference[r['record_id']]['pool']['operational'][z]['candidate'] == 'origin'
                and reference[r['record_id']]['pool']['mathematical'][z]['candidate'] != 'origin'
                for z in m.METALS)), None)
        if marginal and marginal['record_id'] not in {r['record_id'] for r in records}:
            records.append(marginal)
        rows = []; mr = {}; sr = []
        for r in records:
            cid = r['record_id']; rows.append({'case_id': cid})
            for z in m.METALS:
                for q in m.ORDER:
                    cell = r['matrix'][z][q]; native = read(cell['native_result'])
                    assert native['status'] == 'available'
                    close(native['G_qvac_hartree'], cell['G_qvac_hartree'])
                    use = next(u for u in native['input']['uses'] if u['MACE_eV'] == cell['MACE_eV'])
                    old_mace = read(use['MACE'])
                    assert old_mace['status'] in ('complete', 'computed')
                    assert old_mace['energy_eV'] == cell['MACE_eV']
                    tid = cid+'__'+z+'__'+q
                    # Adapter views of actual receipts; no new scientific result is created.
                    mr[tid] = {'status': 'complete', 'energy_eV': cell['MACE_eV'], 'receipt': use['MACE']}
                    sr.append({'task_id': tid, 'status': 'available', 'G_qvac_hartree': native['G_qvac_hartree'],
                               'result': cell['native_result']})
        view = {'protocol_id': 'actual_saved_matrices_implementation_test', 'rows': rows}
        observed = m.summarize(view, mr, sr)
        assert observed['denominator'] == len(records)
        for old, now in zip(records, observed['rows']):
            assert now['case_id'] == old['record_id'] and now['status'] == 'available'
            assert now['classification'] is None and now['classification_reason'] == 'changed_context_has_no_calibration'
            expected = reference[old['record_id']]['pool']
            for z in m.METALS:
                for q in m.ORDER:
                    work = now['works'][z][q]
                    close(work['total_kcal_mol'], old['works'][z][q])
                    close(work['MACE_kcal_mol']+work['solvent_kcal_mol'], work['total_kcal_mol'])
                for policy in ('operational', 'mathematical'):
                    actual = now['selected'][policy][z]; wanted = expected[policy][z]
                    assert actual['candidate'] == wanted['candidate']
                    close(actual['work']['total_kcal_mol'], wanted['work_kcal_mol'])
            for q in m.ORDER:
                wanted = old['paired_candidates'][q]['R']-old['paired_candidates']['origin']['R']
                close(now['same_geometry_differential_works_kcal_mol'][q], wanted)
            for policy in ('operational', 'mathematical'):
                wanted = expected[policy]['R_model_kcal_mol']-expected['origin']['R_model_kcal_mol']
                close(now['selected'][policy]['differential_work_kcal_mol'], wanted)
        receipt['tests'].append({'name': 'actual_nine_target_matrices_and_optional_marginal',
                                 'passed': True, 'records': [r['record_id'] for r in records],
                                 'tolerance_kcal_mol': 1e-8})
        # Intentionally incomplete copies of real fixtures, never scientific observations.
        absent = []
        for r in targets:
            cid = r['record_id']; one = {'protocol_id': view['protocol_id'], 'rows': [{'case_id': cid}]}
            for z in m.METALS:
                for q in m.ORDER:
                    tid = cid+'__'+z+'__'+q
                    reduced = [x for x in sr if x['task_id'] != tid]
                    row = m.summarize(one, mr, reduced)['rows'][0]
                    assert row['status'] == 'unavailable' and row['selected'] is None and row['works'] is None
                    assert len(row['failures']) == 1 and row['classification'] is None
                    absent.append(tid)
        receipt['tests'].append({'name': 'missing_required_real_cell_never_falls_back', 'passed': True,
                                 'corrupted_fixture_count': len(absent), 'removed_cell_ids': absent})
        receipt['tests'].append({'name': 'strict_origin_tolerance_real_marginal_example',
                                 'status': 'passed' if marginal else 'unrun_no_marginal_case_in_declared_saved_panel',
                                 'record_id': marginal['record_id'] if marginal else None,
                                 'note': 'No synthetic boundary energies were manufactured.'})
        receipt.update(status='passed', denominator_targets=9, matrix_records=len(records),
                       scientific_scope='saved-artifact algebra regression only; new endpoints unexecuted')
    except BaseException as exc:
        receipt['reason'] = type(exc).__name__+': '+str(exc)
        raise
    finally:
        receipt['elapsed_seconds'] = time.monotonic()-began
        with (out/'RESULT.json').open('x') as f:
            json.dump(receipt, f, indent=2, allow_nan=False); f.write('\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', required=True); p.add_argument('--output', required=True)
    a = p.parse_args(); main(a.config, a.output)
