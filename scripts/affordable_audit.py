"""Read-only audit of released evidence, live provenance and unopened tasks."""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import subprocess

from affordable_common import read_json, record, verify, digest, energy, contrast, classify_raw, write_new


def audit(root):
    root = Path(root).resolve()
    base = root / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914'
    release = read_json(base / 'result.json')
    holdout = read_json(base / 'reserved_crystal_holdout/result/holdout_result.json')
    rows, costs, issues = [], [], []
    for stratum, payload in [('calibration', release), ('consumed_structural_transfer', holdout)]:
        for row in payload['scores']:
            ident = row.get('panel_id', row.get('target_id', row.get('pdb_id')))
            endpoints, artifact_issues, receipts = {}, [], {}
            for metal, art in row['artifacts'].items():
                for kind in ('input', 'xyz', 'output', 'execution'):
                    rec = art[kind]
                    try:
                        verify(rec)
                    except ValueError as exc:
                        artifact_issues.append({'metal': metal, 'kind': kind,
                                                'published': rec, 'live_sha256': digest(rec['path']) if Path(rec['path']).is_file() else None,
                                                'error': str(exc)})
                op = Path(art['output']['path'])
                try:
                    endpoints[metal] = energy(op)
                except ValueError as exc:
                    endpoints[metal] = None
                    artifact_issues.append({'metal': metal, 'error': str(exc)})
                ep = Path(art['execution']['path'])
                if ep.is_file():
                    receipt = read_json(ep)
                    receipts[metal] = record(ep)
                    seconds = (datetime.fromisoformat(receipt['finished_at_utc']) -
                               datetime.fromisoformat(receipt['started_at_utc'])).total_seconds()
                    costs.append({'case': ident, 'metal': metal, 'wall_seconds': seconds,
                                  'mpi_ranks': receipt['parallelism']['nprocs'],
                                  'active_rank_seconds': seconds * receipt['parallelism']['nprocs'],
                                  'allocated_core_seconds': None,
                                  'allocation': receipt['allocation'], 'receipt': record(ep),
                                  'published_receipt_hash_matches': digest(ep) == art['execution']['sha256'],
                                  'peak_memory_bytes': None, 'gpu_seconds': None})
            archived = {m: row['artifacts'][m]['energy_hartree'] for m in ('La','Ca')}
            score = contrast(archived['Ca'], archived['La'], release['aquo_reporting_gauge']['delta_E_aquo_hartree'])
            parsed_match = all(endpoints[m] is not None and abs(endpoints[m] - archived[m]) < 1e-10 for m in archived)
            rows.append({'case': ident, 'evidence_stratum': stratum,
                         'use': 'threshold_calibration' if stratum == 'calibration' else 'method_development',
                         'baseline_protocol_id': payload['protocol_id'], 'released_endpoints_hartree': archived,
                         'live_parsed_endpoints_hartree': endpoints, 'parser_energy_match': parsed_match,
                         'released_score': score, 'released_decision': classify_raw(score['R_kcal_mol'], release, payload['protocol_id']),
                         'algebra_matches_release': abs(score['R_kcal_mol'] - row['R_kcal_mol']) < 1e-7,
                         'provenance_status': 'verified' if not artifact_issues else 'live_artifacts_differ_from_release',
                         'artifact_issues': artifact_issues, 'source_record': row,
                         'repaired_baseline': {'status': 'not_applicable_pqq_fixed_core', 'score': None},
                         'environmental_challenger': {'status': 'not_run', 'score': None}})
            issues.extend(artifact_issues)
    inv = read_json(root / 'diagnostics/nonpqq_direct_site_benchmark_20260915/prepared_task_inventory.json')
    unopened = [{'task': t['task_key'], 'output_exists': Path(t['expected_output_path']).exists(),
                 'energy_inspected_in_this_audit': False} for t in inv['tasks']]
    external = root / 'diagnostics/pqq_q46444_1kb0_external_validation_20260915/prepared/01_1KB0'
    unopened.extend({'task': str(p), 'output_exists': p.exists(), 'energy_inspected_in_this_audit': False}
                    for p in [external / 'sp_pmdh_fc_holdout_1kb0_La.out', external / 'sp_pmdh_fc_holdout_1kb0_Ca.out'])
    names = ['scripts/carve_generic.py', 'scripts/carve_with_pqq.py', 'scripts/protonate_cif.py',
             'scripts/run_orca_task_manifest.py', 'scripts/render_orca_runtime_input.py',
             'scripts/result_protocol.py', 'scripts/site_mechanics.py',
             'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/fixed_core_carver.py']
    jobs = sorted({c['allocation']['slurm_job_id'] for c in costs})
    accounting = subprocess.run(['sacct', '-j', ','.join(jobs), '--format=JobID,State,AllocCPUS,ElapsedRaw,CPUTimeRAW,MaxRSS,NodeList', '-P'],
                                capture_output=True, text=True)
    return {'schema_version': 'alquemia.affordable_audit.v1',
            'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
            'implementation_inventory': [record(root/n) for n in names],
            'released_calibration': record(base/'result.json'), 'released_holdout': record(base/'reserved_crystal_holdout/result/holdout_result.json'),
            'aquo_gauge': release['aquo_reporting_gauge'], 'bands': release['calibration']['released_supported_bands'],
            'rows': rows, 'costs': costs, 'sacct_raw': accounting.stdout, 'sacct_error': accounting.stderr,
            'unopened_outputs': unopened, 'provenance_issue_count': len(issues),
            'parser_matches': sum(r['parser_energy_match'] for r in rows),
            'algebra_matches': sum(r['algebra_matches_release'] for r in rows),
            'fully_verified_pairs': sum(r['provenance_status']=='verified' for r in rows),
            'high_level_evaluations_launched': 0,
            'cost_caution': 'Active MPI-rank seconds are not allocated core-seconds. Use sacct; do not apportion inconsistent receipts.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    result=audit(a.root)
    write_new(a.output,result)
    print(f"{len(result['rows'])} released pairs; {result['parser_matches']} parse matches; "
          f"{result['fully_verified_pairs']} fully hash-verified; {result['provenance_issue_count']} provenance issues")


if __name__ == '__main__':
    main()
