#!/usr/bin/env python3
"""One saved-only review of the exact 18901 tighter W9 Dy scout; no engine calls."""
import csv
import hashlib
import importlib.util
import json
import re
import time
from decimal import Decimal
from pathlib import Path

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
REL = Path('lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/cluster_handoff_v1')
D = ROOT / 'diagnostics' / REL
W = ROOT / 'workspaces' / REL
CID = 'cluster__tighter__aquo9_central_Dy'
TASK = W / 'run_v1/tasks' / CID
OUT = W / 'independent_actual_review_v1'


def pin(path):
    p = Path(path).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def read(path):
    return json.loads(Path(path).read_text())


def main():
    tick, cpu = time.monotonic(), time.process_time()
    OUT.mkdir(parents=True, exist_ok=False)
    checks = {}

    def check(name, value):
        checks[name] = bool(value)
        if not value:
            raise ValueError(name)

    q = dict(status='failed_cluster_handoff_scout_review', reviewer_implementation=pin(__file__),
             new_native_calls=0, threads=1, checks=checks, blocking_findings=[])
    try:
        source = D / 'REFERENCE.py'
        check('frozen_implementation', pin(source)['sha256'] == '5133c853d89f45af81b5bf85a5eebc06e21a14f432df9618e7739beeefd79810')
        spec = importlib.util.spec_from_file_location('frozen_cluster_handoff', source)
        h = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(h)
        panel = W / 'preparation_v1/PANEL.json'
        p = h.load(panel)
        check('frozen_panel', pin(panel)['sha256'] == 'db44019e47a841d5bb89519502ab5acc0b1553270330eb9204d114ca69bf3bda')
        native_path = TASK / 'execution/RESULT.json'
        collection_path = W / 'scout_collection_v1/result/RESULT.json'
        r, col = read(native_path), read(collection_path)
        check('exact_native_result', pin(native_path)['sha256'] == 'e45638de614418941d0c2897ff4e0c45f2de4f0d042775d93f0850afed9c3946')
        c = next(x for x in p['cells'] if x['cell_id'] == CID)
        replay = h.admit_result(native_path, c, p, panel)
        check('all_31_original_native_gates', len(replay['checks']) == 31 and all(replay['checks'].values()))
        check('complete_collector_exact_raw_replay', col['status'] == 'cluster_handoff_scout_complete' and col['declared'] == col['complete'] == 1 and col['cell_ids'] == [CID] and col['rows'] == {CID: replay})
        check('collection_source_and_panel', col['implementation'] == pin(source) and col['panel'] == pin(panel))
        st, release = read(r['start']['path']), read(r['release']['path'])
        check('one_original_native_invocation', r['job_id'] == st['job_id'] == '18901' and r['native_program_invocations_attempted'] == 1 and release['stage'] == 'scout' and release['cell_ids'] == [CID])
        check('unchanged_runtime_template', Path(r['input']['path']).read_bytes() == Path(c['template_input']['path']).read_bytes())
        check('only_declared_parent_MPI_delta', Path(r['input']['path']).read_text() == Path(p['original_inputs'][CID]['path']).read_text().replace('%pal nprocs 4 end', '%pal nprocs 24 end'))
        check('24_MPI_one_node_default_allocation', st['ranks'] == 24 and st['nodes'] == 1 and st['memory_MiB'] == 96000 and st['maxcore_MiB'] == 3000)
        check('physical_effective_state_distinction', c['charge'] == 3 and c['multiplicity'] == 1 and c['physical_multiplicity_of_parent'] == 6 and not any(c['ghost_mask']) and len(c['symbols']) == 28)
        raw = Path(r['output']['path']).read_text()
        native = raw.split('****END OF INPUT****')[-1]
        check('normal_termination_empty_stderr', raw.count('****ORCA TERMINATED NORMALLY****') == 1 and Path(r['stderr']['path']).stat().st_size == 0)
        check('native_state_and_counts', bool(re.search(r'Total Charge\s+Charge\s+\.+\s+3\s*$', native, re.M)) and bool(re.search(r'Multiplicity\s+Mult\s+\.+\s+1\s*$', native, re.M)) and r['total_electrons'] == [98, 98] and r['correlated_electrons'] == [98] and r['ecp_core_counts'] == [55] and sum(r['basis_functions_by_center']) == 434 and r['aux_C']['totals'] == [1323])
        check('native_no_frozen_core', 'Frozen core treatment                      ... NO frozen core' in native)
        check('both_localizations_complete', native.count('LOCALIZATION HAS CONVERGED.') == 2)
        cc = re.search(r'^\s*18\s+(-723\.291009260)\s+(-2\.754796398)\s+(0\.000000001)\s+(0\.000000005)\s+[0-9.]+\s*$', native, re.M)
        tri = re.search(r'^\s*4\s+(-0\.075397739719)\s+(0\.000000000973)\s+[0-9.]+\s*$', native, re.M)
        check('independent_CC_terminal_row', cc is not None and Decimal(cc[4]) <= Decimal('1e-8') and r['CC_terminal_residuals'] == [5e-9])
        check('independent_iterative_triples_terminal_row', tri is not None and Decimal(tri[2]) <= Decimal('1e-8') and r['triples_terminal_deltas'] == [9.73e-10])
        check('unchanged_printed_CC_and_triples_limits', r['CC_printed_tolerance'] == [1e-8] and r['triples_printed_energy_tolerance'] == [1e-8])
        energies = {}
        for key, pattern in [('CCSD', r'^E\(CCSD\)\s+\.+\s+(\S+)\s*$'), ('CCSD_T', r'^E\(CCSD\(T\)\)\s+\.+\s+(\S+)\s*$'), ('triples', r'^Triples Correction \(T\)\s+\.+\s+(\S+)\s*$'), ('final', r'^FINAL SINGLE POINT ENERGY\s+(\S+)\s*$')]:
            tokens = re.findall(pattern, native, re.M)
            check('one_raw_energy_token_' + key, len(tokens) == 1 and Decimal(tokens[0]).is_finite())
            energies[key] = tokens[0]
        check('raw_final_and_CCSD_T_identity_to_stored_values', float(energies['CCSD_T']) == r['energy_hartree'] and r['CCSD_energies_hartree'] == [float(energies['CCSD'])] and r['triples_hartree'] == [float(energies['triples'])] and r['printed_final_energies_hartree'] == [float(energies['final'])])
        stages = r['stage_specific_local_evidence']
        check('tighter_main_cutoffs', stages['primary_thresholds'] == {'TCutPairs': [1e-6], 'TCutPNO': [1e-8], 'TCutMKN': [1e-4], 'TCutDO': [5e-3]})
        check('full_local_MP2_terminal', stages['stages'][1]['terminal'] == dict(iteration=7, energy=-2.690486852946, delta=3.9e-8, residual=9.8e-8, converged=True))
        check('D4_and_physical_qualification_not_added', c['d4'] is False and col['physical_accuracy_qualified'] is False and col['cluster_CP_or_B_qualified'] is False)
        check('Tight_failures_preserved', col['original_Tight_controls_passed'] is False and col['tighter_controls_passed'] is True)
        event = read(TASK / 'completion_event.json')['records']['18901']
        account_path = TASK / 'ROOT_SACCT.psv'
        account = {x['JobID']: x for x in csv.DictReader(account_path.read_text().splitlines(), delimiter='|')}
        job, batch, collect = account['18901'], account['18901.batch'], account['18902']
        check('native_scheduler_completion', event['state'] == job['State'] == 'COMPLETED' and event['exit_code'] == job['ExitCode'] == '0:0' and event['elapsed_seconds'] == job['ElapsedRaw'] == '6748' and event['allocated_cpus'] == job['AllocCPUS'] == '24')
        check('collector_scheduler_completion', collect['State'] == 'COMPLETED' and collect['ExitCode'] == '0:0' and collect['ElapsedRaw'] == '3' and collect['AllocCPUS'] == '1')
        check('actual_batch_memory_receipt', batch['MaxRSS'] == '98303396K' and job['ReqMem'] == '96000M' and batch['State'] == 'COMPLETED')
        files = {x.name: x.stat().st_size for x in (TASK / 'execution').iterdir() if x.is_file()}
        check('terminal_artifact_size_inventory', all(files.get(k) == n for k, n in r['artifact_sizes'].items()))
        counts = [dict(cell_id=x['cell_id'], electrons=x['expected_electrons'], orbital_functions=x['expected_basis_functions'], aux_C_functions=(369 if x['metal'] == 'Dy' else 351) + (9 if 'aquo9' in x['geometry_id'] else 8) * 106) for x in p['cells']]
        q.update(status='passed_cluster_handoff_scout_review', implementation=pin(source), panel=pin(panel), collection=pin(collection_path), result=pin(native_path), native_checks_passed=True, cost_reviewed=True,
                 declaration=pin(D / 'DECLARATION.md'), source_review=pin(W / 'independent_review_v1/REVIEW.json'), ready=pin(D / 'READY.json'), release=r['release'], input=r['input'], output=r['output'], start=r['start'], stderr=r['stderr'], parser=p['parser'], native_input_adapter=p['native_input_adapter'],
                 qualified_controls=p['qualified_controls'], failed_Tight_controls=p['failed_tight_controls'],
                 native_summary=dict(energy_tokens_hartree=energies, CC_iteration=18, CC_residual=5e-9, iterative_triples_iteration=4, triples_delta_hartree=9.73e-10, T1_diagnostic=r['T1_diagnostic'], full_local_MP2=stages['stages'][1]['terminal'], orbital_functions=434, aux_C_functions=1323, electrons=98, correlated_electrons=98, ECP_core=55, charge=3, effective_multiplicity=1, physical_parent_multiplicity=6, real_centers=28, main_local_thresholds=stages['primary_thresholds'], all_stage_thresholds=stages['map_threshold_reports'], printed_scf_criteria=r['printed_scf_criteria'], scf_cycles=14, scf_mode=r['native_scf_mode'], scf_stop=r['native_scf_stop_signals'], warnings=[x for x in raw.splitlines() if re.match(r'\s*(WARNING|Warning):', x)]),
                 cost=dict(scheduler_receipt=pin(account_path), scheduler_query=pin(TASK / 'ROOT_SACCT_QUERY.json'), completion_receipt=pin(TASK / 'completion_event.json'), collector_completion_receipt=pin(W / 'scout_collection_v1/completion_event.json'), native_elapsed_seconds=6748, native_allocated_CPUs=24, native_allocated_CPU_seconds=161952, collector_allocated_CPU_seconds=3, combined_scheduler_CPU_seconds=161955, native_TotalCPU=job['TotalCPU'], application_elapsed_seconds=r['elapsed_seconds'], application_allocated_CPU_seconds=r['allocated_cpu_seconds'], batch_MaxRSS_KiB=98303396, batch_MaxRSS_MiB=98303396 / 1024, memory_grant_MiB=96000, MaxCore_per_rank_MiB=3000, ORCA_printed_MDCI_maximum_MB=3033.7, ORCA_printed_average_batch_MB=2405, ORCA_printed_all_IJ_pairs_MB=87333, ORCA_printed_triples_available_MB=2742, retained_artifact_bytes=sum(files.values()), retained_artifacts=files),
                 remaining_dimensions=counts, recommendation='Supports root release of the exact remaining seven frozen cells after explicit acceptance of measured cost and near-limit memory risk; no scout repeat or changed recipe.',
                 limitations=[
                     'Only the largest tighter W9 Dy scout is admitted; neither four-cluster contrast nor eight-cell cutoff sensitivity is available.',
                     'The batch MaxRSS is almost the full grant. ORCA and MPI are descendants of one batch step, with no separately reported rank steps. Do not multiply MaxRSS by24 or identify its accounting with per-rank RSS, anonymous memory, or a measured synchronized whole-job peak. No comfortable headroom or memory guarantee is established.',
                     'Terminal scratch sizes are retained sizes, not peak temporary disk demand. Native reports RI integrals on disk and normal cleanup.',
                     'Remaining W9 cells have the same orbital/electron dimensions and no larger AuxC count; W8 is smaller and Tight cutoffs are looser. This supports feasibility but does not bound their runtime or peak memory.',
                     'SCF stopped by Total+1el-Energy/Energy Check; printed MAX/RMS density and DIIS criteria remain false. The DIIS number is inherited across the SOSCF switch, not a final SOSCF error.',
                     'Small-system tighter18 qualification and converged local CC/T1 do not establish canonical cluster accuracy, basis convergence, frozen-f/SOC/multiplet accuracy, physical state, forces, dipoles, CP/B, sampling, or affinity.',
                     'Tight18 retains three failed controls. Tight remains the declared lower-accuracy comparator. D4 is absent from CC and must remain separate in PBE0 comparisons.',
                     'No native calls, author control-suite reruns, scheduler queries or submissions were made in this independent review.'])
    except Exception as exc:
        q['blocking_findings'].append(repr(exc))
    q.update(check_count=len(checks), wall_seconds=time.monotonic() - tick, CPU_seconds=time.process_time() - cpu)
    (OUT / 'REVIEW.json').write_text(json.dumps(q, indent=2) + '\n')
    print(json.dumps(dict(status=q['status'], review=pin(OUT / 'REVIEW.json'), checks=q['check_count'], blocking_findings=q['blocking_findings'])))
    return 0 if q['status'] == 'passed_cluster_handoff_scout_review' else 1


if __name__ == '__main__':
    raise SystemExit(main())
