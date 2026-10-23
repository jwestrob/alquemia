#!/usr/bin/env python3
"""Saved-only review of the single 06:25–06:28 UTC observation; no polling."""
import datetime as dt
import hashlib
import importlib.util
import json
import math
import re
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path.cwd()
OUT = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/tighter_dy_acetate_v1/independent_progress_review_v2'

def read(p):
    return json.loads(Path(p).read_text())

def pin(p):
    p = Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def original_pin(metadata):
    assert pin(metadata['snapshot_path'])['sha256'] == metadata['sha256']
    assert not metadata['changed_during_read']
    return dict(path=metadata['source_path'], sha256=metadata['sha256'])

def main():
    tick, cpu = time.monotonic(), time.process_time()
    snap = read(OUT / 'SNAPSHOT.json')
    panel = read(OUT / 'snapshots/PANEL.json')
    pp = original_pin(snap['files']['panel'])
    rp = original_pin(snap['files']['root_release'])
    assert panel['adapter'] == original_pin(snap['files']['adapter'])
    assert pin(panel['adapter']['path']) == panel['adapter']
    assert pin(panel['parser']['path']) == panel['parser']
    spec = importlib.util.spec_from_file_location('tighter_reference_progress', panel['adapter']['path'])
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    release = read(OUT / 'snapshots/ROOT_RELEASE.json')
    helper.release_contract(panel, pp['path'], release)
    rows = {}
    for cell in panel['cells']:
        role = cell['role']
        files = snap['roles'][role]
        pins = {name: original_pin(m) for name, m in files.items() if m.get('available', True)}
        st = read(files['START.json']['snapshot_path'])
        raw = Path(files['endpoint.out']['snapshot_path']).read_text()
        inp = Path(files['endpoint.inp']['snapshot_path']).read_text()
        assert st['cell'] == cell and st['panel'] == pp and st['release'] == rp
        assert st['input'] == pins['endpoint.inp']
        assert all(st[k] == panel[k] for k in ('adapter', 'native_input_adapter', 'parser'))
        assert st['ranks'] == 8 and st['memory_MiB'] == 32000 and st['maxcore_MiB'] == 3000
        assert st['maxcore_MiB'] == math.floor(.75 * st['memory_MiB'] / st['ranks'])
        assert st['command'] == [panel['orca']['path'], 'endpoint.inp', '--bind-to none']
        assert inp == helper.deck(panel, cell, 8, st['maxcore_MiB'])
        assert st['maximum_native_calls'] == 1
        common = dict(job_id=st['job_id'], captured_pins=pins, resources={k: st[k] for k in ('host','ranks','memory_MiB','maxcore_MiB')}, exact_input_and_producer_joins=True)
        if files['RESULT.json']['available']:
            result = read(files['RESULT.json']['snapshot_path'])
            assert result['cell'] == cell and result['cell_id'] == cell['cell_id']
            assert all(result[k] == st[k] for k in ('panel','release','input','adapter','native_input_adapter','parser','job_id'))
            assert result['start'] == pins['START.json'] and result['output'] == pins['endpoint.out'] and result['stderr'] == pins['endpoint.err']
            assert result['returncode'] == 0 and result['native_program_invocations_attempted'] == 1
            parsed = json.loads(json.dumps(helper.parse(raw, cell, result['returncode'])))
            assert parsed['status'] == 'complete' and all(parsed['checks'].values())
            assert all(result[k] == v for k, v in parsed.items())
            rows[role] = dict(common, status='saved_raw_replay_passed', checks=len(parsed['checks']), all_checks_passed=True,
                all_parsed_fields_exact=True, native_admission=parsed,
                elapsed_seconds=result['elapsed_seconds'], allocated_cpu_seconds=result['allocated_cpu_seconds'])
        else:
            native = raw.split('****END OF INPUT****')[-1]
            fatal = [line for line in native.splitlines() if re.search(r'error termination|aborting|out of memory|SCF NOT CONVERGED|DID NOT CONVERGE|segmentation fault|MPI_ABORT', line, re.I)]
            last_mtime = files['endpoint.out']['source_mtime_ns_after'] / 1e9
            observed = dt.datetime.fromisoformat(snap['observation_utc']).timestamp()
            rows[role] = dict(common, status='running_no_terminal_result', native_output_tail=raw.splitlines()[-35:],
                output_mtime_utc=dt.datetime.fromtimestamp(last_mtime, dt.timezone.utc).isoformat(),
                output_silence_seconds_at_capture=observed-last_mtime, output_bytes=files['endpoint.out']['bytes_read'],
                stderr_bytes=files['endpoint.err']['bytes_read'], explicit_failure_lines=fatal,
                last_native_stage='RI/C V-Matrix + Cholesky decomposition',
                SCF_or_CC_or_triples_residual_available=False,
                terminal_energy_available=False)
    process_files = [OUT / (host + '_processes.json') for host in ('node-64-768g-5','node-64-768g-6')]
    hosts = {}
    all_ranks = []
    for file in process_files:
        record = read(file)
        ranks = [p for p in record['processes'] if Path(p['command'].split()[0]).name == 'orca_startup_mpi']
        assert len(ranks) == 4 and record['uid'] == 601
        assert all(p['cwd'] == str(Path(panel['execution_root']) / 'tasks/tighter__cc__acetate_bi__Dy__ligand_ghost_ion/execution') for p in record['processes'])
        hosts[record['host']] = dict(source=pin(file), observation_utc=record['observation_utc'], ranks=ranks,
            rank_count=len(ranks), rank_CPU_seconds=sum(p['CPU_seconds'] for p in ranks),
            rank_RSS_MiB=sum(p['RSS_MiB'] for p in ranks), parent_processes=[p for p in record['processes'] if p not in ranks])
        all_ranks += ranks
    scheduler = dict(provenance='Transcribed from the single sstat/scontrol tool response between raw and process snapshots; no second scheduler query.',
        commands=['sstat -j18518 --allsteps --parsable2 --format=JobID,NTasks,AveCPU,AveRSS,MaxRSS,MaxRSSNode,MaxRSSTask,AveDiskRead,AveDiskWrite','scontrol show job -o18518'],
        sstat_stdout='JobID|NTasks|AveCPU|AveRSS|MaxRSS|MaxRSSNode|MaxRSSTask|AveDiskRead|AveDiskWrite\n18518.batch|1|07:40:48|302956K|303860K|node-64-768g-5|0|307933379|29027568\n18518.0|1|07:40:48|182840K|187516K|node-64-768g-6|0|52101358|4811339\n',
        scontrol_extracted=dict(JobId=18518,JobState='RUNNING',Reason='None',Restarts=0,ExitCode='0:0',RunTime='01:55:14',TimeLimit='5-00:00:00',
            StartTime='2026-10-21T21:31:24',scheduler_time_zone='America/Los_Angeles',NodeList='node-64-768g-[5-6]',BatchHost='node-64-768g-5',
            NumNodes=2,NumCPUs=8,NumTasks=8,CPUsPerTask=1,ReqTRES='cpu=8,mem=32000M,node=1,billing=8',AllocTRES='cpu=8,mem=32000M,node=2,billing=8',MinMemoryCPU='4000M'),
        allocation_note='Eight ranks are distributed four per allocated host; sstat NTasks=1 per step describes launcher accounting, not one molecular MPI rank. A one-time CPU snapshot cannot distinguish useful Cholesky arithmetic from busy MPI waiting. No evidence of memory pressure in observed RSS.')
    q = dict(status='saved_progress_review_complete_CP_cycle_unavailable',implementation=pin(__file__),snapshot=pin(OUT/'SNAPSHOT.json'),
        panel=pp,release=rp,rows=rows,declared_roles=3,terminal_raw_roles_admitted=2,CP_cycle_qualified=False,CP_comparison=None,
        original_CP_diagnostic_limit_kcal_mol=.25,full_eighteen_qualified=False,cluster_release_authorized=False,
        scheduler_snapshot=scheduler,owned_process_snapshot=dict(hosts=hosts,total_MPI_ranks=len(all_ranks),
            rank_states=[p['state'] for p in all_ranks],total_rank_CPU_seconds=sum(p['CPU_seconds'] for p in all_ranks),
            total_rank_RSS_MiB=sum(p['RSS_MiB'] for p in all_ranks),receipt=pin(OUT/'OWNED_PROCESS_SNAPSHOT_RECEIPTS.json')),
        disposition='CPU-active MPI startup without printed numerical progress; neither productive progress nor a busy-wait hang is established. No explicit native failure. Await root disposition/terminal evidence; no cancellation, retry or additional monitoring performed.',
        limitations=['Both completed roles stopped at Energy Check; printed density and DIIS criteria remain false.','Native admission is not physical-state or method-accuracy qualification.','Complex T1=0.020818158 is retained as a diagnostic, not a unique explanation.','An incomplete CP cycle cannot supply a method comparison or pass the unchanged 0.25 kcal/mol limit.'],
        cost=dict(completed_native_wall_seconds=sum(v.get('elapsed_seconds',0) for v in rows.values()),
            completed_native_allocated_CPU_seconds=sum(v.get('allocated_cpu_seconds',0) for v in rows.values()),
            running_rank_observed_CPU_seconds=sum(p['CPU_seconds'] for p in all_ranks),
            root_reported_completed_scheduler_CPU_seconds=(4242+64)*8,
            running_terminal_cost_available=False,review_wall_seconds=time.monotonic()-tick,review_CPU_seconds=time.process_time()-cpu,
            threads=1,new_native_calls=0,new_scheduler_queries=0,new_remote_process_queries=0))
    dest=OUT/'REVIEW.json'
    with dest.open('x') as f: json.dump(q,f,indent=2,allow_nan=False);f.write('\n')
    report='''# One-time tighter Dy-acetate progress review

At 06:25–06:28 UTC on 22 October, two roles pass exact saved raw admission (32 native checks each), while the ligand-with-ghost-ion role has no terminal result. The three-role CP cycle remains unavailable; the original 0.25 kcal/mol gate and cluster hold are unchanged.

| Role / job | Saved admission | Native wall | Allocated CPU |
| --- | --- | ---: | ---: |
| Complex / 18516 | Complete | 4241.074 s | 33928.593 CPU-s |
| Ion with ghost ligand / 18517 | Complete | 59.203 s | 473.627 CPU-s |
| Ligand with ghost ion / 18518 | Incomplete | Not terminal | Not terminal |

The completed outputs retain Energy Check stopping and failed printed MAX/RMS density and DIIS criteria. CC residuals are 8e-9 and 3e-9; terminal triples energy changes are 5.32e-9 and zero against 1e-8. The complex T1 diagnostic is 0.020818158. All declared cutoff, basis, core, electron and producer checks replay exactly; this does not qualify physical accuracy or electronic-state uniqueness.

Job 18518's captured output ends at `Calculating RI/C V-Matrix + Cholesky decomp....`, before any SCF/CC/triples residual. It had not changed for about 113 min 45 s and stderr was empty. No explicit native failure appeared. The native input reports 189 orbital and 718 auxiliary-C functions and MaxCore 3000 MiB.

The single scheduler/process observation confirms eight live `orca_startup_mpi` ranks, four on each of node-64-768g-5 and node-64-768g-6. Every rank is runnable with about 99.98–99.99% lifetime CPU and 72–82 MiB RSS. The ORCA/mpirun/srun/orted parent chain and exact endpoint working directories match. The allocated two-node layout is real; rank count is correct. Launcher-step `NTasks=1` does not mean only one MPI rank ran.

This excludes idle or dead workers at the observation, but cannot distinguish useful RI/C work from a busy MPI wait. Low observed RSS gives no indication of memory pressure. No cancellation, retry, new call, threshold change or repeated polling was performed. Root retains the execution decision; no incomplete CP contrast is reported.

Completed native allocated cost is 34,402.221 CPU-s (rounded scheduler cost supplied by root: 34,448 CPU-s). The eight running ranks had accumulated 55,933.46 CPU-s at their process snapshots; this is observed process CPU, not a terminal allocation cost. Raw snapshots, source pins, exact process parentage, all printed SCF flags and the transcribed one-time scheduler response are retained in `REVIEW.json` and its linked evidence.
'''
    assert abs(q['owned_process_snapshot']['total_rank_CPU_seconds']-55933.46)<1e-7
    for folder in (OUT,HERE):
        with (folder/'REPORT.md').open('x') as f:f.write(report)
    print(json.dumps(dict(review=pin(dest),report=pin(HERE/'REPORT.md'),summary={k:q[k] for k in ('status','terminal_raw_roles_admitted','CP_cycle_qualified')},cost=q['cost'])))

if __name__=='__main__':main()
