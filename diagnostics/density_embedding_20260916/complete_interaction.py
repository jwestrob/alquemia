"""Read-only completion after the named existing EDA job; launches no science."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'manifest', 'job', 'previous-job', 'output', 'vault'):
        p.add_argument('--' + name, required=True)
    a = p.parse_args()
    root, out = Path(a.root).resolve(), Path(a.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(root / 'scripts'))
    from affordable_common import read_json, record, snapshot_implementation, verify, write_new
    import shutil
    impl = out / 'implementation'
    snapshot_implementation(impl)
    extras = {}
    for name in ('interaction_decomposition.py', 'global_electrostatic.py', 'result_protocol.py'):
        source = root / 'scripts' / name
        shutil.copyfile(source, impl / name)
        extras[name] = {'source': record(source), 'copy': record(impl / name)}
    shutil.copyfile(__file__, out / 'complete_interaction.py')
    write_new(out / 'extra_implementation.json', extras)
    write_new(out / 'started.json', {'arguments': vars(a), 'manifest': record(a.manifest),
                                   'script': record(out / 'complete_interaction.py')})
    collector = impl / 'interaction_decomposition.py'
    old = verify(read_json(a.manifest)['retry_of'])
    # Verify the frozen collector/imports immediately on the real failed attempt.
    subprocess.run([sys.executable, str(collector), 'collect', '--manifest', str(old),
                    '--output', str(out / 'initial_failed_collection.json')], check=True)
    subprocess.run([sys.executable, str(impl / 'affordable_watch.py'), '--job', a.job,
                    '--output', str(out / 'terminal_accounting.json')], check=True)
    subprocess.run([sys.executable, str(collector), 'collect', '--manifest', a.manifest,
                    '--output', str(out / 'collection.json')], check=True)
    subprocess.run([sys.executable, str(collector), 'report', '--result', str(out / 'collection.json'),
                    '--output', str(out / 'REPORT.md')], check=True)
    jobs = [a.previous_job, a.job]
    acct = subprocess.run(['sacct', '-j', ','.join(jobs),
                          '--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,MaxRSS,NodeList', '-P'],
                          text=True, capture_output=True, check=True)
    lines = [line.split('|') for line in acct.stdout.splitlines()]
    headers = lines[0]
    rows = [dict(zip(headers, line)) for line in lines[1:] if line[0] in jobs]
    if len(rows) != len(jobs):
        raise RuntimeError('incomplete scheduler accounting; do not infer zero cost')
    allocated = sum(int(row['CPUTimeRAW']) for row in rows)
    write_new(out / 'cost.json', {'jobs': rows, 'sacct': acct.stdout,
                                'allocated_cpu_seconds': allocated, 'gpu_seconds': 0,
                                'local_preparation_cost': 'not fully timed'})
    result = read_json(out / 'collection.json')
    with (out / 'REPORT.md').open('a') as f:
        f.write(f'\nRecorded allocation including input-order failure: **{allocated} CPU-seconds**, zero GPU.\n')
    with Path(a.vault).open('a') as f:
        f.write(f'\n\n## Native interaction job {a.job}: terminal collection\n\n'
                f"Status: **{result['status']}**. Allocation including initial input-order failure: "
                f'**{allocated} CPU-seconds**, zero GPU. Baseline/default unchanged. '
                f"Exact attribution to the original partition supported: {result['attribution_to_original_partition_supported']}. "
                'No affinity score or classification. Final report and all native execution '
                f'inventories: `{out / "REPORT.md"}` and `collection.json`.\n')
    write_new(out / 'completion.json', {'status': result['status'], 'report': record(out / 'REPORT.md'),
                                      'collection': record(out / 'collection.json'),
                                      'cost': record(out / 'cost.json')})
    print('Native EDA collection, report, allocation receipt and vault entry recorded.', flush=True)


if __name__ == '__main__':
    main()
