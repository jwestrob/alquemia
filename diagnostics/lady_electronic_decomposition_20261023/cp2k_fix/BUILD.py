#!/usr/bin/env python3
"""Bounded isolated CP2K build; performs zero molecular evaluations."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    base = Path(sys.argv[1]).resolve()
    script = Path(__file__).with_name('build_inside.sh').resolve()
    start = time.monotonic()
    result = {
        'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'job_id': os.environ['SLURM_JOB_ID'],
        'allocated_cpus': int(os.environ['SLURM_CPUS_PER_TASK']),
        'scope': 'Software build and version/loader inspection only; zero molecular calls.',
        'process_limit_seconds': 7200,
        'commands': [], 'status': 'running',
    }
    receipt = base / 'BUILD_RESULT.json'

    def save():
        receipt.write_text(json.dumps(result, indent=2) + '\n')

    env = os.environ.copy()
    env.update(APPTAINER_CACHEDIR=str(base / 'cache'),
               APPTAINER_TMPDIR=str(base / 'tmp'),
               APPTAINER_BUILD_PROCS='24',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for name in ('cache', 'tmp'):
        (base / name).mkdir(exist_ok=True)

    def run(argv, log):
        command = {'argv': argv, 'log': log, 'started_elapsed': time.monotonic() - start}
        result['commands'].append(command)
        save()
        with (base / log).open('x') as stream:
            p = subprocess.run(argv, env=env, stdout=stream, stderr=subprocess.STDOUT,
                               timeout=max(1, 7200 - (time.monotonic() - start)))
        command['returncode'] = p.returncode
        save()
        p.check_returncode()

    try:
        assert result['allocated_cpus'] >= 24
        assert sha(base/'cp2k-2026.2/src/atom_types.F') == '9adc98b30524851370c17aad8f4a9b6a95e1ca8b0ae50c21c479029d2e55c60a'
        result['driver_sha256'] = sha(Path(__file__))
        result['inside_script_sha256'] = sha(script)
        image = base / 'cp2k-2026.2-mpich-generic.sif'
        image_ref = 'docker://cp2k/cp2k@sha256:28e35f63a5253b95cc0f50e786370669185ad9f58eacacf33170cb1257216f36'
        result['image_reference'] = image_ref
        assert not image.exists(), 'This attempt must not silently reuse an unrecorded image.'
        run(['/usr/bin/apptainer', 'pull', str(image), image_ref], 'pull.log')
        result['image_sha256'] = sha(image)
        run(['/usr/bin/apptainer', 'inspect', '--json', str(image)], 'image_metadata.json')
        run(['/usr/bin/apptainer', 'exec', '--cleanenv', '--bind', str(base)+':/work',
             '--bind', str(script)+':/build_inside.sh:ro', '--pwd', '/work', str(image),
             '/bin/bash', '/build_inside.sh'], 'configure_build.log')
        binary = base/'build/bin/cp2k.psmp'
        result['binary_sha256'] = sha(binary)
        result['binary_path'] = str(binary)
        result['version'] = (base/'version.txt').read_text()
        result['status'] = 'built_version_and_loader_checks_passed'
        result['scientific_qualification'] = 'pending_exact_input_and_pilot_checks'
    except BaseException as exc:
        result['status'] = 'failed'
        result['error'] = repr(exc)
        raise
    finally:
        result['elapsed_seconds'] = time.monotonic() - start
        result['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    main()
