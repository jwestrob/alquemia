#!/usr/bin/env python3
"""Reuse the saved image/source and supply missing build tools; no chemistry."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import signal
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
    prior = base.with_name('cp2k_build_v1')
    tools = base.with_name('cp2k_tools_v1')
    source = prior / 'cp2k-2026.2'
    image = prior / 'cp2k-2026.2-mpich-generic.sif'
    script = Path(__file__).with_name('build_inside_v2.sh').resolve()
    start = time.monotonic()
    result = {'job_id': os.environ['SLURM_JOB_ID'],
              'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'status': 'running', 'scope': 'Software build only; zero molecular calls.',
              'process_limit_seconds': 7200, 'commands': []}
    receipt = base/'BUILD_RESULT.json'

    def save():
        receipt.write_text(json.dumps(result, indent=2)+'\n')

    try:
        assert int(os.environ['SLURM_CPUS_PER_TASK']) >= 24
        prior_result = json.loads((prior/'BUILD_RESULT.json').read_text())
        result['image_sha256'] = sha(image)
        assert result['image_sha256'] == prior_result['image_sha256'] == '90c3a7a7368e8025689a3ea13ba773ecf832a94ea5378c984de1b861c9d8deab'
        result['atom_types_sha256'] = sha(source/'src/atom_types.F')
        assert result['atom_types_sha256'] == '9adc98b30524851370c17aad8f4a9b6a95e1ca8b0ae50c21c479029d2e55c60a'
        result['driver_sha256'] = sha(Path(__file__))
        result['inside_script_sha256'] = sha(script)
        result['tools_provenance'] = {'path': str(tools/'PROVENANCE.json'),
                                      'sha256': sha(tools/'PROVENANCE.json')}
        assert sha(tools/'cmake-3.28.3-linux-x86_64/bin/cmake') == 'b801b21dfdb444a91fb9f8ad3d39910f5471d4c32328ee8c9aa0352ef50154f0'
        result['pkgconf_provenance'] = json.loads((tools/'PKGCONF_PROVENANCE.json').read_text())
        argv = ['/usr/bin/apptainer','exec','--cleanenv',
                '--bind',str(base)+':/work', '--bind',str(source)+':/source:ro',
                '--bind',str(tools)+':/tools:ro','--bind',str(script)+':/build_inside.sh:ro',
                '--pwd','/work',str(image),'/bin/bash','/build_inside.sh']
        result['commands'].append({'argv': argv, 'log': 'configure_build.log'})
        save()
        with (base/'configure_build.log').open('x') as stream:
            proc = subprocess.Popen(argv, stdout=stream, stderr=subprocess.STDOUT,
                                    start_new_session=True)
            try:
                returncode = proc.wait(timeout=max(1,7200-(time.monotonic()-start)))
            except subprocess.TimeoutExpired:
                # Enforce only this predeclared software attempt's process cap.
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                raise
        result['commands'][-1]['returncode'] = returncode
        if returncode:
            raise RuntimeError(f'Container build returned {returncode}; see configure_build.log')
        result['binary_sha256'] = sha(base/'build/bin/cp2k.psmp')
        result['binary_path'] = str(base/'build/bin/cp2k.psmp')
        result['version'] = (base/'version.txt').read_text()
        result['status'] = 'built_version_and_loader_checks_passed'
        result['scientific_qualification'] = 'pending_exact_input_and_pilot_checks'
    except BaseException as exc:
        result['status'] = 'failed'
        result['error'] = repr(exc)
        raise
    finally:
        result['elapsed_seconds'] = time.monotonic()-start
        result['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    main()
