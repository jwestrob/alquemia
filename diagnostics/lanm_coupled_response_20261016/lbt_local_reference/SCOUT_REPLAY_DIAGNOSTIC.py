#!/usr/bin/env python3
"""Read the executed scout and diagnose serialization differences; zero calls."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def differences(a, b, path='row'):
    if type(a) is not type(b):
        return [{'path': path, 'native_type': type(a).__name__, 'saved_type': type(b).__name__,
                 'native': repr(a), 'saved': repr(b)}]
    if isinstance(a, dict):
        result = []
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                result.append({'path': path + '.' + key, 'missing_key': True})
            else:
                result.extend(differences(a[key], b[key], path + '.' + key))
        return result
    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            return [{'path': path, 'different_lengths': [len(a), len(b)]}]
        return [item for i, (x, y) in enumerate(zip(a, b))
                for item in differences(x, y, path + '[' + str(i) + ']')]
    return [] if a == b else [{'path': path, 'native': repr(a), 'saved': repr(b)}]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    assert os.environ.get('SLURM_JOB_ID', '').isdigit()
    assert int(os.environ['SLURM_CPUS_PER_TASK']) == 1
    out = Path(args.output); out.mkdir(parents=True, exist_ok=False)
    cfg = json.loads(Path(args.config).read_text())
    source = cfg['implementation']
    assert pin(source['path']) == source
    spec = importlib.util.spec_from_file_location('actual_reference', source['path'])
    ref = importlib.util.module_from_spec(spec); spec.loader.exec_module(ref)
    _, model, _, collector, residual, _, _ = ref.inputs(cfg)
    old_path = Path(args.config).parent / 'execution/RESULT.json'
    old = json.loads(old_path.read_text())
    actual = ref.admit(ref.verify(old['manifest']), model, collector, residual)['7CCO_La']
    saved = old['rows']['7CCO_La']
    canonical = json.loads(json.dumps(actual, allow_nan=False))
    result = {'status': 'complete' if canonical == saved else 'unresolved',
              'job_id': os.environ['SLURM_JOB_ID'], 'source': source,
              'scout_result': pin(old_path), 'raw_equality': actual == saved,
              'JSON_roundtrip_equality': canonical == saved,
              'differences': differences(actual, saved),
              'canonical_differences': differences(canonical, saved),
              'new_molecular_calls': 0}
    (out / 'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
