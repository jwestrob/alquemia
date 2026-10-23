#!/usr/bin/env python3
"""Freeze installed metadata and future unbound inputs; no solver/proposal call."""
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

import clarabel

HERE = Path(__file__).resolve().parent
BASE = next(p for p in HERE.parents if (p/'workspaces').is_dir() and (p/'diagnostics').is_dir())
D = BASE/'diagnostics/lanm_selectivity_free_energy_20261021'
W = BASE/'workspaces/lanm_selectivity_free_energy_20261021'
REL = Path('explicit_model_v1/joint_component_feasibility_v1/component_fit_v1')
WORK = W/REL/'conic_saved_comparison_v1'


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def write(name, value):
    with open(HERE/name, 'x') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    assert clarabel.__version__ == '0.11.1'
    settings = clarabel.DefaultSettings()
    settings.max_iter = 200; settings.time_limit = 60.; settings.max_threads = 1
    settings.direct_solve_method = 'qdldl'
    for key in ('tol_gap_abs', 'tol_gap_rel', 'tol_feas', 'tol_infeas_abs', 'tol_infeas_rel'):
        setattr(settings, key, 1e-10)
    effective = {k: getattr(settings, k) for k in dir(settings) if not k.startswith('_') and not callable(getattr(settings, k))}
    write('SETTINGS.json', effective)
    files = [pin(sys.executable), pin(Path(sys.prefix)/'pyvenv.cfg')]
    for name in ('clarabel', 'numpy', 'scipy'):
        dist = metadata.distribution(name)
        for entry in dist.files:
            path = str(entry)
            if path.endswith(('.dist-info/METADATA', '.dist-info/RECORD')) or path == name+'/__init__.py' or (name == 'clarabel' and path.endswith('.so')):
                files.append(pin(dist.locate_file(entry)))
    write('ENVIRONMENT.json', {'prefix': str(Path(sys.prefix).resolve()), 'python_version': sys.version,
          'versions': {name: metadata.version(name) for name in ('clarabel', 'numpy', 'scipy')}, 'files': files,
          'installation': pin(WORK/'environment_v1/INSTALL.json'), 'pip_report': pin(WORK/'environment_v1/PIP_REPORT.json'),
          'system_numpy_scipy_unchanged': True, 'metadata_only_no_solver': True})
    cfgpin = pin(W/REL/'run_v1/CONFIG_RUN.json')
    campaign = json.loads(Path(cfgpin['path']).read_text())
    inv = json.loads(Path(campaign['reuse_inventory']['path']).read_text())
    parent = json.loads(Path(campaign['parent_configuration']['path']).read_text())
    inputs = {'declaration': pin(HERE/'DECLARATION.md'), 'core': campaign['core'],
              'campaign_driver': campaign['implementation'], 'campaign_configuration': cfgpin,
              'campaign_fingerprint': '69d56ebc8c68adfd04603e2942408aa188a6db5eb594b7788e1ed5899c07c185',
              'campaign_declaration': campaign['declaration'], 'reuse_inventory': campaign['reuse_inventory'],
              'manifest': parent['manifest'], 'domain': parent['domain'], 'decomposition': campaign['decomposition'],
              'origin_key': inv['origin_key'], 'derivative_keys': inv['derivative_keys'], 'vectors': inv['vectors'],
              'mapping_review': pin(D/REL/'conic_formulation_review_v1/REVIEW.json'),
              'mapping': pin(D/REL/'conic_formulation_review_v1/MAPPING.json'),
              'expected_terminal_path': str(W/REL/'run_v1/execution/RESULT.json')}
    assert len(inputs['vectors']) == 25 and len(inputs['derivative_keys']) == 24
    write('INPUTS.json', inputs)
    write('CONFIG_UNBOUND.json', {'implementation': pin(HERE/'COMPARE.py'), 'inputs': pin(HERE/'INPUTS.json'),
          'environment': pin(HERE/'ENVIRONMENT.json'), 'settings': pin(HERE/'SETTINGS.json'),
          'terminal_campaign': None, 'source_review': None, 'root_release': None,
          'output_directory': str(WORK/'actual_v1'), 'one_attempt_marker': str(WORK/'ONE_SOLVER_ATTEMPT.json')})
    print(json.dumps({'metadata_prepared': True, 'solver_objects': 0, 'solver_calls': 0, 'native_calls': 0, 'campaign_outcomes_read': False}))


if __name__ == '__main__':
    main()
