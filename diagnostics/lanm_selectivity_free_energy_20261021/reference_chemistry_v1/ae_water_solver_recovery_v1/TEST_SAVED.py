#!/usr/bin/env python3
"""Narrow saved controls; no native engine entry points are called."""
import argparse
import ast
import importlib.util
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); wall, cpu = time.monotonic(), time.process_time()
    source = Path(__file__).with_name('RUN.py')
    ast.parse(source.read_text())
    spec = importlib.util.spec_from_file_location('recovery', source)
    r = importlib.util.module_from_spec(spec); spec.loader.exec_module(r)
    args.out.mkdir(parents=True, exist_ok=False)
    cp = args.out/'CONFIG_CONTROL.json'; r.write(cp, r.default_config())
    cfg, native, panel, cell = r.load(cp)
    checks = []
    def check(ok, label):
        r.check(ok, label); checks.append(label)
    base, changed = r.recovery_deck(native, panel, cell, 4, 3000)
    check(base == r.verify(cfg['original_input']).read_text(), 'Original17260 deck exact')
    check(changed.replace(r.INSERT, r.ORIGINAL) == base and changed.count('QuadRegionStart') == 1,
          'Exactly one solver setting added')
    check('Guess PModel\n MaxIter 300' in changed and 'VeryTightSCF' in changed and '* xyzfile 3 6 geometry.xyz' in changed,
          'Guess cap tolerance and sextet preserved')
    b2, c2 = r.recovery_deck(native, panel, cell, 4, 6000)
    check(c2.replace('%maxcore 6000', '%maxcore 3000') == changed,
          'Actual-memory scaling changes only MaxCore')
    failed_text = r.verify(cfg['original_output']).read_text()
    failed = r.scf_evidence(failed_text)
    check(failed['macro_rows'] == 300 and failed['TRAH_used'], 'Actual300macro failedTRAH parsed')
    check(failed['last_macro']['orbital_gradient_norm'] == 9.875006e-5,
          'Final macro residual retained instead of smaller final micro residual')
    check(not failed['active_scf_convergence'] and not failed['native_scf_converged'], 'Actualfailure rejected')
    check(failed['declared_TolG'] == 2e-6 and failed['printed_NR_entry_threshold'] == 1e-4,
          'Independent TolG and NR-entry thresholds parsed')
    check(failed['printed_TRAH_convergence_threshold'] == 2e-6 and
          failed['active_orbital_gradient_threshold'] == 2e-6, 'ActualactiveTRAHthreshold separately admitted')
    la = r.WORK/'pilot_v1/tasks/ae_pbe0__water_central_La/execution/endpoint.out'
    ordinary = r.scf_evidence(la.read_text())
    check(not ordinary['TRAH_used'] and ordinary['active_scf_convergence'], 'ActualordinarySOSCF convergence admitted')
    check(ordinary['final_printed_orbital_gradient'] == 4.9145e-7,
          'Actualfinal orbitalgradient distinct from staleDIIS')
    check('cannot be attributed' in ordinary['solver_change_causality'], 'NoTRAH means no controlcausality claim')
    # Explicitly synthetic controls exercise only parser/admission arithmetic.
    fixed = ('Orbital Gradient TolG .... 2.000e-06\n'
             'Converg. threshold (grad. norm) .... 2.000e-06\n'
             'NR start threshold (gradient norm) .... 1.000e-06\n')
    good = fixed+' 12 -12520.123456 1.0e-06 0.013 (TRAH MAcro) No\nSCF CONVERGED AFTER 12 CYCLES\n'
    check(r.scf_evidence(good)['active_scf_convergence'], 'Syntheticconverged macro passes')
    bad = good.replace('1.0e-06 0.013', '9.0e-05 0.013')+' 12 dE -1e-10 1e-9 (TRAH MIcro)\n'
    check(not r.scf_evidence(bad)['active_scf_convergence'], 'Synthetic smallmicro cannot rescue largemacro')
    check(not r.scf_evidence(good.replace('SCF CONVERGED AFTER', 'SCF NOT CONVERGED AFTER'))['active_scf_convergence'],
          'Smallresidual without nativeconvergence cannot pass')
    result = {'status': 'saved_controls_passed', 'implementation': r.pin(source), 'test_source': r.pin(__file__),
              'checks': checks, 'count': len(checks), 'actual_failed_output': cfg['original_output'],
              'actual_converged_format_output': r.pin(la), 'synthetic_results_are_controls_only': True,
              'native_calls': 0, 'threads': 1, 'wall_seconds': time.monotonic()-wall,
              'CPU_seconds': time.process_time()-cpu}
    print(r.write(args.out/'RESULT.json', result))


if __name__ == '__main__':
    main()
