"""Declared current-SCF readback; no electronic calculation or native imports."""
from decimal import Decimal
import re

N = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
LIMITS = {'Energy change': '1e-11', 'MAX-Density change': '1e-8',
          'RMS-Density change': '1e-9', 'DIIS Error': '1e-8',
          'Orbital Gradient': '1e-8', 'Orbital Rotation': '1e-8'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def dec(token):
    x = Decimal(token.replace('D', 'E').replace('d', 'e'))
    require(x.is_finite(), 'Nonfinite SCF token')
    return x


def print_join(a, b):
    x, y = dec(a), dec(b)
    allowance = (Decimal(10)**x.as_tuple().exponent + Decimal(10)**y.as_tuple().exponent)/2
    return abs(abs(x)-abs(y)) <= allowance


def native(text):
    return text.split('****END OF INPUT****')[-1]


def settings(text):
    t = native(text)
    keys = ('IntAcc', 'RadialGrid', 'AngularGrid', 'GridPruning', 'BFCut',
            'WCut', 'Thresh', 'TCut', 'TPreSel')
    result = {k: [' '.join(line.split()) for line in t.splitlines()
                  if re.search(r'\b'+k+r'\b', line)] for k in keys}
    require(all(result.values()), 'Missing actual grid/integral setting')
    return result


def occupations(text):
    sections = native(text).split('ORBITAL ENERGIES')
    require(len(sections) == 2, 'Missing/duplicate restricted orbital energies')
    rows = []
    for line in sections[1].splitlines():
        match = re.fullmatch(r'\s*(\d+)\s+('+N+r')\s+('+N+r')\s+('+N+r')\s*', line)
        if match:
            rows.append([int(match[1]), match[2]])
        elif rows:
            break
    require(rows and [x[0] for x in rows] == list(range(len(rows))), 'Orbital occupation inventory')
    require(all(dec(x[1]) in (0, 2) for x in rows), 'Changed restricted occupations')
    return rows


def current(text):
    t = native(text)
    modes = re.findall(r'Convergence Check Mode\s+ConvCheckMode\s*\.{2,}\s*([^\n]+)', t)
    require(modes == ['All-Criteria'], 'Actual mode must be All-Criteria')
    forced = re.findall(r'Convergence forced\s+ConvForced\s*\.{2,}\s*(\d+)', t)
    require(forced == ['1'], 'Actual ConvForced differs')
    iterations = re.findall(r'Maximum # iterations\s+MaxIter\s*\.{2,}\s*(\d+)', t)
    require(iterations == ['300'], 'Actual MaxIter differs')
    for key, value in [('TolE', '1e-11'), ('TolG', '1e-8'), ('TolX', '1e-8'), ('TolErr', '1e-8')]:
        values = re.findall(r'\b'+key+r'\s*\.{2,}\s*('+N+r')', t)
        require(len(values) == 1 and dec(values[0]) == dec(value), 'Actual threshold differs: '+key)
    blocks = t.split('SCF CONVERGENCE')
    require(len(blocks) == 2, 'Missing/duplicate terminal SCF summary')
    summary = blocks[1].split('ORBITAL ENERGIES')[0]
    rows = re.findall(r'^\s*Last\s+([^\n]*?)\s*\.{2,}\s*('+N+r')\s+Tolerance\s*:\s*('+N+r')\s*$', summary, re.M)
    fields = {name.strip(): {'value_token': value, 'tolerance_token': tolerance} for name, value, tolerance in rows}
    require(len(fields) == len(rows) and set(fields) == set(LIMITS), 'Terminal residual denominator differs')
    for name, field in fields.items():
        require(dec(field['tolerance_token']) == dec(LIMITS[name]), 'Summary tolerance differs: '+name)
        require(name == 'Energy change' or dec(field['value_token']) >= 0, 'Negative residual: '+name)
    segments = t.split('ORCA LEAN-SCF')
    require(len(segments) == 2, 'Unrecognized or duplicate solver entry')
    trace = segments[1].split('SCF CONVERGED AFTER')
    require(len(trace) == 2, 'Missing/duplicate SCF convergence marker')
    final_count = re.match(r'\s*(\d+)\s+CYCLES', trace[1])
    require(final_count is not None, 'SCF iteration count missing')
    algorithm, last = None, None
    for line in trace[0].splitlines():
        if 'D-I-I-S' in line:
            algorithm = 'DIIS'
        elif 'S-O-S-C-F' in line:
            algorithm = 'SOSCF'
        elif re.search(r'T-?R-?A-?H', line):
            raise ValueError('Unrecognized terminal solver path: TRAH')
        tokens = line.split()
        if tokens and tokens[0].isdigit():
            require(algorithm in ('DIIS', 'SOSCF') and len(tokens) == (8 if algorithm == 'DIIS' else 7),
                    'Untraceable iteration row')
            for token in tokens[1:]:
                dec(token)
            last = {'algorithm': algorithm, 'iteration': int(tokens[0]), 'delta_E': tokens[2],
                    'RMSDP': tokens[3], 'MaxDP': tokens[4],
                    'MaxGrad' if algorithm == 'SOSCF' else 'DIISErr': tokens[5], 'raw': line}
    require(last is not None and last['iteration'] == int(final_count[1]) and last['algorithm'] == algorithm,
            'Final iteration is not the converged state')
    joins = [('Energy change', 'delta_E'), ('MAX-Density change', 'MaxDP'), ('RMS-Density change', 'RMSDP')]
    joins += [('Orbital Gradient', 'MaxGrad')] if algorithm == 'SOSCF' else [('DIIS Error', 'DIISErr')]
    for name, column in joins:
        require(print_join(fields[name]['value_token'], last[column]), 'Current residual/last-row mismatch: '+name)
    required = [x[0] for x in joins] + (['Orbital Rotation'] if algorithm == 'SOSCF' else [])
    checks = {name: abs(dec(fields[name]['value_token'])) <= dec(LIMITS[name]) for name in required}
    return {'status': 'passed_current_scf_residuals' if all(checks.values()) else 'failed_current_scf_residuals',
            'algorithm': algorithm, 'mode': modes[0], 'last_iteration': last, 'summary': fields,
            'required_checks': checks, 'all_required_passed': all(checks.values()),
            'rotation_source': 'Associated terminal summary; no per-iteration rotation column' if algorithm == 'SOSCF' else None,
            'DIIS_scope': 'Stale after SOSCF; descriptive only' if algorithm == 'SOSCF' else 'Current terminal DIIS row joined',
            'electronic_root_continuity_proved': False}
