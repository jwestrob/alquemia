#!/usr/bin/env python3
"""Saved-only comparison of two native CP2K outputs with the PySCF ALMO test.

Root supplies actual exit codes and the exact inputs used by its native runner.
This script does not execute CP2K or recompute any electronic state. Native
termination and independently printed stages gate all energy differences.
"""
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import sys
import traceback

HERE = Path(__file__).resolve().parent
NUM = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def dec(s):
    d = Decimal(str(s).replace('D', 'E').replace('d', 'e'))
    if not d.is_finite():
        raise ValueError('Nonfinite value')
    return d


def need(ok, label):
    if not ok:
        raise ValueError(label)


def parse_native(path, exit_code):
    need(exit_code == 0, 'Native exit code was not zero')
    lines = Path(path).read_text().splitlines()
    evidence = {'output': pin(path), 'exit_code': exit_code}

    def rows(pattern, lo=0, hi=None):
        return [{'line': i+1, 'raw': line, 'tokens': list(m.groups())}
                for i, line in enumerate(lines) if i >= lo and (hi is None or i < hi)
                if (m := re.fullmatch(pattern, line))]

    def one(name, pattern, lo=0, hi=None):
        hits = rows(pattern, lo, hi)
        evidence[name] = hits
        need(len(hits) == 1, name+': missing/ambiguous record')
        return hits[0]

    start = one('start', r'.*PROGRAM STARTED AT\s+(.+)')
    end = one('end', r'.*PROGRAM ENDED AT\s+(.+)')
    bad = [line for line in lines if re.search(
        r'\bABORT\b|SCF run NOT converged|not converged!|FAILED to converge|'
        r'Leaving inner SCF loop|EPS_ERROR_EARLY is on|MAX_ITER_EARLY is on|'
        r'Read restart ALMO|ASPC order:', line, re.I)]
    need(not bad, 'Native failed/early-stop/restart marker')
    R = one('R', rf'\s*Single-molecule energy:\s+({NUM})\s*')
    F = one('F', rf'\s*Energy of the initial guess:\s+({NUM})\s*')
    block = one('block', r'\s*-+ Optimization of block-diagonal ALMOs -+\s*')
    full = one('full', r'\s*-+ Optimization of fully delocalized MOs -+\s*')
    U = one('U', rf'\s*ENERGY\| Total FORCE_EVAL \(\s*QS\s*\) energy \[a\.u\.\]\s+({NUM})\s*')
    stage = [x['line'] for x in (start, R, F, block, full, U, end)]
    need(all(b > a for a, b in zip(stage, stage[1:])), 'Native stage order')
    need(one('deloc', r'\s*delocalization:\s+(\S+)\s*')['tokens'] == ['FULL_SCF'], 'Not FULL_SCF')
    need(one('domains', r'\s*Total fragments:\s+(\d+)\s*')['tokens'] == ['2'], 'Not two domains')
    targets = rows(rf'\s*target error:\s+({NUM})\s*')
    need(len(targets) == 2 and all(dec(r['tokens'][0]) == Decimal('1e-8') for r in targets),
         'Changed native ALMO convergence targets')
    heads = rows(r'\s*-\s*MOLECULAR GUESS:\s*FRAGMENT\s+(\d+)\s+OUT OF\s+(\d+)\s*-\s*')
    need([r['tokens'] for r in heads] == [['1', '2'], ['2', '2']], 'Fragment sequence')
    fragE = []
    for i, head in enumerate(heads):
        lo = head['line']
        hi = heads[i+1]['line']-1 if i == 0 else R['line']-1
        cm = one(f'fragment_{i}_state', r'\s*-\s*CHARGE\s+(-?\d+)\s+MULTIPLICITY\s+(\d+)\s*-\s*', lo, hi)
        need(cm['tokens'] == ['0', '1'], 'Fragment charge/spin')
        membership = rows(r'\s*-\s+(\d+)\s+([A-Za-z][A-Za-z0-9_]*)\s*-\s*', lo, hi)
        need([[int(r['tokens'][0]), r['tokens'][1]] for r in membership] ==
             [[3*i+1, 'O'], [3*i+2, 'H'], [3*i+3, 'H']], 'Fragment atom membership')
        conv = one(f'fragment_{i}_converged', r'\s*\*\*\* SCF run converged in\s+(\d+) steps \*\*\*\s*', lo, hi)
        e = one(f'fragment_{i}_energy', rf'\s*Total energy:\s+({NUM})\s*', conv['line'], hi)
        fragE.append(dec(e['tokens'][0]))
        trace = rows(rf'\s*(\d+)\s+(OT\s+(?:CG|LS|SD|DIIS))\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s*', lo, conv['line']-1)
        need(trace and trace[-1]['tokens'][0] == conv['tokens'][0], 'Missing final ordinary SCF row')
        need(Decimal(0) <= dec(trace[-1]['tokens'][4]) <= Decimal('1e-10'), 'Fragment residual')
        evidence[f'fragment_{i}_final_row'] = trace[-1]
    table = one('domain_table', r'\s*Fragment\s+Basis Set\s+Occupied\s+Virtual\s+Charge\s+Deloc Neig\s*')
    dom = rows(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+)\s+(ALL|NONE|N/A|\d+)\s*', table['line'], table['line']+5)
    need([r['tokens'] for r in dom] == [['1','18','5','13','0','ALL'], ['2','18','5','13','0','ALL']],
         'AO/occupation/domain identity')
    last = {}
    for name, lo, hi, allowed in [('block', block['line'], full['line']-1, {'DIIS'}),
                                   ('full', full['line'], U['line']-1, {'CG','LS'})]:
        trace = rows(rf'\s*ALMO SCF (DIIS|CG|LS)\s+(\d+)\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s*', lo, hi)
        possible = [j+1 for j in range(lo, hi) if re.match(r'\s*ALMO SCF\s', lines[j])]
        need(trace and [r['line'] for r in trace] == possible, 'Missing/unknown ALMO trace')
        need(all(r['tokens'][0] in allowed for r in trace), 'Changed optimizer')
        steps = [int(r['tokens'][1]) for r in trace]
        need(steps[0] in (0,1) and all(b == a+1 for a,b in zip(steps,steps[1:])) and steps[-1] <= 300,
             'ALMO iteration sequence')
        need(Decimal(0) <= dec(trace[-1]['tokens'][4]) <= Decimal('1e-8'), 'ALMO final residual')
        last[name] = trace[-1]
        evidence[name+'_last'] = trace[-1]
    e = {'fragment_1':fragE[0], 'fragment_2':fragE[1], 'R':dec(R['tokens'][0]),
         'F':dec(F['tokens'][0]), 'P':dec(last['block']['tokens'][2]), 'U':dec(U['tokens'][0])}
    need(abs(sum(fragE)-e['R']) <= Decimal('1e-8'), 'Native fragment reference sum')
    need(abs(dec(last['full']['tokens'][2])-e['U']) <= Decimal('1e-8'), 'Full trace/final energy mismatch')
    need(e['U'] <= e['P']+Decimal('1e-8') and e['P'] <= e['F']+Decimal('1e-8'), 'Native variational order')
    c = {'FRZ':e['F']-e['R'], 'POL':e['P']-e['F'], 'CT':e['U']-e['P'], 'interaction':e['U']-e['R']}
    evidence['closure_hartree'] = str(c['FRZ']+c['POL']+c['CT']-c['interaction'])
    return e, c, evidence


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pyscf', type=Path, required=True)
    for name in ('coarse', 'fine'):
        ap.add_argument('--'+name, type=Path, required=True)
        ap.add_argument('--'+name+'-input', type=Path, required=True)
        ap.add_argument('--'+name+'-exit-code', type=int, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    need(not a.output.exists(), 'Refuse to overwrite comparison')
    r = {'status':'failed_independent_comparison', 'implementation':pin(__file__),
         'gates':{}, 'error':None, 'scope':'Independent implementation test only; no quartet scientific result.'}
    try:
        manifest = json.loads((HERE/'cp2k/MANIFEST.json').read_text())
        r['manifest'] = pin(HERE/'cp2k/MANIFEST.json')
        for name, cut in [('coarse',400), ('fine',600)]:
            actual = pin(getattr(a,name+'_input'))
            expected = manifest['inputs'][f'water_pbe_{cut}.inp']
            need(actual['sha256'] == expected['sha256'], name+' input differs from frozen deck')
            r[name+'_input'] = actual
        py = json.loads(a.pyscf.read_text())
        r['pyscf'] = pin(a.pyscf)
        need(py['status'] == 'passed_pyscf_benchmark' and py['method']['xc'] == 'PBE', 'PySCF benchmark failed/wrong method')
        # Root may use distinct immutable snapshot directories. The scientific
        # fixture must be byte-identical; an absolute source path is not identity.
        need(py['fixture']['sha256'] == pin(HERE/'FIXTURE.json')['sha256'], 'PySCF fixture differs')
        need(all(x['passed'] for x in py['checks']), 'Failed PySCF prerequisite')
        energies, components = {}, {}
        for name in ('coarse','fine'):
            energies[name], components[name], r[name+'_native'] = parse_native(
                getattr(a,name),getattr(a,name+'_exit_code'))
        pe = {k:dec(v) for k,v in py['energies_hartree'].items()}
        pc = {k:dec(v) for k,v in py['components_hartree'].items()}
        r['absolute_energy_differences_fine_minus_pyscf_hartree'] = {k:str(energies['fine'][k]-pe[k]) for k in pe}
        r['absolute_energy_resolution_changes_hartree'] = {k:str(energies['fine'][k]-energies['coarse'][k]) for k in pe}
        r['compartments'] = {}
        for k in ('FRZ','POL','CT','interaction'):
            error = components['fine'][k]-pc[k]
            grid = components['fine'][k]-components['coarse'][k]
            r['compartments'][k] = {'pyscf':str(pc[k]), 'cp2k_400':str(components['coarse'][k]),
                'cp2k_600':str(components['fine'][k]), 'fine_minus_pyscf':str(error),
                'fine_minus_coarse':str(grid), 'threshold_hartree':'1e-5'}
            r['gates'][k+'_agreement'] = abs(error) <= Decimal('1e-5')
            r['gates'][k+'_resolution_change'] = abs(grid) <= Decimal('1e-5')
        r['limits'] = ['The two-resolution difference is a measured sensitivity, not a rigorous error bound.',
                      'Absolute energies are reported as backend diagnostics, not required to be bit-identical.',
                      'Closure is a bookkeeping identity, independent of the native convergence gates.',
                      'Passing PBE/no-ECP water does not qualify PBE0, Dy55, charged fragments or the quartet.',
                      'Actual runner input/output/binary/exit-code linkage remains root runtime provenance.']
        if all(r['gates'].values()):
            r['status'] = 'passed_independent_cp2k_pbe_comparison'
    except Exception:
        r['error'] = traceback.format_exc()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open('x') as f:
        json.dump(r,f,indent=2);f.write('\n')
    return 0 if r['status'].startswith('passed_') else 1


if __name__ == '__main__':
    sys.exit(main())
