"""Known-interaction and energy-origin controls; no molecular calculations."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import time
import ARITHMETIC as a


def main(out):
    start = time.monotonic(); cpu = time.process_time(); checks = []
    # Distinct element energy origins must cancel, leaving pair 3 and cooperative 4.
    ions = {'La': -100., 'Dy': -200.}
    pair = {'La': -2., 'Dy': 1.}
    clusters = {}; pairs = {}
    for g in a.GIDS:
        n = int(g[4]); metal = g.rsplit('_', 1)[1]
        waters = [10. + i / 2 for i in range(n)]
        pairs[g] = {i: ions[metal] + waters[i] + pair[metal] for i in range(n)}
        extra = 4. if n == 9 and metal == 'Dy' else 0.
        clusters[g] = ions[metal] + sum(waters) + n * pair[metal] + extra
    base = a.decompose(clusters, pairs, ions)
    assert base['pair_contrast'] == 3 and base['cooperative_remainder'] == 4 and base['cluster_contrast'] == 7
    checks.append('known pair=3, cooperative=4, total=7 with distinct ion origins')
    shifted_i = dict(ions); shifted_i['Dy'] += 777.
    shifted_c = dict(clusters); shifted_p = copy.deepcopy(pairs)
    for g in a.GIDS:
        if g.endswith('_Dy'):
            shifted_c[g] += 777.
            shifted_p[g] = {i: v + 777. for i, v in shifted_p[g].items()}
    shifted = a.decompose(shifted_c, shifted_p, shifted_i)
    for k in ('pair_contrast', 'cluster_contrast', 'cooperative_remainder'):
        assert shifted[k] == base[k]
    checks.append('arbitrary element-specific ion origin cancels')
    c = dict(clusters); p = copy.deepcopy(pairs)
    for g in a.GIDS:
        c[g] += 55. * int(g[4]); p[g] = {i: v + 55. for i, v in p[g].items()}
    water = a.decompose(c, p, ions)
    assert all(water[k] == base[k] for k in ('pair_contrast', 'cluster_contrast', 'cooperative_remainder'))
    checks.append('water energy origin cancels between metals')
    for kind in ('missing_pair', 'extra_pair', 'missing_ion', 'nonfinite', 'boolean'):
        p = copy.deepcopy(pairs); ii = dict(ions)
        if kind == 'missing_pair': del p[a.GIDS[0]][0]
        if kind == 'extra_pair': p[a.GIDS[0]][8] = 0.
        if kind == 'missing_ion': del ii['Dy']
        if kind == 'nonfinite': p[a.GIDS[0]][0] = float('nan')
        if kind == 'boolean': ii['La'] = True
        try: a.decompose(clusters, p, ii)
        except ValueError: checks.append('rejects ' + kind)
        else: raise AssertionError(kind)
    e = a.compare(base, dict(base, pair_contrast=5., cooperative_remainder=1., cluster_contrast=6.))
    assert e['pair_contrast'] == 2 and e['cooperative_remainder'] == -3 and e['cluster_contrast'] == -1 and e['closure_residual'] == 0
    checks.append('signed error decomposition with opposing contributions')
    def pin(p): return {'path': str(p.resolve()), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
    result = {'status': 'passed', 'checks': checks, 'checks_passed': len(checks),
              'arithmetic': pin(Path(a.__file__)), 'controls': pin(Path(__file__)),
              'new_native_calls': 0, 'wall_seconds': time.monotonic()-start, 'CPU_seconds': time.process_time()-cpu}
    Path(out).mkdir(parents=True, exist_ok=False)
    (Path(out)/'RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'checks_passed': len(checks)}))


if __name__ == '__main__': main(sys.argv[1])
