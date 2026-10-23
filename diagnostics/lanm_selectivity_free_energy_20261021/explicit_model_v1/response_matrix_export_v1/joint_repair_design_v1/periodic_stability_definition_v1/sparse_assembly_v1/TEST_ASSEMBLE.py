"""Finite artificial assembly checks; no molecular fixture or backend call."""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import ASSEMBLE as A


def run():
    checks = []
    def check(name, condition):
        checks.append({'name': name, 'passed': bool(condition)})
        if not condition:
            raise AssertionError(name)
    xyz = np.array([[0., 0., 0.], [2., 1., .3], [-.4, 2.1, .7]])
    types, alpha, thole = [1, 2, 2], np.array([.4, .6, .3]), np.array([.39, .4, .4])
    settings = dict(cutoff=9., ewald=.4, max_assembly_bytes=10**6,
                    pair_overrides={(1, 2): .31}, stiffness={0: .2}, shift=.01)
    lower, info = A.assemble(xyz, types, alpha, thole, [44.]*3, **settings)
    dense = lower.toarray(); dense = dense+dense.T-np.diag(np.diag(dense))
    reference = np.diag(np.repeat(1-alpha*4*.4**3/(3*np.sqrt(np.pi)), 3))
    reference[:3, :3] += np.eye(3)*.2
    reference -= np.eye(9)*.01
    original = A.parent()
    for i in range(3):
        for j in range(i+1, 3):
            gamma = .31 if types[i] != types[j] else .4
            t = original.tensor((xyz[j]-xyz[i])[None, :], alpha[i:i+1], alpha[j:j+1], np.array([gamma]), .4)[0]
            reference[3*i:3*i+3, 3*j:3*j+3] = -t
            reference[3*j:3*j+3, 3*i:3*i+3] = -t.T
    check('direct_dense_entries_and_orientation', np.allclose(dense, reference, rtol=1e-13, atol=1e-13))
    check('triangle_and_pair_counts', lower.nnz == 36 and info['pairs'] == 3 and np.all(np.triu(lower.toarray(), 1) == 0))
    eig = np.linalg.eigvalsh(dense+.01*np.eye(9))
    check('upper_spectral_bound', info['unshifted_site_upper'] >= eig[-1]-1e-13)
    check('lower_spectral_bound', info['unshifted_site_lower'] <= eig[0]+1e-13)
    shifted = xyz+np.array([[44., -88., 132.], [-44., 44., 0.], [0., 88., -44.]])
    image, _ = A.assemble(shifted, types, alpha, thole, [44.]*3, **settings)
    check('periodic_reimage', np.allclose(image.toarray(), lower.toarray(), rtol=1e-13, atol=1e-13))
    perm = np.array([2, 0, 1]); rotation = xyz[:, perm]
    rot, _ = A.assemble(rotation, types, alpha, thole, [44.]*3, **settings)
    fullrot = rot.toarray(); fullrot += fullrot.T-np.diag(np.diag(fullrot))
    idx = np.concatenate([3*i+perm for i in range(3)])
    check('axis_rotation', np.allclose(fullrot, dense[np.ix_(idx,idx)], rtol=1e-13, atol=1e-13))
    for label, changes, aa in [
            ('zero_alpha', {}, [0., .6, .3]),
            ('nonfinite_alpha', {}, [np.nan, .6, .3]),
            ('negative_shift', {'shift': -.1}, alpha),
            ('budget', {'max_assembly_bytes': 10}, alpha),
            ('scaling_exceptions', {'mutual_scales_all_one': False}, alpha),
            ('duplicate_override', {'pair_overrides': {(1,2): .31, (2,1): .31}}, alpha),
            ('invalid_beta', {'stiffness': {0: -1.}}, alpha),
            ('cutoff_halfbox', {'cutoff': 22.}, alpha)]:
        try:
            A.assemble(xyz, types, aa, thole, [44.]*3, **(settings | changes))
        except (ValueError, MemoryError):
            check(label+'_rejected', True)
        else:
            check(label+'_rejected', False)
    return checks


if __name__ == '__main__':
    out = Path(sys.argv[1]); wall, cpu = time.perf_counter(), time.process_time()
    if out.exists():
        raise ValueError('immutable output already exists')
    result = {'status': 'failed', 'source_sha256': hashlib.sha256(Path(A.__file__).read_bytes()).hexdigest(),
              'controls_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'new_native_calls': 0, 'physical_box_calls': 0, 'backend_factorizations': 0}
    try:
        result['checks'] = run(); result['status'] = 'passed'
    finally:
        result['cost'] = {'wall_seconds': time.perf_counter()-wall, 'CPU_seconds': time.process_time()-cpu, 'numerical_threads': 1}
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open('x') as f:
            json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps(result))
