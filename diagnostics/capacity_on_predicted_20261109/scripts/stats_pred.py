"""Statistics and frozen reading (PLAN.md sections 2-5). Usage: python -I stats_pred.py DIAG"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

DIAG = sys.argv[1]
RES = os.path.join(DIAG, 'results')
L = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    L.append(s)


def rd(f):
    return list(csv.DictReader(open(os.path.join(RES, f)), delimiter='\t'))


def num(v):
    return None if v in ('', 'None') else float(v)


def protein_vals(rs, f):
    """site -> uid -> entry -> protein(c95) (benchmark stats.py)."""
    a = defaultdict(list)
    for r in rs:
        a[(r['c95'], r['entry'], r['uid'])].append(f(r))
    b = defaultdict(list)
    for (c, e, u), v in a.items():
        b[(c, e)].append(np.mean(v))
    c_ = defaultdict(list)
    for (c, e), v in b.items():
        c_[c].append(np.mean(v))
    c30 = {r['c95']: r['c30'] for r in rs}
    return {c: (c30[c], float(np.mean(v))) for c, v in c_.items()}


def boot_est(pv_list):
    """Benchmark boot(): c30 cluster bootstrap, 2000 resamples, default_rng(42). Returns point, est (2000 x m)."""
    keys = sorted(pv_list[0])
    cl = sorted({pv_list[0][k][0] for k in keys})
    ci = {c: i for i, c in enumerate(cl)}
    idx = np.array([ci[pv_list[0][k][0]] for k in keys])
    V = np.array([[pv[k][1] for k in keys] for pv in pv_list])
    sums = np.zeros((len(pv_list), len(cl)))
    for j in range(len(pv_list)):
        np.add.at(sums[j], idx, V[j])
    cnt = np.bincount(idx, minlength=len(cl)).astype(float)
    rng = np.random.default_rng(42)
    draws = rng.integers(0, len(cl), size=(2000, len(cl)))
    W = np.stack([np.bincount(d, minlength=len(cl)) for d in draws])
    est = (W @ sums.T) / (W @ cnt)[:, None]
    return V.mean(1), est, len(keys), len(cl)


def fmt(pt, e):
    return f'{pt:.3f}[{np.percentile(e, 2.5):.3f},{np.percentile(e, 97.5):.3f}]'


S = rd('bench_sites.tsv')
for r in S:
    for k in ('r_v0', 'r_b2', 'an_true', 'site_plddt'):
        r[k] = num(r[k])
tiers = ('a_stripped', 'b_partner', 'c_pred', 'c2_pred_trim')
T = {t: [r for r in S if r['tier'] == t] for t in tiers}
p('sites per tier', {t: len(v) for t, v in T.items()}, 'proteins', {t: len({r['c95'] for r in v}) for t, v in T.items()})
p('stripped re-score matches benchmark site_ranks:', sum(r['bench_match'] == 'True' for r in T['a_stripped']), 'of', len(T['a_stripped']))
p('transfer rmsd (pred) median', np.median([float(r['rmsd']) for r in T['c_pred']]),
  '>2A', sum(float(r['rmsd']) > 2 for r in T['c_pred']), ' seq identity crystal->model median',
  np.median([float(r['ident']) for r in T['c_pred']]))


def hit(k, m):
    return lambda r: float(r[f'r_{m}'] is not None and r[f'r_{m}'] < k)


def table(rs, label):
    out = {}
    if not rs:
        p(label, 'no sites'); return out
    for k in (1, 3, 5):
        pt, est, npro, ncl = boot_est([protein_vals(rs, hit(k, 'v0')), protein_vals(rs, hit(k, 'b2'))])
        out[('v0', k)] = (pt[0], est[:, 0]); out[('b2', k)] = (pt[1], est[:, 1])
        p(f'{label} k={k} sites={len(rs)} proteins={npro} c30={ncl} v0={fmt(pt[0], est[:, 0])} b2={fmt(pt[1], est[:, 1])}')
    pt, est, *_ = boot_est([protein_vals(rs, lambda r: r['an_true']), protein_vals(rs, lambda r: float(r['an_true'] >= 4)),
                            protein_vals(rs, lambda r: float(r['an_true'] >= 3))])
    out['an'] = (pt[0], est[:, 0]); out['F4'] = (pt[1], est[:, 1]); out['F3'] = (pt[2], est[:, 2])
    p(f'{label} mean an_true={fmt(pt[0], est[:, 0])} F4={fmt(pt[1], est[:, 1])} F3={fmt(pt[2], est[:, 2])} '
      f'| site-level median an_true {np.median([r["an_true"] for r in rs])} F4 sites {sum(r["an_true"] >= 4 for r in rs)}/{len(rs)}')
    return out


p('\n== Three-way table (same sites; partner = subset with partner transfers)')
R = {t: table(T[t], t) for t in tiers}
pp = {r['c95'] for r in T['b_partner']}
p('\n== Restricted to proteins with a partner')
for t in ('a_stripped', 'c_pred'):
    table([r for r in T[t] if r['c95'] in pp], t + '|partner_proteins')

# frozen reading with paired ratio
pv_a = protein_vals(T['a_stripped'], lambda r: float(r['an_true'] >= 4))
pv_c = protein_vals(T['c_pred'], lambda r: float(r['an_true'] >= 4))
pv_r5 = protein_vals(T['c_pred'], hit(5, 'v0'))
pt, est, *_ = boot_est([pv_r5, pv_c, pv_a])
ratio = est[:, 1] / est[:, 2]
p('\n== FROZEN READING')
p(f'v0 recall@5 predicted = {fmt(pt[0], est[:, 0])} (threshold 0.45)')
p(f'F4 predicted = {fmt(pt[1], est[:, 1])}; F4 stripped = {fmt(pt[2], est[:, 2])}; ratio = {pt[1] / pt[2]:.3f}'
  f'[{np.percentile(ratio, 2.5):.3f},{np.percentile(ratio, 97.5):.3f}] (threshold 0.8)')
verdict = pt[0] >= 0.45 and pt[1] >= 0.8 * pt[2]
p('READING:', 'the count rule transfers to predicted structures' if verdict else 'does not transfer')
pt2, est2, *_ = boot_est([protein_vals(T['c_pred'], hit(5, 'v0')), protein_vals(T['a_stripped'], hit(5, 'v0'))])
d = est2[:, 0] - est2[:, 1]
p(f'paired v0 recall@5 predicted - stripped = {pt2[0] - pt2[1]:+.3f}[{np.percentile(d, 2.5):+.3f},{np.percentile(d, 97.5):+.3f}]')

# site-level paired an_true
pa = {(r['name'], r['sid']): r for r in T['a_stripped']}
dd = [r['an_true'] - pa[(r['name'], r['sid'])]['an_true'] for r in T['c_pred']]
p('site-level an_true(pred) - an_true(crystal): lower', sum(x < 0 for x in dd), 'equal', sum(x == 0 for x in dd),
  'higher', sum(x > 0 for x in dd), 'median', np.median(dd), 'mean', round(float(np.mean(dd)), 3))
dct = [(r['r_v0'] is not None and r['r_v0'] < 5, pa[(r['name'], r['sid'])]['r_v0'] is not None and pa[(r['name'], r['sid'])]['r_v0'] < 5) for r in T['c_pred']]
p('site-level hit@5 v0 (pred, crystal): both', sum(a and b for a, b in dct), 'pred only', sum(a and not b for a, b in dct),
  'crystal only', sum(b and not a for a, b in dct), 'neither', sum(not a and not b for a, b in dct))

p('\n== pLDDT strata (predicted; crystal shown for the same sites)')
for lab, f in (('<70', lambda r: r['site_plddt'] is not None and r['site_plddt'] < 70),
               ('>=70', lambda r: r['site_plddt'] is not None and r['site_plddt'] >= 70)):
    keys = {(r['name'], r['sid']) for r in T['c_pred'] if f(r)}
    table([r for r in T['c_pred'] if (r['name'], r['sid']) in keys], f'c_pred|pLDDT{lab}')
    table([r for r in T['a_stripped'] if (r['name'], r['sid']) in keys], f'a_stripped|pLDDT{lab}')
p('site pLDDT missing', sum(r['site_plddt'] is None for r in T['c_pred']), 'median site pLDDT',
  np.median([r['site_plddt'] for r in T['c_pred'] if r['site_plddt'] is not None]))

p('\n== False-positive burden (merged centres >8 A from every observed Ln)')
F = rd('bench_fp.tsv')
FT = defaultdict(dict)
for r in F:
    FT[r['tier']][r['name']] = r
for t in ('a_stripped', 'c_pred', 'c2_pred_trim'):
    rows = list(FT[t].values())
    for th in (4, 5):
        per = [float(r[f'off{th}_per100']) for r in rows]
        cnt = [int(r[f'off{th}']) for r in rows]
        p(f'{t} t>={th}: proteins {len(rows)} off-site total {sum(cnt)} per-protein median {np.median(cnt)} '
          f'per100 median {np.median(per):.3f} mean {np.mean(per):.3f} IQR [{np.percentile(per, 25):.3f},{np.percentile(per, 75):.3f}] '
          f'| all merged {sum(int(r[f"n{th}"]) for r in rows)} unmappable {sum(int(r[f"unmap{th}"]) for r in rows)} '
          f'| residues {sum(int(r["nres"]) for r in rows)}')
for th in (4, 5):
    dd = [float(FT['c_pred'][n][f'off{th}_per100']) - float(FT['a_stripped'][n][f'off{th}_per100']) for n in FT['a_stripped']]
    p(f'paired per100 pred - crystal t>={th}: median {np.median(dd):+.3f}, pred higher {sum(x > 0 for x in dd)}, '
      f'equal {sum(x == 0 for x in dd)}, lower {sum(x < 0 for x in dd)}')
    dd = [float(FT['c2_pred_trim'][n][f'off{th}_per100']) - float(FT['a_stripped'][n][f'off{th}_per100']) for n in FT['a_stripped']]
    p(f'paired per100 trimmed-pred - crystal t>={th}: median {np.median(dd):+.3f}, higher {sum(x > 0 for x in dd)}, '
      f'equal {sum(x == 0 for x in dd)}, lower {sum(x < 0 for x in dd)}')

p('\n== Capacity curve (v0 centres; Ln within 4 A; unmappable model centres excluded)')
cap = defaultdict(lambda: [0, 0, 0])
for r in rd('capacity.tsv'):
    c = cap[(r['tier'], int(r['n_anion_bin']))]
    lab = int(r['ln_within4'])
    if lab == -1:
        c[2] += 1
    else:
        c[0] += 1; c[1] += lab
for t in ('a_stripped', 'c_pred', 'c2_pred_trim'):
    p(t, ' '.join(f'{b}:{cap[(t, b)][1]}/{cap[(t, b)][0]}={cap[(t, b)][1] / max(1, cap[(t, b)][0]):.2f}' for b in range(9)))

p('\n== Reference proteins')
for r in rd('ref_sites.tsv'):
    p('\t'.join(r[k] for k in ('name', 'pdb', 'site', 'cls', 'n_don', 'types', 'ident', 'rmsd', 'site_plddt', 'model_plddt',
                               'x_an_true', 'm_an_true', 'x_r_v0', 'm_r_v0', 'x_r_b2', 'm_r_b2')))

p('\n== Candidates')
M = rd('cand_models.tsv')
for g in ('unannotated', 'comparator', 'TPR'):
    q = [r for r in M if r['group'] == g]
    med = lambda k: np.median([float(r[k]) for r in q if r[k] not in ('', 'None')]) if q else None  # noqa: E731
    p(f'{g}: n={len(q)} pLDDT med {med("mean_plddt"):.1f} | n4/100 med {med("n4_per100"):.2f} n5/100 med {med("n5_per100"):.2f} '
      f'| in-region>=4 med {med("inreg4")} | period-matched>=4 med {med("pm4")} proteins with >=1 {sum(int(r["pm4"]) > 0 for r in q)} '
      f'total {sum(int(r["pm4"]) for r in q)} | period-matched>=5 total {sum(int(r["pm5"]) for r in q)} '
      f'| units/centre(>=4, in-region) med of medians {med("units4_median")} max {max([int(r["units4_max"]) for r in q if r["units4_max"] not in ("", "None")] or [0])} '
      f'| pLDDT>=70 models {sum(float(r["mean_plddt"]) >= 70 for r in q)}')
top = sorted(M, key=lambda r: (-int(r['pm4']), -float(r['n4_per100']), r['record_id']))[:10]
p('top10 by period-matched >=4 centres:')
for r in top:
    p('\t'.join(r[k] for k in ('group', 'family', 'record_id', 'protein_id', 'fold_len', 'seq_length', 'mean_plddt', 'n4', 'n4_per100',
                               'inreg4', 'pm4', 'units4_max', 'n5', 'pm5', 'period', 'consensus')) + '\tUNVERIFIED')
open(os.path.join(RES, 'stats.log'), 'w').write('\n'.join(L) + '\n')
