"""Score crystal chains and ESMFold2 models with the frozen v0 / B2 (PLAN.md sections 2-5).
Usage: python -I score_pred.py WS DIAG BENCH_WS BENCH_RESULTS AUDIT_CIF
WS/frozen holds the unchanged benchmark modules; WS/models holds <name>.pdb.
Writes DIAG/results/{bench_sites.tsv, bench_fp.tsv, capacity.tsv, ref_sites.tsv, cand_models.tsv,
cand_centres.tsv} and WS/out/*.json."""
import csv
import json
import os
import sys
import time
from collections import defaultdict

import gemmi
import numpy as np

WS, DIAG, BWS, BRES, ACIF = sys.argv[1:6]
sys.path.insert(0, os.path.join(WS, 'frozen'))
import common as C  # noqa: E402
import detect as Dt  # noqa: E402

C.CIF_DIRS[:] = [os.path.join(BWS, 'cif'), ACIF, os.path.join(WS, 'cif')]
RES = os.path.join(DIAG, 'results')
MOD = os.path.join(WS, 'models')
HIT = 4.0
LNSET = set(C.LN_EL)


# ---------- verbatim from benchmark score.py
def first_hit(centres, x):
    if len(centres) == 0:
        return None
    d = np.linalg.norm(np.asarray(centres) - np.asarray(x), axis=1)
    w = np.where(d <= HIT)[0]
    return int(w[0]) if len(w) else None


def residue_reach(P, x):
    """Per-residue v0 weight at point x (no clash filter)."""
    S, an, _ = Dt.v0_scores(P, x[None, :])
    out = {}
    for i, rn in enumerate(P['names']):
        w = 0.0
        if rn in Dt.REACH and not P['is_gly'][i]:
            cb, ca = P['anchor'][i], P['ca'][i]
            v = x - cb
            D = np.linalg.norm(v)
            ax = cb - ca
            ok_ang = np.dot(ax, v) / (np.linalg.norm(ax) * max(D, 1e-9)) >= np.cos(np.radians(Dt.MAX_ANGLE))
            for (rmin, rmax, wt) in Dt.REACH[rn]:
                if ok_ang and max(0, rmin - Dt.BAND[1], Dt.BAND[0] - rmax) <= D <= rmax + Dt.BAND[1]:
                    w = max(w, wt)
        O, Cc = P['bbO'][i], P['bbC'][i]
        if np.isfinite(O).all() and np.isfinite(Cc).all():
            d = np.linalg.norm(x - O)
            u1, u2 = Cc - O, x - O
            if Dt.BAND[0] <= d <= Dt.BAND[1] and np.dot(u1, u2) / (np.linalg.norm(u1) * np.linalg.norm(u2)) <= np.cos(np.radians(120)):
                w = max(w, Dt.BB_W)
        out[i] = w
    return float(S[0]), float(an[0]), out


def clash(P, x):
    okO = np.isfinite(P['bbO'][:, 0])
    return bool((np.linalg.norm(P['bbNCAC'] - x, axis=1) < Dt.CLASH['bbNCAC']).any()
                or (np.linalg.norm(P['bbO'][okO] - x, axis=1) < Dt.CLASH['bbO']).any()
                or (np.linalg.norm(P['anchor'][~P['is_gly']] - x, axis=1) < Dt.CLASH['CB']).any())
# ---------- end verbatim


def load_model(name):
    st = gemmi.read_structure(os.path.join(MOD, name + '.pdb'))
    st.setup_entities()
    st.remove_hydrogens()
    ch = st[0][0].name
    res = C.chain_residues(st[0], ch)
    pl = []
    for r in st[0][ch]:
        if r.name in C.AA:
            a = r.find_atom('CA', '*')
            pl.append(a.b_iso if a else np.nan)
    assert len(pl) == len(res)
    return res, np.array(pl)


def score_chain(res):
    P = Dt.prepare(res)
    pts = Dt.candidate_points(P)
    S, an, tb = Dt.v0_scores(P, pts)
    k0 = Dt.rank(S, tb, pts)
    n2, tb2 = Dt.b2_scores(P, pts)
    k2 = Dt.rank(n2, tb2, pts)
    return {'P': P, 'pts': pts, 'S': S, 'an': an, 'tb': tb, 'k0': k0, 'k2': k2}


def site_eval(sc, x):
    Sx, anx, contrib = residue_reach(sc['P'], x)
    return {'r_v0': first_hit(sc['pts'][sc['k0']], x), 'r_b2': first_hit(sc['pts'][sc['k2']], x),
            'S_true': Sx, 'an_true': anx, 'clash': clash(sc['P'], x)}


def merged(sc, t):
    m = np.where(sc['an'] >= t)[0]
    if not len(m):
        return np.zeros(0, int)
    order = m[np.lexsort((sc['tb'][m], -sc['S'][m], -sc['an'][m]))]
    return Dt.nms(sc['pts'], order)


def transfer(x, src_ca, dst_ca, mp):
    """Local superposition src->dst (Calpha within 12 A of x in src frame), as units.py."""
    near = np.where(np.linalg.norm(src_ca - x, axis=1) <= 12.0)[0]
    pairs = [(i, mp[i]) for i in near if i in mp]
    if len(pairs) < 8:
        return None, None
    R, t, rmsd = C.kabsch(dst_ca[[j for _, j in pairs]], src_ca[[i for i, _ in pairs]])
    return R @ x + t, rmsd


class Crystal:
    def __init__(self, entry, chain):
        self.st = C.read(entry)
        self.res = C.chain_residues(self.st[0], chain)
        self.ca = np.array([r[3]['CA'] for r in self.res])
        self.seq = C.seq_of(self.res)
        self.ns = gemmi.NeighborSearch(self.st[0], self.st.cell, 8).populate()

    def ln_within(self, x, r):
        for mk in self.ns.find_atoms(gemmi.Position(*x), '\0', radius=r):
            a = mk.to_cra(self.st[0]).atom
            if a.element.name in LNSET and C.alt_ok(a):
                return True
        return False

    def near_res(self, x):
        return [i for i, r in enumerate(self.res)
                if any(np.linalg.norm(p - x) <= 3.0 for n, p in r[3].items() if n[0] in 'ON')]


def ln_label_points(cr, pts_model, model_ca, inv_mp, r):
    """For model points: map to crystal frame; return array of 1 (Ln within r), 0 (none), -1 (unmappable)."""
    out = np.zeros(len(pts_model), int)
    for j, y in enumerate(pts_model):
        x, _ = transfer(y, model_ca, cr.ca, inv_mp)
        out[j] = -1 if x is None else int(cr.ln_within(x, r))
    return out


def fp_counts(sc, labeller, nres):
    o = {'nres': nres}
    for t in (4, 5):
        k = merged(sc, t)
        lab = labeller(sc['pts'][k]) if len(k) else np.zeros(0, int)
        o[f'n{t}'] = int(len(k))
        o[f'off{t}'] = int((lab != 1).sum())
        o[f'unmap{t}'] = int((lab == -1).sum())
        o[f'off{t}_per100'] = round(100 * o[f'off{t}'] / nres, 3)
    return o


def cap_rows(sc, labeller, tag, name):
    k = sc['k0']
    lab = labeller(sc['pts'][k])
    return [(tag, name, int(min(sc['an'][i], 8)), int(l)) for i, l in zip(k, lab)]


def wtsv(path, rows, cols):
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in cols) + '\n')


t0 = time.time()
units = {u['uid']: u for u in json.load(open(os.path.join(BWS, 'out', 'units.json')))}
bench_ranks = defaultdict(dict)
for r in csv.DictReader(open(os.path.join(BRES, 'site_ranks.tsv')), delimiter='\t'):
    bench_ranks[(r['tier'], r['sid'])][r['uid']] = r
summ = {r['name']: r for r in csv.DictReader(open(os.path.join(MOD, 'summary.tsv')), delimiter='\t')}

# ======================= Part 2 / 3: benchmark sample
site_rows, fp_rows, capr, fails = [], [], [], []
for s in csv.DictReader(open(os.path.join(RES, 'bench_sample.tsv')), delimiter='\t'):
    name = s['name']
    if not os.path.exists(os.path.join(MOD, name + '.pdb')):
        fails.append((name, 'no model')); continue
    u = units[s['uid']]
    cr = Crystal(s['entry'], s['chain'])
    mres, pl = load_model(name)
    ident, mp = C.align_map(cr.seq, C.seq_of(mres))
    inv = {v: k for k, v in mp.items()}
    mca = np.array([r[3]['CA'] for r in mres])
    trim_idx = sorted(mp.values())
    tres = [mres[j] for j in trim_idx]
    tpos = {j: n for n, j in enumerate(trim_idx)}
    sc_a = score_chain(cr.res)
    sc_c = score_chain(mres)
    sc_t = score_chain(tres)
    mp_t = {i: tpos[j] for i, j in mp.items()}
    for site in u['sites']:
        if not (site['cls'] == 'single-chain' and site['n_don'] >= 3):
            continue
        x = np.array(site['xyz'])
        y, rmsd = transfer(x, cr.ca, mca, mp)
        if y is None:
            fails.append((name, site['sid'] + ' transfer<8')); continue
        yt, _ = transfer(x, cr.ca, mca[trim_idx], mp_t)
        near = cr.near_res(x)
        sp = [pl[mp[i]] for i in near if i in mp]
        base = {'name': name, 'uid': s['uid'], 'entry': s['entry'], 'c95': s['c95'], 'c30': s['c30'], 'sid': site['sid'],
                'n_don': site['n_don'], 'rmsd': round(rmsd, 3), 'ident': round(ident, 3),
                'site_plddt': round(float(np.mean(sp)), 2) if sp else None, 'n_near': len(near),
                'model_plddt': round(float(np.nanmean(pl)), 2), 'model_len': len(mres), 'cryst_len': len(cr.res)}
        ea = site_eval(sc_a, x)
        br = bench_ranks[('stripped', site['sid'])].get(s['uid'])
        ok = br is not None and (br['r_v0'] or None) == (None if ea['r_v0'] is None else str(ea['r_v0'])) \
            and float(br['an_true']) == ea['an_true']
        site_rows.append({**base, 'tier': 'a_stripped', **ea, 'bench_match': ok})
        site_rows.append({**base, 'tier': 'c_pred', **site_eval(sc_c, y)})
        site_rows.append({**base, 'tier': 'c2_pred_trim', **site_eval(sc_t, yt)})
        for puid, pr in bench_ranks[('partner', site['sid'])].items():
            pu = units[puid]
            site_rows.append({**base, 'tier': 'b_partner', 'uid': puid, 'entry': pu['entry'], 'kind': pu['kind'],
                              'r_v0': None if pr['r_v0'] == '' else int(pr['r_v0']),
                              'r_b2': None if pr['r_b2'] == '' else int(pr['r_b2']),
                              'an_true': float(pr['an_true']), 'S_true': float(pr['S_true']), 'rmsd': pr['rmsd']})
    lab_a = lambda P, r=8.0: np.array([int(cr.ln_within(x, r)) for x in P], int)  # noqa: E731
    lab_c = lambda P, r=8.0: ln_label_points(cr, P, mca, inv, r)  # noqa: E731
    inv_t = {tpos[j]: i for i, j in mp.items()}
    lab_t = lambda P, r=8.0: ln_label_points(cr, P, mca[trim_idx], inv_t, r)  # noqa: E731
    for tag, sc, lab, n in (('a_stripped', sc_a, lab_a, len(cr.res)), ('c_pred', sc_c, lab_c, len(mres)),
                            ('c2_pred_trim', sc_t, lab_t, len(tres))):
        fp_rows.append({'name': name, 'c95': s['c95'], 'c30': s['c30'], 'tier': tag, **fp_counts(sc, lab, n)})
        capr += cap_rows(sc, lambda P, lab=lab: lab(P, HIT), tag, name)
    print('bench', name, len(cr.res), len(mres), round(time.time() - t0), flush=True)

cols = ['name', 'tier', 'uid', 'entry', 'kind', 'c95', 'c30', 'sid', 'n_don', 'r_v0', 'r_b2', 'S_true', 'an_true', 'clash',
        'rmsd', 'ident', 'site_plddt', 'n_near', 'model_plddt', 'model_len', 'cryst_len', 'bench_match']
wtsv(os.path.join(RES, 'bench_sites.tsv'), site_rows, cols)
wtsv(os.path.join(RES, 'bench_fp.tsv'), fp_rows, ['name', 'c95', 'c30', 'tier', 'nres', 'n4', 'off4', 'unmap4', 'off4_per100',
                                                   'n5', 'off5', 'unmap5', 'off5_per100'])
with open(os.path.join(RES, 'capacity.tsv'), 'w') as fh:
    fh.write('tier\tname\tn_anion_bin\tln_within4\n')
    for r in capr:
        fh.write('\t'.join(map(str, r)) + '\n')

# ======================= Part 4: reference proteins
lnsites = [r for r in csv.DictReader(open(os.path.join(BRES, 'ln_sites.tsv')), delimiter='\t')]
lnsites += [r for r in csv.DictReader(open(os.path.join(WS, 'ref_label', 'out', 'ln_sites.tsv')), delimiter='\t')]
ref_rows = []
for s in csv.DictReader(open(os.path.join(RES, 'ref_sample.tsv')), delimiter='\t'):
    name = s['name']
    if not os.path.exists(os.path.join(MOD, name + '.pdb')):
        fails.append((name, 'no model')); continue
    cr = Crystal(s['pdb'], 'A')
    mres, pl = load_model(name)
    ident, mp = C.align_map(cr.seq, C.seq_of(mres))
    mca = np.array([r[3]['CA'] for r in mres])
    sc_a, sc_c = score_chain(cr.res), score_chain(mres)
    for r in lnsites:
        if r['pdb'] != s['pdb'] or r['chain'] != 'A' or int(r['n_prot_don']) < 2:
            continue
        x = np.array([float(r['qx']), float(r['qy']), float(r['qz'])])
        y, rmsd = transfer(x, cr.ca, mca, mp)
        near = cr.near_res(x)
        sp = [pl[mp[i]] for i in near if i in mp]
        base = {'name': name, 'pdb': s['pdb'], 'site': r['site'], 'cls': r['cls'], 'n_don': r['n_prot_don'],
                'types': r['types'], 'ident': round(ident, 3), 'rmsd': None if rmsd is None else round(rmsd, 3),
                'site_plddt': round(float(np.mean(sp)), 2) if sp else None, 'model_plddt': round(float(np.nanmean(pl)), 2)}
        ea = site_eval(sc_a, x)
        ec = site_eval(sc_c, y) if y is not None else {}
        ref_rows.append({**base, **{f'x_{k}': v for k, v in ea.items()}, **{f'm_{k}': v for k, v in ec.items()}})
    print('ref', name, round(time.time() - t0), flush=True)
wtsv(os.path.join(RES, 'ref_sites.tsv'), ref_rows,
     ['name', 'pdb', 'site', 'cls', 'n_don', 'types', 'ident', 'rmsd', 'site_plddt', 'model_plddt',
      'x_an_true', 'm_an_true', 'x_S_true', 'm_S_true', 'x_r_v0', 'm_r_v0', 'x_r_b2', 'm_r_b2', 'x_clash', 'm_clash'])

# ======================= Part 5: candidates
cand = json.load(open(os.path.join(RES, 'cand_sample.json')))
mrows, crows = [], []
for c in cand:
    name = c['name']
    if not os.path.exists(os.path.join(MOD, name + '.pdb')):
        fails.append((name, 'no model')); continue
    mres, pl = load_model(name)
    assert len(mres) == c['fold_len'], name
    sc = score_chain(mres)
    per, rs, re_ = c['period'], c['reg_start'], c['reg_end']
    tol = max(2.0, 0.15 * per)
    row = {'name': name, 'family': c['family'], 'group': c['group'], 'record_id': c['record_id'],
           'protein_id': c['protein_id'], 'seq_length': c['seq_length'], 'fold_start': c['fold_start'],
           'fold_len': c['fold_len'], 'reg_start': rs, 'reg_end': re_, 'period': per,
           'mean_plddt': round(float(np.nanmean(pl)), 2), 'consensus': c['consensus']}
    for t in (4, 5):
        k = merged(sc, t)
        cents = []
        for i in k:
            x = sc['pts'][i]
            _, _, contrib = residue_reach(sc['P'], x)
            ac = [j for j, w in contrib.items() if w >= 1.0 and sc['P']['names'][j] in ('ASP', 'GLU')]
            pos = [c['fold_start'] + j for j in ac]  # model residue j (0-based) -> protein position (1-based)
            mpos = float(np.median(pos)) if pos else np.nan
            inreg = bool(pos) and rs <= mpos <= re_
            units_ = sorted({int((p - rs) // per) for p in pos if rs <= p <= re_})
            cents.append({'name': name, 't': t, 'an': int(sc['an'][i]), 'S': float(sc['S'][i]), 'mpos': mpos,
                          'inreg': inreg, 'n_units': len(units_), 'units': ','.join(map(str, units_)),
                          'pos': ','.join(map(str, pos)),
                          'plddt': round(float(np.mean([pl[j] for j in ac])), 2) if ac else None})
        ir = sorted([e for e in cents if e['inreg']], key=lambda e: e['mpos'])
        for a, b in zip(ir, ir[1:]):
            sp = b['mpos'] - a['mpos']
            if any(abs(sp - m * per) <= tol for m in (1, 2)):
                a['pm'] = b['pm'] = True
        for e in cents:
            e['pm'] = bool(e.get('pm', False))
        crows += cents
        row[f'n{t}'] = len(cents)
        row[f'n{t}_per100'] = round(100 * len(cents) / len(mres), 3)
        row[f'inreg{t}'] = len(ir)
        row[f'pm{t}'] = sum(e['pm'] for e in cents)
        row[f'units{t}_median'] = float(np.median([e['n_units'] for e in ir])) if ir else None
        row[f'units{t}_max'] = max([e['n_units'] for e in ir]) if ir else None
    mrows.append(row)
    print('cand', name, len(mres), round(time.time() - t0), flush=True)
wtsv(os.path.join(RES, 'cand_models.tsv'), mrows,
     ['name', 'family', 'group', 'record_id', 'protein_id', 'seq_length', 'fold_start', 'fold_len', 'reg_start', 'reg_end',
      'period', 'mean_plddt', 'n4', 'n4_per100', 'inreg4', 'pm4', 'units4_median', 'units4_max', 'n5', 'n5_per100',
      'inreg5', 'pm5', 'units5_median', 'units5_max', 'consensus'])
wtsv(os.path.join(RES, 'cand_centres.tsv'), crows, ['name', 't', 'an', 'S', 'mpos', 'inreg', 'pm', 'n_units', 'units', 'plddt', 'pos'])
json.dump(fails, open(os.path.join(RES, 'failures.json'), 'w'), indent=1)
print('fails', fails, 'sec', round(time.time() - t0))
