"""Smoke test on non-benchmark Ca structures (code check only; no tuning).
Usage: python -I smoke.py WS SCRIPTS AUDIT_SCRIPTS SCANNER_ROOT"""
import os
import sys
import time

WS, SCR, AUD, SROOT = sys.argv[1:5]
sys.path.insert(0, AUD)
import rcsb as r  # noqa: E402
sm = os.path.join(WS, 'smoke_cif')
os.makedirs(sm, exist_ok=True)
for e in ('4ICB', '1EXR'):
    r.download_cif(e, sm)
sys.argv = ['score.py', WS, SCR, sm, SROOT, '1']
sys.path.insert(0, SCR)
import score as S  # noqa: E402
S.C.CIF_DIRS[:] = [sm]
for e, ch in (('4ICB', 'A'), ('1EXR', 'A')):
    st = S.C.read(e)
    cas = [a.pos.tolist() for c in st[0] for res in c for a in res if a.element.name == 'Ca']
    u = {'uid': f'T|{e}|{ch}', 'tier': 'smoke', 'entry': e, 'chain': ch,
         'sites': [{'sid': str(i), 'xyz': x, 'cls': 'single-chain'} for i, x in enumerate(cas)]}
    t = time.time()
    o = S.work(u)
    if 'error' in o:
        print(o['error']); continue
    print(e, 'res', o['n_res'], 'pts', o['n_pts'], 'cent', o['n_cent'], 'b1err', o['b1_err'], 'sec', round(time.time() - t, 1))
    for s in o['sites']:
        print('  ', s['rank'], 'S_true', s['S_true'], 'an', s['an_true'], 'clash', s['clash_true'],
              'near', s['n_near_res'], 'reach', s['n_near_reach'], s['near_res'], 'Stop1', s['S_top1'])
