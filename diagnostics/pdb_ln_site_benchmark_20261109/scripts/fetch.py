"""Download mmCIFs listed one-per-line in LIST into WS/cif unless present in WS/cif or AUDIT_CIF.
Usage: python -I fetch.py WS AUDIT_SCRIPTS_DIR AUDIT_CIF LIST"""
import os
import sys

WS, AUD, ACIF, LIST = sys.argv[1:5]
sys.path.insert(0, AUD)
import rcsb as r  # noqa: E402

ids = [x.strip().upper() for x in open(LIST) if x.strip()]
n_new = n_aud = n_have = n_fail = 0
for e in ids:
    if os.path.exists(os.path.join(ACIF, e + '.cif')):
        n_aud += 1; continue
    if os.path.exists(os.path.join(WS, 'cif', e + '.cif')):
        n_have += 1; continue
    try:
        r.download_cif(e, os.path.join(WS, 'cif')); n_new += 1
    except Exception as ex:  # noqa: BLE001
        n_fail += 1; print('FAIL', e, ex, flush=True)
print('audit', n_aud, 'already', n_have, 'new', n_new, 'failed', n_fail)
