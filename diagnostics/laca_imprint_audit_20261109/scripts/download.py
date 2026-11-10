"""Step 1b: download the mmCIF files listed in WS/out/candidates.json. Usage: python -I download.py WS"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rcsb as r  # noqa: E402
WS = sys.argv[1]
c = json.load(open(os.path.join(WS, 'out', 'candidates.json')))
ids = sorted({e for v in c['clusters'].values() for e in v['ln'] + v['ca']})
for e in ids:
    r.download_cif(e, os.path.join(WS, 'cif'))
print('downloaded', len(ids))
