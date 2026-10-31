"""Retrieve only the pinned data, two required best-model ensembles and metadata."""
import json, hashlib, urllib.request, concurrent.futures
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
W=ROOT/'workspaces/ln_structure_increment_20261030/upstream'
a=json.loads((HERE/'assets.json').read_text())
def fetch(f):
 p=W/f['path']; p.parent.mkdir(parents=True,exist_ok=True)
 url='https://raw.githubusercontent.com/Critical-Materials-Institute/stabilityconstant-ml-models/'+a['sha']+'/'+urllib.parse.quote(f['path'])
 if not p.exists():
  with urllib.request.urlopen(url,timeout=180) as r: b=r.read()
  assert len(b)==f['size'], (p,len(b),f['size'])
  p.write_bytes(b)
 b=p.read_bytes()
 assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==f['sha']
 return {'path':str(p.relative_to(W)), 'url':url,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: out=list(ex.map(fetch,a['files']))
(W/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
print('verified assets',len(out),'bytes',sum(x['bytes'] for x in out))
