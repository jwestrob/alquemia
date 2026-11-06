import json,hashlib,urllib.request,urllib.parse,concurrent.futures
from pathlib import Path
D=Path(__file__).resolve().parent;W=D.parents[1]/'workspaces/so3lr_tauri_screen_20261106'
m=json.load(open(D/'UPSTREAM.json'));root=W/'software/so3lr'
def fetch(f):
 p=root/f['path'];p.parent.mkdir(parents=True,exist_ok=True)
 url='https://raw.githubusercontent.com/'+m['repository']+'/'+m['commit']+'/'+urllib.parse.quote(f['path'])
 if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=90).read())
 b=p.read_bytes();assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==f['sha'],p
 return dict(path=f['path'],sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),url=url)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool: out=list(pool.map(fetch,m['files']))
(W/'software/RECEIPTS.json').write_text(json.dumps(out,indent=2)+'\n');print('verified',len(out),'source/assets',sum(r['bytes'] for r in out),'bytes')
