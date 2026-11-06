import json,hashlib,os
from pathlib import Path
D=Path(__file__).resolve().parent;ROOT=D.parents[1];W=ROOT/'workspaces/so3lr_tauri_screen_20261106'
DATA=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/external_review_20261023/data'
RESULTS_DIR=Path(os.environ.get('TAURI_RESULTS_DIR',str(D)));RESULTS_DIR.mkdir(parents=True,exist_ok=True)
EV_KCAL=23.060547830619;EH_KCAL=627.5094740631;BOHR_A=.529177210903
SYMBOLS={1:'H',6:'C',7:'N',8:'O',57:'La',66:'Dy'}
def read(p):return json.loads(Path(p).read_text())
def save(p,j):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
def pin(p):return {'path':str(Path(p).resolve()),'sha256':hashlib.sha256(Path(p).read_bytes()).hexdigest()}
def xyzread(p):
 l=Path(p).read_text().splitlines();n=int(l[0]);a=[s.split() for s in l[2:2+n]]
 return [next(z for z,s in SYMBOLS.items() if s==r[0]) for r in a],[[float(v) for v in r[1:4]] for r in a]
def xyzwrite(p,z,x,comment=''):
 Path(p).write_text(str(len(z))+'\n'+comment+'\n'+'\n'.join(SYMBOLS[a]+' '+' '.join(format(v,'.17g') for v in r) for a,r in zip(z,x))+'\n')
