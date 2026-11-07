import json,hashlib
from pathlib import Path
D=Path(__file__).resolve().parent;ROOT=D.parents[1];W=ROOT/'workspaces/so3lr_neutral_pair_20261107'
EV_KCAL=23.060547830619
def read(p):return json.loads(Path(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def pin(p):return dict(path=str(Path(p).resolve()),sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest())
def xyz(p,z,x,q=0,m=6):
 from ase.data import chemical_symbols
 Path(p).write_text(f'{len(z)}\nProperties=species:S:1:pos:R:3 charge={q} multiplicity={m} pbc="F F F"\n'+'\n'.join(chemical_symbols[int(a)]+' '+' '.join(format(v,'.17g') for v in row) for a,row in zip(z,x))+'\n')
