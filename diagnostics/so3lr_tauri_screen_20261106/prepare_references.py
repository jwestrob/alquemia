import ast,shutil
from common import *
assert read(W/'OLD_REFERENCE_AUDIT.json')['compatible_reference_available'] is False
src=ROOT/'diagnostics/lanm_coupled_mechanics_20261017/lbt_omol_target_reference_v1/assets';basis=src/'def2-tzvpd.bas'
assert pin(basis)['sha256']=='7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6'
code=ast.parse((src/'calc.py').read_text());vals={}
for node in code.body:
 if isinstance(node,ast.Assign):
  for t in node.targets:
   if isinstance(t,ast.Name) and t.id in ('ORCA_SIMPLE_INPUT','ORCA_BLOCKS'):
    vals[t.id]=ast.literal_eval(node.value)
rows=read(W/'PANEL.json')['rows'];wanted=['aquo8_central_La','aquo9_central_La','aquo8_central_Dy','aquo9_central_Dy','hydrated_acetate_La','hydrated_acetate_Dy'];out=[]
for cid in wanted:
 r=next(r for r in rows if r['id']==cid);d=W/'references'/cid;d.mkdir(parents=True,exist_ok=False)
 shutil.copyfile(r['xyz'],d/'start.xyz');shutil.copyfile(basis,d/'def2-tzvpd.bas')
 core=28 if 66 in r['Z'] else 46;ne=sum(r['Z'])-core-r['charge'];rot=''
 if r['multiplicity']==1:rot=f'%scf rotate {{{ne//2-1}, {ne//2}, 20, 1, 1}} end end\n'
 inp='! UKS wB97M-V def2-TZVPD '+' '.join(vals['ORCA_SIMPLE_INPUT'])+' NONBO NONPA\n'+'\n'.join(vals['ORCA_BLOCKS'])+'\n'
 inp+='%scf Guess PModel STABPerform true STABRestartUHFifUnstable false end\n'+rot+'%pal nprocs @NPROCS@ end\n%maxcore @MAXCORE@\n'+f'* xyzfile {r["charge"]} {r["multiplicity"]} start.xyz\n'
 assert not any(s in inp for s in (' D4',' SCNL','CPCM','Opt\n'))
 (d/'prepared.inp').write_text(inp)
 out.append(dict(id=cid,directory=str(d),charge=r['charge'],multiplicity=r['multiplicity'],expected_electrons=ne,expected_ecp_core=core,
                 input=pin(d/'prepared.inp'),xyz=pin(d/'start.xyz'),basis=pin(d/'def2-tzvpd.bas')))
save(W/'REFERENCE_MANIFEST.json',dict(rows=out,source=pin(src/'calc.py'),source_basis=pin(basis),plan=pin(D/'PLAN.md'),native_calls_cap=6,
 port='ORCA6.1.1; PModel initialization and diagnostic non-restarting orbital stability; no appended SCNL/D4'))
print('six inputs prepared; no molecular calls executed')
