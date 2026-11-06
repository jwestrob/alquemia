from common import *
rows=[]
def add(cid,z,x,source,kind,old=None,mode=None):
 metal=next(a for a in z if a in (57,66));q=3 if kind=='aquo' else 2;mult=6 if metal==66 else 1
 r=dict(id=cid,Z=z,positions_A=x,charge=q,multiplicity=mult,kind=kind,source=pin(source),historical_reference=old,coordinate=mode)
 r['xyz']=str(W/'inputs'/f'{cid}.xyz');xyzwrite(r['xyz'],z,x,f'{cid} Q={q} multiplicity={mult}');rows.append(r)
(W/'inputs').mkdir(exist_ok=True)
forces=read(DATA/'water/FORCES.json')
for n in (8,9):
 for met in ('La','Dy'):
  cid=f'aquo{n}_central_{met}';p=DATA/'water/fixtures'/f'{cid}.json';j=read(p);g=j['geometry']
  add(cid,g['atomic_numbers'],g['xyz_A'],p,'aquo',forces[cid]['reference'])
for met,fn in [('La','donor_center.json'),('Dy','donor_Dy_report_only.json')]:
 p=DATA/'input_pins'/fn;j=read(p);g=j['geometry'];add(f'dry_acetate_{met}',g['atomic_numbers'],g['xyz_A'],p,'dry_acetate',j['reference_targets'])
p=DATA/'input_pins/donor_displacements.json'
for cid,j in read(p)['rows'].items():
 g=j['geometry'];add(cid,g['atomic_numbers'],g['xyz_A'],p,'dry_acetate',j['reference'],j['coordinate'])
for met in ('La','Dy'):
 p=ROOT/'diagnostics/lady_electronic_decomposition_20261023/inputs/xyz'/f'hydrated_acetate_mono_central_{met}.xyz';z,x=xyzread(p)
 add(f'hydrated_acetate_{met}',z,x,p,'hydrated_acetate')
assert len(rows)==24
save(W/'PANEL.json',dict(rows=rows,source_plan=pin(D/'PLAN.md'),all_coordinates_archived=True))
save(D/'INPUT_PINS.json',dict(plan=pin(D/'PLAN.md'),panel=pin(W/'PANEL.json'),sources=list({r['source']['path']:r['source'] for r in rows}.values()),software=pin(W/'software/RECEIPTS.json'),environment=pin(D/'environment-lock.txt')))
print('prepared exact panel',len(rows))
