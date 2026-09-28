"""Read-only, predeclared metal-load projections from actual compact gradients."""
import argparse,json
from pathlib import Path
import numpy as np
from affordable_common import read_json,record,verify,xyz,write_new,InvalidArtifact,HA_TO_KCAL
from affordable_response import read_engrad
from affordable_common import BOHR_TO_A

SOURCES=('Hans8DQ2','Hans8FNR','Mex8FNS')

def geometry(atoms):
    metals=[i for i,a in enumerate(atoms) if a[0] in ('La','Dy')]
    if len(metals)!=1:raise InvalidArtifact('Expected exactly one target metal')
    mi=metals[0];q=np.array([a[1:] for a in atoms],float)
    oxy=sorted(((float(np.linalg.norm(q[i]-q[mi])),i) for i,a in enumerate(atoms) if a[0]=='O'))
    if not oxy or oxy[0][0]<=0:raise InvalidArtifact('Real oxygen direction unavailable')
    distance,oi=oxy[0]
    return dict(metal_index=mi,oxygen_index=oi,distance_A=distance,unit_vector=((q[oi]-q[mi])/distance).tolist(),oxygen_distances_A=[dict(index=i,distance_A=d) for d,i in oxy])

def freeze(manifest,agreement):
    m=read_json(manifest);tasks={t['task_id']:t for t in m['tasks']}
    om=read_json(verify(m['reused_origin']['manifest']))
    configs={}
    for s in SOURCES:
        refs={metal:(om['tasks'][0]['xyz'] if s=='Hans8DQ2' and metal=='Dy' else tasks[s+'_'+metal]['xyz']) for metal in ('La','Dy')}
        ats={k:xyz(verify(v)) for k,v in refs.items()}
        if not np.array_equal([a[1:] for a in ats['La']],[a[1:] for a in ats['Dy']]):raise InvalidArtifact('Paired coordinates differ')
        if [a[0].replace('La','Ln') for a in ats['La']] != [a[0].replace('Dy','Ln') for a in ats['Dy']]:raise InvalidArtifact('Paired atoms differ')
        configs[s]=dict(xyz=refs,**geometry(ats['La']))
    return dict(protocol='compact_metal_force_diagnostic_v1',manifest=record(manifest),agreement=record(agreement),implementation=record(__file__),sources=configs)

def endpoint(row,coordinates,direction):
    if row.get('status')!='complete':return dict(status='unavailable',reason=row.get('reason',row.get('status','missing')))
    try:
        atoms=xyz(verify(coordinates));raw=read_engrad(verify(row['engrad']))
        nums={'H':1,'C':6,'N':7,'O':8,'La':57,'Dy':66}
        if raw['atomic_numbers'].tolist()!=[nums[a[0]] for a in atoms]:raise InvalidArtifact('Gradient atom order differs')
        if not np.allclose(raw['coordinates_bohr']*BOHR_TO_A,[a[1:] for a in atoms],rtol=0,atol=1e-6):raise InvalidArtifact('Gradient coordinates differ')
        if abs(raw['energy_Ha']-row['energy_hartree'])>1e-8:raise InvalidArtifact('Gradient energy differs')
        g=raw['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A
        if not np.isfinite(g).all() or not np.allclose(g,row['gradient_kcal_mol_per_A'],atol=1e-10,rtol=1e-12):raise InvalidArtifact('Collected gradient differs from artifact')
        gm=g[direction['metal_index']];load=float(gm@np.array(direction['unit_vector']));total=g.sum(axis=0)
        return dict(status='available',gradient_load_kcal_mol_A=load,force_load_kcal_mol_A=-load,metal_gradient_kcal_mol_A=gm.tolist(),metal_gradient_norm_kcal_mol_A=float(np.linalg.norm(gm)),translation_gradient_sum_kcal_mol_A=total.tolist(),translation_gradient_sum_norm_kcal_mol_A=float(np.linalg.norm(total)),engrad=row['engrad'],electronic_state=row.get('electronic_state'),force_consistency_qualified=row.get('force_consistency_qualified',False),physical_state_qualified=row.get('physical_state_qualified',False))
    except (ValueError,OSError,KeyError) as e:return dict(status='invalid',reason=str(e))

def analyze(design,collection):
    d=read_json(design);verify(d['manifest']);verify(d['agreement']);verify(d['implementation']);c=read_json(collection);rows={}
    for s,conf in d['sources'].items():
        ep={k:endpoint(c.get('rows',{}).get(s+'_'+k,{}),conf['xyz'][k],conf) for k in ('La','Dy')}
        delta=ep['Dy']['gradient_load_kcal_mol_A']-ep['La']['gradient_load_kcal_mol_A'] if all(v['status']=='available' for v in ep.values()) else None
        rows[s]=dict(direction=conf,endpoints=ep,Dy_minus_La_gradient_load_kcal_mol_A=delta,Dy_minus_La_force_load_kcal_mol_A=-delta if delta is not None else None)
    return dict(design=record(design),collection=record(collection),rows=rows,interpretation='Conditional fixed-core projected gradients, not affinity or relaxation work',protein_coordinate_mapping_qualified=False)

def report(r):
    lines=['# Compact metal-load diagnostic','',r['interpretation'],'','Positive gradient load means moving the metal toward its preselected nearest oxygen raises energy. Force has the opposite sign. Units: kcal/mol/angstrom.','', '| Source | nearest O index / distance Å | La gradient | Dy gradient | Dy−La |','|---|---|---:|---:|---:|']
    def fmt(v):return 'unavailable' if v is None else f'{v:.6g}'
    for s,row in r['rows'].items():
        d=row['direction']; ep=row['endpoints'];lines.append(f"| {s} | {d['oxygen_index']} / {d['distance_A']:.6f} | {fmt(ep['La'].get('gradient_load_kcal_mol_A'))} | {fmt(ep['Dy'].get('gradient_load_kcal_mol_A'))} | {fmt(row['Dy_minus_La_gradient_load_kcal_mol_A'])} |")
    lines+=['','Full endpoint status, gradient norms, oxygen distances and capped-system translation residuals are retained in JSON. No force is subtracted to impose translational invariance. Directions are source-specific, so their load differences are not a common-coordinate finite difference between different sources. No physical protein/link-force projection or state/affinity qualification is claimed.']
    return '\n'.join(lines)+'\n'

def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='op',required=True)
    a=sub.add_parser('freeze');a.add_argument('--manifest',required=True);a.add_argument('--agreement',required=True);a.add_argument('--output',required=True)
    a=sub.add_parser('analyze');a.add_argument('--design',required=True);a.add_argument('--collection',required=True);a.add_argument('--output',required=True)
    a=p.parse_args();r=freeze(a.manifest,a.agreement) if a.op=='freeze' else analyze(a.design,a.collection);write_new(a.output,r)
    if a.op=='analyze':
        md=Path(a.output).with_suffix('.md')
        with md.open('x') as f:f.write(report(r))
if __name__=='__main__':main()
