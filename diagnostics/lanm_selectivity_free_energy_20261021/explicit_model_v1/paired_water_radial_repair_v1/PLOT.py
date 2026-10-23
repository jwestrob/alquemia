#!/usr/bin/env python3
"""Plot saved water-distance responses; no new molecular evaluations."""
import argparse, hashlib, json, os
from pathlib import Path
os.environ['MPLBACKEND']='Agg'
import matplotlib.pyplot as plt
import numpy as np

def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    a=argparse.ArgumentParser(); a.add_argument('--training',required=True); a.add_argument('--result',required=True); a.add_argument('--out',required=True); x=a.parse_args()
    out=Path(x.out); out.mkdir(parents=True,exist_ok=False)
    data=json.loads(Path(x.training).read_text()); result=json.loads(Path(x.result).read_text())
    assert data['complete']==8 and result['status']=='bounded_radial_refinement_complete'
    models={'PBE0-D4':{k:r['reference'] for k,r in data['rows'].items()},'Original parameters':{k:r['original_baseline'] for k,r in data['rows'].items()},'Center fit':{k:r['baseline'] for k,r in data['rows'].items()},'Distance-inclusive fit':result['selected']['model']}
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex='col',layout='constrained'); saved=[]
    styles=[('black','o'),('0.55','x'),('#e69f00','s'),('#0072b2','D')]
    for col,metal in enumerate(('La','Dy')):
        ids=['water_'+s+'_'+metal for s in ('minus','central','plus')]
        rs=[]; directions=[]
        for gid in ids:
            g=data['rows'][gid]['geometry']; pos=np.array(g['positions_A']); d=np.array(g['radial_direction']); assert np.isclose(np.linalg.norm(d),1)
            assert np.dot(pos[0]-pos[1],d)>0
            rs.append(float(np.linalg.norm(pos[0]-pos[1]))); directions.append(d)
        for (name,model),(color,marker) in zip(models.items(),styles):
            e0=model[ids[1]]['energy']; energies=[model[k]['energy']-e0 for k in ids]
            deriv=[float(-np.dot(model[k]['forces'][0],d)) for k,d in zip(ids,directions)]
            axes[0,col].plot(rs,energies,color=color,marker=marker,label=name)
            axes[1,col].plot(rs,deriv,color=color,marker=marker,label=name)
            saved.append({'metal':metal,'model':name,'metal_O_distance_A':rs,'relative_energy_kcal_mol':energies,'outward_dE_dR_kcal_mol_A':deriv})
        axes[0,col].set_title(metal); axes[1,col].set_xlabel('Metal–oxygen distance (Å)')
    for ax in axes.flat: ax.axhline(0,color='0.8',lw=.7); ax.grid(alpha=.15)
    axes[0,0].set_ylabel('Energy relative to center (kcal/mol)'); axes[1,0].set_ylabel('Outward energy derivative (kcal/mol/Å)')
    axes[0,0].legend(fontsize=8); fig.suptitle('Consumed water configurations: energy and mechanical response')
    fig.savefig(out/'water_response.png',dpi=180); fig.savefig(out/'water_response.pdf'); plt.close(fig)
    (out/'DATA.json').write_text(json.dumps({'inputs':[pin(x.training),pin(x.result)],'implementation':pin(__file__),'series':saved,'interpretation':'Three fixed configurations joined for readability; not continuous curves or independent validation. Derivative displaces only the metal outward; water and internal geometry fixed.'},indent=2)+'\n')

if __name__=='__main__':main()
