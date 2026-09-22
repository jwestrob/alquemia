"""Plot all six archived fixed-Asp differential responses; no model or fit."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CASES = [('4MAE', '4MAE', 'known reference'),
         ('PQQSEQ_83440678cbbd658047c9', 'PLM8344 · sample0', 'biological preference unknown'),
         ('PQQSEQ_07ab500e3df76b30d71c', 'PLM07ab · sample0', 'biological preference unknown')]
METHODS = [('DFT_delta_R', 'Native r2SCAN-3c/CPCM', '#263B53'),
           ('MACE_delta_R', 'Native MACE', '#D88428'),
           ('composite_delta_R', 'MACE + GFN2 solvent transfer', '#8273AA')]


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def render(source,report,output):
    source=Path(source).resolve();report=Path(report).resolve();data=json.loads(source.read_text())
    expected={(c,a) for c,_,_ in CASES for a in (-.2,.2)}
    actual={(r['case_id'],r['angle_radian']) for r in data['rows']}
    if actual!=expected or len(data['rows'])!=6:
        raise ValueError('Expected exactly all three contexts and both predeclared displacements')
    if data['metrics']['differential_denominator']!=6 or data['metrics']['endpoint_denominator']!=12:
        raise ValueError('Archived endpoint denominator differs')
    if pin(data['analysis']['path'])!=data['analysis']:
        raise ValueError('Archived final analysis changed')
    for r in data['rows']:
        for field,name,_ in METHODS:
            method={'DFT_delta_R':'DFT','MACE_delta_R':'MACE','composite_delta_R':'composite'}[field]
            direct=r['endpoints']['Ca'][method]-r['endpoints']['La'][method]
            if not np.isfinite(r[field]) or abs(direct-r[field])>1e-10:
                raise ValueError('Archived differential does not match endpoint work algebra')
        if abs(r['MACE_delta_R']+r['solvent_delta_R']-r['composite_delta_R'])>1e-10:
            raise ValueError('Archived component sum differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    with (out/'figure_data.json').open('x') as f:
        json.dump({'source':pin(source),'report':pin(report),'analysis':data['analysis'],
                   'coordinate':'intact terminal extra-Asp chi2 rotation, radian',
                   'quantity':'[E_Ca(q)-E_Ca(0)]-[E_La(q)-E_La(0)]',
                   'units':'kcal/mol','rows':data['rows'],'new_scientific_calls':0},f,indent=2);f.write('\n')
    with (out/'figure_data.csv').open('x',newline='') as f:
        writer=csv.writer(f);writer.writerow(['case_id','angle_radian','DFT_delta_R_kcal_mol',
                                             'MACE_delta_R_kcal_mol','composite_delta_R_kcal_mol'])
        for c,_,_ in CASES:
            for a in (-.2,.2):
                r=next(r for r in data['rows'] if (r['case_id'],r['angle_radian'])==(c,a))
                writer.writerow([c,a,*[r[k] for k,_,_ in METHODS]])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none',
                         'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.65,
                         'xtick.major.width':.65,'ytick.major.width':.65,
                         'savefig.facecolor':'white','svg.hashsalt':'fixed-donor-response-v1'})
    fig,axes=plt.subplots(1,3,figsize=(8.0,3.65))
    fig.subplots_adjust(left=.095,right=.985,bottom=.235,top=.69,wspace=.40)
    x=np.arange(2);width=.23
    for i,(ax,(case,title,subtitle)) in enumerate(zip(axes,CASES)):
        rows=[next(r for r in data['rows'] if (r['case_id'],r['angle_radian'])==(case,a)) for a in (-.2,.2)]
        ax.axhline(0,color='#6B7280',linewidth=.65,zorder=1)
        ax.set_axisbelow(True);ax.grid(axis='y',color='#E3E6E9',linewidth=.5)
        for j,(field,label,color) in enumerate(METHODS):
            ax.bar(x+(j-1)*width,[r[field] for r in rows],width=width*.88,
                   color=color,label=label,zorder=3)
        ax.set_xticks(x,['−0.2','+0.2']);ax.set_xlabel('Asp χ₂ displacement (rad)',labelpad=5)
        ax.set_xlim(-.52,1.52)
        if i==0:
            ax.set_ylim(-2.05,2.05);ax.set_yticks([-2,-1,0,1,2])
            ax.text(.5,1.005,'expanded y scale',transform=ax.transAxes,ha='center',va='bottom',
                    fontsize=7.5,color='#626B76')
            ax.set_ylabel('Change in Ca − La contrast\nΔR (kcal mol⁻¹)',labelpad=8)
        else:
            ax.set_ylim(-36,36);ax.set_yticks([-30,-15,0,15,30])
        ax.text(0,1.23,f'{chr(65+i)}  {title}',transform=ax.transAxes,fontweight='bold',fontsize=10)
        ax.text(0,1.115,subtitle,transform=ax.transAxes,fontsize=8,color='#465364')
        ax.spines[['top','right']].set_visible(False)
        ax.spines[['left','bottom']].set_color('#6B7280')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.53,.985),ncol=3,
               frameon=False,fontsize=8.5,handlelength=1.3,columnspacing=1.7)
    fig.text(.095,.052,'Positive ΔR: shift toward La  •  Fixed displacements only; no relaxed-candidate validation',
             fontsize=8.3,color='#394656')
    stem='fixed_displacement_response'
    fig.savefig(out/(stem+'.svg'),metadata={'Date':None,'Description':'All six archived fixed extra-Asp differential responses; no fitted curves.'})
    fig.savefig(out/(stem+'.pdf'),metadata={'Title':'Fixed extra-Asp displacement response',
                'Subject':'Native DFT, native MACE and composite; three consumed contexts, both directions',
                'Creator':'render_fixed_displacement_response.py','CreationDate':None,'ModDate':None})
    fig.savefig(out/(stem+'.png'),dpi=220)
    plt.close(fig)
    provenance={'renderer':pin(__file__),'source':pin(source),'report':pin(report),'analysis':data['analysis'],
                'matplotlib_version':matplotlib.__version__,'numpy_version':np.__version__,
                'source_rows':6,'plotted_model_values':18,'endpoint_work_denominator':12,
                'fits':0,'new_scientific_calls':0,'error_bars':'not estimated; single fixed-geometry evaluations',
                'axis_policy':'4MAE expanded y scale; both PLM panels share -36..36 kcal/mol',
                'editable_text':'SVG text retained; PDF TrueType fonts embedded',
                'artifacts':{p.name:pin(p) for p in sorted(out.iterdir()) if p.is_file()}}
    with (out/'PROVENANCE.json').open('x') as f:json.dump(provenance,f,indent=2);f.write('\n')
    return provenance


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('source','report','output'):p.add_argument('--'+k,required=True)
    a=vars(p.parse_args());r=render(**a)
    print(json.dumps({'source_rows':r['source_rows'],'plotted_model_values':r['plotted_model_values'],
                      'output':str(Path(a['output']).resolve()),'new_scientific_calls':0},indent=2))
