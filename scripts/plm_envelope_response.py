"""Read-only geometry and response decomposition of the exact six PLM sources."""
from __future__ import annotations
import argparse,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from coordination_policy import COORDINATION_CUTOFF_A,donor_type
from accommodation_reference_geometry import CUTOFF_A as LA_WARNING_A
from mace_site_kinematics import Kinematics

EXPECTED_SHA='06e2c0b7727ff0aa53e6085b144a307c2eb6ca71cd0cc8f06112e3b8ec9e4be7'
ROLES=('anchor_glutamate','anchor_asparagine','catalytic_aspartate','extra_acidic_ligand_homolog')


def load(pin):return read_json(verify(pin))


def geometry(kin,q,selectors):
    physical,context,*_=kin.evaluate(q);meta=kin.data['source_atom_metadata'];roles={};moving=defaultdict(list)
    for i,a in enumerate(meta):
        if not a or a['element'] in ('H','D'):continue
        extent=float(np.linalg.norm(physical[i]-kin.positions[i]))
        if extent>1e-8:moving[(a['chain'],a['resnum'],a['resname'])].append({'atom':a['atom'],'displacement_A':extent,'source_id':kin.data['physical_ids'][i]})
    for role in ROLES:
        sel=selectors[role];points=[]
        for i,a in enumerate(meta):
            if not a or (a['chain'],a['resnum'],a.get('insertion_code',''))!=(sel['chain'],sel['resnum'],sel.get('icode','')):continue
            if a['atom'] not in ('OD1','OD2','OE1','OE2') or donor_type(a['resname'],a['atom'],a['element'])!='sidechain_O':continue
            points.append({'atom':a['atom'],'physical_index':i,'source_id':kin.data['physical_ids'][i],
                           'distance_A':float(np.linalg.norm(physical[i]-physical[kin.metal])),
                           'displacement_A':float(np.linalg.norm(physical[i]-kin.positions[i]))})
        if len(points)!=(1 if role=='anchor_asparagine' else 2):raise InvalidArtifact('complete typed role donor atoms unavailable')
        contacts=sum(p['distance_A']<=COORDINATION_CUTOFF_A for p in points)
        roles[role]={'source_selector':sel,'atoms':points,'typed_contacts_within_3p1_A':contacts,
                     'nearest_A':min(p['distance_A'] for p in points),'farthest_A':max(p['distance_A'] for p in points)}
    pqq=[i for i,a in enumerate(meta) if a and a['resname']=='PQQ']
    if len(pqq)!=24 or not np.array_equal(physical[pqq],kin.positions[pqq]):raise InvalidArtifact('PQQ heavy identity/fixed coordinates differ')
    if not np.array_equal(physical[kin.metal],kin.positions[kin.metal]):raise InvalidArtifact('metal moved')
    return {'roles':roles,'moving_residues':[{'chain':c,'resnum':n,'resname':r,'heavy_atoms':v} for (c,n,r),v in sorted(moving.items())],
            'maximum_heavy_displacement_A':float(np.linalg.norm(physical[kin.heavy]-kin.positions[kin.heavy],axis=1).max()),
            'PQQ_heavy_atoms_fixed':24,'metal_fixed':True},context


def analyze(result,output):
    if record(result)['sha256']!=EXPECTED_SHA:raise InvalidArtifact('authoritative six-source result differs')
    r=read_json(result);pc=load(r['source_collections']['pool']);pm=load(pc['manifest']);mp=pm['proposal_manifest'];m=load(mp)
    if r['source_denominator']!=6 or r['source_available']!=6 or len(r['rows'])!=6:raise InvalidArtifact('all six required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);rows=[]
    for row in r['rows']:
        case=next(c for c in pc['cases'] if c['case_id']==row['case_id']);source=case['source'];parent=load(source['original_core']['parent']);selectors=parent['fixed_core']['requested_roles']
        fragments={f['role']:f for f in parent['qm_fragments']}
        for role in ('anchor_glutamate','catalytic_aspartate','extra_acidic_ligand_homolog'):
            if fragments[role]['formal_charge']!=-1 or fragments[role].get('carboxyl_hydrogens'):raise InvalidArtifact('unexpected acidic state; do not assume deprotonation')
        ts={z:next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(row['case_id'],z)) for z in ('Ca','La')};maps={z:load(t['mapping'])['context'] for z,t in ts.items()}
        if maps['Ca']!=maps['La'] or ts['La']['charge']-ts['Ca']['charge']!=1:raise InvalidArtifact('paired physical/state mapping differs')
        kin=Kinematics(maps['Ca']);points={};receipts={}
        for name in ('origin','adaptive_Ca','adaptive_La'):
            q=np.zeros(len(kin.modes))
            if name!='origin':
                candidate=next(c for c in row['candidates'] if c['id']==name);receipt=load(candidate['proposal_receipt'])
                if receipt['manifest']!=mp or receipt['status']!='proposal_available':raise InvalidArtifact('actual proposal receipt mismatch')
                q=np.asarray(receipt['proposal']['full_q']);receipts[name]=candidate
            geom,coords=geometry(kin,q,selectors)
            for z in ('Ca','La'):
                atoms=xyz(verify(row['matrix'][z][name]['xyz']))
                if not np.allclose(coords,[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('cross-metal score cell geometry differs')
            geom.update(full_q=q.tolist(),coordinate=row['matrix']['Ca'][name]['xyz'],boundary_flag=False if name=='origin' else receipts[name]['boundary_flag'])
            points[name]=geom
        selections={z:row['pool']['rows'][z]['operational_candidate'] for z in ('Ca','La')}
        works={z:row['pool']['rows'][z]['work_from_origin_kcal_mol'][selections[z]] for z in ('Ca','La')}
        differential={key:works['Ca'][key]-works['La'][key] for key in ('native_MACE_kcal_mol','solvent_transfer_kcal_mol','composite_kcal_mol')}
        delta=row['pool']['operational']['composite_R_model_kcal_mol']-row['union_origin_R_model_kcal_mol']
        if abs(delta-differential['composite_kcal_mol'])>1e-7:raise InvalidArtifact('component score algebra differs')
        rows.append({'case_id':row['case_id'],'protein_id':row['protein_id'],'source_evidence':row['source_evidence'],
                     'source':source['source'],'parent':source['original_core']['parent'],'context_preparation':source['representations']['context']['preparation'],
                     'mapping':ts['Ca']['mapping'],'state_signature':source['state_signature'],'atom_count':source['atom_count'],'cap_count':source['cap_count'],
                     'charges':{z:t['charge'] for z,t in ts.items()},'selected_modes':{z:t['active_mode_ids'] for z,t in ts.items()},
                     'points':points,'proposal_receipts':receipts,'selected_candidates':selections,'selected_work':works,'differential_work':differential,
                     'origin_R':row['union_origin_R_model_kcal_mol'],'pool_R':row['pool']['operational']['composite_R_model_kcal_mol']})
    result={'authoritative_result':record(result),'proposal_manifest':mp,'rows':rows,'group_spreads':[{k:g[k] for k in ('protein_id','members','union_origin','variants')} for g in r['groups']],
            'policy':{'typed_contact_cutoff_A':COORDINATION_CUTOFF_A,'previous_La_short_contact_warning_A':LA_WARNING_A,
                      'contact_policy':record(Path(__file__).with_name('coordination_policy.py')),'warning_policy':record(Path(__file__).with_name('accommodation_reference_geometry.py')),
                      'denticity':'count actual same-residue carboxylate oxygen donors within existing typed3.1A cutoff; not an experimental bond assertion'},
            'implementation':record(__file__),'new_molecular_calls':0,'new_proposals':0,'new_thresholds':0,'biological_accuracy_evaluated':False}
    write_new(out/'GEOMETRY.json',result);return {'geometry':record(out/'GEOMETRY.json'),'sources':len(rows)}


def plot(geometry,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.size':10})
    a=read_json(geometry);out=Path(output);out.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(2,2,figsize=(11,7),gridspec_kw={'width_ratios':[1.1,1]},layout='constrained')
    colors=['#0072B2','#D55E00','#009E73'];groups=sorted({r['protein_id'] for r in a['rows']})
    for index,protein in enumerate(groups):
        rows=[r for r in a['rows'] if r['protein_id']==protein];left,right=axes[index];short=protein.split('_')[1][:4]
        for i,r in enumerate(rows):
            for atom,style in (('OD2','-'),('OD1','--')):
                y=[next(v['distance_A'] for v in r['points'][name]['roles']['extra_acidic_ligand_homolog']['atoms'] if v['atom']==atom) for name in ('origin',r['selected_candidates']['Ca'],r['selected_candidates']['La'])]
                left.plot([0,1,2],y,style,marker='o',color=colors[i],lw=1.7,label=f'Sample {i}' if atom=='OD2' else None)
                if atom=='OD2':
                    for pos,z in ((1,'Ca'),(2,'La')):
                        if r['points'][r['selected_candidates'][z]]['boundary_flag']:left.scatter(pos,y[pos],s=95,marker='s',facecolors='none',edgecolors=colors[i],linewidths=1.4)
            d=r['differential_work'];right.bar(i-.16,d['native_MACE_kcal_mol'],.30,color='#4477AA',label='Native MACE' if i==0 else None)
            right.bar(i+.16,d['solvent_transfer_kcal_mol'],.30,color='#EE9944',label='ALPB − vacuum' if i==0 else None)
            right.scatter(i,d['composite_kcal_mol'],color='black',marker='D',s=26,label='Composite total' if i==0 else None,zorder=4)
        left.axhline(COORDINATION_CUTOFF_A,color='.45',lw=.8,ls=':',label='Typed contact 3.1 Å' if index==0 else None)
        left.axhline(LA_WARNING_A,color='.6',lw=.8,ls='-.',label='Prior La warning 2.2 Å' if index==0 else None)
        left.set(xticks=[0,1,2],xticklabels=['Source q0','Ca selected','La selected'],ylabel='Extra-Asp metal–O distance (Å)',ylim=(1.55,4.65),title=f'{short}: intact extra-Asp response')
        left.text(.02,.98,'Solid OD2; dashed OD1\nOpen square: 0.8 Å boundary',transform=left.transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
        right.axhline(0,color='.4',lw=.8);right.set(xticks=[0,1,2],xticklabels=['Sample 0','Sample 1','Sample 2'],ylabel='Ca work − La work (kcal/mol)',title=f'{short}: actual selected contrast change')
        group=next(g for g in a['group_spreads'] if g['protein_id']==protein)
        right.text(.02,.98,f"Source score range {group['union_origin']['range_kcal_mol']:.2f} → {group['variants']['operational']['range_kcal_mol']:.2f}",transform=right.transAxes,va='top',fontsize=8)
        right.set_ylim(-9,100)
    axes[0,0].legend(loc='center left',bbox_to_anchor=(0,.40),fontsize=8,ncol=2);axes[0,1].legend(loc='upper right',fontsize=8,bbox_to_anchor=(1,.88))
    fig.suptitle('Six consumed PLM structures: compression relief and model response\nExperimental metal-use labels unknown',fontsize=13)
    for ext in ('svg','pdf'):fig.savefig(out/f'plm_geometry_response.{ext}')
    plt.close(fig);write_new(out/'PINS.json',{'geometry':record(geometry),'implementation':record(__file__),'figures':{ext:record(out/f'plm_geometry_response.{ext}') for ext in ('svg','pdf')}})
    return read_json(out/'PINS.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'analyze':('result','output'),'plot':('geometry','output')}.items():
        q=sub.add_parser(op)
        for field in fields:q.add_argument('--'+field,required=True,type=Path)
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
