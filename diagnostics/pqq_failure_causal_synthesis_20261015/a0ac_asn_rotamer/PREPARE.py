#!/usr/bin/env python3
"""Freeze metadata or prepare the declared ten-source amide flips inside Slurm."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_rotamer'
PROTOCOL='A0AC_Asn290_discrete_terminal_amide_180_degree_preparation_v1'
INVENTORY_SHA='c91b152b1fb8b548395c41aac72fb2370b5ba73d015abbfce14767ac2787d4d8'
METALS=('Ca','La');MOVED=('OD1','ND2','HD21','HD22')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def require(v,reason):
    if not v:raise ValueError(reason)
def load(name,p):
    s=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def freeze(output):
    ip=pin(ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/a0ac_inventory_v1/EVIDENCE.json')
    require(ip['sha256']==INVENTORY_SHA,'Wrong complete ten-source inventory')
    inv=read(ip);require(len(inv['rows'])==10,'Wrong source denominator')
    old=read(ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/pose_preparation_v1/manifest.json')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir();(out/'inputs').mkdir()
    snapshots={};sources=[]
    def snapshot(p):
        source=check(p);key=p['sha256']
        if key not in snapshots:
            dest=out/'inputs'/(key+source.suffix);shutil.copyfile(source,dest)
            snapshots[key]={'original':p,'snapshot':pin(dest)}
        return snapshots[key]['snapshot']
    for r in inv['rows']:
        inputs={k:snapshot(r[k]) for k in ('source_preparation','state_signature','protonated_source')}
        inputs['maps']={z:snapshot(r['maps'][z]) for z in METALS}
        inputs['xyz']={z:snapshot(r['matrix'][z]['origin']['xyz']) for z in METALS}
        inputs['protonation']=snapshot(pin(Path(r['protonated_source']['path']).with_name('protonated_protonation_manifest.json')))
        sources.append({'source_case_id':r['source_case_id'],'record':r,'inputs':inputs})
    implementation={}
    for name in ('PREPARE.py','RUN.sbatch'):
        dest=out/'implementation'/name;shutil.copyfile(HERE/name,dest);implementation[name]=pin(dest)
    m={'protocol':PROTOCOL,'implementation':implementation,'inventory':ip,'sources':sources,'input_snapshots':snapshots,
       'python':old['python'],'source_snapshot':old['source_snapshot'],'source_modules':old['source_modules'],
       'source_parser_reference':pin(ROOT/'diagnostics/native_cavity_shielding_20261014/prepare/PREPARE.py'),
       'contact_policy_reference':pin(ROOT/'diagnostics/native_cavity_accommodation_20261014/prepare/PREPARE.py'),
       'plan':pin(HERE.parent/'ASN_ROTAMER_PLAN.md'),
       'settings':{'angle_degrees':180,'moved_atom_names':list(MOVED),'source_coordinate_tolerance_A':1e-12,
         'kinematic_replay_A':1e-8,'covalent_rigid_tolerance_A':1e-10,'new_overlap_heavy_A':1.0,'new_overlap_H_A':0.55,
         'old_local_search_domain_applied':False,'serialization':'.17g; every nonterminal coordinate line unchanged'},
       'declared_sources':10,'declared_new_cells':20,'reused_baseline_cells':60,'new_energy_calls':0,'output':str(out/'execution')}
    write(out/'manifest.json',m);print(json.dumps(pin(out/'manifest.json')))

def physical_source(path,base):
    """Existing cavity PDB parser policy, extended explicitly to native La sources."""
    atoms=[];seen=set();models=0
    for line in Path(path).read_text().splitlines():
        if line.startswith('MODEL '):models+=1
        if not line.startswith(('ATOM  ','HETATM')):continue
        a={'chain':line[21:22].strip(),'resnum':int(line[22:26]),'insertion_code':line[26:27].strip(),
           'resname':line[17:20].strip(),'atom':line[12:16].strip(),'altloc':line[16:17].strip(),
           'element':line[76:78].strip().capitalize()}
        require(a['element'] in ('H','C','N','O','S','Ca','La') and not a['altloc'],'Unsupported source element/alternate location')
        k=base.key(a);require(k not in seen,'Duplicate source identity');seen.add(k)
        atoms.append({'source':a,'element':a['element'],'xyz':[float(line[i:j]) for i,j in ((30,38),(38,46),(46,54))]})
    require(models<=1 and atoms,'Missing or multiple source models');return atoms

def contact_check(np,before,after,symbols,moved,identities,caps):
    """Same strict relative thresholds as prior local and source-contact admission."""
    bad=[];inherited=[];nearest={'physical_heavy':None,'physical_H_involving':None,'synthetic_cap':None}
    omitted_nearest={'heavy_heavy':None,'H_involving':None}
    for i in moved:
        other=[j for j in range(len(before)) if j!=i and not (j in moved and j<i)]
        old=np.linalg.norm(before[other]-before[i],axis=1);new=np.linalg.norm(after[other]-after[i],axis=1)
        limits=np.array([.55 if 'H' in (symbols[i],symbols[j]) else 1.0 for j in other])
        for k,j in enumerate(other):
            category='synthetic_cap' if j in caps else ('physical_H_involving' if 'H' in (symbols[i],symbols[j]) else 'physical_heavy')
            record={'moving_local_index':i,'other_union_index':j,'moving_identity':identities[i],'other_identity':identities[j],
                    'before_A':float(old[k]),'after_A':float(new[k]),'threshold_A':float(limits[k]),'category':category}
            if nearest[category] is None or new[k]<nearest[category]['after_A']:nearest[category]=record
            if j>=204:
                category_omitted='H_involving' if 'H' in (symbols[i],symbols[j]) else 'heavy_heavy'
                if omitted_nearest[category_omitted] is None or new[k]<omitted_nearest[category_omitted]['after_A']:
                    omitted_nearest[category_omitted]=record
            if old[k]<limits[k]:inherited.append(record)
            if new[k]<limits[k] and old[k]>=limits[k]:bad.append(record)
    return {'pass':not bad,'new_severe_contacts':bad,'inherited_short_contacts':inherited,'nearest_contacts':nearest,
            'nearest_omitted_source_contacts':omitted_nearest,'omitted_membership':'Exact full-source identities absent from the fixed204-atom scoring context; no distance-based omission',
            'policy':'new < threshold AND origin >= threshold; heavy1.0A / anyH0.55A; no invented covalent exclusions',
            'synthetic_caps_separate_from_physical_source':True}

def prepare(s,directory,m,np,base,Kinematics):
    r=s['record'];inp=s['inputs'];cp=read(inp['source_preparation']);state=read(inp['state_signature']);proton=read(inp['protonation'])
    maps={z:read(inp['maps'][z]) for z in METALS}
    require(maps['Ca']==maps['La'],'Paired physical maps differ')
    require(cp['source']==r['protonated_source']==proton['output'],'Source/protonation pin differs')
    require(not proton['detected_missing_residues'] and proton['repaired_missing_atom_count']==0 and proton['repaired_missing_terminal_atom_count']==0,'Unqualified source completeness')
    require(cp['new_atom_count']==204 and state['endpoint_charges']=={'Ca':-4,'La':-3} and state['multiplicities']=={'Ca':1,'La':1} and not state['water_inventory'],'Wrong204-atom dry paired state')
    require(['A:ASN290',0,None] in state['core_microstates'] and ['C:PQQ1',-3,'pqq_ox_3minus_v1'] in state['core_microstates'],'Asn/PQQ microstate differs')
    pair={z:base.xyz_atoms(check(inp['xyz'][z])) for z in METALS}
    require(all(len(v)==204 and v[0]['element']==z for z,v in pair.items()),'Wrong paired origins')
    require(pair['Ca'][0]['xyz']==pair['La'][0]['xyz'] and [(v['element'],v['xyz']) for v in pair['Ca'][1:]]==[(v['element'],v['xyz']) for v in pair['La'][1:]],'Paired origin coordinates differ')
    before=np.asarray([a['xyz'] for a in pair['Ca']]);symbols=[a['element'] for a in pair['Ca']]
    entries={x['qm_index']:x for x in cp['mapping']['source_to_qm']};require(set(entries)==set(range(1,204)),'Incomplete local identity map')
    named={x['source']['atom']:i for i,x in entries.items() if x.get('source',{}).get('chain')=='A' and x.get('source',{}).get('resnum')==290}
    require(set(named)=={'N','CA','C','O','CB','CG','OD1','ND2','H','HA','HB2','HB3','HD21','HD22'},'Incomplete or unexpected Asn atom inventory')
    require({name:symbols[named[name]] for name in MOVED}=={'OD1':'O','ND2':'N','HD21':'H','HD22':'H'},'Terminal amide identity differs')
    kin=Kinematics(maps['Ca']['context']);meta=kin.data['source_atom_metadata'];phys={x['atom']:i for i,x in enumerate(meta) if x and x['chain']=='A' and x['resnum']==290}
    mode,=[x for x in kin.modes if x['id']=='A/290/chi2'];mode_index=kin.modes.index(mode)
    require(mode['axis_indices']==[phys['CB'],phys['CG']] and set(mode['moving_indices'])=={phys[n] for n in (*MOVED,'CG')},'Actual source chi2 moves an unexpected group')
    bonds={frozenset((i,j)) for i,j in kin.data['bonds']}
    for a,b in (('CB','CG'),('CG','OD1'),('CG','ND2'),('ND2','HD21'),('ND2','HD22')):
        require(frozenset((phys[a],phys[b])) in bonds,'Missing physical covalent bond '+a+'-'+b)
    moving_phys={phys[n] for n in MOVED}
    require(all(not ((i in moving_phys)!=(j in moving_phys)) or (i==phys['CG'] or j==phys['CG']) for i,j in kin.data['bonds']),'Unsupported external bond from moving amide')
    hparents={tuple(h):tuple(p) for h,p in state['H_covalent_parents']}
    for h in ('HD21','HD22'):
        require(hparents[('A',290,'','ASN',h,'H')]==('A',290,'','ASN','ND2','N'),'Asn physical H parent differs')
    q=np.zeros(len(kin.modes));origin_physical,rendered,*_=kin.evaluate(q)
    require(np.max(np.abs(rendered-before))<=1e-8,'Original physical map does not replay')
    q[mode_index]=np.pi;physical,rendered,*_=kin.evaluate(q)
    ids=[named[n] for n in MOVED];fixed=[i for i in range(204) if i not in ids]
    anchor=before[named['CG']];axis=anchor-before[named['CB']];length=np.linalg.norm(axis);require(length>1e-10,'Degenerate CB-CG axis');axis/=length
    rotation=2*np.outer(axis,axis)-np.eye(3);after=before.copy();after[ids]=(before[ids]-anchor)@rotation.T+anchor
    require(np.max(np.abs(rotation@rotation.T-np.eye(3)))<=1e-10 and abs(np.linalg.det(rotation)-1)<=1e-10,'Improper amide rotation')
    require(np.max(np.abs(rendered-after))<=1e-8 and np.array_equal(after[fixed],before[fixed]),'Physical-map flip or fixed-coordinate agreement failed')
    physical_fixed=[i for i in range(len(physical)) if i not in mode['moving_indices']]
    require(np.array_equal(physical[physical_fixed],origin_physical[physical_fixed]),'Physical exterior moved')
    axis_roundoff=float(np.max(np.abs(physical[phys['CG']]-origin_physical[phys['CG']])))
    require(axis_roundoff<=1e-10 and np.array_equal(after[named['CG']],before[named['CG']]),'CG axis invariant failed')
    bond_error=max(abs(float(np.linalg.norm(physical[i]-physical[j])-np.linalg.norm(origin_physical[i]-origin_physical[j]))) for i,j in kin.data['bonds'])
    group=[named[n] for n in ('CB','CG',*MOVED)];pair_error=max(abs(float(np.linalg.norm(after[i]-after[j])-np.linalg.norm(before[i]-before[j]))) for i in group for j in group)
    require(bond_error<=1e-10 and pair_error<=1e-10,'Covalent or rigid-amide invariant failed')
    physical_atoms=physical_source(check(inp['protonated_source']),base);lookup={base.key(a['source']):a for a in physical_atoms}
    metal_atoms=[a for a in physical_atoms if a['element'] in METALS];require(len(metal_atoms)==1,'Ambiguous source-metal occupancy')
    active=base.key(metal_atoms[0]['source']);require(base.close(metal_atoms[0]['xyz'],before[0],1e-12),'Source active metal coordinate differs')
    used={active};identities=[{'kind':'exchanged_source_metal','source':metal_atoms[0]['source']}];caps=[];carver=[]
    for i in range(1,204):
        x=entries[i];require(base.close(x['xyz_A'],before[i],1e-12),'Local preparation coordinate differs')
        if x['kind']=='sigma_link_H':
            require(symbols[i]=='H' and all(base.key(x[k]) in lookup for k in ('retained','omitted')),'Cap source identity missing')
            require(all(not (x[k]['chain']=='A' and x[k]['resnum']==290 and x[k]['atom'] in MOVED) for k in ('retained','omitted')),'Moving cap anchor')
            caps.append(i);identities.append({'kind':'synthetic_cap','mapping':x})
        elif 'source' in x:
            k=base.key(x['source']);require(k in lookup and k not in used,'Duplicate or missing local/full-source identity')
            require(lookup[k]['element']==symbols[i] and base.close(lookup[k]['xyz'],before[i],1e-12),'Local/full-source element or frame differs')
            used.add(k);identities.append({'kind':'physical_source','source':x['source']})
        elif x['kind']=='opaque_cofactor' and x['element']=='H' and x['fragment']=='C:PQQ1' and x['name'] in ('H_N1','H_C3','H_C8'):
            carver.append(x['name']);identities.append({'kind':'unchanged_carver_PQQ_H','mapping':x})
        else:raise ValueError('Unsupported local physical identity')
    require(len(caps)==14 and set(carver)=={'H_N1','H_C3','H_C8'} and np.array_equal(after[caps],before[caps]),'Cap/cofactor H inventory or coordinates changed')
    exterior=[a for a in physical_atoms if base.key(a['source']) not in used]
    union_before=np.vstack((before,np.asarray([a['xyz'] for a in exterior])));union_after=union_before.copy();union_after[:204]=after
    identities.extend({'kind':'physical_source','source':a['source']} for a in exterior);union_symbols=symbols+[a['element'] for a in exterior]
    require(len(used)+len(exterior)==len(physical_atoms),'Full source was omitted or duplicated')
    source_index={base.key(x['source']):i for i,x in enumerate(identities) if 'source' in x}
    for i,x in enumerate(meta):
        if not x:continue
        k=base.key(x);require(k in source_index,'Missing physical kinematic atom from full source')
        j=source_index[k];require(np.max(np.abs(union_after[j]-physical[i]))<=1e-8,'Full-source physical graft mismatch')
    contacts=contact_check(np,union_before,union_after,union_symbols,ids,identities,set(caps))
    admission={'status':'prepared' if contacts['pass'] else 'unavailable','pass':contacts['pass'],'source_case_id':r['source_case_id'],
       'rotation_degrees':180,'axis_source_names':['CB','CG'],'rotation_matrix':rotation.tolist(),'axis_unit_vector':axis.tolist(),
       'moved_indices_zero_based':ids,'named_Asn_indices_zero_based':named,'covalent_max_error_A':bond_error,'rigid_amide_max_error_A':pair_error,
       'physical_map_max_error_A':float(np.max(np.abs(rendered-after))),'maximum_heavy_displacement_A':float(np.linalg.norm(after[[named['OD1'],named['ND2']]]-before[[named['OD1'],named['ND2']]],axis=1).max()),
       'kinematic_CG_axis_roundoff_A':axis_roundoff,'emitted_CG_fixed_exact':True,
       'all_other_coordinate_lines_fixed':True,'CA_backbone_stereocentre_and_substituents_fixed':True,'cap_indices':caps,
       'all_cap_anchors_fixed':True,'full_source_atom_count':len(physical_atoms),'contact_union_atom_count':len(union_before),
       'full_source_identity_count':len(source_index),'appended_source_atoms':len(exterior),'all_physical_H_included':True,
       'source':r['protonated_source'],'contacts':contacts,'old_0p8_local_domain_applied':False,'continuous_path_admission':False}
    write(directory/'GEOMETRY_ADMISSION.json',admission);ap=pin(directory/'GEOMETRY_ADMISSION.json')
    write(directory/'CONTACT_IDENTITIES.json',{'zero_based':True,'centers':identities,'source':r['protonated_source'],'moving_indices':ids})
    require(contacts['pass'],'New severe contact; see GEOMETRY_ADMISSION.json')
    matrix={z:{} for z in METALS}
    for z in METALS:
        original=check(r['matrix'][z]['origin']['xyz']);lines=original.read_text().splitlines(keepends=True)
        for i in ids:lines[i+2]=symbols[i]+' '+' '.join(format(v,'.17g') for v in after[i])+'\n'
        xp=directory/(z+'_amide_flip.xyz')
        with xp.open('x') as f:f.write(''.join(lines))
        actual=base.xyz_atoms(xp);require(all(actual[i]['line']==pair[z][i]['line'] for i in fixed),'A fixed coordinate line changed')
        require(np.array_equal(np.asarray([a['xyz'] for a in actual]),after),'Serialized coordinates differ')
        matrix[z]['amide_flip']={'status':'prepared','xyz':pin(xp),'charge':state['endpoint_charges'][z],'multiplicity':1,'atom_count':204,
          'geometry_admission':ap,'original_xyz':r['matrix'][z]['origin']['xyz'],'original_MACE':r['matrix'][z]['origin']['MACE']}
    return {'status':'prepared','reason':None,'case_id':r['source_case_id'],'source_case_id':r['source_case_id'],'case_key':r['case_key'],'record_id':r['record_id'],
       'selection_id':r['selection_id'],'conditioning':r['conditioning'],'sample':int(r['source_case_id'].rsplit('sample-',1)[1]),
       **{k:r[k] for k in ('maps','state_signature','source_preparation','protonated_source')},'geometry_admission':ap,
       'moved_indices_zero_based':ids,'named_Asn_indices_zero_based':named,'matrix':matrix,'baseline_matrix':r['matrix']}

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports')
    mp={'path':str(Path(path).resolve()),'sha256':sha};m=read(mp);base=load('_asn_source_helpers',m['source_parser_reference']);place=base.allocation()
    require(m['protocol']==PROTOCOL and m['implementation']['PREPARE.py']==pin(__file__),'Wrong frozen preparation')
    require(m['settings']['angle_degrees']==180 and m['settings']['moved_atom_names']==list(MOVED) and not m['settings']['old_local_search_domain_applied'],'Discrete protocol changed')
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time();rows=[]
    receipt={'status':'failed','manifest':mp,'placement':place,'new_energy_calls':0};write(out/'START.json',receipt)
    try:
        for p in (m['inventory'],m['plan'],m['python'],m['contact_policy_reference'],*m['implementation'].values(),*m['source_modules'].values()):check(p)
        require(pin(sys.executable)==m['python'],'Wrong numerical interpreter')
        for p in m['input_snapshots'].values():check(p['original']);check(p['snapshot'])
        sys.path.insert(0,m['source_snapshot'])
        import numpy as np
        import mace_site_kinematics as kin
        require(pin(kin.__file__)==m['source_modules']['mace_site_kinematics.py'],'Wrong kinematics helper')
        inv=read(m['inventory']);require(len(m['sources'])==10 and [s['record'] for s in m['sources']]==inv['rows'],'Changed complete input population')
        for s in m['sources']:
            r=s['record'];d=out/r['source_case_id'];d.mkdir()
            row={'status':'unavailable','reason':'not prepared','source_case_id':r['source_case_id'],'case_id':r['source_case_id'],
                 'record_id':r['record_id'],'case_key':r['case_key'],'selection_id':r['selection_id'],'conditioning':r['conditioning'],
                 **{k:r[k] for k in ('maps','state_signature','source_preparation','protonated_source')},'baseline_matrix':r['matrix'],
                 'matrix':{z:{'amide_flip':{'status':'unavailable','xyz':None}} for z in METALS}}
            try:row=prepare(s,d,m,np,base,kin.Kinematics)
            except Exception as e:
                row.update(reason=f'{type(e).__name__}: {e}',traceback=traceback.format_exc())
                if (d/'GEOMETRY_ADMISSION.json').is_file():row['geometry_admission']=pin(d/'GEOMETRY_ADMISSION.json')
            rows.append(row);write(d/'PREPARATION.json',row)
        receipt['status']='complete' if all(r['status']=='prepared' for r in rows) else 'complete_with_unavailable'
    except BaseException as e:receipt.update(error=str(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'manifest':mp,'inventory':m['inventory'],'plan':m['plan'],'terminal':True,
          'status':receipt['status'],'declared_sources':10,'declared_new_cells':20,'reused_baseline_cells':60,
          'prepared_cells':2*sum(r['status']=='prepared' for r in rows),'rows':rows,'execution':pin(out/'EXECUTION.json'),'new_energy_calls':0})

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('freeze');a.add_argument('--output',required=True)
    a=s.add_parser('execute');a.add_argument('--manifest',required=True);a.add_argument('--manifest-sha',required=True)
    a=p.parse_args()
    if a.action=='freeze':freeze(a.output)
    else:execute(a.manifest,a.manifest_sha)
