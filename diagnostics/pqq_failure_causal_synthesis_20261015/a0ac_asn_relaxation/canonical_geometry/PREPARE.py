#!/usr/bin/env python3
"""Canonical25 metadata freeze and Slurm-only fixed terminal-amide geometry."""
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
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_geometry'
PROTOCOL='canonical25_anchor_Asn_fixed_terminal_amide_180_geometry_v1'
SCOPE_SHA='e0f7c3d11e7e92f618dd2f10b83eca3616b5e4a4b9b4840077124d701d1c2ccf'
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

def same_res(x,a):
    return all(x.get(k,'')==a[k] for k in ('chain','resnum','insertion_code','resname'))

def identity(x):
    return [x['chain'],x['resnum'],x.get('insertion_code',''),x['resname'],x['atom'],x['element']]

def base_row(s):
    r=s['record'];g=s['graph']
    return {**{k:r[k] for k in ('source_case_id','case_id','record_id','selection_id','protein','maps','state_signature','source_preparation','protonated_source')},
      'role':g['role'],'atom_count':g['atom_count'],'endpoint_charges':g['endpoint_charges'],'water_inventory':g['water_inventory'],
      'terminal_mode':g['terminal_mode'],'terminal_mode_selected_by_metal':g['terminal_mode_selected_by_metal'],
      'original_four_modes':g['original_four_modes'],'status':'unavailable','reason':'not prepared',
      'matrix':{z:{'amide_flip':{'status':'unavailable','xyz':None}} for z in METALS}}

def freeze(output):
    """Read/hash/copy metadata only; no coordinate parsing or scientific imports."""
    scope=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_scope_v1'
    rp=pin(scope/'execution/RESULT.json');require(rp['sha256']==SCOPE_SHA,'Wrong completed graph inventory')
    result=read(rp);sm=read(result['manifest']);execution=pin(scope/'execution/EXECUTION.json')
    require(read(execution)['status']=='complete' and read(execution)['manifest']==result['manifest'],'Graph execution incomplete')
    require(result['status_counts']=={'identity_and_graph_supported':25} and len(sm['rows'])==len(result['rows'])==25,'Wrong canonical denominator/support')
    require(sm['canonical_collection']==result['canonical_collection'],'Canonical authority differs')
    require([r['record_id'] for r in result['rows']]==[r['record_id'] for r in read(sm['canonical_collection'])['rows']],'Canonical order differs')
    oldroot=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_preparation_v1'
    oldmp=pin(oldroot/'manifest.json');old=read(oldmp);oldrp=pin(oldroot/'execution/PREPARATION.json');oldresult=read(oldrp)
    require(oldresult['manifest']==oldmp and oldresult['status']=='complete','Original flip preparation incomplete')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir();(out/'inputs').mkdir()
    snapshots={};sources=[];regressions=[]
    def snapshot(p):
        source=check(p);key=p['sha256']
        if key not in snapshots:
            dest=out/'inputs'/(key+source.suffix);shutil.copyfile(source,dest)
            snapshots[key]={'original':p,'snapshot':pin(dest)}
        return snapshots[key]['snapshot']
    for r,g in zip(sm['rows'],result['rows']):
        require(r['record_id']==g['record_id'] and r['source_case_id']==g['source_case_id'] and g['role']['resname']=='ASN','Source/neutral Asn join differs')
        inputs={k:snapshot(r[k]) for k in ('source_preparation','state_signature','protonated_source','protonation')}
        inputs['maps']={z:snapshot(r['maps'][z]) for z in METALS};inputs['xyz']={z:snapshot(r['original_xyz'][z]) for z in METALS}
        sources.append({'record':r,'graph':g,'inputs':inputs})
        for prior in oldresult['rows']:
            if prior['source_case_id']!=r['source_case_id']:continue
            require(prior['status']=='prepared' and all(prior[k]==r[k] for k in ('maps','state_signature','source_preparation','protonated_source')),'Same-source prior flip is incompatible')
            require(all(prior['matrix'][z]['amide_flip']['original_xyz']==r['original_xyz'][z] for z in METALS),'Prior flip origins differ')
            regressions.append({'source_case_id':r['source_case_id'],'prior':prior})
    require(len(regressions)==1 and regressions[0]['source_case_id']=='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-3','Exact real canonical A0AC regression unavailable')
    implementation={}
    for name in ('PREPARE.py','RUN.sbatch','PLAN.md'):
        dest=out/'implementation'/name;shutil.copyfile(HERE/name,dest);implementation[name]=pin(dest)
    m={'protocol':PROTOCOL,'implementation':implementation,'scope_result':rp,'scope_manifest':result['manifest'],'scope_execution':execution,
      'canonical_collection':sm['canonical_collection'],'old_preparation_manifest':oldmp,'old_preparation':oldrp,
      'old_helper':old['implementation']['PREPARE.py'],'source_parser_reference':old['source_parser_reference'],
      'contact_policy_reference':old['contact_policy_reference'],'python':old['python'],'source_snapshot':old['source_snapshot'],
      'source_modules':{k:old['source_modules'][k] for k in ('mace_site_kinematics.py','affordable_common.py','affordable_response.py')},
      'sources':sources,'input_snapshots':snapshots,'regression':regressions[0],'settings':old['settings'],
      'declared_sources':25,'declared_new_geometries':50,'new_energy_calls':0,'output':str(out/'execution')}
    write(out/'manifest.json',m);print(json.dumps(pin(out/'manifest.json')))

def contact_check(np,before,after,symbols,moved,identities,caps,ncore):
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
            if j>=ncore:
                category_omitted='H_involving' if 'H' in (symbols[i],symbols[j]) else 'heavy_heavy'
                if omitted_nearest[category_omitted] is None or new[k]<omitted_nearest[category_omitted]['after_A']:
                    omitted_nearest[category_omitted]=record
            if old[k]<limits[k]:inherited.append(record)
            if new[k]<limits[k] and old[k]>=limits[k]:bad.append(record)
    return {'pass':not bad,'new_severe_contacts':bad,'inherited_short_contacts':inherited,'nearest_contacts':nearest,
            'nearest_omitted_source_contacts':omitted_nearest,'omitted_membership':'Exact full-source identities absent from the pinned source-specific scoring context; no distance-based omission',
            'policy':'new < threshold AND origin >= threshold; heavy1.0A / anyH0.55A; no invented covalent exclusions',
            'synthetic_caps_separate_from_physical_source':True}

def prepare(s,directory,m,np,base,Kinematics,old):
    r=s['record'];g=s['graph'];a=g['role'];n=g['atom_count'];inp=s['inputs'];cp=read(inp['source_preparation']);state=read(inp['state_signature']);proton=read(inp['protonation'])
    maps={z:read(inp['maps'][z]) for z in METALS}
    require(maps['Ca']==maps['La'],'Paired physical maps differ')
    require(cp['source']==r['protonated_source']==proton['output'],'Source/protonation pin differs')
    require(not proton['detected_missing_residues'] and proton['repaired_missing_atom_count']==0 and proton['repaired_missing_terminal_atom_count']==0,'Unqualified source completeness')
    require(cp['new_atom_count']==n and state['endpoint_charges']==g['endpoint_charges'] and state['multiplicities']==g['multiplicities']=={'Ca':1,'La':1} and state['water_inventory']==g['water_inventory'],'Pinned paired state/count/waters differ')
    require(g['status']=='identity_and_graph_supported' and a['resname']=='ASN' and g['amide_microstate'] in state['core_microstates'] and g['amide_microstate'][1:]==[0,None],'Neutral declared anchor Asn differs')
    pair={z:base.xyz_atoms(check(inp['xyz'][z])) for z in METALS}
    require(all(len(v)==n and v[0]['element']==z for z,v in pair.items()),'Wrong paired origins')
    require(pair['Ca'][0]['xyz']==pair['La'][0]['xyz'] and [(v['element'],v['xyz']) for v in pair['Ca'][1:]]==[(v['element'],v['xyz']) for v in pair['La'][1:]],'Paired origin coordinates differ')
    before=np.asarray([a['xyz'] for a in pair['Ca']]);symbols=[a['element'] for a in pair['Ca']]
    entries={x['qm_index']:x for x in cp['mapping']['source_to_qm']};require(set(entries)==set(range(1,n)),'Incomplete local identity map')
    named={x['source']['atom']:i for i,x in entries.items() if same_res(x.get('source',{}),a)}
    require(set(named)=={'N','CA','C','O','CB','CG','OD1','ND2','H','HA','HB2','HB3','HD21','HD22'},'Incomplete or unexpected Asn atom inventory')
    require({name:symbols[named[name]] for name in MOVED}=={'OD1':'O','ND2':'N','HD21':'H','HD22':'H'},'Terminal amide identity differs')
    kin=Kinematics(maps['Ca']['context']);meta=kin.data['source_atom_metadata'];phys={x['atom']:i for i,x in enumerate(meta) if x and same_res(x,a)}
    require(named==g['qm_indices'] and phys==g['physical_indices'],'Declared named mapping differs')
    mode,=[x for x in kin.modes if x['id']==g['terminal_mode']['id']];mode_index=kin.modes.index(mode)
    require(mode==g['terminal_mode'],'Existing terminal mode differs')
    require(mode['axis_indices']==[phys['CB'],phys['CG']] and set(mode['moving_indices'])=={phys[n] for n in (*MOVED,'CG')},'Actual source chi2 moves an unexpected group')
    bonds={frozenset((i,j)) for i,j in kin.data['bonds']}
    for left,right in (('CB','CG'),('CG','OD1'),('CG','ND2'),('ND2','HD21'),('ND2','HD22')):
        require(frozenset((phys[left],phys[right])) in bonds,'Missing physical covalent bond '+left+'-'+right)
    moving_phys={phys[n] for n in MOVED}
    require(all(not ((i in moving_phys)!=(j in moving_phys)) or (i==phys['CG'] or j==phys['CG']) for i,j in kin.data['bonds']),'Unsupported external bond from moving amide')
    hparents={tuple(h):tuple(p) for h,p in state['H_covalent_parents']}
    for h in ('HD21','HD22'):
        require(hparents[tuple([a[k] for k in ('chain','resnum','insertion_code','resname')]+[h,'H'])]==tuple([a[k] for k in ('chain','resnum','insertion_code','resname')]+['ND2','N']),'Asn physical H parent differs')
    q=np.zeros(len(kin.modes));origin_physical,rendered,*_=kin.evaluate(q)
    require(np.max(np.abs(rendered-before))<=1e-8,'Original physical map does not replay')
    q[mode_index]=np.pi;physical,rendered,*_=kin.evaluate(q)
    ids=[named[n] for n in MOVED];fixed=[i for i in range(n) if i not in ids]
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
    physical_atoms=old.physical_source(check(inp['protonated_source']),base);lookup={base.key(a['source']):a for a in physical_atoms}
    metal_atoms=[a for a in physical_atoms if a['element'] in METALS];require(len(metal_atoms)==1,'Ambiguous source-metal occupancy')
    active=base.key(metal_atoms[0]['source']);require(base.close(metal_atoms[0]['xyz'],before[0],1e-12),'Source active metal coordinate differs')
    used={active};identities=[{'kind':'exchanged_source_metal','source':metal_atoms[0]['source']}];caps=[];carver=[]
    for i in range(1,n):
        x=entries[i];require(base.close(x['xyz_A'],before[i],1e-12),'Local preparation coordinate differs')
        if x['kind']=='sigma_link_H':
            require(symbols[i]=='H' and all(base.key(x[k]) in lookup for k in ('retained','omitted')),'Cap source identity missing')
            require(all(not (same_res(x[k],a) and x[k]['atom'] in MOVED) for k in ('retained','omitted')),'Moving cap anchor')
            caps.append(i);identities.append({'kind':'synthetic_cap','mapping':x})
        elif 'source' in x:
            k=base.key(x['source']);require(k in lookup and k not in used,'Duplicate or missing local/full-source identity')
            require(lookup[k]['element']==symbols[i] and base.close(lookup[k]['xyz'],before[i],1e-12),'Local/full-source element or frame differs')
            used.add(k);identities.append({'kind':'physical_source','source':x['source']})
        elif x['kind']=='opaque_cofactor' and x['element']=='H' and ['cofactor',x['fragment'],x['name'],'H'] in state['source_atoms_and_cofactor']:
            carver.append([x['fragment'],x['name']]);identities.append({'kind':'unchanged_carver_cofactor_H','mapping':x})
        else:raise ValueError('Unsupported local physical identity')
    actual_caps=[[identity(entries[i]['retained']),identity(entries[i]['omitted']),entries[i]['length_A']] for i in caps]
    require(sorted(map(json.dumps,actual_caps))==sorted(map(json.dumps,state['caps'])) and np.array_equal(after[caps],before[caps]),'Cap identities/lengths or coordinates differ')
    expected_carver=[[x['fragment'],x['name']] for x in cp['mapping']['source_to_qm']
       if x['kind']=='opaque_cofactor' and 'source' not in x and x['element']=='H'
       and ['cofactor',x['fragment'],x['name'],'H'] in state['source_atoms_and_cofactor']]
    require(carver==expected_carver,'Cofactor H inventory differs')
    exterior=[a for a in physical_atoms if base.key(a['source']) not in used]
    union_before=np.vstack((before,np.asarray([a['xyz'] for a in exterior])));union_after=union_before.copy();union_after[:n]=after
    identities.extend({'kind':'physical_source','source':a['source']} for a in exterior);union_symbols=symbols+[a['element'] for a in exterior]
    require(len(used)+len(exterior)==len(physical_atoms),'Full source was omitted or duplicated')
    source_index={base.key(x['source']):i for i,x in enumerate(identities) if 'source' in x}
    for i,x in enumerate(meta):
        if not x:continue
        k=base.key(x);require(k in source_index,'Missing physical kinematic atom from full source')
        j=source_index[k];require(np.max(np.abs(union_after[j]-physical[i]))<=1e-8,'Full-source physical graft mismatch')
    contacts=contact_check(np,union_before,union_after,union_symbols,ids,identities,set(caps),n)
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
        original=check(r['original_xyz'][z]);lines=original.read_text().splitlines(keepends=True)
        for i in ids:lines[i+2]=symbols[i]+' '+' '.join(format(v,'.17g') for v in after[i])+'\n'
        xp=directory/(z+'_amide_flip.xyz')
        with xp.open('x') as f:f.write(''.join(lines))
        actual=base.xyz_atoms(xp);require(all(actual[i]['line']==pair[z][i]['line'] for i in fixed),'A fixed coordinate line changed')
        require(np.array_equal(np.asarray([a['xyz'] for a in actual]),after),'Serialized coordinates differ')
        matrix[z]['amide_flip']={'status':'prepared','xyz':pin(xp),'charge':state['endpoint_charges'][z],'multiplicity':1,'atom_count':n,
          'geometry_admission':ap,'original_xyz':r['original_xyz'][z]}
    return {**base_row(s),'status':'prepared','reason':None,'geometry_admission':ap,
       'moved_indices_zero_based':ids,'named_Asn_indices_zero_based':named,'matrix':matrix}

def regression(row,prior):
    require(row['status']==prior['status']=='prepared','Real A0AC regression did not prepare')
    for z in METALS:
        require(check(row['matrix'][z]['amide_flip']['xyz']).read_bytes()==check(prior['matrix'][z]['amide_flip']['xyz']).read_bytes(),'Real A0AC emitted XYZ differs')
    actual=read(row['geometry_admission']);old=read(prior['geometry_admission'])
    fields=('rotation_matrix','axis_unit_vector','moved_indices_zero_based','named_Asn_indices_zero_based','covalent_max_error_A',
       'rigid_amide_max_error_A','physical_map_max_error_A','maximum_heavy_displacement_A','kinematic_CG_axis_roundoff_A',
       'cap_indices','full_source_atom_count','contact_union_atom_count','full_source_identity_count','appended_source_atoms')
    require(all(actual[k]==old[k] for k in fields),'Real A0AC geometric admission differs')
    def contacts(x):
        if isinstance(x,dict):return {k:contacts(v) for k,v in x.items() if k!='omitted_membership'}
        if isinstance(x,list):return [contacts(v) for v in x]
        return 'unchanged_carver_cofactor_H' if x=='unchanged_carver_PQQ_H' else x
    require(contacts(actual['contacts'])==contacts(old['contacts']),'Real A0AC contact admission differs')
    return {'status':'passed','source_case_id':row['source_case_id'],'prior_admission':prior['geometry_admission'],
       'new_admission':row['geometry_admission'],'exact_paired_XYZ_bytes':True,'exact_geometric_and_contact_evidence':True,
       'only_report_label_changes':['source-specific context size wording','carver cofactor H identity label']}

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports')
    mp={'path':str(Path(path).resolve()),'sha256':sha};m=read(mp)
    require(m['protocol']==PROTOCOL and m['implementation']['PREPARE.py']==pin(__file__),'Wrong frozen implementation')
    require(m['declared_sources']==25 and m['declared_new_geometries']==50 and len(m['sources'])==25,'Wrong declared geometry denominator')
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    rows=[base_row(s) for s in m['sources']];receipt={'status':'failed','manifest':mp,'new_energy_calls':0,'job_id':os.environ['SLURM_JOB_ID']}
    try:
        for k in ('scope_result','scope_manifest','scope_execution','canonical_collection','old_preparation_manifest','old_preparation','old_helper','source_parser_reference','contact_policy_reference','python'):check(m[k])
        for p in (*m['implementation'].values(),*m['source_modules'].values()):check(p)
        for v in m['input_snapshots'].values():check(v['original']);check(v['snapshot'])
        require(m['scope_result']['sha256']==SCOPE_SHA and m['settings']==read(m['old_preparation_manifest'])['settings'],'Frozen graph/settings changed')
        require([s['record'] for s in m['sources']]==read(m['scope_manifest'])['rows'] and [s['graph'] for s in m['sources']]==read(m['scope_result'])['rows'],'Frozen source graph population differs')
        require(pin(sys.executable)==m['python'],'Wrong numerical interpreter')
        base=load('_canonical_source_parser',m['source_parser_reference']);receipt['placement']=base.allocation()
        old=load('_canonical_existing_amide_helpers',m['old_helper'])
        sys.path.insert(0,m['source_snapshot'])
        import numpy as np
        import mace_site_kinematics as kin
        for name,p in m['source_modules'].items():require(pin(sys.modules[Path(name).stem].__file__)==p,'Imported source differs: '+name)
        first=next(i for i,s in enumerate(m['sources']) if s['record']['source_case_id']==m['regression']['source_case_id'])
        for i in (first,*[j for j in range(25) if j!=first]):
            s=m['sources'][i];d=out/s['record']['source_case_id'];d.mkdir()
            try:rows[i]=prepare(s,d,m,np,base,kin.Kinematics,old)
            except Exception as e:
                rows[i].update(reason=f'{type(e).__name__}: {e}',traceback=traceback.format_exc())
                if (d/'GEOMETRY_ADMISSION.json').is_file():rows[i]['geometry_admission']=pin(d/'GEOMETRY_ADMISSION.json')
            write(d/'PREPARATION.json',rows[i])
            if i==first:
                write(out/'A0AC_REGRESSION.json',regression(rows[i],m['regression']['prior']))
                receipt['actual_A0AC_regression']=pin(out/'A0AC_REGRESSION.json')
        receipt['status']='complete' if all(r['status']=='prepared' for r in rows) else 'complete_with_unavailable'
    except BaseException as e:
        receipt.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'manifest':mp,'scope_result':m['scope_result'],'canonical_collection':m['canonical_collection'],
          'terminal':True,'status':receipt['status'],'declared_sources':25,'declared_new_geometries':50,
          'prepared_cells':2*sum(r['status']=='prepared' for r in rows),'rows':rows,'execution':pin(out/'EXECUTION.json'),
          'new_energy_calls':0,'new_optimization':False,'active_mode_reselection':False,'calibration_or_scoring':False})

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('freeze');a.add_argument('--output',required=True)
    a=sub.add_parser('execute');a.add_argument('--manifest',required=True);a.add_argument('--manifest-sha',required=True)
    args=p.parse_args()
    if args.action=='freeze':freeze(args.output)
    else:execute(args.manifest,args.manifest_sha)
