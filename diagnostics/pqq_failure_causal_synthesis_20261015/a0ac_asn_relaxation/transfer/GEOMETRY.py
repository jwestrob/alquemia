"""Versioned copy of canonical fixed-flip prepare; only active-ion routing differs.
Source: a430adea58947507019e9d3ef9456f44b2d9633caf8ed1cb7d4767905bac9d11.
Canonical geometry/contact functions are supplied unchanged by the caller.
"""
from COMMON import *
METALS=('Ca','La');MOVED=('OD1','ND2','HD21','HD22')
# The caller supplies frozen same_res/identity/base_row/contact_check functions.
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
    selector=r['active_selector']
    def matches(src,sel):
        return all(src.get(k,'')==sel.get(k,'') for k in ('chain','resnum','resname','atom')) and src.get('insertion_code','')==sel.get('insertion_code',sel.get('icode','')) and ('element' not in sel or src['element']==sel['element'])
    chosen=[x for x in physical_atoms if matches(x['source'],selector)]
    require(len(chosen)==1 and chosen[0]['element'] in METALS,'Recorded active source ion is absent or ambiguous')
    active_atom=chosen[0];active=base.key(active_atom['source'])
    require(base.close(active_atom['xyz'],before[0],1e-12),'Recorded source active metal coordinate differs')
    spectator=None
    if r.get('spectator_admission'):
        sp=read(r['spectator_admission'])
        require(sp['status']=='passed' and sp['passed'] and sp['full_protonated_source']==r['protonated_source'] and selector==sp['active_selector'],'Spectator admission/source differs')
        spectators=[x for x in physical_atoms if matches(x['source'],sp['spectator_selector'])]
        require(len(spectators)==1 and base.key(spectators[0]['source'])!=active,'Recorded spectator missing or aliases active ion')
        spectator=spectators[0]
    used={active};identities=[{'kind':'exchanged_source_metal','source':active_atom['source']}];caps=[];carver=[]
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
    spectator_evidence=None
    if spectator is not None:
        si=source_index[base.key(spectator['source'])]
        require(si!=0 and np.array_equal(union_before[si],union_after[si]) and base.close(union_after[si],spectator['xyz'],1e-12),'Spectator coordinate/identity changed')
        require(union_symbols[si]==spectator['element'],'Spectator element changed')
        spectator_evidence={'admission':r['spectator_admission'],'source':spectator['source'],'union_index':si,'membership':'core' if si<n else 'exterior','coordinates_unchanged':True,'element_unchanged':True}
    contacts=contact_check(np,union_before,union_after,union_symbols,ids,identities,set(caps),n)
    admission={'status':'prepared' if contacts['pass'] else 'unavailable','pass':contacts['pass'],'source_case_id':r['source_case_id'],
       'rotation_degrees':180,'axis_source_names':['CB','CG'],'rotation_matrix':rotation.tolist(),'axis_unit_vector':axis.tolist(),
       'moved_indices_zero_based':ids,'named_Asn_indices_zero_based':named,'covalent_max_error_A':bond_error,'rigid_amide_max_error_A':pair_error,
       'physical_map_max_error_A':float(np.max(np.abs(rendered-after))),'maximum_heavy_displacement_A':float(np.linalg.norm(after[[named['OD1'],named['ND2']]]-before[[named['OD1'],named['ND2']]],axis=1).max()),
       'kinematic_CG_axis_roundoff_A':axis_roundoff,'emitted_CG_fixed_exact':True,
       'all_other_coordinate_lines_fixed':True,'CA_backbone_stereocentre_and_substituents_fixed':True,'cap_indices':caps,
       'all_cap_anchors_fixed':True,'full_source_atom_count':len(physical_atoms),'contact_union_atom_count':len(union_before),
       'full_source_identity_count':len(source_index),'appended_source_atoms':len(exterior),'all_physical_H_included':True,
       'source':r['protonated_source'],'active_selector':selector,'spectator_preservation':spectator_evidence,'contacts':contacts,'old_0p8_local_domain_applied':False,'continuous_path_admission':False}
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
