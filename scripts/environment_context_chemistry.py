"""Versioned complete thiol/disulfide neighbors; original supported contexts exact."""
from __future__ import annotations
from collections import Counter
import numpy as np
from affordable_common import InvalidArtifact,record
from second_shell_context import atom_key,cap_key,discover,sidechain,expansion

POLICY_ID='source_graph_complete_polar_neighbors_with_source_disulfides_v2'


def cysteine_units(state,fragments):
    """Use the installed source-PDB topology reader, never residue-number+1."""
    import openmm.app.topology
    from openmm.app import PDBFile
    g=state['graph'];wanted=[r for kind,r in fragments if kind=='sidechain' and g.residues[r].canonical_resname=='CYS']
    if not wanted:return {},[]
    top=PDBFile(str(g.source)).topology;ss={};metadata=[]
    def key(atom):
        r=atom.residue
        source=dict(chain=r.chain.id,resnum=int(r.id),resname=r.name,insertion_code=r.insertionCode.strip())
        return g.key(g.locate(source).key,atom.name)
    for a,b in top.bonds():
        if a.name!='SG' or b.name!='SG':continue
        x,y=key(a),key(b)
        if x in ss or y in ss:raise InvalidArtifact('multiply bonded source disulfide')
        ss[x]=y;ss[y]=x
    units={}
    for rkey in wanted:
        observed=Counter(p[2] for h,p in g.hparents.items() if p[:2]==rkey and p[2] in ('CB','SG'))
        sg=g.key(rkey,'SG')
        if observed=={'CB':2,'SG':1}:
            if sg in ss:raise InvalidArtifact('thiol hydrogen conflicts with source disulfide')
            units[rkey]={g.key(rkey,'CB'),sg};continue
        if observed!={'CB':2} or sg not in ss:raise InvalidArtifact('unsupported cysteine protonation/connectivity')
        other=ss[sg];partner=other[:2]
        partner_h=Counter(p[2] for h,p in g.hparents.items() if p[:2]==partner and p[2] in ('CB','SG'))
        if g.residues[partner].canonical_resname!='CYS' or partner_h!={'CB':2}:raise InvalidArtifact('unsupported disulfide partner')
        distance=float(g.atoms[sg].pos.dist(g.atoms[other].pos))
        if not 1.8<=distance<=2.3:raise InvalidArtifact('source disulfide outside declared covalent-distance range')
        if not any(a['atoms']==[g.meta[k] for k in sorted((sg,other))] for a in metadata):
            metadata.append({'atoms':[g.meta[k] for k in sorted((sg,other))],'distance_A':distance,
                             'formal_charge':0,'source_topology_reader':record(openmm.app.topology.__file__),
                             'rule':'existing_OpenMM_PDBFile_disulfide_bond_plus_source_H_and_distance_checks'})
        g.add_edge(sg,other)
        units[rkey]={g.key(rkey,'CB'),sg,g.key(partner,'CB'),other}
    return units,metadata


def complete_expansion(state):
    g=state['graph'];anchors,fragments,contacts,excluded=discover(state)
    units,disulfides=cysteine_units(state,fragments)
    if not units:
        # All established supported chemistry is literally the existing function.
        return expansion(state)
    selected=set(state['selected']);charges=dict(state['charges']);added=[]
    for kind,rkey in sorted(fragments):
        if kind=='peptide':selected|=g.peptide(rkey);charge=0
        elif rkey in units:
            selected|=units[rkey];charge=0
            for k in units[rkey]:charges[k[:2]]=0
        else:
            nodes,charge=sidechain(g,rkey);selected|=nodes;charges[rkey]=charge
        added.append({'kind':kind,'source':g.residues[rkey].source_dict(),'formal_charge':charge})
    atoms,mapping=g.materialize(selected)
    for a in state['opaque']:
        mapping['source_to_qm'].append(dict(a,qm_index=len(atoms)+1));atoms.append(state['original']['La'][a['qm_index']])
    for a in state['maps']:
        if a['kind']=='source' and atom_key(a['source']) in state['water_keys']:
            mapping['source_to_qm'].append(dict(a,qm_index=len(atoms)+1));atoms.append(state['original']['La'][a['qm_index']])
    caps={cap_key(a):a['qm_index'] for a in state['maps'] if a['kind']=='sigma_link_H'}
    opaque={(a['fragment'],a['name']):a['qm_index'] for a in state['opaque']};index={0:0};new_sources=[]
    for a in mapping['source_to_qm']:
        if a['kind']=='opaque_cofactor':i=opaque[(a['fragment'],a['name'])]
        elif a['kind']=='source':
            k=atom_key(a['source']);i=state['old_sources'].get(k)
            if i is None:new_sources.append(a['source'])
        else:i=caps.get(cap_key(a))
        if i is not None:index[i]=a['qm_index']
    if set(state['old_sources'].values())-set(index):raise InvalidArtifact('source context lost original core atoms')
    paired={}
    for metal,original in state['original'].items():
        rows=[original[0]]+list(atoms)
        for i,j in index.items():rows[j]=original[i]
        for a in mapping['source_to_qm']:a['xyz_A']=list(rows[a['qm_index']][1:])
        coords=np.array([a[1:] for a in rows]);d=np.linalg.norm(coords[:,None]-coords[None,:],axis=2);np.fill_diagonal(d,np.inf)
        if d.min()<.45:raise InvalidArtifact('overlapping complete context atoms/caps')
        paired[metal]=rows
    audit={'anchors':[g.meta[k] for k in sorted(anchors)],'added_fragments':added,'contacts':contacts,'excluded_contacts':excluded,
           'mapping':mapping,'core_to_context':index,'removed_caps':[a for a in state['maps'] if a['kind']=='sigma_link_H' and a['qm_index'] not in index],
           'added_source_atoms':new_sources,'added_formal_charge':sum(charges.values())-sum(state['charges'].values()),
           'source':state['source'],'new_atom_count':len(paired['La']),'original_atom_count':len(state['original']['La']),
           'core_source_coordinates_unchanged':True,'water_inventory_unchanged':True,'mapping_coordinate_endpoint':'La',
           'paired_differing_nonmetal_indices':[i for i in range(1,len(paired['La'])) if paired['Ca'][i]!=paired['La'][i]],
           'coordinate_note':'core source atoms and surviving caps exact; source-complete thiol/disulfide sidechain units',
           'disulfide_closures':disulfides,'cysteine_policy':POLICY_ID}
    return paired,audit
