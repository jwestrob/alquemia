"""Source-backed contextual water preparation using the established shell policy."""
from __future__ import annotations
from collections import defaultdict
import numpy as np
import hydration_square as square
from affordable_common import InvalidArtifact, paired, read_json, record, verify, xyz


def residue_id(meta):return tuple(meta[k] for k in ('chain','resnum','insertion_code','resname'))
def source_id(meta):return residue_id(meta)+(meta['atom'],)


def policy(path):
    p=read_json(path);source=read_json(verify(p['source_config']))
    if p['contact_cutoff_A']!=3.5 or p['orientation_seeds']!=['source','radial_away']:
        raise InvalidArtifact('unvalidated contextual water policy; new protocol required')
    return {'source_policy':record(path),'contact_cutoff_A':p['contact_cutoff_A'],
            'starts':p['orientation_seeds'],'water_reference_geometry':source['water_reference_geometry']}


def source_parent(preparation):
    """Preparation requires real chemistry/coordinates, not a precomputed DFT score."""
    p=read_json(verify(preparation['source_preparation']))
    if p.get('protocol_id')!=square.PARENT_PROTOCOL:
        raise InvalidArtifact('water-bearing context requires supported source-backed amide-v3 core; cofactor/PQQ extensions are unsupported')
    outputs=p['outputs']
    paired(verify(outputs['La']['xyz']),verify(outputs['Ca']['xyz']),outputs['La']['charge'],outputs['Ca']['charge'])
    if any(outputs[z]['multiplicity']!=1 for z in ('Ca','La')):
        raise InvalidArtifact('closed-shell paired preparation required')
    atoms=xyz(verify(outputs['La']['xyz']))
    waters=square.water_groups(p,atoms)
    if not waters:raise InvalidArtifact('context requested for a zero-water core')
    return p,atoms,waters


def discover(cases,settings):
    """Same 3.5Å polar shell and homolog union; group membership is explicit input."""
    import affordable_peptide as peptide
    from hydration_network import materialize,radial_seed
    groups=defaultdict(list)
    for case in cases:
        prep=read_json(verify(case['preparation']))
        if not prep['explicit_waters']:continue
        parent,atoms,waters=source_parent(prep)
        graph=peptide.SourceGraph(verify(parent['source_structure']),verify(parent['topology_definition']))
        site_chains={w['source']['chain'] for w in waters}
        sequence={(r.chain,r.resnum,r.insertion_code):r.canonical_resname for r in graph.residues.values()
                  if r.key in graph.supported_topology_residues and r.chain in site_chains}
        selected={source_id(a['source']) for a in parent['atom_graph']['source_to_qm'] if a['kind']=='source'}
        assembly_chains={a['id'].split('/')[0] for a in prep['physical_atoms'] if a['kind'].startswith('protein')}
        contacts=[];fragments=set();outer=set();variable={residue_id(w['source']) for w in waters}
        for water in waters:
            origin=np.array(atoms[water['oxygen_index']][1:])
            for key,atom in graph.atoms.items():
                if atom.element.name not in ('N','O','S'):continue
                meta=graph.meta[key]
                if residue_id(meta)==residue_id(water['source']):continue
                distance=float(np.linalg.norm(np.array(tuple(atom.pos))-origin))
                if distance>settings['contact_cutoff_A']:continue
                if meta['canonical_resname']!='HOH' and meta['chain'] not in assembly_chains:
                    raise InvalidArtifact('polar context crosses the explicitly prepared protein assembly')
                contacts.append({'water':water['source'],'neighbor':meta,'distance_A':distance,
                                 'in_original_core':source_id(meta) in selected})
                if meta['canonical_resname']=='HOH':
                    if residue_id(meta) not in variable:outer.add(residue_id(meta))
                elif key[:2] not in graph.supported_topology_residues:
                    raise InvalidArtifact(f'unsupported nearby polar atom: {meta}')
                elif meta['atom'] in ('C','O','OXT'):fragments.add(('peptide',residue_id(meta)))
                elif meta['atom']=='N':
                    previous=graph.amide_prev.get(key[:2])
                    if previous is None:raise InvalidArtifact(f'unsupported terminal amine contact: {meta}')
                    fragments.add(('peptide',residue_id(graph.meta[previous])))
                else:fragments.add(('sidechain',residue_id(meta)))
        group=case.get('context_group',case['case_id'])
        groups[group].append({'case':case,'item':{'preparation':prep['source_preparation']},'parent':parent,'graph':graph,
            'contacts':contacts,'outer':outer,'variable':variable,'waters':waters,'sequence':sequence,'fragments':fragments})
    contexts={};reference=xyz(verify(settings['water_reference_geometry']))
    for group,members in groups.items():
        seqs=[d['sequence'] for d in members];common=set.intersection(*(set(s) for s in seqs))
        if not common or any(any(s[k]!=seqs[0][k] for k in common) for s in seqs[1:]):
            raise InvalidArtifact('context group requires matching indexed protein sequences')
        union=sorted(set().union(*(d['fragments'] for d in members)))
        for d in members:
            atoms,mapping,ledger,waters=materialize(d,union,reference)
            charge=sum(v['formal_charge'] for v in ledger)
            contexts[d['case']['case_id']]={'case_id':d['case']['case_id'],'parent':d['item']['preparation'],
                'source_structure':d['parent']['source_structure'],'context_group':group,'protein_fragments_union':union,
                'contacts':d['contacts'],'atom_graph':mapping,'charge_ledger':ledger,'water_groups':waters,
                'atoms':atoms,'charges':{'Ca':charge+2,'La':charge+3},
                'sequence_audit':{'common_residue_count':len(common),'common_sequence_identity':True},
                'starts':{'source':atoms,'radial_away':radial_seed(atoms,waters,mapping)}}
    return contexts


def transfer(preparation,context,proposed,metal):
    """Reuse original-core transfer/serialization, then the proven whole-chain map."""
    from hydration_core_transfer import transfer_water_hydrogens
    from hydration_scanner import transfer as transfer_whole
    parent=read_json(verify(preparation['source_preparation']))
    if context['parent']!=preparation['source_preparation']:raise InvalidArtifact('context and whole chain use different source cores')
    original=xyz(verify(parent['outputs'][metal]['xyz']))
    core,mobile=transfer_water_hydrogens(parent,original,context,proposed)
    # This is the established core-transfer XYZ precision, not new geometry tuning.
    serialized=[(a[0],*(float(f'{v:.10f}') for v in a[1:])) for a in core]
    rows,mapping=transfer_whole(preparation,parent,serialized,metal)
    return rows,{'whole_H_mapping':mapping,'source_core_water_H_indices':mobile,
                 'water_H_transfer_serialization_decimal_places':10,'water_inventory_changed':False}


def verify_proposal_geometry(task,result):
    """No failed/incomplete optimizer is a prepared-water success."""
    if (not result.get('optimization_eligible') or result.get('selected_seed') not in ('source','radial_away')
            or not result.get('proposed_xyz') or len(result.get('optimization_runs',[]))!=2
            or not all(r.get('converged') for r in result['optimization_runs'])):
        raise InvalidArtifact('contextual water optimization not converged/eligible')
    old=xyz(verify(task['xyz']));new=xyz(verify(result['proposed_xyz']));mobile=set(task['mobile_indices'])
    if len(old)!=len(new) or any(old[i]!=new[i] for i in range(len(old)) if i not in mobile):
        raise InvalidArtifact('proposal changed fixed coordinates, oxygen or atom inventory')
    if [a[0] for a in old]!=[a[0] for a in new]:raise InvalidArtifact('proposal changed element ordering')
    for w in task['water_groups']:
        if w['role']!='variable':continue
        a=np.array([old[i][1:] for i in w['indices']]);b=np.array([new[i][1:] for i in w['indices']])
        if not np.allclose(np.linalg.norm(a[:,None]-a[None,:],axis=2),np.linalg.norm(b[:,None]-b[None,:],axis=2),atol=1e-10,rtol=0):
            raise InvalidArtifact('proposal changed rigid water shape')
    return new
