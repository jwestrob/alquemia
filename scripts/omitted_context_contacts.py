"""Read-only actual donor-proposal contacts with omitted source protein atoms."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from functools import lru_cache
import json
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from affordable_peptide import SourceGraph
from second_shell_context import parent_state,atom_key
from environment_context_chemistry import cysteine_units
from mace_site_kinematics import Kinematics

POLICY={'severe_distance_A':2.0,'descriptive_distance_A':3.5,'bond_exclusion_max_edges':3,
        'motion_identity_tolerance_A':1e-8,'coordinate_replay_tolerance_A':1e-12,
        'outside':'raw_source_protein_heavy_not_in_scored_context','outside_coordinates':'frozen_raw_source',
        'cofactor':'fully_mapped_fixed_PQQ_excluded_from_protein_contact_sets',
        'hydrogens_and_caps':'excluded_no_reconstruction','classification_changed':False}

@lru_cache(None)
def load(path):return read_json(path)

def pinned(pin):return load(str(verify(pin)))

def coords(atoms):return np.asarray([a[1:] for a in atoms],dtype=float)

def label(meta):return f"{meta['chain']}:{meta['canonical_resname']}{meta['resnum']}{meta.get('insertion_code','')}:{meta['atom']}"

@lru_cache(None)
def source_state(parent_path,raw_path,topology_path,core_serialized):
    core=json.loads(core_serialized);state=parent_state(core,record(topology_path),require_endpoint_receipts=False)
    g=state['graph'];raw=SourceGraph(raw_path,topology_path)
    if g.structure.connections or raw.structure.connections:raise InvalidArtifact('explicit source covalent connections require separate supported graph mapping')
    fragments=[('sidechain',key) for key,r in g.residues.items() if r.canonical_resname=='CYS']
    _,disulfides=cysteine_units(state,fragments)
    protein={k for k in g.atoms if k[:2] in g.supported_topology_residues and g.meta[k]['element'] not in ('H','D')}
    raw_heavy={k for k in raw.atoms if raw.meta[k]['element'] not in ('H','D')}
    additions=protein-raw_heavy;observed=protein&raw_heavy;maximum=0.
    for key in observed:
        a,b=g.meta[key],raw.meta[key]
        if (a['chain'],a['resnum'],a['insertion_code'],a['atom'],a['canonical_resname'],a['element'])!=(b['chain'],b['resnum'],b['insertion_code'],b['atom'],b['canonical_resname'],b['element']):
            raise InvalidArtifact('raw/protonated protein source identity differs')
        delta=g.atoms[key].pos.dist(raw.atoms[key].pos);maximum=max(maximum,delta)
        if delta>1e-8:raise InvalidArtifact('source preparation moved an observed protein heavy atom')
    unexpected={g.meta[k]['canonical_resname'] for k in g.atoms if k[:2] not in g.supported_topology_residues and g.meta[k]['element'] not in ('H','D')}-{ 'PQQ','LA','CA' }
    if unexpected:raise InvalidArtifact('unsupported nonprotein source species: '+','.join(sorted(unexpected)))
    adjacency={k:set() for k in g.atoms}
    for a,b in g.edges:adjacency[a].add(b);adjacency[b].add(a)
    return state,raw,observed,adjacency,{'prepared_source':state['source'],'raw_source':record(raw_path),
        'source_protein_heavy_atoms':len(observed),'prepared_only_heavy_atoms':[g.meta[k] for k in sorted(additions)],
        'raw_heavy_max_difference_A':maximum,'disulfides':disulfides,'explicit_source_connections':0}


def within_three(key,adjacency,g):
    seen={key};front={key}
    for depth in range(3):
        nxt=set()
        for current in front:
            absent={k for k in g.required[current] if k not in g.atoms}
            if absent:raise InvalidArtifact('missing relevant source-graph neighbor: '+str(sorted(absent)))
            nxt.update(adjacency[current])
        front=nxt-seen;seen.update(front)
    return seen


def logical_case(collection,cid,exact=False):
    cases=pinned(collection)['cases']
    if exact:matches=[c for c in cases if c['case_id']==cid]
    else:matches=[c for c in cases if c['source'].get('actual_union_case_id',c['source'].get('source_case_id',c['case_id']))==cid]
    if len(matches)!=1:raise InvalidArtifact('one actual pool counterpart required')
    return matches[0]


def metrics(distances,origin,keep,query_keys,outside_keys,g):
    d=np.where(keep,distances,np.inf);old=np.where(keep,origin,np.inf)
    if not np.isfinite(d).any():return {'minimum_A':None,'below_3_5_A':0,'below_2_0_A':0,'new_below_2_0_A':0,'contacts_below_3_5_A':[],'new_severe_contacts':[]}
    close=np.argwhere(d<POLICY['descriptive_distance_A']);new=(d<2.)&(old>=2.)
    def contact(i,j):return {'query':label(g.meta[query_keys[i]]),'outside':label(g.meta[outside_keys[j]]),
                           'distance_A':float(d[i,j]),'origin_distance_A':float(old[i,j])}
    nearest=np.unravel_index(np.argmin(d),d.shape)
    return {'minimum_A':float(d[nearest]),'nearest':contact(*nearest),'below_3_5_A':len(close),
            'below_2_0_A':int(np.sum(d<2.)),'new_below_2_0_A':int(np.sum(new)),
            'origin_below_2_0_A':int(np.sum(old<2.)),
            'contacts_below_3_5_A':[contact(i,j) for i,j in close],
            'new_severe_contacts':[contact(i,j) for i,j in np.argwhere(new)]}


def context_contacts(source,collection,cid,topology,exact=False):
    case=logical_case(collection,cid,exact);pm=pinned(pinned(collection)['manifest']);sp=verify(pm['source_manifest']);sm=load(str(sp))
    if case['pool']['status']!='available':raise InvalidArtifact('actual common pool unavailable')
    tasks={z:next(t for t in sm['tasks'] if (t['case_id'],t['metal'])==(case['case_id'],z)) for z in ('Ca','La')}
    maps={z:pinned(t['mapping'])['context'] for z,t in tasks.items()}
    if maps['Ca']!=maps['La']:raise InvalidArtifact('actual paired physical maps differ')
    kin=Kinematics(maps['Ca']);q0=np.zeros(len(kin.modes));core=source['original_core']
    state,raw,protein,adjacency,graph_audit=source_state(core['parent']['path'],str(verify(source['source']['source_structure'])),str(verify(topology)),json.dumps(core,sort_keys=True))
    g=state['graph'];audit=pinned(source['representations']['context']['preparation'])
    if audit['source']!=state['source']:raise InvalidArtifact('context and full source identity differ')
    # Check every real retained heavy atom against the exact source graph; caps never enter.
    mappings={atom_key(a['source']):a['qm_index'] for a in audit['mapping']['source_to_qm'] if 'source' in a and a['source']['element'] not in ('H','D')}
    pqq={k for k in g.atoms if g.meta[k]['canonical_resname']=='PQQ' and g.meta[k]['element'] not in ('H','D')}
    if not pqq or not pqq<=mappings.keys():raise InvalidArtifact('complete mapped PQQ heavy inventory required')
    query=sorted(set(mappings)&protein);outside=sorted(protein-set(mappings))
    if not query or not outside:raise InvalidArtifact('empty observed source protein query/outside set')
    omitted_observed=set(k for k in mappings if k[:2] in g.supported_topology_residues)-protein
    if omitted_observed:raise InvalidArtifact('query contains preparation-added protein heavy atoms')
    kin_sources={tuple(v):i for i,v in enumerate(maps['Ca']['physical_ids']) if isinstance(v,list)}
    for key,index in mappings.items():
        if key not in kin_sources:raise InvalidArtifact('missing physical source identity')
        if [index,'source',kin_sources[key]] not in maps['Ca']['core_links']:raise InvalidArtifact('source-to-context physical link differs')
    points={};receipts={};full_q={};origin_atoms=xyz(verify(case['matrix']['Ca']['origin']['xyz']))
    for name in ('origin','adaptive_Ca','adaptive_La'):
        if name=='origin':q=q0;pin=case['matrix']['Ca']['origin']['xyz'];receipt=None
        else:
            z=name.split('_')[1];t=tasks[z];rp=sp.parent/'proposals'/t['task_id']/'result.json';result=load(str(rp))
            if result['status']!='proposal_available' or result['manifest']!=record(sp):raise InvalidArtifact('actual admitted proposal receipt differs')
            q=np.asarray(result['proposal']['full_q']);pin=result['proposal']['coordinate'];receipt=record(rp)
        atoms=xyz(verify(pin));physical,actual,*_=kin.evaluate(q)
        if not np.allclose(actual,coords(atoms),atol=1e-12,rtol=0):raise InvalidArtifact('saved proposal does not replay physical mapping')
        for z in ('Ca','La'):
            paired_atoms=xyz(verify(case['matrix'][z][name]['xyz']))
            if not np.allclose(coords(paired_atoms),coords(atoms),atol=1e-12,rtol=0):raise InvalidArtifact('cross-metal candidate coordinates differ')
        for key in pqq:
            if not np.array_equal(actual[mappings[key]],kin.core[mappings[key]]):raise InvalidArtifact('PQQ heavy atom moved')
        points[name]=coords(atoms)[[mappings[k] for k in query]];receipts[name]={'coordinate':pin,'proposal_receipt':receipt};full_q[name]=q.tolist()
    rawq=np.asarray([[raw.atoms[k].pos.x,raw.atoms[k].pos.y,raw.atoms[k].pos.z] for k in query])
    if not np.allclose(points['origin'],rawq,rtol=0,atol=1e-8):raise InvalidArtifact('query q0 differs from raw source-heavy coordinates')
    outside_xyz=np.asarray([[raw.atoms[k].pos.x,raw.atoms[k].pos.y,raw.atoms[k].pos.z] for k in outside])
    forbidden=[within_three(k,adjacency,g) for k in query]
    keep=np.asarray([[k not in forbidden[i] for k in outside] for i in range(len(query))])
    mobile=np.zeros(len(query),dtype=bool)
    for name in ('adaptive_Ca','adaptive_La'):mobile|=np.linalg.norm(points[name]-points['origin'],axis=1)>1e-8
    distances={name:np.linalg.norm(x[:,None,:]-outside_xyz[None,:,:],axis=2) for name,x in points.items()}
    summary={name:{'all_context_protein':metrics(d,distances['origin'],keep,query,outside,g),
                   'moving_donor_subset':metrics(d,distances['origin'],keep&mobile[:,None],query,outside,g)} for name,d in distances.items()}
    return {'status':'available','collection':collection,'source_case_id_in_collection':case['case_id'],'source_manifest':record(sp),
            'mapping':tasks['Ca']['mapping'],'source_preparation':source['representations']['context']['preparation'],
            'graph':graph_audit,'query_heavy_count':len(query),'omitted_protein_heavy_count':len(outside),
            'moving_donor_heavy_count':int(mobile.sum()),'moving_donor_atoms':[label(g.meta[k]) for k,m in zip(query,mobile) if m],
            'bond_excluded_pairs':int(np.sum(~keep)),'PQQ_heavy_atoms_fully_mapped_and_fixed':len(pqq),
            'cofactor_graph_scope':'internal PQQ bonds unused; no PQQ atom is a protein query/outside atom',
            'coordinates':receipts,'full_q':full_q,'selected_mode_ids':{z:t['active_mode_ids'] for z,t in tasks.items()},
            'operational_selected_candidates':{z:case['pool']['rows'][z]['operational_candidate'] for z in ('Ca','La')},
            'contacts':summary}


def run(inputs,graph_policy,output):
    start=time.monotonic();m=read_json(inputs);a=pinned(m['audit']);expected={r['pair_id']:r for r in a['rows'] if not r['pool_reuse']}
    if len(m['pairs'])!=55 or set(expected)!={r['pair_id'] for r in m['pairs']}:raise InvalidArtifact('fixed55 population differs')
    verify(m['agreement']);verify(record(graph_policy));rows=[];old_cache={}
    for r in m['pairs']:
        src=expected[r['pair_id']]['source'];row={**r,'representations':{}}
        for mode in ('threefold','tenfold'):
            try:
                if mode=='threefold':value=context_contacts(src,r['threefold_collection'],r['pair_id'],a['config']['topology'],exact=True)
                else:
                    if r['case_id'] not in old_cache:
                        case=logical_case(r['tenfold_collection'],r['case_id']);oldsrc=case['source']['union']
                        if oldsrc['original_core']!=src['original_core']:raise InvalidArtifact('tenfold original source/state differs')
                        old_cache[r['case_id']]=context_contacts(oldsrc,r['tenfold_collection'],r['case_id'],a['config']['topology'])
                    value=old_cache[r['case_id']]
                row['representations'][mode]=value
            except Exception as exc:row['representations'][mode]={'status':'unsupported','reason':str(exc)}
        rows.append(row)
        print(json.dumps({'pair_id':r['pair_id'],'statuses':{k:v['status'] for k,v in row['representations'].items()}}),flush=True)
    summary={}
    for mode in ('threefold','tenfold'):
        values=[r['representations'][mode] for r in rows if r['representations'][mode]['status']=='available']
        summary[mode]={'pair_denominator':55,'available':len(values),
            'pairs_with_any_new_severe_contact':sum(any(v['contacts'][q]['all_context_protein']['new_below_2_0_A'] for q in ('adaptive_Ca','adaptive_La')) for v in values),
            'proposal_denominator':110,'proposals_with_new_severe_contact':sum(v['contacts'][q]['all_context_protein']['new_below_2_0_A']>0 for v in values for q in ('adaptive_Ca','adaptive_La')),
            'pairs_with_preexisting_severe_contact':sum(v['contacts']['origin']['all_context_protein']['below_2_0_A']>0 for v in values)}
    result={'protocol_id':'Nikasha_omitted_protein_contact_audit_v1','inputs':record(inputs),'graph_policy':record(graph_policy),'policy':POLICY,
            'pair_denominator':55,'distinct_sources':len({r['case_id'] for r in rows}),'protein_groups':len({r['protein_id'] for r in rows}),
            'rows':rows,'summary':summary,'wall_seconds':time.monotonic()-start,'new_molecular_calls':0,'new_H':0,
            'scores_or_classifications_changed':False,'implementation':record(__file__)}
    write_new(output,result);return {'output':record(output),'summary':summary,'wall_seconds':result['wall_seconds']}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('inputs','graph-policy','output'):p.add_argument('--'+k,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2))
