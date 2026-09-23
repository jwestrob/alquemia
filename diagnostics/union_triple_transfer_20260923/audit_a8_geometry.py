"""Read saved A8 geometry: fixed source Ser352 OG to actual frozen donor anchors."""
import argparse
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from consistent_context import source_id

@lru_cache(None)
def load(path):return read_json(path)
def pinned(pin):return load(str(verify(pin)))
def case(summary):return next(c for c in pinned(summary['collection'])['cases'] if c['case_id']==summary['case_id'])
def mapping(c):return pinned(c['source']['union']['representations']['context']['preparation'])['mapping']['source_to_qm']
def indices(entries):return {source_id(a['source']):a['qm_index'] for a in entries if 'source' in a}

def audit(audit,output):
    data=read_json(audit);rows=[]
    for row in data['rows']:
        old=case(row['before']);new=case(row['after']);old_indices=indices(mapping(old));og=('A',352,'','SER','OG','O')
        oldcoords=np.array([a[1:] for a in xyz(verify(old['matrix']['Ca']['origin']['xyz']))]);og0=oldcoords[old_indices[og]]
        local=pinned(row['source_preparation']['old_context']['preparation']);anchors=local['anchors']
        stages={}
        for label,c,summary in [('tenfold',old,row['before']),('threefold',new,row['after'])]:
            ix=indices(mapping(c));m=pinned(summary['source_manifest']);tasks=[t for t in m['tasks'] if t['case_id']==c['case_id']]
            for t in tasks:
                g=pinned(t['mapping'])['context']
                for mi in t['active_indices']:
                    moved=[g['source_atom_metadata'][i] for i in g['modes'][mi]['moving_indices']]
                    if any(a and source_id(a)==og for a in moved):raise InvalidArtifact('Ser OG is not fixed in selected subspace')
            points=[]
            for candidate in c['candidates']:
                coords=np.array([a[1:] for a in xyz(verify(candidate['xyz']))]);ogpos=coords[ix[og]] if og in ix else og0
                if np.max(np.abs(ogpos-og0))>1e-12:raise InvalidArtifact('saved Ser OG moved')
                ds=[(float(np.linalg.norm(ogpos-coords[ix[source_id(a)]])),a) for a in anchors]
                distance,nearest=min(ds,key=lambda x:x[0])
                points.append({'candidate':candidate['id'],'xyz':candidate['xyz'],'nearest_anchor':nearest,'distance_A':distance,
                    'within_original_3_5_A':distance<=3.5,'within_proposed_q0_4_3_A':distance<=4.3,
                    'OG_in_scored_context':og in ix,'OG_fixed_source_xyz_A':og0.tolist(),
                    'all_anchor_distances':[{'anchor':a,'distance_A':d} for d,a in ds]})
            stages[label]=points
        rows.append({'pair_id':row['pair_id'],'source_case_id':row['case_id'],'atom_count':row['after_atoms'],
            'anchor_preparation':row['source_preparation']['old_context']['preparation'],'frozen_anchor_count':len(anchors),
            'stages':stages,'source_OG_is_fixed_under_both_actual_selected_subspaces':True})
    result={'energy_audit':record(audit),'implementation':record(__file__),'rows':rows,
        'distance_policy':'distance to original source-specific donor anchors (3.3 A direct plus actual functional-group partners), source OG fixed; no re-selection at moved geometry',
        'proposal_only':{'original_cutoff_A':3.5,'maximum_anchor_motion_A':0.8,'candidate_envelope_A':4.3,
          'fixed_neighbor_bound':'triangle inequality for fixed omitted neighbor and anchors moving at most0.8 A',
          'both_atoms_mobile_bound_A':5.1},'new_molecular_calls':0}
    write_new(output,result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit',required=True);p.add_argument('--output',required=True);audit(**vars(p.parse_args()))
