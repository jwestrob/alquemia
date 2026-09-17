"""One shared learned charge-zero feature; physical charge inputs stay explicit."""
from __future__ import annotations
import hashlib
from types import MethodType
from affordable_common import InvalidArtifact

ADAPTER='omol_shared_neutral_charge_feature_before_joint_projection_v1'
PROTOCOL='omol_shared_neutral_charge_feature_descriptor_v1'
COMPONENT='MACE_OMOL_shared_neutral_charge_feature_descriptor'
RAW_SHA='27552d5a0a4049402bf2daf39e9fbfa587fb77c410656bcb5176649eba201c46'
INDEX=100


def install(model):
    import torch
    joint=getattr(model,'joint_embedding',None)
    if (type(model).__name__!='ScaleShiftMACE' or type(joint).__name__!='GenericJointEmbedding'
            or set(joint.embedders)!={'total_charge','total_spin'}):
        raise InvalidArtifact('unsupported neutral-feature architecture')
    spec=dict(model.embedding_specs)['total_charge']
    if spec!={'type':'categorical','per':'graph','in_dim':1,'emb_dim':1024,'num_classes':201,'offset':100}:
        raise InvalidArtifact('native charge category mapping differs')
    charge=joint.embedders['total_charge']
    if (not isinstance(charge,torch.nn.Embedding) or tuple(charge.weight.shape)!=(201,1024)
            or hasattr(model,'_alquemia_charge_ablation') or hasattr(model,'_alquemia_neutral_feature')):
        raise InvalidArtifact('unexpected or already modified charge feature')
    raw=charge.weight[INDEX].detach().cpu().numpy()
    if hashlib.sha256(raw.tobytes()).hexdigest()!=RAW_SHA:
        raise InvalidArtifact('learned neutral vector differs from pinned checkpoint inspection')
    original=charge.forward
    state={'id':ADAPTER,'calls':0,'rows':0,'physical_charge_input_preserved':True,
           'spin_embedding_modified':False,'conditioned_charge_category':0,'conditioned_table_index':INDEX,
           'raw_feature_sha256':RAW_SHA,'all_rows_match_learned_feature':None,
           'original_index_min':None,'original_index_max':None,
           'output_semantics':'energy_like_descriptor_not_quantum_endpoint'}

    def forward(this,indices):
        # This fixed categorical input has no nuclear-coordinate dependence.
        with torch.no_grad():
            out=original(torch.full_like(indices,INDEX))
            same=torch.equal(out,this.weight[INDEX].expand_as(out))
        state.update(calls=state['calls']+1,rows=state['rows']+indices.numel(),
                     original_index_min=int(indices.min().item()),original_index_max=int(indices.max().item()),
                     all_rows_match_learned_feature=same)
        return out
    charge.forward=MethodType(forward,charge);model._alquemia_neutral_feature=state
    return state


def receipt(model):
    r=dict(model._alquemia_neutral_feature)
    if r['calls']!=1 or r['all_rows_match_learned_feature'] is not True:
        raise InvalidArtifact('shared neutral feature did not execute exactly once')
    return r


def accepted_receipt(r,task,atoms):
    return r=={'id':ADAPTER,'calls':1,'rows':atoms,'physical_charge_input_preserved':True,
        'spin_embedding_modified':False,'conditioned_charge_category':0,'conditioned_table_index':INDEX,
        'raw_feature_sha256':RAW_SHA,'all_rows_match_learned_feature':True,
        'original_index_min':task['charge']+INDEX,'original_index_max':task['charge']+INDEX,
        'output_semantics':'energy_like_descriptor_not_quantum_endpoint'}
