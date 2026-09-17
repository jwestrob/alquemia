"""Explicit charge-feature ablation; outputs are descriptors, not QM energies."""
from __future__ import annotations
from types import MethodType

ADAPTER='omol_zero_charge_embedding_before_joint_projection_v1'
COMPONENT='MACE_OMOL_charge_feature_ablated_descriptor'


def install(model):
    import torch
    joint=getattr(model,'joint_embedding',None)
    if (type(model).__name__!='ScaleShiftMACE' or joint is None
            or type(joint).__name__!='GenericJointEmbedding'
            or set(joint.embedders)!={'total_charge','total_spin'}):
        raise ValueError('charge ablation supports only the declared native OMOL architecture')
    charge=joint.embedders['total_charge']
    if not isinstance(charge,torch.nn.Embedding) or charge.weight.shape!=(201,1024):
        raise ValueError('unexpected charge embedding dimensions')
    if hasattr(model,'_alquemia_charge_ablation'):
        raise ValueError('charge feature is already modified')
    state={'id':ADAPTER,'calls':0,'rows':0,'output_max_abs':None,
           'physical_charge_input_preserved':True,'spin_embedding_modified':False,
           'output_semantics':'energy_like_descriptor_not_quantum_endpoint'}

    def zero_feature(this,indices):
        if torch.is_grad_enabled():
            raise RuntimeError('this research adapter is qualified for energy-like outputs only')
        output=this.weight.new_zeros((*indices.shape,this.embedding_dim))
        state['calls']+=1;state['rows']+=indices.numel();state['output_max_abs']=0.
        return output

    charge.forward=MethodType(zero_feature,charge)
    model._alquemia_charge_ablation=state
    return state


def receipt(model):
    state=dict(model._alquemia_charge_ablation)
    if state['calls']!=1 or state['output_max_abs']!=0.:
        raise ValueError('the declared single zero-charge-feature forward did not occur')
    return state
