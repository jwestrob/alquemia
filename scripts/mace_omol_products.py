"""Exact atom batching of the native OMOL product block, following edge batching."""
from __future__ import annotations
import types
from affordable_common import InvalidArtifact
from mace_omol_edges import install as install_edges, receipt as edge_receipt

ADAPTER = 'omol_exact_edge_and_product_batches_v1'
BLOCK = 'EquivariantProductBasisBlock'


def forward(self, node_feats, sc, node_attrs):
    import torch
    if torch.is_grad_enabled():
        raise InvalidArtifact('product batching supports no-grad energy inference only')
    count = node_feats.shape[0]
    if count != node_attrs.shape[0] or (sc is not None and count != sc.shape[0]):
        raise InvalidArtifact('product input atom counts disagree')
    chunk = self._omol_product_chunk_size
    outputs = []
    processed = 0
    for start in range(0, count, chunk):
        stop = min(start + chunk, count)
        output = self._omol_product_original_forward(
            node_feats=node_feats[start:stop], sc=None if sc is None else sc[start:stop],
            node_attrs=node_attrs[start:stop])
        if output.shape[0] != stop - start:
            raise InvalidArtifact('native product changed atom count')
        outputs.append(output)
        processed += stop - start
    if not outputs or processed != count:
        raise InvalidArtifact('product atom inventory was not fully processed')
    self._omol_product_calls.append({'atoms':count, 'processed_atoms':processed,
                                     'chunk_size':chunk, 'batches':len(outputs)})
    return torch.cat(outputs, dim=0)


def install(model, edge_chunk_size, product_chunk_size):
    if product_chunk_size not in (32, 1024) or edge_chunk_size != 1024:
        raise InvalidArtifact('undeclared product/edge batch size')
    if len(model.products) != 3:
        raise InvalidArtifact('unsupported native product count')
    for block in model.products:
        if (type(block).__name__ != BLOCK or getattr(block, 'cueq_config', None)
                or getattr(block, 'oeq_config', None)
                or hasattr(block, '_omol_product_chunk_size')):
            raise InvalidArtifact('unsupported or already adapted native product')
    install_edges(model, edge_chunk_size)
    for block in model.products:
        block._omol_product_chunk_size = product_chunk_size
        block._omol_product_calls = []
        block._omol_product_original_forward = block.forward
        block.forward = types.MethodType(forward, block)


def receipt(model, edge_chunk_size, product_chunk_size):
    result = edge_receipt(model, edge_chunk_size)
    calls = [block._omol_product_calls for block in model.products]
    if any(len(rows) != 1 for rows in calls):
        raise InvalidArtifact('each native product must execute exactly once')
    layers = [rows[0] for rows in calls]
    count = result['layers'][0]['atoms']
    if any(r['atoms'] != count or r['processed_atoms'] != count
           or r['chunk_size'] != product_chunk_size for r in layers):
        raise InvalidArtifact('product adapter changed atom inventory')
    result.update(id=ADAPTER, product_chunk_size=product_chunk_size, product_layers=layers)
    return result
