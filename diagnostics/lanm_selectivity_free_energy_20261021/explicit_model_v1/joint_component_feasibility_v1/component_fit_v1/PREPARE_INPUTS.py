#!/usr/bin/env python3
"""Freeze only existing admitted vector/cell inventory; no new point evaluation."""
import copy
import importlib.util
import json
from pathlib import Path
import socket
import sys
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('component_inputs',HERE/'ADAPTER.py');A=importlib.util.module_from_spec(s);s.loader.exec_module(A)
ROOT=HERE.parents[4]
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
J='explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1'
F='explicit_model_v1/joint_component_feasibility_v1'

def main():
    out=W/F/'component_fit_v1/inputs_v1';out.mkdir(parents=True,exist_ok=False)
    declaration=A.pin(HERE/'ROOT_DECLARATION.json');decl=A.read(declaration)
    rootpin=A.pin(W/J/'fit_v1/run_v2/ROOT_FINAL_REVIEW.json');root=A.read(rootpin);fit=A.read(root['fit_result']);it=fit['iterations'][23]
    parent=A.pin(W/J/'fit_v1/CONFIG_EVALUATOR.json');parentcfg=A.read(parent)
    inv={'status':'immutable_admitted_component_reuse_inventory','declaration':declaration,'execution_context':A.pin(HERE/'ROOT_EXECUTION_CONTEXT_v2.json'),
         'parent_configuration':parent,'root_fit_review':rootpin,'fit_result':root['fit_result'],'origin_key':it['origin_key'],'derivative_keys':it['derivative_keys'],
         'vectors':{k:fit['vectors'][k]['artifact'] for k in [it['origin_key']]+it['derivative_keys']},'pilot_root_review':decl['prior_native_pilot_admission'],'cells':{},
         'inherited_fit_cost':fit['cost'],'no_old_ledger_import':True,'new_native_calls':0}
    vectors={k:A.read(p) for k,p in inv['vectors'].items()}
    pilotroot=A.read(inv['pilot_root_review']);pilot=A.read(pilotroot['result']);pv=copy.deepcopy(pilot['answer']);pv.update(status='complete',parameters=pilot['candidate']['parameters'])
    inv.update(pilot_result=pilotroot['result'],pilot_key=pilot['candidate']['key'],pilot_answer_sha256=A.digest(pilot['answer']),inherited_pilot_cost=pilot['cost'])
    vectors[inv['pilot_key']]=pv
    for key,v in vectors.items():
        A.require(A.digest(v['parameters'])==key,'Exact saved vector token')
        for cid,r in v['cells'].items():
            token=r['cache_key'];entry={'vector_key':key,'cell_id':cid,'record_sha256':A.cell_digest(r)}
            if token in inv['cells']:A.require(inv['cells'][token]['record_sha256']==entry['record_sha256'],'Duplicate immutable cell equality')
            else:inv['cells'][token]=entry
    inv['counts']={'origin_vectors':1,'derivative_vectors':24,'additional_pilot_vectors':1,'unique_cells':len(inv['cells']),'vector_cell_references':26*46}
    A.save(out/'REUSE_INVENTORY.json',inv)
    cfg={'implementation':A.pin(HERE/'DRIVER.py'),'adapter':A.pin(HERE/'ADAPTER.py'),'core':A.pin(HERE/'CORE.py'),'declaration':declaration,'execution_context':inv['execution_context'],
         'parent_evaluator':parentcfg['implementation'],'parent_configuration':parent,'reuse_inventory':A.pin(out/'REUSE_INVENTORY.json'),
         'decomposition':A.pin(W/J/'interaction_decomposition_v1/actual_v1/RESULT.json'),'interpreter':A.pin(sys.executable),'hostname':socket.gethostname(),
         'limits':A.LIMITS,'reserved_outcomes_authorized':False,'output_directory':str(W/F/'component_fit_v1/run_v1/execution'),'root_release':None,'source_review':None,'preparation':None}
    A.save(out/'CONFIG_PREPARE.json',cfg)
    print(json.dumps({'inventory':A.pin(out/'REUSE_INVENTORY.json'),'configuration':A.pin(out/'CONFIG_PREPARE.json'),'counts':inv['counts'],'new_native_calls':0}))
if __name__=='__main__':main()
