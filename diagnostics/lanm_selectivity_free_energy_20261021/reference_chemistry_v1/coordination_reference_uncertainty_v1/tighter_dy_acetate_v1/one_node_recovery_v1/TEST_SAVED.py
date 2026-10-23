#!/usr/bin/env python3
"""Focused one-node/resource and preserved-attempt release controls; zero calls."""
import copy, importlib.util, json, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('recovery_control',HERE/'RECOVER.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)

def main():
    tick,cpu=time.monotonic(),time.process_time()
    w=Path('workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/tighter_dy_acetate_v1/one_node_recovery_v1').resolve()
    mp=w/'preparation_v1/MANIFEST.json';m,p,c=r.load(mp);checks={}
    def check(name,fn,accept):
        try:fn();passed=accept
        except Exception:passed=not accept
        checks[name]=passed
        assert passed,name
    env=dict(SLURM_JOB_ID='99999',SLURM_NTASKS='8',SLURM_CPUS_PER_TASK='1',SLURM_MEM_PER_CPU='4000',SLURM_JOB_NUM_NODES='1')
    check('one_node_default_grant',lambda:r.resources(env),True)
    for name,change in [('two_nodes',dict(SLURM_JOB_NUM_NODES='2')),('four_ranks',dict(SLURM_NTASKS='4')),('changed_grant',dict(SLURM_MEM_PER_CPU='2000'))]:
        check(name,lambda change=change:r.resources(dict(env,**change)),False)
    rel=dict(status='approved_single_node_tighter_attribution',implementation=m['implementation'],manifest=r.b.pin(mp),root_decision=m['root_decision'],cell_id=r.CID,
        maximum_new_native_calls=1,nodes=1,ranks=8,execution_out=m['execution_out'],preparation='prep',source_review='review',original_terminal_receipt='terminal',original_collection='collection',
        full_eighteen_qualified=False,cluster_release_authorized=False)
    values=dict(prep=dict(status='passed_single_node_recovery_preparation',implementation=m['implementation'],manifest=r.b.pin(mp),new_native_calls=0),
        review=dict(status='passed_single_node_recovery_source_review',implementation=m['implementation'],manifest=r.b.pin(mp),preparation='prep'),
        terminal=dict(records={'18518':dict(state='CANCELLED',allocated_cpus='8')}),
        collection=dict(status='unavailable',declared=3,complete=2,panel=m['original_panel'],implementation=m['original_adapter'],
            rows={**{cid:dict(status='complete',result=pin) for cid,pin in m['original_results'].items()},r.CID:dict(status='unavailable')}))
    original_read=r.old.read
    r.old.read=lambda key:values[key]
    try:
        check('canceled_old_incomplete_collection_passes',lambda:r.release_contract(m,mp,rel),True)
        for name,container,key,value in [('running_old',values['terminal']['records']['18518'],'state','RUNNING'),
            ('completed_old',values['terminal']['records']['18518'],'state','COMPLETED'),
            ('old_complete_cycle',values['collection'],'complete',3),('missing_old_role',values['collection']['rows'][r.CID],'status','complete'),
            ('wrong_source_review',values['review'],'implementation',None),('extra_calls',rel,'maximum_new_native_calls',2),('cluster_release',rel,'cluster_release_authorized',True)]:
            previous=container[key];container[key]=value
            check(name,lambda:r.release_contract(m,mp,rel),False);container[key]=previous
        cid=next(iter(m['original_results']));previous=values['collection']['rows'][cid]['result'];values['collection']['rows'][cid]['result']=None
        check('changed_reused_result',lambda:r.release_contract(m,mp,rel),False);values['collection']['rows'][cid]['result']=previous
    finally:r.old.read=original_read
    check('literal_input_bytes',lambda:r.b.check(r.b.verify(m['input']).read_bytes()==r.b.verify(m['original_input']).read_bytes(),'changed'),True)
    out=w/'source_controls_v1';out.mkdir(exist_ok=False)
    q=dict(status='passed',implementation=r.b.pin(HERE/'RECOVER.py'),test=r.b.pin(__file__),manifest=r.b.pin(mp),
        checks=checks,passed=len(checks),new_native_calls=0,threads=1,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu)
    print(json.dumps(r.b.write(out/'RESULT.json',q)))

if __name__=='__main__':main()
