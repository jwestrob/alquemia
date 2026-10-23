#!/usr/bin/env python3
"""Bounded saved input/parser/arithmetic/failure controls; no molecular calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import time
import resource

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('pair_model_controls',HERE/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main(config,preparation,out):
    start=time.monotonic();cpu=time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False);checks={}
    def check(name,value):
        if not value:raise AssertionError(name)
        checks[name]=True
    def rejects(name,fn):
        try:fn()
        except Exception:checks[name]=True
        else:raise AssertionError(name)
    cfg=json.loads(Path(config).read_text());prep=json.loads(Path(preparation).read_text());manifest=m.read(prep['manifest']);data=m.read(cfg['data']);panel=m.read(cfg['panel'])
    m.bound_inputs(cfg,manifest);check('all34_exact_shared_literal_inputs',len(manifest['rows'])==34)
    check('four_selected_clusters_and_two_explicit_zero_ions',len(manifest['saved_cluster_models'])==4 and all(v['status']=='analytically_zero' and v['missing_data_substitution'] is False for v in manifest['analytic_bare_ions'].values()))
    for gid,v in manifest['saved_cluster_models'].items():
        raw=m.parse_analyze(m.verify(v['native']['analyze_output']));old=data['baseline_model'][gid]
        check('saved_native_energy_components_'+gid,raw['total_energy_kcal_mol']==old['energy'] and raw['components']==old['native']['components'])
    sample=panel['cells'][0];gid=sample['parent_geometry_id'];source=manifest['rows'][sample['cell_id']]['selected_source_input'];parent=data['rows'][gid]['geometry'];text=m.verify(source['xyz']).read_text()
    for name,edit in [('map',lambda c:c['parent_atom_indices'].__setitem__(1,4)),('water_index',lambda c:c.__setitem__('water_index',2)),('coordinate',lambda c:c['geometry']['positions_A'][1].__setitem__(0,99.))]:
        v=copy.deepcopy(sample);edit(v);rejects('reject_changed_'+name,lambda:m.extract(v,parent,text))
    first=next(iter(manifest['saved_cluster_models'].values()))['native']['analyze_output'];actual=m.verify(first).read_text()
    for name,bad in [('nonfinite',actual+'\nNaN\n'),('nonconvergence',actual+'\nNOT CONVERGED\n'),('version',actual.replace('Version 26.2','Version wrong')),('missing_energy',actual.replace('Total Potential Energy','Absent potential'))]:
        p=out/(name+'.out');p.write_text(bad);rejects('reject_native_'+name,lambda:m.parse_analyze(p))
    c=m.Counter()
    for i in range(34):check('reserve_'+str(i+1),c.reserve()==i+1)
    rejects('hard35th_program_rejected',c.reserve);c=m.Counter();c.stop.set();rejects('global_stop_prevents_start',c.reserve)
    ions=manifest['analytic_bare_ions'];clusters={g:{'energy_kcal_mol':float(i+1)} for i,g in enumerate(m.GIDS)}
    rows={cid:{'status':'complete','parent_geometry_id':v['parent_geometry_id'],'energy_kcal_mol':float(1 if v['parent_geometry_id'].endswith('La') else 3)} for cid,v in manifest['rows'].items()}
    cycle=m.model_cycle(rows,clusters,ions);check('cycle_signs_and_units',cycle['C']==0. and cycle['P']==2. and cycle['B']==-2. and cycle['units']=='kcal/mol')
    broken=copy.deepcopy(rows);broken[m.IDS[0]]['status']='not_attempted';check('missing_slot_makes_cycle_unavailable',m.model_cycle(broken,clusters,ions) is None)
    broken=copy.deepcopy(rows);del broken[m.IDS[0]];rejects('missing_denominator_rejected',lambda:m.model_cycle(broken,clusters,ions))
    broken=copy.deepcopy(rows);broken[m.IDS[0]]['energy_kcal_mol']=float('inf');rejects('nonfinite_cycle_rejected',lambda:m.model_cycle(broken,clusters,ions))
    badion=copy.deepcopy(ions);badion['Dy']['status']='missing';rejects('missing_ion_not_zero',lambda:m.model_cycle(rows,clusters,badion))
    resolver=m.module(cfg['resolver'],'pair_model_control_resolver')
    rejects('external_field_invalidates_zero_ion',lambda:m.bare_ion_proof(m.verify(source['key']).read_text()+'\nexternal-field 1 0 0\n',m.verify(source['parameters']).read_text(),'La',resolver,cfg['physics_sources']))
    rejects('no_release_cannot_run',lambda:m.release_gate(cfg,prep,manifest))
    # Exercise the actual scheduler-free driver failure path with an explicitly
    # substituted native invocation; no executable or chemistry is called.
    original_module=m.module;original_gate=m.release_gate;starts=[]
    class FailureNative:
        @staticmethod
        def invoke(argv,cwd,output,receipts):
            starts.append(str(cwd));receipts.append({'fixture_only':True,'real_native_calls':0,'returncode':1})
            raise RuntimeError('FIXTURE_ONLY_NATIVE_FAILURE_NO_EXECUTION')
    m.module=lambda p,n:FailureNative if n=='literal_pair_native' else original_module(p,n)
    m.release_gate=lambda *args,**kwargs:None
    import os
    previous=os.environ.get('SLURM_JOB_ID');os.environ['SLURM_JOB_ID']='0'
    try:
        fake=copy.deepcopy(cfg);fake['preparation_result']=m.pin(preparation);directory=out/'explicit_mock_failure';directory.mkdir();failure=m.run(fake,directory)
    finally:
        m.module=original_module;m.release_gate=original_gate
        if previous is None:os.environ.pop('SLURM_JOB_ID',None)
        else:os.environ['SLURM_JOB_ID']=previous
    check('mock_failure_retains_all34_slots',set(failure['rows'])==set(m.IDS))
    check('mock_failure_stops_within_active_four',1<=len(starts)<=4 and failure['native_program_invocations']==len(starts) and failure['global_stop'] and failure['status']=='failed_or_unavailable')
    check('mock_failure_no_partial_cycle',failure['model_cycle'] is None and failure['error'] is not None)
    result={'status':'passed','implementation':m.pin(HERE/'RUN.py'),'control_source':m.pin(__file__),'configuration':m.pin(config),'preparation':m.pin(preparation),
        'checks':checks,'check_count':len(checks),'new_native_calls':0,'mock_invocations_not_native':len(starts),'threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    m.save(out/'RESULT.json',result);print(json.dumps({'status':'passed','checks':len(checks),'result':m.pin(out/'RESULT.json')}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--preparation',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(a.config,a.preparation,a.output)
