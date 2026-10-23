#!/usr/bin/env python3
"""Focused saved controls only; synthetic future reviews are never authorities."""
import argparse,copy,importlib.util,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('cluster_handoff',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def main(panel,preparation,out):
    tick,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False);checks={};p=m.load(panel)
    def ok(name,value):
        checks[name]=bool(value)
        if not value:raise AssertionError(name)
    def rejects(name,fn):
        try:fn()
        except (RuntimeError,ValueError,KeyError,AssertionError):ok(name,True)
        else:ok(name,False)
    try:
        ok('original_deck_and_recovered_parser_reused',m.deck is m.old.deck and m.parse is m.rec.parse)
        ok('exact_eight_original_ids',len(p['cells'])==8 and len(set(c['cell_id'] for c in p['cells']))==8)
        ok('tighter_scout_explicit_old_scout_preserved',m.subsets(p)['scout']==['cluster__tighter__aquo9_central_Dy'] and p['original_scout_cell_id']=='cluster__tight__aquo9_central_Dy' and len(m.subsets(p)['remaining'])==7)
        q=m.read(p['qualified_controls']);bad=m.read(p['failed_tight_controls'])
        ok('actual_new35_original32_kept_separate',all(q['checks'].values()) and sum(bad['checks'].values())==32)
        env={'SLURM_JOB_ID':'123','SLURM_NTASKS':'24','SLURM_JOB_NUM_NODES':'1','SLURM_CPUS_PER_TASK':'1','SLURM_MEM_PER_NODE':'96000'}
        ok('per_node_default_memory',m.resources(env)==(24,1,96000,3000))
        cpuenv=dict(env);cpuenv.pop('SLURM_MEM_PER_NODE');cpuenv['SLURM_MEM_PER_CPU']='4000'
        ok('per_cpu_default_memory',m.resources(cpuenv)==(24,1,96000,3000))
        for key,value in [('SLURM_JOB_ID',''),('SLURM_NTASKS','8'),('SLURM_JOB_NUM_NODES','2'),('SLURM_CPUS_PER_TASK','2'),('SLURM_MEM_PER_NODE','0')]:
            rejects('reject_resource_'+key,lambda key=key,value=value:m.resources(dict(env,**{key:value})))
        fake={'path':'SYNTHETIC_ONLY_SOURCE_REVIEW','sha256':'0'*64}
        source_review=dict(status='passed_cluster_handoff_source_review',implementation=p['adapter'],panel=m.b.pin(panel),preparation=m.b.pin(preparation))
        r=dict(status='released_consumed_cluster_subset',implementation=p['adapter'],panel=m.b.pin(panel),parser=p['parser'],declaration=p['declaration'],maximum_new_native_calls=8,ranks=24,nodes=1,execution_root=p['execution_root'],**{k:p[k] for k in m.AUTH},stage='scout',cell_ids=[m.SCOUT],preparation=m.b.pin(preparation),source_review=fake,scout_actual_review=None,scout_root_review=None,scout_collection=None)
        original_read=m.read
        def read(pin):return source_review if pin==fake else original_read(pin)
        m.read=read
        try:
            ok('exact_one_scout_release',m.release_gate(p,panel,r)==r)
            for name,change in [('old_tight_scout',{'cell_ids':[p['original_scout_cell_id']]}),('all8_initial',{'cell_ids':[c['cell_id'] for c in p['cells']]}),('extra_call',{'maximum_new_native_calls':9}),('wrong_ranks',{'ranks':8}),('wrong_actual_authority',{'qualified_controls':p['failed_tight_controls']})]:
                rejects('reject_release_'+name,lambda change=change:m.release_gate(p,panel,dict(r,**change)))
            # Deliberately missing future native/cost authorities; never create fake actual files.
            rr=dict(r,stage='remaining',cell_ids=m.subsets(p)['remaining'])
            try:m.release_gate(p,panel,rr)
            except (TypeError,KeyError,RuntimeError):ok('remaining7_blocked_without_actual_scout',True)
            else:ok('remaining7_blocked_without_actual_scout',False)
        finally:m.read=original_read
        energies={'aquo8_central_La':10.,'aquo8_central_Dy':13.,'aquo9_central_La':20.,'aquo9_central_Dy':27.}
        rows={c['cell_id']:dict(status='complete',energy_hartree=energies[c['geometry_id']]/m.b.HA) for c in p['cells']}
        rows[m.SCOUT]['energy_hartree']+=.125/m.b.HA;s=m.score(p,rows)
        ok('signed_four_term_contrast',abs(s['raw_CC_C_kcal_mol']['tight']-4)<1e-12 and abs(s['raw_CC_C_kcal_mol']['tighter']-4.125)<1e-12)
        ok('all_four_endpoint_tier_shifts',len(s['endpoint_tighter_minus_tight_kcal_mol'])==4 and abs(s['endpoint_tighter_minus_tight_kcal_mol']['aquo9_central_Dy']-.125)<1e-12)
        ok('unchanged_quarter_kcal_gate_pass',s['target_threshold_sensitivity_passed'] and s['cutoff_sensitivity_limit_kcal_mol']==.25)
        rows[m.SCOUT]['energy_hartree']+=.135/m.b.HA;ok('unchanged_quarter_kcal_gate_fail',not m.score(p,rows)['target_threshold_sensitivity_passed'])
        missing=copy.deepcopy(rows);missing.pop(m.SCOUT);rejects('missing_endpoint_rejected',lambda:m.score(p,missing))
        failed=copy.deepcopy(rows);failed[m.SCOUT]['status']='failed';rejects('failed_endpoint_rejected',lambda:m.score(p,failed))
        ok('PBE0_electronic_D4_separate',abs(s['PBE0_D4_C_kcal_mol']-s['PBE0_electronic_C_kcal_mol']-s['original_D4_C_kcal_mol'])<1e-12 and s['original_D4_C_kcal_mol']>0)
        result=dict(status='passed_cluster_handoff_saved_controls',implementation=m.b.pin(m.__file__),test_source=m.b.pin(__file__),panel=m.b.pin(panel),preparation=m.b.pin(preparation),checks=checks,check_count=len(checks),new_native_calls=0,threads=1,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu)
    except Exception as exc:
        result=dict(status='failed_saved_controls',error=repr(exc),checks=checks,new_native_calls=0,threads=1,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu)
        m.b.write(out/'RESULT.json',result);raise
    pin=m.b.write(out/'RESULT.json',result);print(json.dumps(dict(result=pin,status=result['status'],checks=len(checks))))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--preparation',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();main(a.panel,a.preparation,a.out)
