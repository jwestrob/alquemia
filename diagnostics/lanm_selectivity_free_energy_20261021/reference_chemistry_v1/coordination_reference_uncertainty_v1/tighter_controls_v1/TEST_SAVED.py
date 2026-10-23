#!/usr/bin/env python3
"""Bounded deck/resource/unchanged35-score/collector controls; zero native calls."""
import copy,importlib.util,json,sys,time
from pathlib import Path
H=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('tighter_controls_under_test',H/'REFERENCE.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
pp=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();p=a.load(pp);checks={}
def check(k,v):checks[k]=bool(v);assert v,k
def rejects(k,f):
    try:f()
    except Exception:check(k,True)
    else:check(k,False)
check('exact15new3reuse18controls',len(p['new_cell_ids'])==15 and len(p['reused'])==3 and len(p['cells'])==18 and not set(p['new_cell_ids'])&set(p['reused']))
check('new15onlysame_tighter_recipe',all(c['tier']=='tighter' and c['stage']=='controls' and 'TriTolE 1e-8' in a.b.verify(c['template_input']).read_text() and 'NoFrozenCore' in a.b.verify(c['template_input']).read_text() for c in p['cells']))
for metal in ('La','Dy'):
    cells=[c for c in p['cells'] if c['metal']==metal];check(metal+'_ghost_ECP_inventory',all(('DelECP '+metal in a.deck(p,c,8,3000))==c['ghost_mask'][0] and ('NewECP' in a.deck(p,c,8,3000))==bool(c['expected_ecp_core']) for c in cells))
env={'SLURM_JOB_ID':'999002','SLURM_NTASKS':'8','SLURM_CPUS_PER_TASK':'1','SLURM_JOB_NUM_NODES':'1','SLURM_MEM_PER_NODE':'32000'};check('8ranks1node_defaultgrant_MaxCore',a.resources(env)==(8,32000,3000,1));e=dict(env,SLURM_JOB_NUM_NODES='2');rejects('reject_splitnode',lambda:a.resources(e));e=dict(env,SLURM_NTASKS='4');rejects('reject_wrongranks',lambda:a.resources(e))
synthetic={}
for c in p['cells']:
    old=p['controls'][c['canonical_id']];synthetic[c['cell_id']]={'status':'complete','energy_hartree':old['energy_hartree'],'RHF_energies_hartree':[old['RHF_energy_hartree']],'aux_C':{'synthetic_same_parent':c['donor']+'_'+c['metal']},'allocated_cpu_seconds':0.0}
q=a.score(p,synthetic);check('synthetic_canonical_identity_all35',len(q['checks'])==35 and q['all_controls_passed'] and all(v==0 for v in q['CP_errors_kcal_mol'].values()))
check('exact_gate_categories',sum(k.startswith('RHF__') for k in q['checks'])==18 and sum(k.startswith('AuxC__') for k in q['checks'])==6 and sum(k.startswith('CP__') for k in q['checks'])==6 and sum(k.startswith('paired_CP__') for k in q['checks'])==3 and sum(k.startswith('paired_exchange__') for k in q['checks'])==2)
cid=next(c['cell_id'] for c in p['cells'] if c['metal']=='Dy' and c['donor']=='water' and c['role']=='complex');bad=copy.deepcopy(synthetic);bad[cid]['energy_hartree']+=.3/a.b.HA;s=a.score(p,bad)
check('CP_error_sign_and_original_point25',abs(s['CP_errors_kcal_mol']['water__Dy']-.3)<1e-8 and not s['checks']['CP__water__Dy'])
check('paired_water_reference_exchange_sign',all(abs(v[0]-v[1]+.3)<1e-8 for v in s['paired_exchange_local_canonical'].values()) and all(not v for k,v in s['checks'].items() if k.startswith('paired_exchange__')))
bad=copy.deepcopy(synthetic);bad[cid]['RHF_energies_hartree'][0]+=2e-6;check('original_RHF_gate',not a.score(p,bad)['checks']['RHF__'+cid]);bad=copy.deepcopy(synthetic);bad[cid]['aux_C']={'mismatch':True};check('original_AuxC_gate',not a.score(p,bad)['checks']['AuxC__water__Dy']);bad=copy.deepcopy(synthetic);bad.pop(cid);rejects('missing18_noqualification',lambda:a.score(p,bad));bad=copy.deepcopy(synthetic);bad[cid]['status']='unavailable';rejects('failed18_noqualification',lambda:a.score(p,bad))
# Collector orchestration only. The real three-role admission was performed in actual preparation.
original_load,original_reuse,original_admit=a.load,a.reused_rows,a.admit_result;a.load=lambda _:p;a.reused_rows=lambda *args:copy.deepcopy(p['reused'])
try:
    inventory={cid:None for cid in p['new_cell_ids']};ip=out/'SYNTHETIC_MISSING15.json';a.b.write(ip,inventory);q=a.collect(pp,ip,out/'missing15')
    check('missing15_preserves18_and_reused3',q['status']=='unavailable' and q['complete']==q['reused_complete']==3 and q['new_complete']==0 and len(q['rows'])==18 and q['all_controls_passed'] is False and 'checks' not in q)
    a.admit_result=lambda path,c,*args:dict(copy.deepcopy(synthetic[c['cell_id']]),reused_from_tighter_scout=False)
    ip=out/'SYNTHETIC_NEW15.json';a.b.write(ip,{cid:'SYNTHETIC_ONLY_'+cid for cid in p['new_cell_ids']});q=a.collect(pp,ip,out/'mock15_complete')
    check('mock15_plus_real3_full35',q['complete']==18 and q['new_complete']==15 and q['reused_complete']==3 and len(q['checks'])==35 and q['all_controls_passed'])
    check('passing35_never_releases_clusters',q['cluster_release_authorized'] is False and q['physical_accuracy_qualified'] is False)
finally:a.load,a.reused_rows,a.admit_result=original_load,original_reuse,original_admit
prep=a.b.pin(pp.with_name('RESULT.json'));rev=out/'SYNTHETIC_SOURCE_REVIEW.json';a.b.write(rev,{'status':'passed_tighter_controls_source_review','implementation':a.b.pin(H/'REFERENCE.py'),'panel':a.b.pin(pp),'preparation':prep})
r={'status':'released_tighter_controls15','implementation':a.b.pin(H/'REFERENCE.py'),'panel':a.b.pin(pp),'parser':p['parser'],'declaration':p['declaration'],'cell_ids':p['new_cell_ids'],'maximum_new_native_calls':15,'ranks':8,'nodes':1,'execution_root':p['execution_root'],'cluster_release_authorized':False,'scout_actual_review':p['scout_actual_review'],'scout_collection':p['scout_collection'],'preparation':prep,'source_review':a.b.pin(rev)}
check('concrete_oneway_release',a.release_gate(p,pp,r)==r);bad=dict(r,cell_ids=p['new_cell_ids']+p['reused_cell_ids']);rejects('reject_reexecution_of_reused3',lambda:a.release_gate(p,pp,bad));bad=dict(r,scout_actual_review=None);rejects('reject_unbound_actual_scoutreview',lambda:a.release_gate(p,pp,bad));bad=dict(r,cluster_release_authorized=True);rejects('reject_automatic_clusterrelease',lambda:a.release_gate(p,pp,bad))
result={'status':'passed_saved_controls','implementation':a.b.pin(H/'REFERENCE.py'),'controls_source':a.b.pin(__file__),'panel':a.b.pin(pp),'checks':checks,'check_count':len(checks),'all_passed':all(checks.values()),'new_native_calls':0,'synthetic_new_native_calls':0,'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1};a.b.write(out/'RESULT.json',result);print(json.dumps(result,indent=2))
