"""Focused saved controls for the three-role change; never starts a process."""
import copy, importlib.util, json, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('tighter_three',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main(panel,out):
    wall,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False)
    checks={}
    def ok(name,value):
        checks[name]=bool(value);assert value,name
    def rejects(name,fn):
        try:fn()
        except (AssertionError,ValueError,RuntimeError,KeyError,StopIteration):checks[name]=True
        else:raise AssertionError(name)
    m.subprocess.run=lambda *a,**k:(_ for _ in ()).throw(AssertionError('FORBIDDEN_NATIVE_CALL'))
    p=m.load(panel);parent=m.rec.load(m.b.verify(p['parent_panel']));failed=m.read(p['failed_controls'])
    ok('exact_three_roles',[c['role'] for c in p['cells']]==m.ROLES)
    ok('all_tighter_Dy_acetate',all(c['tier']=='tighter' and c['metal']=='Dy' and c['donor']=='acetate_bi' for c in p['cells']))
    ok('default_32GB_grant',m.resources(dict(SLURM_JOB_ID='123',SLURM_NTASKS='8',SLURM_MEM_PER_NODE='32000'))==(8,32000,3000))
    ok('actual_grant_changes_MaxCore',m.resources(dict(SLURM_JOB_ID='123',SLURM_NTASKS='8',SLURM_MEM_PER_CPU='5000'))==(8,40000,3750))
    rejects('reject_four_ranks',lambda:m.resources(dict(SLURM_JOB_ID='123',SLURM_NTASKS='4',SLURM_MEM_PER_NODE='32000')))
    rejects('reject_missing_grant',lambda:m.resources(dict(SLURM_JOB_ID='123',SLURM_NTASKS='8')))
    rejects('reject_login_native',lambda:m.resources(dict(SLURM_NTASKS='8',SLURM_MEM_PER_NODE='32000')))
    simulated={};baseline_settings={}
    for c in p['cells']:
        cid=c['cell_id'];q=failed['rows'][c['source_tight_cell_id']];raw=m.b.verify(q['output']).read_text()
        oldc=next(x for x in parent['cells'] if x['cell_id']==c['source_tight_cell_id'])
        parsed=m.parse(raw,oldc,0)
        ok(cid+'_old_native_admits',parsed['status']=='complete')
        baseline_settings[cid]=dict(output=q['output'],settings=parsed['actual_triples_settings'])
        mislabeled=m.parse(raw,c,0)
        ok(cid+'_old_output_cannot_be_tighter',mislabeled['status']=='failed_admission' and all(not mislabeled['checks']['local_'+k] for k in ('TCutPairs','TCutPNO','TCutMKN')))
        ok(cid+'_missing_triples_readback_fails',not m.triples_settings(raw.replace('TNOScaleW','MissingScaleW'))['complete_and_finite'])
        ok(cid+'_malformed_triples_readback_fails',not m.triples_settings(raw.replace('TNOScaleW                   ... 1.000e+02','TNOScaleW                   ... NaN'))['complete_and_finite'])
        simulated[cid]=q
    sc=m.comparison(p,simulated)
    ok('saved_Tight_error_reproduced',abs(sc['CP_errors_vs_canonical_kcal_mol']['tighter']['total']-.44369751639204386)<1e-7)
    ok('saved_Tight_remains_failed',not sc['diagnostic_limit_passed'])
    exact=copy.deepcopy(simulated)
    for c in p['cells']:
        q=exact[c['cell_id']];v=p['baseline_components'][c['source_tight_cell_id']]['canonical']['values']
        for key,term in [('CCSD_energies_hartree','CCSD'),('triples_hartree','triples'),('RHF_energies_hartree','RHF')]:q[key]=[float(v[term])]
        q['energy_hartree']=float(v['total'])
    sc=m.comparison(p,exact);ok('synthetic_canonical_self_target',sc['all_scout_diagnostics_passed'] and abs(sc['CP_errors_vs_canonical_kcal_mol']['tighter']['total'])<1e-10)
    bad=copy.deepcopy(exact);bad[p['cells'][0]['cell_id']]['energy_hartree']+=.251/m.b.HA
    ok('original_point25_limit_unchanged',not m.comparison(p,bad)['diagnostic_limit_passed'])
    bad=copy.deepcopy(exact);bad[p['cells'][1]['cell_id']]['aux_C']['counts'][0]+=1
    ok('mismatched_AuxC_fails',not m.comparison(p,bad)['diagnostic_checks']['same_parent_AuxC'])
    bad=copy.deepcopy(exact);bad.pop(p['cells'][0]['cell_id']);rejects('partial_CP_rejected',lambda:m.comparison(p,bad))
    prep=m.read(m.b.pin(Path(panel).parent/'RESULT.json'));original_read=m.read
    synthetic_review=dict(status='passed_tighter_dy_acetate_source_review',implementation=p['adapter'],panel=m.b.pin(panel),preparation=m.b.pin(Path(panel).parent/'RESULT.json'))
    release=dict(status='approved_tighter_dy_acetate_attribution',implementation=p['adapter'],parser=p['parser'],panel=m.b.pin(panel),declaration=p['declaration'],
        cell_ids=[c['cell_id'] for c in p['cells']],maximum_native_calls=3,ranks=8,execution_root=p['execution_root'],preparation=m.b.pin(Path(panel).parent/'RESULT.json'),
        source_review={'path':'synthetic-control-only'},full_eighteen_qualified=False,cluster_release_authorized=False)
    m.read=lambda pin:synthetic_review if pin==release['source_review'] else original_read(pin)
    try:
        ok('synthetic_release_shape_accepted',m.release_contract(p,panel,release)==release)
        for key,value in [('maximum_native_calls',4),('ranks',4),('cell_ids',release['cell_ids'][:-1]),('cluster_release_authorized',True),('full_eighteen_qualified',True)]:
            bad=dict(release);bad[key]=value;rejects('release_rejects_'+key,lambda bad=bad:m.release_contract(p,panel,bad))
    finally:m.read=original_read
    inventory={c['cell_id']:None for c in p['cells']};f=out/'UNAVAILABLE_INVENTORY.json';m.b.write(f,inventory)
    absent=m.collect(panel,f,out/'unavailable_collection')
    ok('full_three_missing_denominator',absent['declared']==3 and absent['complete']==0 and len(absent['rows'])==3 and absent['status']=='unavailable')
    ok('unavailable_cannot_release',absent['full_eighteen_qualified'] is False and absent['cluster_release_authorized'] is False and 'CP_components_kcal_mol' not in absent)
    q=dict(status='passed_saved_controls',implementation=m.b.pin(m.__file__),test_implementation=m.b.pin(__file__),panel=m.b.pin(panel),checks=checks,
        baseline_actual_triples_settings=baseline_settings,new_native_calls=0,synthetic_values_are_controls_only=True,threads=1,
        wall_seconds=time.monotonic()-wall,parent_cpu_seconds=time.process_time()-cpu)
    print(json.dumps(dict(result=m.b.write(out/'RESULT.json',q),passed=len(checks))))

if __name__=='__main__':main(sys.argv[1],Path(sys.argv[2]))
