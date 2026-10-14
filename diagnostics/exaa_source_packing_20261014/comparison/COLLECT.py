#!/usr/bin/env python3
"""Pinned saved-output comparison only; execute under Slurm, without calculators."""
import argparse, hashlib, json, os, resource, socket, time, traceback
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
POLICIES=('origin','operational','mathematical')
METALS=('Ca','La')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,('changed input',p);return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def prepare(a):
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    review=read(a.review);scores=read(a.scores);assert scores['input']==pin(a.review)
    canonical=ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/CANONICAL_COLLECTION.json'
    reference=read(canonical);assert reference['calibrations']==review['calibrations']
    selected=[r for r in review['cases'] if r['record']['biological_group']=='q9z4j7']
    expected={f'q9z4j7-pqq-la_model__conditioned_{condition}__seed-1_sample-{i}' for condition in ('La','2Ca_v1') for i in range(5)}
    assert len(selected)==10 and {r['record']['source_case_id'] for r in selected}==expected
    scored={r['record_id']:r for r in scores['rows']}
    for r in selected:
        assert scored[r['record_id']]['case_key']==r['case_key']
        assert scored[r['record_id']]['source_case_id']==r['record']['source_case_id']
    missing=[rid for rid in review['original_unavailable'] if '::q9z4j7-' in rid]
    inventory=read(review['inventory']);records={r['record_id']:r for r in inventory['records']}
    excluded=[records[rid] for rid in missing];assert len(excluded)==5
    assert {r['source_case_id'] for r in excluded}=={f'q9z4j7-pqq-la_model__conditioned_Ca__seed-1_sample-{i}' for i in range(5)}
    native=read(a.native_config);prep=read(native['preparation'])
    assert prep['denominator']==2 and [r['chain'] for r in prep['rows']]==['A','B']
    assert [r['chain'] for r in native['declared']]==['A','B']
    c={'protocol':'ExaA_conditional_crystal_saved_source_comparison_v1','implementation':pin(__file__),
       'plan':pin(HERE/'PLAN.md'),'wrapper':pin(HERE/'RUN.sbatch'),'review':pin(a.review),'scores':pin(a.scores),
       'inventory':review['inventory'],'canonical':pin(canonical),'native_config':pin(a.native_config),
       'crystal_preparation':native['preparation'],'saved_record_ids':[r['record_id'] for r in selected],
       'saved_source_ids':[r['record']['source_case_id'] for r in selected],'original_exclusions':excluded,
       'crystal_declaration':native['declared'],'frozen_calibrations':reference['calibrations'],
       'crystal_collection_path':str(Path(a.crystal_collection).resolve()),'output':str(out/'execution'),
       'denominators':{'current_modeled_sources':10,'crystal_chains':2,'current_comparison':12,
                       'historical_single_Ca_exclusions':5,'all_visible_records':17,'biological_proteins':1},
       'upstream_jobs':{'native':'4095','solvent':'4096'},'new_molecular_calls':0,'new_calibration':False}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def policy_view(pool,decisions,calibrations):
    rows={}
    for p in POLICIES:
        v=pool.get(p) if pool else None
        rows[p]={'status':'available' if v else 'unavailable','bands':calibrations[p]['bands'],
          'R_model_kcal_mol':v['R_model_kcal_mol'] if v else None,
          'decision':decisions.get(p,{'decision':None,'outcome':'unavailable'}),
          'metals':{z:{'candidate':v[z]['candidate'],'work_kcal_mol':v[z]['work_kcal_mol'],
                      'all_works_kcal_mol':v[z]['all_works_kcal_mol']} for z in METALS} if v else None}
    return rows

def run(path):
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_NODELIST'),'actual Slurm required'
    assert int(os.environ['SLURM_CPUS_ON_NODE'])==1
    c=read(path);assert pin(__file__)==c['implementation'];verify(c['plan']);verify(c['wrapper'])
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    receipt={'status':'failed','job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'config':pin(path),'new_molecular_calls':0}
    try:
        review=read(c['review']);old=read(c['scores']);canonical=read(c['canonical']);native=read(c['native_config']);prep=read(c['crystal_preparation'])
        assert canonical['calibrations']==review['calibrations']==c['frozen_calibrations']
        assert old['input']==c['review'] and native['preparation']==c['crystal_preparation']
        refs={r['record_id']:r for r in review['cases']};saved={r['record_id']:r for r in old['rows']};rows=[]
        for rid in c['saved_record_ids']:
            ref=refs[rid];r=ref['record'];s=saved[rid];assert s['case_key']==ref['case_key'] and s['source_case_id']==r['source_case_id']
            assert s['selected_call']==ref['score']['decisions']['operational']['decision']
            assert s['separately_calibrated_origin_call']==ref['score']['decisions']['origin']['decision']
            rows.append({'source_case_id':r['source_case_id'],'record_id':rid,'case_key':ref['case_key'],
              'group':r['source_conditioning_metal'],'source_kind':'saved_predicted_source','historical_exclusion':False,'status':s['status'],'reason':s['reason'],
              'original_source_id':r.get('original_source_id'),'family':r['family'],'role':r['role'],'aliases':ref['aliases'],
              'policies':policy_view(ref['score']['pool'],ref['score']['decisions'],c['frozen_calibrations']),
              'saved_component_works':s['works'],'saved_accommodation_R_shift':s['accommodation_R_shift'],
              'origin_under_accommodated_band_diagnostic':s['same_band_origin_call'],
              'state':s['state'],'state_signature':r['state_signature'],'source_preparation':r['source_preparation'],
              'saved_score_evidence':c['scores'],'prepared_admission':r.get('prepared_admission')})
        for r in c['original_exclusions']:
            rows.append({'source_case_id':r['source_case_id'],'record_id':r['record_id'],'case_key':r['source_case_id']+'::'+r['selection_id'],
              'group':'Ca_historical_exclusion','source_kind':'original_single_Ca_preparation','historical_exclusion':True,'status':'unavailable',
              'reason':r['preparation_reason'],'preparation_status':r['preparation_status'],'input_status':r['input_status'],
              'policies':policy_view(None,{},c['frozen_calibrations']),'inventory_evidence':c['inventory']})
        cp=Path(c['crystal_collection_path']);upstream={'collection':pin(cp) if cp.is_file() else None,'receipts':{},'reason':None};crystal={}
        for label,rp in [('solvent_execution',cp.parent/'EXECUTION.json'),('native_collection',Path(c['native_config']['path']).parent/'FINAL_COLLECTION.json')]:
            if rp.is_file():upstream['receipts'][label]={'pin':pin(rp),'record':read(rp)}
        if cp.is_file():
            try:
                data=read(cp);nc=read(data['config']);assert nc['native_config']==c['native_config']
                assert data['frozen_calibration']==c['canonical'] and data['frozen_calibrations']==c['frozen_calibrations']
                assert data['denominator']==2 and {r['chain'] for r in data['rows']}=={'A','B'}
                assert {r['case_id'] for r in data['rows']}=={r['source_case_id'] for r in c['crystal_declaration']}
                crystal={r['case_id']:r for r in data['rows']};upstream['config']=data['config']
            except Exception as exc:upstream['reason']='unaccepted_collection: '+type(exc).__name__+': '+str(exc)
        else:upstream['reason']='upstream_collection_missing_after_dependency_terminal'
        for declared in c['crystal_declaration']:
            admission=read(declared['admission']);assert admission==next(r for r in prep['rows'] if r['chain']==declared['chain'])
            r=crystal.get(declared['source_case_id']);reason=(r.get('reason') if r else None) or declared['reason'] or upstream['reason']
            rows.append({'source_case_id':declared['source_case_id'],'record_id':None,'case_key':None,'chain':declared['chain'],
              'group':'1FLG_crystal','source_kind':'conditional_oxidized_PQQ_crystal','historical_exclusion':False,'status':r['status'] if r else 'unavailable',
              'reason':reason,'admission':declared['admission'],'preparation_status':declared['preparation_status'],
              'policies':policy_view(r['pool'] if r else None,r['conditional_decisions'] if r else {},c['frozen_calibrations']),
              'saved_component_works':r['works'] if r else None,'selected_differential_works':r['selected'] if r else None,
              'failures':r.get('failures',[]) if r else [],'native_pool_complete':r.get('native_pool_complete') if r else None,
              'state_signature':admission.get('state_signature'),'reference_state_signature':admission.get('reference_state_signature'),
              'local_state_identical_to_Q9':admission.get('local_state_identical_to_Q9'),'construct_compatibility':admission.get('construct_compatibility'),
              'terminal_completion':admission.get('terminal_completion'),'spectator_admission':admission.get('spectator_admission'),
              'source_redox_experimentally_qualified':False,'conditional_source_state':True,'source':admission.get('source')})
        assert len(rows)==17 and len({r['source_case_id'] for r in rows})==17
        result={'protocol':c['protocol'],'status':'collected','denominators':c['denominators'],'rows':rows,'upstream':upstream,
          'frozen_calibrations':c['frozen_calibrations'],'calibration':c['canonical'],'config':pin(path),
          'new_molecular_calls':0,'new_calibration':False,'new_class_decisions':False,'production_changed':False,
          'limitations':['Consumed development evidence for one protein; folds/chains are not independent proteins.',
             'Crystal PQQ oxidation state is assumed, not experimentally qualified.',
             'Crystal mature subset maps exactly where admitted; full constructs differ and terminal completion is documented.',
             'Conditioning and original preparation failures remain separate; no physiological occupancy inference.',
             '249-atom packing-context outputs are not assigned these calibrated bands.']}
        write(out/'RESULT.json',result)
        lines=['# Crystal and saved ExaA sources','','Existing scores only; no new fit or class decisions. Origin/static uses its own frozen bands.','',
          '| Policy | Ca maximum | La minimum |','|---|---:|---:|']
        for p in POLICIES:
            b=c['frozen_calibrations'][p]['bands'];lines.append(f'| {p} | {b["Ca_max"]:.6f} | {b["La_min"]:.6f} |')
        lines+=['','| Group | Source | Static | Operational | Mathematical | Ca selected work | La selected work |','|---|---|---|---|---|---:|---:|']
        for r in rows:
            pol=r['policies'];name=r.get('chain') or r['source_case_id'].rsplit('sample-',1)[-1]
            calls=[pol[p]['decision']['decision'] or 'unavailable' for p in POLICIES];m=pol['operational']['metals']
            w=[f'{m[z]["work_kcal_mol"]:.6f}' if m else 'unavailable' for z in METALS]
            lines.append('| '+' | '.join([r['group'],str(name),*calls,*w])+' |')
        lines+=['','## Availability and interpretation','','Current comparison: ten modeled sources plus both crystal chains. The five historical single-Ca exclusions remain visible as recovery provenance, not additional current observations. All records concern one biological protein.']
        for r in rows:
            if r['status']!='available':lines.append('- '+r['source_case_id']+': '+str(r.get('reason') or 'required endpoint unavailable; see RESULT.json'))
        lines+=['','Crystal redox remains unverified; the mature construct differs from the full predicted sequence. These are correlated development structures, not a new blind test. Full per-policy scores, all candidate works, source states and failure receipts are in RESULT.json.']
        (out/'REPORT.md').write_text('\n'.join(lines)+'\n');receipt.update(status='collected',result=pin(out/'RESULT.json'))
    except BaseException as exc:receipt.update(reason=type(exc).__name__+': '+str(exc),traceback=traceback.format_exc());raise
    finally:
        u=resource.getrusage(resource.RUSAGE_SELF);receipt.update(elapsed_seconds=time.monotonic()-started,CPU_seconds=u.ru_utime+u.ru_stime,peak_RSS_KiB=u.ru_maxrss)
        write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('prepare');a.add_argument('--review',required=True);a.add_argument('--scores',required=True);a.add_argument('--native-config',required=True);a.add_argument('--crystal-collection',required=True);a.add_argument('--output',required=True)
    a=s.add_parser('run');a.add_argument('--config',required=True);a=p.parse_args()
    prepare(a) if a.action=='prepare' else run(a.config)
