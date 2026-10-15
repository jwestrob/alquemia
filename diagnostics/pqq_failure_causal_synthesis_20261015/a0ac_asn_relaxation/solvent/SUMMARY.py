#!/usr/bin/env python3
"""Finite saved-result aggregation only; numerical execution requires Slurm."""
import argparse, csv, hashlib, importlib.util, json, math, os, shutil, subprocess, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/solvent'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
POOLS=('common3','common4','common6'); POLICIES=('operational','mathematical','origin')
METALS=('Ca','La'); COMPONENTS={'Gelec','Ghb','Gsasa','Gshift'}
def require(v,s):
    if not v: raise ValueError(s)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):
    require(pin(p['path'])==p,'Changed artifact '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def freeze(output):
    base=WORK/'asn_relaxation_solvent_v1'; old=ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1'
    coll=pin(base/'execution/FINAL_COLLECTION.json');require(coll['sha256']=='64668c7ec758da42d34b039855bdc370b39d06815db47a00dfaf69c41a69368d','Wrong completed4971 collection')
    j=read(coll);prior=read(old/'FINAL_COLLECTION.json');native=[v['native_result'] for row in j['rows'] for z in METALS for v in row['evidence'][z].values()]
    for p in native:check(p)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir()
    impl={}
    for name in ('SUMMARY.py','SUMMARY.sbatch','SUMMARY_PLAN.md'):
        shutil.copyfile(HERE/name,out/'implementation'/name);impl[name]=pin(out/'implementation'/name)
    c={'implementation':impl,'collection':coll,'execution':pin(base/'execution/EXECUTION.json'),'preflight':pin(base/'execution/PREFLIGHT.json'),
       'canonical':pin(old/'CANONICAL_COLLECTION.json'),'old_collection':pin(old/'FINAL_COLLECTION.json'),'scorer':prior['scorer'],
       'contact_result':pin(WORK/'asn_relaxation_contact_review_v1/execution/RESULT.json'),'native_results':native,
       'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11'),'output':str(out/'execution'),
       'new_molecular_calls':0,'new_calibration':False}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))
def table(path,rows):
    with Path(path).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Real Slurm required before numeric aggregation')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp)
    require(c['implementation']['SUMMARY.py']==pin(__file__) and c['python']==pin(sys.executable),'Implementation/interpreter changed')
    raw=subprocess.check_output(['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']],text=True)
    fields=dict(x.split('=',1) for x in raw.split() if '=' in x)
    require((fields.get('JobState'),fields.get('NumNodes'),fields.get('NumCPUs'))==('RUNNING','1','1'),'Expected one running CPU/node')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();ex={'status':'failed','config':cp,'job_id':os.environ['SLURM_JOB_ID'],'placement':raw,'new_molecular_calls':0}
    try:
        for p in c['implementation'].values():check(p)
        j=read(c['collection']);old=read(c['old_collection']);cal=read(c['canonical']);e=read(c['execution']);pf=read(c['preflight']);contact=read(c['contact_result'])
        require(j['terminal'] and j['declared_sources']==len(j['rows'])==10 and j['available_new_cells']==40 and j['reused_original_cells']==80,'Incomplete declared denominator')
        require(e['status']=='complete' and e['SCC_terminal']==160 and e['native_component_calls']==80 and pf['status']=='passed','Actual execution/admission incomplete')
        require(old['canonical_calibration_receipt']==c['canonical'] and cal['scorer']==old['scorer']==c['scorer'],'Old-band provenance differs')
        spec=importlib.util.spec_from_file_location('_saved_scorer',check(c['scorer']));s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
        oldrows={r['record_id']:r for r in old['rows']};source_rows=[];flat=[];all_checks=[];parity=[];bounds=[];absolute=[];nativepins=[]
        for row in j['rows']:
            previous=oldrows[row['record_id']];require(previous['source_case_id']==row['source_case_id'],'Old source join differs')
            saved={'source_case_id':row['source_case_id'],'conditioning':row['conditioning'],'sample':row['sample'],'expected_class':previous['expected_class'],'pools':{},'changes':{k:row[k] for k in ('common3_to4','common4_to6','common3_to6')}}
            for z in METALS:
                require(set(row['evidence'][z])==set(j['pool_orders']['common6']),'Incomplete six-candidate evidence')
                for q,v in row['evidence'][z].items():
                    require(v['status']=='available' and v['component_v2']['accepted'],'Saved component admission failed')
                    n=read(v['native_result']);nativepins.append(v['native_result']);require(n['status']=='available' and n['G_qvac_hartree']==v['G_qvac_hartree'],'Native result mismatch')
                    require({a['component'] for a in n['checks']}==COMPONENTS and all(a['passed'] and a['informative'] for a in n['checks']),'Native diagonal unqualified')
                    require({a['component'] for a in n['vacuum_precision']}==COMPONENTS and all(a['informative'] for a in n['vacuum_precision']),'Vacuum precision unqualified')
                    require(all(n['water'][k]['sha256']==n['vacuum'][k]['sha256'] for k in ('cavity','born_matrix_with_HB')),'Paired native cavity/matrix differ')
                    all_checks.extend(n['checks']);bounds.extend(n['vacuum_precision']);parity.append(v['component_v2']);absolute.extend(v['historical_absolute_point_accepted'].values())
            for p in POLICIES:
                require(abs(row['pools']['common3'][p]['R_model_kcal_mol']-previous['pool'][p]['R_model_kcal_mol'])<=j['saved_algebra_checks']['tolerance_kcal_mol'],'Historical common3 score mismatch')
            for pool in POOLS:
                saved['pools'][pool]={}
                for policy in POLICIES:
                    x=row['pools'][pool][policy];v=x['R_model_kcal_mol'];require(math.isfinite(v),'Nonfinite saved score')
                    d=s.decision(v,cal['calibrations'][policy]['bands'],previous['expected_class'])
                    saved['pools'][pool][policy]={'R_kcal_mol':v,'selected':{z:x[z]['candidate'] for z in METALS},'old_band_sensitivity':d}
                    flat.append({'source':row['source_case_id'],'conditioning':row['conditioning'],'sample':row['sample'],'pool':pool,'policy':policy,'R_kcal_mol':v,'Ca_candidate':x['Ca']['candidate'],'La_candidate':x['La']['candidate'],'OLD_band_call':d['decision'],'OLD_band_outcome':d['outcome']})
            source_rows.append(saved)
        require(nativepins==c['native_results'] and len(nativepins)==120,'Exact120 native result joins differ')
        ranges=[]
        for group in ('Ca','La','all10'):
            members=[r for r in source_rows if group=='all10' or r['conditioning']==group];require(len(members)==(10 if group=='all10' else 5),'Group denominator changed')
            for pool in POOLS:
                for policy in POLICIES:
                    values=[r['pools'][pool][policy]['R_kcal_mol'] for r in members]
                    ranges.append({'conditioning_group':group,'pool':pool,'policy':policy,'count':len(values),'minimum_R_kcal_mol':min(values),'maximum_R_kcal_mol':max(values),'range_kcal_mol':max(values)-min(values)})
        result={'config':cp,'status':'complete','rows':source_rows,'ranges':ranges,'old_bands':cal['calibrations'],'new_protocol_calibration':'unavailable_not_evaluated',
          'saved_admission':{'native_cells':120,'all_component_v2_accepted':True,'all_diagonal_checks_passed':True,'all_vacuum_bounds_informative':True,'all_cavity_matrix_pairs_identical':True,
             'maximum_water_diagonal_absolute_residual_Eh':max(abs(x['residual_Eh']) for x in all_checks),'maximum_vacuum_rounding_bound_Eh':max(x['rounding_bound_Eh'] for x in bounds),
             'maximum_component_parity_error_Eh':max(abs(x['component_error_hartree']) for x in parity),'maximum_population_error_e':max(x['maximum_population_error_e'] for x in parity),
             'historical_absolute_point_accepted':sum(absolute),'historical_absolute_point_failed':len(absolute)-sum(absolute)},
          'contact_all_geometries_available':contact['all_candidate_geometries_available'],'actual_solvent_execution':c['execution'],'new_molecular_calls':0,'production_changed':False}
        table(out/'SCORES.tsv',flat);table(out/'RANGES.tsv',ranges)
        try:
            import matplotlib
            matplotlib.use('Agg');matplotlib.rcParams['svg.fonttype']='none';import matplotlib.pyplot as plt
            anchor,=[r for r in cal['rows'] if r['source_case_id']=='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-3']
            offset=anchor['pool']['mathematical']['R_model_kcal_mol'];fig,ax=plt.subplots(figsize=(8,5))
            for row in source_rows:
                vals=[row['pools'][p]['mathematical']['R_kcal_mol']-offset for p in POOLS]
                ax.plot(range(3),vals,marker='o',label=f"{row['conditioning']}{row['sample']}",color=plt.get_cmap('tab10')(row['sample']),linestyle='-' if row['conditioning']=='Ca' else '--')
            ax.set_xticks(range(3),['Original common3','+ fixed amide flip\ncommon4','+ relaxed flips\ncommon6']);ax.set_ylabel('R − fixed common3 canonical A0AC La3 R (kcal/mol)');ax.set_title('A0AC alternative-amide response: all10 sources');ax.grid(alpha=.2);ax.legend(ncol=2,title='Source conditioning / sample',fontsize=8);fig.tight_layout()
            for ext in ('svg','png'):fig.savefig(out/('COMMON_POOLS.'+ext),dpi=180)
            plt.close(fig);result['plot']={'status':'complete','matplotlib_version':matplotlib.__version__,'fixed_reference':anchor['record_id'],'offset_R_kcal_mol':offset}
        except Exception as err:result['plot']={'status':'unavailable','reason':str(err)}
        write(out/'RESULT.json',result);ex.update(status='complete',result=pin(out/'RESULT.json'))
    except BaseException as err:ex.update(reason=str(err),traceback=traceback.format_exc());raise
    finally:ex['wall_seconds']=time.monotonic()-start;write(out/'EXECUTION.json',ex)
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);a=s.add_parser('freeze');a.add_argument('output');a=s.add_parser('execute');a.add_argument('config');a.add_argument('sha');a=p.parse_args()
    if a.action=='freeze':freeze(a.output)
    else:execute(a.config,a.sha)
