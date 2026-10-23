#!/usr/bin/env python3
"""Saved-only complete 68+12 common-parent-basis energy decomposition."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('common_basis_collection_adapter',HERE.parent/'REFERENCE.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
b=a.b
SIGNS={'aquo8_central_La':1,'aquo8_central_Dy':-1,'aquo9_central_La':-1,'aquo9_central_Dy':1}
QUANTITIES=('cluster_contrast','pair_contrast','cooperative_remainder')
def delta(v):return math.fsum(SIGNS[g]*v[g] for g in a.GIDS)
def finite(v):return isinstance(v,(float,int)) and not isinstance(v,bool) and math.isfinite(v)
def cid(g,i,role):return f'pbe0__aquo{g[4]}__{g.rsplit("_",1)[1]}__water_{i:02d}__{role}'
def reused_id(g,role):return f'pbe0__aquo{g[4]}__{g.rsplit("_",1)[1]}__{role}'
def unavailable(path,exc=None):return a.prior.unavailable(path,exc)

def calculate(p,new,reused,cp_cycle,pair):
    b.check(set(new)=={c['cell_id'] for c in p['cells']} and len(new)==68 and set(reused)==set(p['saved_reuse']) and len(reused)==12,'Full68/new12/reused denominator required')
    for row in list(new.values())+list(reused.values()):b.check(row.get('status')=='complete' and finite(row.get('electronic_energy_hartree')),'Unavailable/nonfinite term cannot enter decomposition')
    b.check(all(v.get('D4_energy_hartree')==0. for v in new.values()),'Ghost D4 forbidden')
    parents={}
    for g in a.GIDS:
        n=p['geometries'][g]['water_count'];E=reused[g]['electronic_energy_hartree']
        ion=reused[reused_id(g,'ion_ghost_shell')]['electronic_energy_hartree'];shell=reused[reused_id(g,'shell_ghost_ion')]['electronic_energy_hartree']
        b.check(reused[reused_id(g,'ion_ghost_shell')]['D4_energy_hartree']==reused[reused_id(g,'shell_ghost_ion')]['D4_energy_hartree']==0.,'Reused ghost D4 forbidden')
        d=[new[cid(g,i,'metal_water')]['electronic_energy_hartree'] for i in range(1,n+1)]
        w=[new[cid(g,i,'water_only')]['electronic_energy_hartree'] for i in range(1,n+1)]
        interactions=[math.fsum((di,-ion,-wi)) for di,wi in zip(d,w)]
        C=math.fsum((E,-ion,-shell));P=math.fsum(interactions);B=C-P
        direct=math.fsum([E,*[-v for v in d],(n-1)*ion,-shell,*w])
        b.check(abs(B-direct)<1e-9,'Per-parent N-fold ion closure failed')
        parents[g]=dict(water_count=n,complex_electronic_hartree=E,ion_electronic_hartree=ion,shell_electronic_hartree=shell,
            metal_water_electronic_hartree=d,water_electronic_hartree=w,pair_interactions_hartree=interactions,
            N_fold_ion_subtraction_hartree=n*ion,cluster_interaction_hartree=C,pair_sum_hartree=P,cooperative_remainder_hartree=B,
            expanded_remainder_hartree=direct,closure_residual_hartree=B-direct)
    electronic={'cluster_contrast':delta({g:v['cluster_interaction_hartree'] for g,v in parents.items()}),
        'pair_contrast':delta({g:v['pair_sum_hartree'] for g,v in parents.items()}),
        'cooperative_remainder':delta({g:v['cooperative_remainder_hartree'] for g,v in parents.items()})}
    b.check(abs(electronic['cluster_contrast']-cp_cycle['hartree']['CP_electronic'])<1e-9,'C does not replay prior CP8 electronic contrast')
    b.check(abs(electronic['cluster_contrast']-electronic['pair_contrast']-electronic['cooperative_remainder'])<1e-9,'Paired electronic C-P-B closure failed')
    d4={k:pair['reference_decomposition']['D4_energy_hartree']['hartree'][k] for k in QUANTITIES}
    b.check(all(finite(v) for v in d4.values()) and abs(d4['cluster_contrast']-d4['pair_contrast']-d4['cooperative_remainder'])<1e-12,'Saved D4 closure differs')
    actual_d4=delta({g:reused[g]['D4_energy_hartree'] for g in a.GIDS})
    b.check(abs(d4['cluster_contrast']-actual_d4)<1e-12 and abs(actual_d4-cp_cycle['hartree']['original_D4_unchanged'])<1e-12,'Saved full-complex D4 source differs')
    total={k:electronic[k]+d4[k] for k in QUANTITIES};kcal={k:total[k]*b.HA for k in QUANTITIES}
    model={k:pair['model_decomposition'][k] for k in QUANTITIES};own={k:pair['reference_decomposition']['final_energy_hartree']['kcal_mol'][k] for k in QUANTITIES}
    b.check(all(finite(v) for v in list(model.values())+list(own.values())),'Saved comparison unavailable')
    return dict(delta_convention='Dy9 - La9 - Dy8 + La8',per_parent_electronic=parents,
        electronic_hartree=electronic,saved_original_D4_hartree=d4,electronic_plus_original_D4_hartree=total,
        electronic_kcal_mol={k:v*b.HA for k,v in electronic.items()},saved_original_D4_kcal_mol={k:v*b.HA for k,v in d4.items()},
        electronic_plus_original_D4_kcal_mol=kcal,unchanged_model_kcal_mol=model,old_own_basis_reference_kcal_mol=own,
        model_minus_common_basis_reference_kcal_mol={k:model[k]-kcal[k] for k in QUANTITIES},
        model_minus_old_own_basis_reference_kcal_mol={k:model[k]-own[k] for k in QUANTITIES},
        common_minus_old_own_basis_reference_kcal_mol={k:kcal[k]-own[k] for k in QUANTITIES},
        CP8_electronic_replay_residual_hartree=electronic['cluster_contrast']-cp_cycle['hartree']['CP_electronic'],
        N_fold_ion_subtraction=True,ghost_terms_assumed_to_cancel=False,D4_ghost_evaluations=0,D4_additions_per_quantity=1)

def admit_reused(key,p,old,cp_result,data):
    if key in a.GIDS:
        q=a.prior.admit_original(key,old,data,a.helper);expected=a.AUX[key]
        native=a.read(q['result']);q['printed_scf_criteria']=native['printed_scf_criteria'];q['spin_squared_values']=native['spin_squared_values']
        q['physical_state_qualified']=native['physical_state_qualified'];q['method_accuracy_qualified']=native['method_accuracy_qualified']
    else:
        c=next(c for c in old['cells'] if c['cell_id']==key);expected=a.AUX[c['geometry_id']]
        source=cp_result['fragments'][key]['result'];q=a.prior.admit_fragment(b.verify(source),c,old,p['cp_panel'],a.helper)
    counts,passed=a.auxiliary_admission(b.verify(q['output']).read_text(),expected);b.check(passed,'Reused parent auxiliary basis changed')
    b.check(q['result']==p['saved_reuse'][key]['result'] and q['electronic_energy_hartree']==p['saved_reuse'][key]['electronic_energy_hartree'],'Reused preparation authority/energy differs')
    q.update(native_aux_J_basis_function_counts=counts,expected_aux_J_basis_functions=expected);return q

def report(r):
    lines=['# Common-parent-basis decomposition','',f"Status: **{r['status']}**. New fragments {r['new_complete']}/68; reused terms {r['reused_complete']}/12.",'']
    x=r['decomposition']
    if x is None:lines+=['The complete decomposition is unavailable. All missing or failed terms remain in their original slots.','']
    else:
        lines+=['| Quantity | Electronic | Saved D4 | Sum | Fixed model | Model − sum |','|---|---:|---:|---:|---:|---:|']
        for k in QUANTITIES:lines.append('| '+k+' | '+' | '.join(f'{x[f][k]:+.9f}' for f in ('electronic_kcal_mol','saved_original_D4_kcal_mol','electronic_plus_original_D4_kcal_mol','unchanged_model_kcal_mol','model_minus_common_basis_reference_kcal_mol'))+' |')
        lines+=['','All values are kcal/mol; paired order is Dy9 − La9 − Dy8 + La8.',
            f"Common-basis B minus the old own-basis B: {x['common_minus_old_own_basis_reference_kcal_mol']['cooperative_remainder']:+.9f} kcal/mol.",
            'Each parent pair sum subtracts its own full-basis ion N times. No ghost terms cancel by assumption. Original D4 contributions are displayed separately and added once.','']
    failures=[k for k,v in r['new_fragments'].items() if v['status']!='complete']+[k for k,v in r['reused_terms'].items() if v['status']!='complete']
    if failures:lines+=['Unavailable terms: '+', '.join(failures),'']
    if r['error']:lines+=['Collection error: '+r['error'],'']
    lines+=['Native output, allocation, charge/electron/ECP and orbital/Aux-J admissions, stopping signals and printed SCF flags remain in RESULT.json. Native energy-check convergence does not qualify every printed density/DIIS criterion, physical state or basis-limit accuracy.',
        '', 'B is a conditional nonadditive remainder, not a unique polarization or charge-transfer energy. This finite basis sensitivity does not replace fitting targets or release held validation/sampling.',
        '',f"New native allocated CPU seconds recorded: {r['cost']['recorded_new_allocated_CPU_seconds']}. Reused cost is separate. Collector made zero molecular calls."]
    return '\n'.join(lines)+'\n'

def collect(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text());p=a.load(b.verify(cfg['panel']));out=Path(out);out.mkdir(parents=True,exist_ok=False)
    ids=[c['cell_id'] for c in p['cells']]
    r=dict(status='common_basis_collection_unavailable',implementation=b.pin(__file__),configuration=b.pin(config),panel=cfg['panel'],
        new_fragments={cid:unavailable(cfg.get('result_paths',{}).get(cid,'MISSING_PATH')) for cid in ids},
        reused_terms={key:dict(status='not_admitted',result=row['result']) for key,row in p['saved_reuse'].items()},
        new_complete=0,reused_complete=0,decomposition=None,error=None,new_native_calls=0,fit_targets_changed=False,physical_accuracy_qualified=False,fresh_validation_released=False)
    try:
        b.check(cfg['implementation']==b.pin(__file__) and cfg['adapter']==p['adapter'],'Collection implementation differs')
        review=a.read(cfg['independent_review']);b.check(review['status']=='passed_collection_review' and review['implementation']==b.pin(__file__) and review['panel']==cfg['panel'],'Collection review/source/panel differs')
        b.check(set(cfg['result_paths'])==set(ids),'68-result configured denominator differs')
        old=a.helper.load(b.verify(p['cp_panel']));cp_result=a.read(p['cp_collection']);pair=a.read(p['pair_collection']);data=a.read(p['consumed_input_data'])
        b.check(cp_result['status']=='coordination_cp_collection_complete' and pair['status']=='matched_pair_collection_complete','Prior complete collections unavailable')
        r.update(cp_collection=p['cp_collection'],pair_collection=p['pair_collection'],declaration=p['declaration'],independent_review=cfg['independent_review'])
        for key in r['reused_terms']:
            try:r['reused_terms'][key]=admit_reused(key,p,old,cp_result,data)
            except Exception as exc:r['reused_terms'][key]=unavailable(p['saved_reuse'][key]['result']['path'],exc)
        for c in p['cells']:
            key=c['cell_id'];path=cfg['result_paths'][key]
            try:r['new_fragments'][key]=a.admit_result(path,c,p,b.verify(cfg['panel']))
            except Exception as exc:r['new_fragments'][key]=unavailable(path,exc)
        r['new_complete']=sum(v['status']=='complete' for v in r['new_fragments'].values());r['reused_complete']=sum(v['status']=='complete' for v in r['reused_terms'].values())
        if r['new_complete']==68 and r['reused_complete']==12:
            r['decomposition']=calculate(p,r['new_fragments'],r['reused_terms'],cp_result['cycle'],pair);r['status']='common_basis_collection_complete'
        else:r['status']='common_basis_collection_incomplete'
    except Exception as exc:r['error']=repr(exc)
    recorded_cost=0.;missing_cost=[]
    for key,row in r['new_fragments'].items():
        value=row.get('allocated_CPU_seconds')
        if value is None and row.get('result'):
            try:value=a.read(row['result']).get('allocated_cpu_seconds')
            except Exception:pass
        if finite(value):recorded_cost+=value
        else:missing_cost.append(key)
    r['cost']=dict(recorded_new_allocated_CPU_seconds=recorded_cost,new_rows_without_cost=missing_cost,
        reused_allocated_CPU_seconds=sum(v.get('allocated_CPU_seconds') or 0 for v in r['reused_terms'].values() if v['status']=='complete'),
        collector_wall_seconds=time.monotonic()-tick,collector_CPU_seconds=time.process_time()-cpu,local_threads=1,collector_native_calls=0)
    b.write(out/'RESULT.json',r);(out/'REPORT.md').write_text(report(r));print(json.dumps({k:r[k] for k in ('status','new_complete','reused_complete','error')}));return r

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True,type=Path);args=ap.parse_args();raise SystemExit(0 if collect(args.config,args.out)['status']=='common_basis_collection_complete' else 1)
