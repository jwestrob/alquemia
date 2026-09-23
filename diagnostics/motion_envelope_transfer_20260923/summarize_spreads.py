"""Predeclared structural-spread comparison on exactly91 complete matched triples."""
from pathlib import Path
from statistics import median
import json,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new
base=ROOT/'workspaces/motion_envelope_transfer_20260923/COMPARISON_v1.json';x=read_json(base)
rows=[t for t in x['triples']if t['methods']['envelope_operational']['complete']]
if len(rows)!=91:raise RuntimeError('frozen primary matched denominator differs')
comparisons={}
for method in ('envelope_origin','released_static','tenfold_strict'):
 pairs=[]
 for t in rows:
  a=t['methods'][method];b=t['methods']['envelope_operational']
  if not a['complete']:raise RuntimeError('baseline missing on declared matched rows')
  pairs.append({'triple_id':t['triple_id'],'protein_id':t['protein_id'],'baseline_range_model_kcal_mol':a['range'],
    'envelope_adaptive_range_model_kcal_mol':b['range'],'delta_range_model_kcal_mol':b['range']-a['range']})
 comparisons[method]={'matched_triples':91,'range_decreases':sum(r['delta_range_model_kcal_mol']<0 for r in pairs),
 'range_increases':sum(r['delta_range_model_kcal_mol']>0 for r in pairs),'unchanged':sum(r['delta_range_model_kcal_mol']==0 for r in pairs),
 'median_baseline_range_model_kcal_mol':median(r['baseline_range_model_kcal_mol']for r in pairs),
 'median_envelope_adaptive_range_model_kcal_mol':median(r['envelope_adaptive_range_model_kcal_mol']for r in pairs),
 'median_paired_range_change_model_kcal_mol':median(r['delta_range_model_kcal_mol']for r in pairs),'rows':pairs}
result={'comparison':record(base),'implementation':record(__file__),'scope':'same91 complete overlapping triples,25 consumed proteins; no independent-sample claim',
 'comparisons':comparisons,'new_molecular_calls':0,'new_thresholds':False}
d=Path(__file__).parent;write_new(d/'SPREAD.json',result)
lines=['# Matched triple ranges','', 'All comparisons use exactly the same91 complete triples. A range is maximum minus minimum raw R among its three declared members. These overlapping subsets are correlated. No member is removed; no score/reference is changed.','', '| Baseline → envelope adaptive | Range decreases | Increases | Equal | Median baseline range | Median envelope range | Median paired change |','|---|---:|---:|---:|---:|---:|---:|']
for name,r in comparisons.items():lines.append(f"| {name} | {r['range_decreases']} | {r['range_increases']} | {r['unchanged']} | {r['median_baseline_range_model_kcal_mol']:.6f} | {r['median_envelope_adaptive_range_model_kcal_mol']:.6f} | {r['median_paired_range_change_model_kcal_mol']:.6f} |")
lines+=['','Units are model kcal/mol. Median paired change is the median of within-triple differences, not the difference of medians. These range summaries describe structural variation and do not by themselves establish improved accuracy. Exact unrounded values and all91 pairs are in [SPREAD.json](SPREAD.json).']
(d/'SPREAD.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:{a:b for a,b in v.items()if a!='rows'}for k,v in comparisons.items()},indent=2))
