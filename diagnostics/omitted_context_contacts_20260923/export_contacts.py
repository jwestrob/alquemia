"""Compact table/summary of the frozen geometry result; no new molecular calls."""
import argparse
import csv
import json
from pathlib import Path
import statistics
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from affordable_common import read_json,record,write_new

def export(result,output):
 d=read_json(result);out=Path(output);out.mkdir(parents=True,exist_ok=False)
 fields=['pair_id','source_id','protein_id','representation','geometry','query_atoms','outside_atoms','moving_atoms','minimum_all_A','count_all_below_3_5_A','count_all_below_2_A','new_all_below_2_A','minimum_moving_A','count_moving_below_3_5_A']
 rows=[];summary={}
 for r in d['rows']:
  for mode,v in r['representations'].items():
   if v['status']!='available':continue
   for geom,c in v['contacts'].items():
    a,b=c['all_context_protein'],c['moving_donor_subset']
    rows.append(dict(zip(fields,[r['pair_id'],r['case_id'],r['protein_id'],mode,geom,v['query_heavy_count'],v['omitted_protein_heavy_count'],v['moving_donor_heavy_count'],a['minimum_A'],a['below_3_5_A'],a['below_2_0_A'],a['new_below_2_0_A'],b['minimum_A'],b['below_3_5_A']])))
 with (out/'CONTACTS.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
 for mode in ('threefold','tenfold'):
  summary[mode]={}
  for geom in ('origin','adaptive_Ca','adaptive_La'):
   vals=[r for r in rows if (r['representation'],r['geometry'])==(mode,geom)]
   summary[mode][geom]={'rows':len(vals),'minimum_all_A':min(r['minimum_all_A'] for r in vals),
    'minimum_moving_A':min(r['minimum_moving_A'] for r in vals),'median_moving_minimum_A':statistics.median(r['minimum_moving_A'] for r in vals),
    'total_contacts_all_below_3_5_A':sum(r['count_all_below_3_5_A'] for r in vals),'total_contacts_moving_below_3_5_A':sum(r['count_moving_below_3_5_A'] for r in vals)}
 value={'result':record(result),'rows':record(out/'CONTACTS.tsv'),'counts':d['summary'],'statistics':summary,
        'pair_denominator':55,'distinct_sources':41,'protein_groups':11,'new_molecular_calls':0,
        'wall_seconds':d['wall_seconds'],'outside_fixed_check':record(Path(result).parent/'OUTSIDE_FIXED_CHECK.json'),
        'interpretation':'correlated descriptive geometry flags only; no correction, rejection or classification change',
        'implementation':record(__file__)}
 write_new(out/'SUMMARY.json',value);return {'summary':record(out/'SUMMARY.json'),'table':record(out/'CONTACTS.tsv')}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--output',required=True)
 print(json.dumps(export(**vars(p.parse_args())),indent=2))
