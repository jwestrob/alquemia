"""Collect the four unchanged native origin cells; no accommodation claim."""
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[2]
W=ROOT/'workspaces/lanm_global_occupancy_20260923'
sys.path.insert(0,str(W/'scoring_v2/implementation'))
from affordable_common import read_json,record,verify,write_new
from compact_solvation import completed,diagnostics
from structure_informed_starts import scf_details

def main():
    mp=W/'native_feasibility_retry_v1/manifest.json';m=read_json(mp);rows=[]
    for t in m['tasks']:
        row={'task_id':t['task_id'],'metal':t['metal'],'medium':t['medium'],
             'status':'unavailable','energy_hartree':None}
        try:
            p=completed(mp,t['task_id'])
            if p is None:raise ValueError('no successful native execution receipt')
            audit=diagnostics(p,t);scf=scf_details(verify(p['output']).read_text())
            if audit['charge_sanity_status']!='pass' or scf['energy']['tolerance']!=1e-10 or not scf['native_mixer_observed']:
                raise ValueError('native state/numerical audit failed')
            row.update(status='complete',**p,audit=audit,scf=scf)
        except Exception as exc:
            row['reason']=str(exc)
            row['retained_output_path']=t['output_path']
        rows.append(row)
    result={'manifest':record(mp),'rows':rows,'complete_cells':sum(r['status']=='complete' for r in rows),
        'cell_denominator':4,'MACE_reused':record(W/'mace_feasibility_v2/result.json'),
        'accommodated_score':None,'accommodation_status':'unavailable_both_nonorigin_proposals_rejected_geometry',
        'new_MACE_calls':0,'new_DFT_calls':0,'continuation_to_eight_other_systems':False}
    write_new(W/'native_feasibility_retry_v1/COLLECTION.json',result)
    print(json.dumps({k:result[k] for k in ('complete_cells','cell_denominator','accommodation_status')}))

if __name__=='__main__':main()
