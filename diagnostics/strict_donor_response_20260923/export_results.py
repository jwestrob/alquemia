"""Write compact tables from the completed, unmodified strict donor collection."""
import csv,datetime,json,statistics,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from affordable_common import HA_TO_KCAL,read_json,record,verify,write_new
p=Path(sys.argv[1]);out=Path(sys.argv[2]);a=read_json(p);m=read_json(verify(a['manifest']));old={t['task_id']:t for t in m['tasks']+m['reused']};cells=[];times=[]
for r in a['rows']:
    src=old[r['task_id']]['source'];d=r.get('details',{});receipt=read_json(verify(r['actual']['receipt'])) if r.get('actual') else {}
    wall=(datetime.datetime.fromisoformat(receipt['finished_at_utc'])-datetime.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds() if receipt else None
    if not r['reused'] and wall is not None:times.append(wall)
    cells.append(dict(task_id=r['task_id'],status=r['status'],reused=r['reused'],old_hartree=src['source_energy_hartree'],strict_hartree=r['energy_hartree'],
        change_kcal_mol=(r['energy_hartree']-src['source_energy_hartree'])*HA_TO_KCAL if r['status']=='complete' else None,
        SCF_cycles=d.get('cycles'),wall_seconds=wall,printed_max_density_within=d.get('max_density',{}).get('within_printed_tolerance'),printed_rms_density_within=d.get('rms_density',{}).get('within_printed_tolerance'),
        input_sha=receipt.get('artifacts',{}).get('template_input',{}).get('sha256'),output_sha=r.get('actual',{}).get('output',{}).get('sha256')))
diff=[{k:v for k,v in r.items() if k!='old'}|{'old_solvent':r['old']['solvent_delta_R'],'old_composite':r['old']['composite_delta_R'],'old_composite_error':r['old']['composite_error']} for r in a['differentials']]
for name,rows in [('ENDPOINT_WORKS.csv',a['endpoint_works']),('DIFFERENTIALS.csv',diff),('CELL_SHIFTS.csv',cells)]:
    with (out/name).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
execinfo=read_json(p.parent/'EXECUTION.json')
raw=subprocess.check_output(['sacct','-j',execinfo['slurm_job_id'],'-P','-n','--format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,MaxRSS,NodeList'],text=True)
main=next(line.split('|') for line in raw.splitlines() if line.split('|')[0]==execinfo['slurm_job_id'])
cost={'job_id':execinfo['slurm_job_id'],'sacct_raw':raw,'allocated_wall_seconds':int(main[2]),'allocated_CPUs':int(main[3]),'allocated_core_seconds':int(main[4]),'GPU_seconds':0,
      'new_calls':len(times),'molecular_task_wall_sum_seconds':sum(times),'molecular_task_wall_median_seconds':statistics.median(times),'molecular_task_wall_max_seconds':max(times),
      'engine_wall_seconds':execinfo['wall_seconds'],'reused_cells':4,'execution':record(p.parent/'EXECUTION.json'),
      'cost_scope':'new allocation and molecular calls; old reused calls and local preparation/report CPU separate'}
write_new(p.parent/'COSTS.json',cost)
metrics={'collection':record(p),'complete_cells':a['complete_cells'],'endpoint_works':len(a['endpoint_works']),'differentials':len(diff),
         'max_cell_change_kcal_mol':max(abs(c['change_kcal_mol']) for c in cells),
         'max_differential_change_kcal_mol':max(abs(r['strict_composite']-r['old_composite']) for r in diff),
         'max_endpoint_absolute_error':{k:max(abs(r[k]-r['DFT']) for r in a['endpoint_works']) for k in ('MACE','old_composite','strict_composite')},
         'max_differential_absolute_error':{k:max(abs(r[k]-r['DFT']) for r in diff) for k in ('MACE','old_composite','strict_composite')},
         'endpoint_signs_preserved':sum(r['strict_composite']*r['DFT']>0 for r in a['endpoint_works']),
         'differential_signs_preserved':sum(r['strict_composite']*r['DFT']>0 for r in diff),
         'SCF_cycle_min_max':[min(c['SCF_cycles'] for c in cells),max(c['SCF_cycles'] for c in cells)],
         'printed_density_exceptions':[c['task_id'] for c in cells if c['printed_max_density_within'] is False or c['printed_rms_density_within'] is False],
         'response_accuracy_pass_threshold':None,'costs':record(p.parent/'COSTS.json')}
write_new(out/'SUMMARY.json',metrics);print(json.dumps(metrics,indent=2));print(json.dumps(cost,indent=2))
