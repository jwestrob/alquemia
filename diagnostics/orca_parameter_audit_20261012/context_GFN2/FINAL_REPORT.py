"""Summarize actual finite factorial receipts; no scientific executable invoked."""
from pathlib import Path
import csv, json, sys, hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import COMPONENT_V2 as driver

def main():
    base=driver.BASE
    result=driver.f.read(base/'FINAL_COLLECTION.json')
    manifest=driver.f.read(base/'manifest.json')
    raw=driver.raw_results(manifest)
    out=Path(__file__).resolve().parent
    rows=result['rows']
    fields=['id','protein','context','geometry','metal','charge','uhf','status','accepted',
            'primary_vacuum_energy_hartree','primary_water_energy_hartree','primary_solvation_hartree',
            'independent_solvation_hartree','component_error_hartree','max_population_error_e',
            'original_absolute_energy_gate_pass','failure_reason','xyz_sha256']
    with (out/'COMPONENTS.tsv').open('x') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:
            d={k:r.get(k) for k in fields};d['xyz_sha256']=r['xyz']['sha256'];w.writerow(d)
    costs=[];inputrows=[]
    for task,r,pin in raw.values():
        assert r is not None
        ckey='__'.join((task['protein'],task['context'],task['geometry'],task['metal']))
        comp=next(v for v in rows if v['id']==ckey)
        row={k:task[k] for k in ('id','protein','context','geometry','metal','medium','charge','uhf')}
        row.update(xyz_sha256=task['xyz']['sha256'],component_status=comp['status'],
                   original_absolute_pass=r['accepted'],result_path=pin['path'],reused_scout=task['scout'],
                   archive_check_available=task['reference'] is not None,
                   archive_error_Eh=r.get('reference_energy_difference_hartree'))
        for b in ('xtb','tblite'):
            br=r['backends'][b];row[b+'_status']=br['status'];row[b+'_energy_hartree']=br['energy_hartree']
            row[b+'_charge_sum_residual_e']=br.get('charge_sum_residual_e')
            rec=driver.f.read(br['execution']);rec.update(backend=b,id=task['id'],kind='factorial',reused=task['scout']);costs.append(rec)
        if r['parity']:
            row['absolute_energy_error_Eh']=r['parity']['energy_difference_hartree']
            row['maximum_population_error_e']=r['parity']['maximum_population_difference_e']
        else:row.update(absolute_energy_error_Eh=None,maximum_population_error_e=None)
        inputrows.append(row)
    with (out/'INPUTS.tsv').open('x') as h:
        w=csv.DictWriter(h,list(inputrows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(inputrows)
    diag=ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_diagnosis_v1'
    for path in sorted(diag.glob('points/*/*/EXECUTION.json')):
        r=driver.f.read(path);r.update(backend=path.parent.name,kind='dispersion_diagnostic',reused=False);costs.append(r)
    summaries={}
    for scope in ['factorial','dispersion_diagnostic','all']:
        recs=[r for r in costs if scope=='all' or r['kind']==scope]
        summaries[scope]={'process_calls':len(recs),'worker_wall_seconds_sum':sum(r['wall_seconds'] for r in recs),
            'measured_process_cpu_seconds':sum(r['cpu_user_seconds']+r['cpu_system_seconds'] for r in recs),
            'peak_single_process_RSS_KiB':max(r['peak_RSS_KiB'] for r in recs)}
    errors=[abs(r['component_error_hartree']) for r in rows if r['component_error_hartree'] is not None]
    m={'protocol_id':result['protocol_id'],'source_collection':driver.f.pin(base/'FINAL_COLLECTION.json'),
       'input_count':len(inputrows),'qualified_components':sum(r['accepted'] for r in rows),'declared_components':len(rows),
       'maximum_component_error_Eh':max(errors) if errors else None,
       'maximum_population_error_e':max(r['maximum_population_error_e'] for r in inputrows if r['maximum_population_error_e'] is not None),
       'absolute_gate_failed_inputs':sum(not r['original_absolute_pass'] for r in inputrows),
       'archive_checked_inputs':sum(r['archive_check_available'] for r in inputrows),
       'costs':summaries,'per_process_receipts':costs,
       'artifact_pins':[driver.f.pin(out/'COMPONENTS.tsv'),driver.f.pin(out/'INPUTS.tsv')]}
    driver.f.write(out/'RESULT.json',m)
    print(json.dumps({k:v for k,v in m.items() if k not in ('per_process_receipts','artifact_pins')},indent=2))

if __name__=='__main__':main()
