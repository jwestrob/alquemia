"""Fresh execution-only memory recovery of a pinned four-cell origin manifest."""
import argparse
from copy import deepcopy
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
W=ROOT/'workspaces/lanm_global_occupancy_20260923'
sys.path.insert(0,str(W/'scoring_v2/implementation'))
from affordable_common import cache_key,read_json,record,verify,write_new
from affordable_workflow import dry_run
from lanm_global_occupancy import PROTOCOL
from affordable_common import xyz

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source-manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--maxcore-mb',type=int,required=True)
    p.add_argument('--ranks',type=int,required=True)
    p.add_argument('--memory-mib',type=int,required=True)
    a=p.parse_args()
    if a.maxcore_mb<6000:raise ValueError('below observed whole-source SCF memory requirement')
    if a.ranks<1 or 4*a.ranks*a.maxcore_mb>a.memory_mib*1.048576*.76:
        raise ValueError('rank memory exceeds declared allocation with headroom')
    m=read_json(a.source_manifest);parent=read_json(verify(m['parent_manifest']))
    assert m['protocol_id']==PROTOCOL
    assert len(m['tasks'])==4 and not m['reused']
    assert {(t['state_id'],t['candidate'],t['metal'],t['medium']) for t in m['tasks']}=={
        ('Hans_8DQ2__EF12','origin',z,s) for z in ('La','Dy') for s in ('vacuum','alpb')}
    states={read_json(verify(pin))['state_id']:pin for pin in parent['states']}
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    result=deepcopy(m)
    for t in result['tasks']:
        old=verify(t['input']).read_text()
        assert old.count('%maxcore 2000')==1
        body=old.replace('%maxcore 2000',f'%maxcore {a.maxcore_mb}')
        assert body.replace(f'%maxcore {a.maxcore_mb}','%maxcore 2000')==old
        coords=verify(t['xyz']).read_bytes()
        d=out/'tasks'/t['task_id'];d.mkdir(parents=True)
        (d/'core.xyz').write_bytes(coords);(d/'endpoint.inp').write_text(body)
        t['recovery_from_scientific_key']=t['scientific_key']
        t.update(input=record(d/'endpoint.inp'),xyz=record(d/'core.xyz'),output_path=str(d/'endpoint.out'))
        t['scientific_key']=cache_key({'protocol':PROTOCOL,'state':states[t['state_id']],
            'atoms':xyz(d/'core.xyz'),'input_text':body,'ranks':a.ranks})
    result['execution_resources'].update(mpi_ranks=a.ranks,concurrent_tasks=4,
        maxcore_mb_per_rank=a.maxcore_mb,scheduler_memory_mib=a.memory_mib)
    result['technical_recovery']={'source_manifest':record(a.source_manifest),
        'change':'Allocation-derived MPI ranks and MaxCore; four concurrent unchanged molecular tasks',
        'previous_mpi_ranks':m['execution_resources']['mpi_ranks'],'new_mpi_ranks':a.ranks,
        'old_maxcore_mb':2000,'new_maxcore_mb':a.maxcore_mb,
        'source_coordinates_unchanged':True,'Hamiltonian_and_SCF_tolerances_unchanged':True,
        'preparation_script':record(Path(__file__)), 'new_MACE_calls':0,'new_DFT_calls':0}
    write_new(out/'manifest.json',result)
    r=dry_run(out/'manifest.json');write_new(out/'PREFLIGHT.json',r);print(r)

if __name__=='__main__':main()
