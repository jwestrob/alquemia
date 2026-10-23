#!/usr/bin/env python3
"""Compare captured startup to the already completed exact-role Tight control."""
import difflib, hashlib, json, time
from pathlib import Path

def pin(p):
    p=Path(p).resolve(); return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def main():
    tick,cpu=time.monotonic(),time.process_time()
    base=Path('workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1').resolve()
    out=base/'tighter_dy_acetate_v1/independent_progress_review_v2'
    p=json.loads((out/'snapshots/PANEL.json').read_text())
    controls=json.loads(Path(p['failed_controls']['path']).read_text())
    row=controls['rows']['control__cc__acetate_bi__Dy__ligand_ghost_ion']
    evidence={}
    for key in ('native_result','output','input','start'):
        f=out/('previous_tight_'+Path(row[key]['path']).name)
        assert pin(f)['sha256']==row[key]['sha256']
        evidence[key]=dict(authority=row[key],saved_copy=pin(f))
    old_result=json.loads((out/'previous_tight_RESULT.json').read_text())
    old_start=json.loads((out/'previous_tight_START.json').read_text())
    old_input=(out/'previous_tight_endpoint.inp').read_text()
    new_input=(out/'snapshots/ligand_ghost_ion/endpoint.inp').read_text()
    expected=old_input.replace('%pal nprocs 4 end','%pal nprocs 8 end').replace(' TCutPairs 1e-05\n',' TCutPairs 1e-06\n').replace(' TCutPNO 1e-07\n',' TCutPNO 1e-08\n').replace(' TCutMKN 0.001\n',' TCutMKN 0.0001\n')
    assert new_input==expected
    old_raw=(out/'previous_tight_endpoint.out').read_text()
    new_raw=(out/'snapshots/ligand_ghost_ion/endpoint.out').read_text()
    def shark(t):return t[t.index('SHARK INTEGRAL PACKAGE'):t.index('Calculating RI/C V-Matrix')]
    assert shark(old_raw)==shark(new_raw)
    old_receipt=base/'runtime_recovery_v1/remaining14_submission_v1/completion_event.json'
    accounting=json.loads(old_receipt.read_text())['records']['18382']
    assert accounting['state']=='COMPLETED' and accounting['nodes']=='node-64-768g-14'
    evidence['original_scheduler_receipt']=pin(old_receipt)
    q=dict(status='saved_same_role_startup_comparison_complete',implementation=pin(__file__),evidence=evidence,
        prior_job=18382,current_job=18518,prior_native_wall_seconds=old_result['elapsed_seconds'],
        prior_native_allocated_CPU_seconds=old_result['allocated_cpu_seconds'],prior_scheduler=accounting,
        prior_native_startup_seconds=2.594,prior_RI_C_metric_cholesky_seconds_printed=.1,
        native_timing_evidence=[dict(line=i+1,text=line) for i,line in enumerate(old_raw.splitlines()) if 'Calculating RI/C V-Matrix' in line or 'Startup calculation ' in line or 'TOTAL RUN TIME:' in line],
        exact_shark_startup_settings_block=True,shark_startup_block=shark(old_raw),
        input_only_changes_verified=True,input_diff=''.join(difflib.unified_diff(old_input.splitlines(True),new_input.splitlines(True),fromfile='saved Tight control18382',tofile='captured tighter18518')),
        prior_resources={k:old_start[k] for k in ('ranks','memory_MiB','maxcore_MiB','host')},
        current_resources=dict(ranks=8,memory_MiB=32000,maxcore_MiB=3000,nodes=['node-64-768g-5','node-64-768g-6'],observed_ranks_per_node=[4,4]),
        matched_dimensions=dict(atoms=8,orbital_basis_functions=189,orbital_shells=71,AuxC_functions=718,AuxC_shells=176),
        interpretation='The declared tighter cutoffs change downstream local-correlation work, not the global pre-SCF Coulomb fitting metric. Identical printed startup settings and a prior 0.1-second metric operation do not explain the observed >113-minute silence. CPU-active MPI startup is anomalous; a busy communication wait or runtime problem is plausible, not proven.',
        recommendation='A root-declared single-node resource-only recovery is justified after preserving and disposing of the stalled attempt. Retain eight MPI ranks, default grant-derived MaxCore, exact tighter chemistry/input, all stopping and admission gates; a one-node constraint isolates placement without reverting scientific cutoffs. No cancellation or submission performed here.',
        limits=['Rank count and placement both differ from the old four-rank run; these observations alone do not isolate their individual contributions.','No explicit native failure or stack trace establishes a particular MPI or transport defect.','No current live files or scheduler/process status were queried in this comparison.','The three-role CP cycle and original 0.25 kcal/mol test remain unavailable.'],
        wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1,new_native_calls=0,new_live_queries=0)
    dest=out/'STARTUP_COMPARISON.json'
    with dest.open('x') as f:json.dump(q,f,indent=2,allow_nan=False);f.write('\n')
    note='''# Exact-role startup comparison

The completed Tight control **18382** used four MPI ranks on one node (`node-64-768g-14`). Native total time was **1330.206 s**, the whole startup **2.594 s**, and the RI/C V-Matrix plus Cholesky operation printed **0.1 s**.

The captured current input differs only in four-to-eight MPI ranks and the three declared downstream local-correlation cutoffs. The entire printed SHARK startup-settings block is byte-identical: 8 centers, 189 orbital functions/71 shells, 718 auxiliary-C functions/176 shells, and 3000 MiB MaxCore. Geometry, basis assets, SCF settings, ghost treatment and effective state are unchanged. Current placement is four ranks on each of two allocated nodes.

The known workload does **not** explain more than 113 minutes at this pre-SCF operation. The active CPU snapshot could represent a busy MPI wait; it does not prove useful progress or identify a particular transport defect. Eight versus four ranks also differs, so placement is not uniquely established as the cause.

A separately preserved **single-node, eight-MPI resource-only recovery is justified** if root closes the anomalous attempt. Preserve the exact tighter input, default grant-derived MaxCore, original admission and CP threshold. This is a prospective recommendation only: no cancellation, new job, live polling or method change occurred. The incomplete CP cycle remains unavailable.
'''
    for folder in (out,Path(__file__).resolve().parent):
        with (folder/'STARTUP_COMPARISON.md').open('x') as f:f.write(note)
    print(json.dumps(dict(comparison=pin(dest),report=pin(Path(__file__).resolve().parent/'STARTUP_COMPARISON.md'))))

if __name__=='__main__':main()
