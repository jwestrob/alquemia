"""Export compact records from the completed immutable34-source pilot."""
from pathlib import Path
from collections import Counter
import csv,json,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
B=ROOT/'workspaces/motion_envelope_20260923';D=ROOT/'diagnostics/motion_envelope_20260923'

def main():
    ref=read_json(B/'REFERENCE_PINNED_v1.json');d=read_json(B/'COMPARISON_v1.json');cost=read_json(D/'COSTS_v1.json');inv=read_json(B/'INVENTORY_v1.json');pool=read_json(B/'pilot34_pool_v1/COLLECTION_FINAL.json')
    proposals=[read_json(p)for p in (B/'pilot34_searches_v1/proposals').glob('*/result.json')]
    probes=[]
    for r in d['rows']:
        if r['role']!='separate_probe':continue
        bands=ref['variants']['operational']['bands'];R=r['scores']['operational'];z=r['expected_class'];v=r['variants']['operational']
        probes.append({'source_id':r['source_case_id'],'expected_class':z,'R_model_kcal_mol':R,'margin_to_own_class_edge_kcal_mol':R-bands['La_min']if z=='La'else bands['Ca_max']-R,
            'adaptive_outcome':r['own_reference']['operational']['outcome'],'origin_outcome':r['origin_own_reference']['outcome'],
            'strict_tenfold_outcome':v['strict_tenfold_own_reference']['outcome'],'final_shift_vs_strict_tenfold':v['pool_shift']['composite_R_model_kcal_mol'],
            'origin_shift_vs_strict_tenfold':v['origin_shift']['composite_R_model_kcal_mol'],'differential_accommodation_shift':v['differential_accommodation_shift']['composite_R_model_kcal_mol']})
    ca3=next(r for r in d['rows']if r['source_case_id']=='a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-3');c=next(c for c in pool['cases']if c['case_id']==ca3['case_id'])
    bound={z:{k:read_json(verify(c['aliases']['adaptive_'+z]['proposal_receipt']))[k]for k in ('boundary_flag','final_geometry','optimizer')}for z in ('Ca','La')}
    result={'protocol_id':d['protocol_id'],'reference':record(B/'REFERENCE_PINNED_v1.json'),'comparison':record(B/'COMPARISON_v1.json'),'collection':record(B/'pilot34_pool_v1/COLLECTION_FINAL.json'),
        'preparation':inv['preparation'],'inventory':record(B/'INVENTORY_v1.json'),'costs':record(D/'COSTS_v1.json'),'counts':d['counts'],
        'adaptive_reference':{k:v for k,v in ref['variants']['operational'].items()if k!='rows'},'origin_only_reference':{k:v for k,v in ref['origin_only_reference'].items()if k!='rows'},
        'probes':probes,'pair_spreads':d['prescribed_probe_pair_spreads'],'Ca3':{'row':{k:ca3[k]for k in ('case_id','source_case_id','same_original_core_max_A','preparation_changes','origin_components','variants')},'new_endpoint_geometry':bound},
        'calls':{'MACE_origins':68,'MACE_search_evaluations':998,'MACE_cross':68,'MACE_total':1134,'GFN2_origins':136,'GFN2_candidates':272,'GFN2_total':408,'DFT':0,'searches':68,'valid_proposals':68,'boundary_proposals':25,
            'evaluated_trials_outside_final_admissibility':sum(r['infeasible_completed_MACE_requests']for r in proposals),'maximum_trial_heavy_extent_A':max(r['maximum_trial_heavy_extent_A']for r in proposals)},
        'qualification':'canonical_and_crystal_fidelity_retained_with_one_material_Ca_conditioned_probe_abstention','full100_molecular_transfer_launched':False,'production_changed':False,'new_molecular_calls_in_export':0,
        'implementation':record(__file__)}
    write_new(D/'RESULT.json',result)
    with (D/'PROBES.csv').open('x')as f:
        w=csv.DictWriter(f,fieldnames=list(probes[0]));w.writeheader();w.writerows(probes)
    lines=['# Three-source motion envelope: mixed pilot result','',
        '**The4.3 Å envelope retains all25 canonical and three crystal calls and repairs both earlier A8 probe abstentions, but introduces one A0A3 Ca-sample3 abstention.** Its six probes are5correct/0wrong/1inconclusive, versus6/0/0 for matched strict ten-fold scoring. This supports wider qualification of the practical three-source preparation; it does not establish equivalent fidelity or justify promotion.','',
        'All34 pools are available. Expected classes and source membership are unchanged. These are already consumed development references and fold probes, not independent proteins or direct affinity validation.','',
        '## Separate preparation and search results','',
        '| Method, with its own canonical-only reference | Canonical | Crystals | Six probes |',
        '|---|---:|---:|---|','| Strict ten-fold origins |25/25|3/3|3correct,1wrong,2inconclusive|',
        '|4.3 Å envelope origins|25/25|3/3|2correct,1wrong,3inconclusive|',
        '| Strict ten-fold adaptive pool |25/25|3/3|6correct|',
        '|4.3 Å envelope adaptive pool|25/25|3/3|5correct,1inconclusive|','',
        'The envelope search repairs three of its own origin probe calls. Preparation alone does not solve this panel. Its adaptive class gap is7.071396249113604 modelkcal/mol; the separate origin-only gap is10.117019388068002. A larger origin gap is not an accuracy gain. Mathematical and operational adaptive choices give the same calls here. Through the unchanged strict ten-fold adaptive bands, envelope scores give31correct/3inconclusive, including two canonical abstentions: the new representation requires its own reference.','',
        'Only the25 designated canonical geometries fit either reference, using the existing extrema/minimum-gap rule. Crystals and probes are excluded. The separately authorized origin ablation was added after origin completion and remains distinct from the adaptive calibration.','',
        '## Prespecified probes','',
        '| Source | Expected | Envelope adaptive | Margin to own class edge, kcal/mol | Origin call |',
        '|---|---|---|---:|---|']
    for r in probes:lines.append(f"|{r['source_id'].replace('-pqq-la_model__conditioned_',' ').replace('__seed-1_sample-',' sample')}|{r['expected_class']}|{r['adaptive_outcome']}|{r['margin_to_own_class_edge_kcal_mol']:.6f}|{r['origin_outcome']}|")
    lines +=['','## Why A0A3 Ca-sample3 becomes inconclusive','',
        'The matched original core coordinates agree exactly, and both contexts retain Ca/La charges−3/−2. The envelope has202 atoms versus193. Its origin contrast shifts−2.960946848kcal/mol: nativeMACE−5.287715162, solvent subtraction+2.326768314. Differential accommodation improves+1.399701620, leaving a final−1.561245228 shift from the strict ten-fold result. The resulting score is0.569342826 below its new La edge. This is mainly a starting-context offset, partly repaired by search; it is not an excess accommodation penalty. Composition and continuum cavity change together.','',
        'The unchanged four-mode selector keeps197χ1,197χ3 and329χ2, but replaces197χ2 with327χ1. Neither endpoint is boundary limited: new maximum heavy displacement is0.741058Å(Ca)/0.638445Å(La), versus0.765690/0.674366Å before. No Ca3-specific optimization, radius change, threshold adjustment or member replacement was performed.','',
        'Probe-pair spreads are mixed. A0A3 changes4.752626→4.874830kcal/mol; A0AC changes15.032615→13.996623; A8 changes0.421140→0.572308. This is not a uniform structural-robustness improvement.','',
        '## Actual scope, integrity and cost','',
        'Preparation retains94complete/6unavailable declared triples and prepares135/135 unique source-context pairs. Only the frozen34-source pilot was scored. All68 searches succeed;25 final proposals touch a bound, with no unconstrained-minimum claim. The existing SLSQP kernel also records175 evaluated trial points outside final admissibility; none is selected as a final proposal. The4.3 Å reachability argument applies to admitted anchor movements against omitted fixed atoms, not every optimizer trial or continuum-solvent convergence.','',
        'Actual new work:68 MACE origins +998 search evaluations +68 cross evaluations =1,134 MACE calls;136 origin +272 candidate nativeGFN2 cells =408, all accepted. Scalar inputs use the separately qualified nativeTolE1e−10,300K,MaxIter500,freshNoAutostart,one-rank profile. No loose scalar result is reused. No DFT, folding, new protonation/waters or production change.','',
        'The five outer allocations used67,040 CPU-seconds and1,250 requestedGPU-seconds. Their elapsed times are20s(originMACE),1,178s(search),629s(originGFN2),52s(crossMACE),216s(candidateGFN2). They overlap; their sum is not end-to-end latency.','',
        'Search executor wall is1,090.311s; summed per-search wall862.759s; summed actual MACE-evaluation wall246.640s; one-time model-load time1.556s. The executor minus recorded GPU-worker interval is210.446s, retained as unattributed overhead. Model construction occurs once per batch; remaining initialization/collection/setup work is not assigned an unmeasured cause or assumed fully reusable. The origin/candidate scalar runners report440.028/170.007s inside their allocations. Local preparation took73.108s; local preflight/report work is additional. This34-source batch does not measure isolated three-source latency.','',
        '## Reproducibility and limitations','',
        'The reference writer was archived byte-for-byte before correcting a report-only component-key alias. REFERENCE_v1 remains intact; REFERENCE_PINNED_v1 preserves every scientific field and its frozen UTC while pointing to that exact archived implementation. Singleton crystals explicitly have no declared triple ID. An earlier local terminal-file watcher exited before the completed collection became visible; direct preparation used the final artifact. Neither reporting repair nor watcher retry made molecular calls.','',
        'All18 focused actual-fixture tests pass across preparation, execution preflight, source joining, pool integrity and final algebra. Raw matrices, gradients, traces, scalar diagnostics and failed preliminary report logs remain in the workspace. The full100 transfer has not been launched; a separate finite inventory is the next step. Keep the Ca-conditioned probe regression visible in any later La-conditioned transfer result.','',
        'Compact artifacts: [RESULT](RESULT.json), [probe table](PROBES.csv), [costs](COSTS_v1.json), [commands](COMMANDS.md). Full numerical artifacts are under `workspaces/motion_envelope_20260923/`.']
    (D/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'result':record(D/'RESULT.json'),'report':record(D/'REPORT.md'),'counts':d['counts']['all34']},indent=2))
if __name__=='__main__':main()
