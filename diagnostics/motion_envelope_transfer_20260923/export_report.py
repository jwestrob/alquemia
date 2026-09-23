"""Compact report from the completed fixed full100 comparison and actual costs."""
from pathlib import Path
import csv,json,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new
B=ROOT/'workspaces/motion_envelope_transfer_20260923';D=Path(__file__).parent
x=read_json(B/'COMPARISON_v1.json');cost=read_json(B/'costs_v1/COSTS.json')
searches=[];timing=[]
for i in (0,1):
 s=B/f'run_v2/shard_{i}/searches';summary=read_json(next(s.glob('gpu_summary*.json')));executor=read_json(next(s.glob('EXECUTION*.json')))
 rows=[read_json(p)for p in (s/'proposals').glob('*/result.json')];searches.extend(rows)
 timing.append({'shard':i,'searches':len(rows),'actual_MACE_calls':summary['new_MACE_calls'],'worker_wall_seconds':summary['wall_seconds'],
  'model_load_seconds':summary['model_load_seconds'],'executor_wall_seconds':executor['wall_seconds'],
  'sum_per_search_wall_seconds':sum(r['wall_seconds']for r in rows)})
noncorrect=[{'case_id':r['case_id'],'selection_id':r['selection_id'],'expected_class':r['expected_class'],
 'envelope':r['methods']['envelope_operational'],'strict_tenfold':r['methods']['tenfold_strict']}for r in x['rows']if r['methods']['envelope_operational']['outcome']!='correct']
repairs=[{'triple_id':t['triple_id'],'protein_id':t['protein_id'],'old':t['methods']['threefold3p5'],'new':t['methods']['envelope_operational']}
 for t in x['triples']if t['methods']['threefold3p5']['outcome']=='inconclusive'and t['methods']['envelope_operational']['outcome']=='correct']
result={'comparison':record(B/'COMPARISON_v1.json'),'reference':x['reference'],'counts':x['counts'],'matched':x['matched_triples'],
 'denominators':x['denominators'],'unique_physical_sources':98,'counts_are_correlated_developmental_samples':True,'noncorrect_pairs':noncorrect,
 'old_threefold_abstention_repairs':repairs,'actual_new_calls':{'MACE_origins':202,'MACE_searches':sum(t['actual_MACE_calls']for t in timing),
 'MACE_cross':202,'MACE_total':404+sum(t['actual_MACE_calls']for t in timing),'GFN2_origin_attempts':404,'GFN2_candidate_attempts':808,
 'GFN2_attempts':1212,'GFN2_accepted':1211,'DFT':0,'folds':0},'searches':{'attempted':202,'admissible':sum(r['status']=='proposal_available'for r in searches),
 'boundary_flags':sum(bool(r['boundary_flag'])for r in searches),'infeasible_evaluated_MACE_trials':sum(r['infeasible_completed_MACE_requests']for r in searches)},
 'timing':timing,'costs':record(B/'costs_v1/COSTS.json'),'allocated_CPU_seconds':cost['allocated_CPU_seconds'],'requested_GPU_seconds':cost['requested_GPU_seconds'],
 'failure':record(D/'FAILURE.json'),'pilot_regression':'A0A3 Ca-conditioned sample3 remains inconclusive; separate from this La-only transfer',
 'new_calibration':False,'production_changed':False,'implementation':record(__file__)}
write_new(D/'RESULT.json',result)
for name,rows,key in [('TRIPLES.csv',x['triples'],'triple_id'),('PAIRS.csv',x['rows'],'case_id')]:
 with (D/name).open('x',newline='')as f:
  fields=[key,'protein_id','selection_id','expected_class','method','R_model_kcal_mol','decision','outcome'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for r in rows:
   for method,c in r['methods'].items():w.writerow({key:r[key],'protein_id':r.get('protein_id',r.get('root_case_id')),'selection_id':r['selection_id'],
    'expected_class':r['expected_class'],'method':method,'R_model_kcal_mol':c['R'],'decision':c['decision'],'outcome':c['outcome']})
text='''# Motion-envelope full100 transfer: fidelity retained where complete, coverage reduced

**All 91 complete three-source summaries are correct.** The new envelope recovers
both earlier A8 threefold abstentions and matches released/strict-tenfold decisions
on the same 91 triples. One real SCF failure makes three additional triples
unavailable: the primary result is **91 correct / 0 wrong / 0 inconclusive / 9
unavailable**, versus 94 correct / 6 unavailable for strict tenfold scoring.
This is a useful practical three-source result with a numerical coverage limitation;
it does not establish equivalent coverage or justify default promotion.

## Frozen comparisons

| Method, each with its own canonical-only reference | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Released static | 94 | 0 | 0 | 6 |
| Tenfold strict adaptive | 94 | 0 | 0 | 6 |
| Original 3.5 Å threefold adaptive | 92 | 0 | 2 | 6 |
| 4.3 Å envelope origin only | 92 | 0 | 2 | 6 |
| 4.3 Å envelope adaptive | 91 | 0 | 0 | 9 |

On exactly the common 91 complete triples, envelope adaptive and strict tenfold
are both 91 correct; envelope origins and old threefold each give 89 correct and
two inconclusive. Search repairs origin abstentions in A0A3 triples000/002. Both old
A8 abstentions (008/010) are correct. The old3.5 method used the historical looser
scalar policy, so its difference is not a pure context-radius effect. The strict
tenfold comparator uses the same current scalar profile. Mathematical and
operational envelope variants have the same classification counts.

The reference was frozen from the original25 canonical geometries before these
transfer energies. No new thresholds, member selection, label fitting or production
change occurred. All100 original memberships remain. These are correlated subsets
of25 consumed proteins, not100 independent biological validations or affinity tests.

## Individual sources and coverage

There are104 distinct source/context pairs from98 physical source structures and
300 declared triple-member occurrences. These denominators are reported separately.
Unique pairs give99 correct / 4 inconclusive / 1 unavailable, versus103 correct /
1 inconclusive for strict tenfold evaluated on those same104 occurrences. The
inconclusive pairs are C5AX La-sample3 (already inconclusive), Q60AR6 La-sample3 under
two different envelope selections, and Q88JH5 La-sample3. Their complete triple
medians are correct without removing or replacing members. Three newly inconclusive
pairs are a limitation even though aggregation resolves their calls.

The six old preparation exclusions remain: triples037–039 (MMOL1770) and
076/077/079 (Q88JH5). The new numerical failure invalidates A0AC triples004–006.
The separate Ca-conditioned A0A3 sample3 pilot regression remains inconclusive;
a favorable La-conditioned transfer does not erase it.

## Actual numerical failure

La evaluated at the Ca-proposed geometry of A0ACD6B9F2 La-sample0 vacuum fails
native SCF after500 cycles. Its last energy changes remain around10^-4–10^-3
Hartree, well above TolE1e-10. The native process returned0 but explicitly reported
error termination; its receipt correctly records no normal termination or SCF
convergence. No energy is accepted and the required whole pool remains unavailable.
The cell took201.459371s. It used one rank, maxcore2000MB; the outer allocation had
64GiB/32CPU and batch peakRSS approximately6.37GiB, which is not a per-cell memory
measurement. There is no OOM evidence. See FAILURE.json for exact pins and accounting.
Any separately owned restart diagnostic is excluded from this primary result.

## Execution, costs and integrity

Three exact completed pilot pools were reused;101 new pools were executed in51/50
source shards. All202 searches yielded admissible proposals. There were3,408 new
MACE calls (202 origins,3,004 search evaluations,202 cross evaluations) and1,212
strict nativeGFN2 attempts (1,211 accepted,one failed). No DFT, new folds, states,
waters or protonation changes occurred. Both metals score the same candidate pool.
97 final proposals touch a bound, with no unconstrained-minimum claim. The unchanged
SLSQP engine records641 evaluated trials outside final admissibility; none is an
admitted final proposal. The motion-envelope argument applies to admitted original
anchor displacement against omitted fixed atoms, not every optimizer trial.

All21 outer allocations, including empty/failed technical stages, cost82,389 allocated
CPU-seconds and1,571 requested GPU-seconds. GPU time is requested allocation time,
not measured utilization; this cluster omits GPU from AllocTRES, so pinned wrapper
requests supply that count. Local preparation and read-only collection/recovery time
are additional and not included as metered allocation cost. No nested runner costs
are added twice. Search executors took503.780/466.984s for102/100 searches, with
worker times497.090/459.265s and model loads1.174/1.185s. This is prepared-batch timing,
not folding-inclusive or isolated three-source latency.

Two metadata defects were recovered without repeating molecules. First, exact-byte
runner/renderer copies at different paths caused receipt rejection; recovery retains
all actual paths and checks the original receipts. Second, a missing scalar cell join
field was reconstructed from exact case/metal/candidate identity plus physical state.
Original inputs, receipts and rejected collections remain intact. Separate JOIN files
record the recovery. Ten failed allocation statuses comprise nine technical/reporting
failures plus the one allocation containing the real SCF failure.

Nine focused real-artifact tests pass: five scope/recovery checks and four final
matrix/median/denominator/old-context-join checks. Raw scores, components, works,
spreads and every member occurrence remain in the full comparison. The original3.5
context is joined separately for each triple where larger envelopes merge older
selections. No favorable old context is chosen.

[Compact result](RESULT.json) · [triples](TRIPLES.csv) · [pairs](PAIRS.csv) ·
[commands](COMMANDS.md) · [failure](FAILURE.json).
Full comparison: `workspaces/motion_envelope_transfer_20260923/COMPARISON_v1.json`.
'''
(D/'REPORT.md').write_text(text)
print(json.dumps({'result':record(D/'RESULT.json'),'report':record(D/'REPORT.md')},indent=2))
