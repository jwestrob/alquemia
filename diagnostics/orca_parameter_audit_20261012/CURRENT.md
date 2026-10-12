# Corrected Nikasha parameter/context work — 12 October 2026

## Result and current decision

The live Linux native-ORCA parameter correction is reproduced. Sixty actual
controls are terminal:52 admitted, eight original/default-roundtrip density
failures preserved. Corrected La and14 corrected transition-metal aquo cells
match the archived independent controls. Zn/water population parity remains
unqualified. Read runtime/REPORT.md and the independent runtime_review report.

The strict208 corrected calibration remains unsupported (gap−4.556602kcal/mol).
The complete prior2796-point correction was imported, not rerun; see
recovered_strict208/IMPORT.json. Corrected origin-only203/3/2 and frozen geometry
208/208 remain separately scoped results over correlated folds from25proteins.
No threshold tuning or production promotion. Current operating guide and PLM SOP
now identify the earlier native energetic performance claims as historical.

## Completed context experiment

All40 composite endpoints are now available. The independent review agrees with
all energies, selections, interactions and costs. Read context_comparison/REPORT.md
and RESULT.json. Root acknowledged the actual context-factorial-solvent-v2 wake.
No factorial workers remain; do not rerun completed80input pairs.

Two-anchor ordering (La reference c5 minus Ca reference q9) in original/expanded
contexts, kcal/mol:

- Origin: +8.867782 / +12.175484.
- Original three proposals: −4.556602 / −1.212842.
- Expanded three proposals: +0.437421 / +3.472171.
- Common five: −3.714195 / +1.265264.

Expanded context alone does not restore the original proposal pool's ordering.
q9's original La proposal still gives~12.1kcal/mol La-favoring motion response
in both contexts. Context adds a mostly static offset; geometry selection also
matters. A positive two-anchor ordering is not a full canonical calibration.
No production promotion or claim that accommodation is repaired.

MACE:24 reused energies+16 new cross-context+four exact CPU/GPU replays. Actual
20calls used113.884s elapsed/590.696 processCPU-s/4.25GiB peakRSS, zeroGPU.
GFN2:80primary+80independent+four causal diagnostic calls,2490.296 measuredCPU-s.
All40 solvent components pass v2; every original absolute-energy comparison fails
v1 and stays recorded. Density-independent nonsc dispersion cancels exactly
within same-coordinate water−vacuum terms; no fitted offset or changed parameters.
Reports: context_MACE/, context_GFN2/, context_comparison_review/.

## Corrected expanded transfer complete: limited useful gain

Read `expanded_transfer_review/REPORT.md` and RESULT.json (independent raw audit,
commitf397fac). All135 pools and810 solvent components are available; calibration
is exactly unchanged from the canonical collection. Completed1620input pairs
include1596new pairs plus24reuses. Both actual completion wakes were acknowledged.
No canonical/transfer worker or collector remains active; do not rerun them.

| Population | Corrected origin C/W/I/U | Corrected accommodation C/W/I/U |
|---|---|---|
| Canonical calibration25 |25/0/0/0|25/0/0/0|
| Crystal controls3 |3/0/0/0|3/0/0/0|
| All noncalibration pools110 |106/2/2/0|109/1/0/0|
| Primary unique contextpairs104 |103/1/0/0|104/0/0/0|
| Declared triples100 |94/0/0/6|94/0/0/6|

Rows overlap and are not additive. Primary104pairs are98structures from25proteins;
all examples are consumed development evidence. Accommodation fixes two
A0A3F2YLY8 Ca-conditioned inconclusives and one A0ACD6B9F2 La-conditioned wrong
call; its Ca-conditioned sample4 remains wrong. Both methods use their own fixed
canonical-only bands. The expanded accommodated gap is3.472171047kcal/mol;
origin gap10.297678912. Smaller-context strict208 calibration is still invalid.
The separate two-anchor common-five sensitivity is not erased or folded into
this own-three method. No production promotion or unknown PLM inference.

Spread is mixed:51/94 triples narrow,43widen; median range5.104313→5.522710.
Only3completed sources are Ca-conditioned. The six old preparation exclusions
remain. Historical native nonconvergence stays failed; its newly complete
standalone result is a separate corrected-backend outcome.

All810 component checks pass (max2.124440925e-7Eh); max population error
5.411258489e-7e. All1620 absolute backend-energy comparisons still fail the
originalv1 gate, preserved explicitly. Complete raw audit and rows are available.
Actual3192fresh backend calls cost59339.275963childCPU-s; numerical-stage times
744.688672s canonical and3563.979321s transfer on14 one-thread workers. These
exclude preparation/collection and reuse all MACE work; no whole-scanner speedup
is claimed. Peak single-childRSS129020KiB. No GPU or Slurm allocation.

Final collection:
`workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/FINAL_COLLECTION.json`
SHA256 `49bcb0040f31d7f19864591c2392364cb847123ed113095476cdd56c1aa39fab`.
Manifest SHA256 `4f842a2d92364e1559dd5d266ae86e6b682a8ed81181a9e58f0a17cbd59ccfe1`.
Canonical SHA256 `929dedab7ba46c7088de7a4a695a8f19fc5af60066ecd3599645f687bd4e6896`.
Email update to Jacob accepted by local relay; delivery is not independently
verified. Message/receipt are in that workspace as TRANSFER_UPDATE_EMAIL*.

## Next owned task: complete Ca-conditioned structural coverage

The real Ca inventory125 contains110 original supported sources and15 retained
PQQ exclusions. Only3had exact existing expanded pools. Curie's same-policy
preparation of107missing sources admitted106; MMOL2048 Ca sample0 additionally
fails unsupported-disulfide-partner chemistry. Total potential coverage109/125,
with16 explicit exclusions. No repair, new fold or new energetic result implied.
Uses each protein's frozen canonical La-union4.3Å selection, as the old probes
did. Current preparation:
`workspaces/orca_parameter_audit_20261012/ca_conditioned_preparation_v3/`.
First sorted new admitted source A0A3F2YLY8 Ca sample0 is pinned before energies.
Three old-probe replays preserve exact discrete state/maps and coordinates within
existing1e-12Å tolerance (cap roundoff≤8.88e-16Å); retain original old-probe files
for exact cache reuse. Curie owns final report/vault and preparation receipts.

CPU qualification passed: five archived energy/force replays and both original
optimizer starts reproduce GPU energies exactly, force error at most1.47e-13eV/Å,
final coordinate differences at most1.95e-14Å. Same modes, optimizer statistics
and boundary admissions.33 actual forwards cost443.418s elapsed and2766.676
combined processCPU-s at eight active numerical threads; zero failures/GPU.
Both processes are terminal. Read cpu_proposal_feasibility/REPORT.md and
SCOUT_RESULT.json (commit2cc4bf2), independently checked by Gauss.

Curie and Heisenberg are preparing one finite106-source extension: unchanged
four-mode policy, original optimizer parent/kernel, maxiter200, common
origin/Ca/La pool and cross-scoring. Three completed Ca probes will be reused;
all125 original sources stay in the denominator. Frozen expanded canonical bands
are not refitted. No new folds, chemistry, starts, time/call ceiling or production
change. Planned two warm CPU model workers use seven active numerical threads
each, with their one-thread optimizer parents active only while their respective
model waits; aggregate at most14, leaving two headroom. Batch not launched yet.
Root owns final manifest review, launch/completion wake and the separate corrected
GFN2 stage. Agent preparation/code review performs no additional molecular calls.

The q9 structural diagnosis is complete in `q9_response_diagnosis/REPORT.md`.
CN7/denticity stay unchanged; two included hydrogen-bond geometries change in
opposite directions; both La proposals reach0.8Å displacement boundary. Omitted
packing has no new heavy contact within the fixed3.5Å bin. Archived73-atom DFT
origins do not validate the145/192-atom donor motions. Zero new molecular calls.

Exact read-only completed comparison reproduction (choose a fresh output path):

```bash
python scripts/compare_context_factorial.py \
  --mace diagnostics/orca_parameter_audit_20261012/context_MACE/RESULT.json \
  --solvent workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json \
  --output /tmp/nikasha_context_factorial_recheck_20261012
```

## Ownership and safeguards

Jacob cleared the pause and explicitly authorizes<=16 TOTAL numerical login
threads. Recorded verbatim atop AGENTS.md. CPU scout is terminal. Planned finite Ca extension owns at most14 after launch;
all other numerical work must share the remaining two.
Earlier runtime/factorial/reference workers and all independent audits are
complete. Ca-conditioned preparation and CPU qualification are complete. Independent guard/runtime/factorial/expanded-launch reviews
are complete. New guarded dispatch is committed/pushed asdd758db; its molecular branch
remains explicitly unexecuted. Root owns scoring/reporting, current-facing docs and integration.

Parameter guard commitcd97a66 implements exact task-local full-file imports,
shell-aware audits, runtime/profile/state admission and pre/postflight hashes.
Historical original/default controls remain diagnostic-only. The guard's Ca/La
registry uses actual current receipts; it does not globally patch old production
executors or qualify new-metal mixed profiles automatically. The new guarded
dispatch wrapper is separately tested; inspect its report before use.

Two headless Claude attempts produced zero science: first OAuth refresh failed;
one technical retry refused the explicitly authorized login scope. Their receipts
remain under workspaces/claude_parameter_audit_20261012/. Codex completed runtime
work without changing authentication or asking Jacob to repeat authorization.

No Slurm science jobs were submitted. The30minute readiness monitor remains
activePID669874; latest probe122 still fails account/partition admission. On recovery
inspect actual live jobs; do not duplicate login work or automatically restart
old LanM campaigns. Historical outputs, existing production and concurrent files
are preserved. Scoped commits/pushes remain authorized.
