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

## Next owned tasks (preparation only at this checkpoint)

Curie: inventory/prepare corrected evaluation of the *existing* expanded-context
three-candidate method on its canonical25 and known-reference transfer population.
Use current actual pool/reference/group manifests, existing MACE, and24 compatible
GFNinput reuses. No repeat of2796original-panel points, no newoptimization, no
unknown PLM labels in calibration. Canonical charges differ; derive each state
from exact manifests rather than copying the two-boundary charge.

Preliminary finite inventory:135pools (25canonical,3crystals,6pilot,101additional
transfer); full100transfer uses104source/context pairs from98structures. This is
La-conditioned coverage plus three Ca-conditioned probes, not full two-condition
validation. Maximum1596new primary+1596independent points; first canonical stage
276+276. Integrator must inspect final pinned manifest/dryrun and freeze the same
existing reference rule before launch. This is scope review under already-given
autonomy, not a new user approval gate. Own expanded_transfer_inventory/.

Heisenberg: source-derived q9 donor/denticity/contact diagnosis and exact compatible
archivedDFT inventory, no molecularcalls. Own q9_response_diagnosis/.

Root: scientific comparison/report, next finite execution decision, aggregate
thread allocation and durable wake. Until declared otherwise these follow-ups
use <=2 and <=1 threads, respectively. All previous molecular workers are done.

Exact read-only completed comparison reproduction (choose a fresh output path):

```bash
python scripts/compare_context_factorial.py \
  --mace diagnostics/orca_parameter_audit_20261012/context_MACE/RESULT.json \
  --solvent workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json \
  --output /tmp/nikasha_context_factorial_recheck_20261012
```

## Ownership and safeguards

Jacob cleared the pause and explicitly authorizes<=16 TOTAL numerical login
threads. Recorded verbatim atop AGENTS.md. Runtime4/rootMACE8 allocations are
released after completion. Curie owns<=2 preparation threads; Heisenberg owns<=1
coordinate-analysis thread. Dalton completed independent guard/runtime/factorial
review. New guarded dispatch is committed/pushed asdd758db; its molecular branch
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
activePID669874; latest probes still fail account/partition admission. On recovery
inspect actual live jobs; do not duplicate login work or automatically restart
old LanM campaigns. Historical outputs, existing production and concurrent files
are preserved. Scoped commits/pushes remain authorized.
