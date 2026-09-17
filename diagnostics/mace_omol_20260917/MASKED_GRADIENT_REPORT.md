# Intact-protein masked-score gradients work on one A5000

**Numerical qualification passed; physical response and predictive improvement
remain unvalidated.** The production ORCA baseline, references, masked scalar
score and all biological classifications are unchanged. The new adapter is
`omol_exact_checkpointed_edge_product_autograd_v3`; it supplies analytic
derivatives of the existing masked descriptor, not validated nuclear forces.

## Actual results

The final 73-atom 1H4I core qualification passes23/23 checks. Maximum component
error versus native gradients is9.02e-14 eV-equivalent/A; rotated-gradient error
is at most1.32e-7, below the frozen0.001 tolerance. Native/adapted energies,
readout closure and real signed-displacement checks pass.

The real4698-atom GGR1GLG test passes11/11 checks. Both center energies reproduce
the archived energy-only result exactly. All atoms retain source mappings;
endpoint gradients and grad(R=Ca-La) are exported as arrays and a TSV.
The prescribed metal displacement is0.01 A toward source GLN140/O:

| Quantity | Analytic projected derivative | Actual odd signed energy change | Predicted odd change |
|---|---:|---:|---:|
| Ca | +0.270638768 | +0.002725241 | +0.002706388 |
| La | -8.566307436 | -0.085614377 | -0.085663074 |
| R=Ca-La | +8.836946204 | +0.088339618 | +0.088369462 |

Derivatives use model-kcal/A; changes use model-kcal. These units express the
same model energy scale as the masked score, not measured binding energies.
Maximum odd-change error is0.00004870 versus the frozen0.01 absolute floor.
The disconnected reference constants have zero coordinate derivative.

## Memory recovery and preserved failures

Initial checkpointing passed the core but the whole-GGR forward exhausted the
24GB A5000 at24,969,669,632 peak allocated bytes. Installed PyTorch inspection
showed that `index_add` retained each source edge-message tensor outside the
checkpoint. Replacing that sum with `scatter_add` and an expanded receiver-index
view preserves messages, receivers and edge order while avoiding source retention.
The replacement was requalified on all ten core tasks before the whole retry.
The final whole gradient peaks at12,011,144,704 bytes (11.19GiB).

No custom backward, reduced precision, model retraining or graph detachment was
introduced. CUDA summation order is not promised bit-exact on every device;
actual native/rotated comparisons enforce the declared tolerances. The
[PyTorch2.8 API](https://docs.pytorch.org/docs/2.8/generated/torch.Tensor.scatter_add_.html)
documents the sum and equal source/index shapes required for backward.

Other failures remain visible: exact coordinate/report replay rejected
compute-host roundoff; a TorchScript checkpoint recomputation failed until
checkpoint early stopping was disabled; an initial collector could not serialize
NumPy booleans. Coordinate files remained byte-identical across recovery.
Scientific criteria and decisions did not change. Corrected reports replay
actual receipts; they do not fabricate or rerun successful calculations.
Report-field tolerance1e-10 is limited to derived arithmetic; source artifacts,
criteria and pass/fail decisions remain exact. Actual host discrepancies were
1.42e-14 A for rotated coordinates and4.44e-16 in report reductions.

## Measured cost and tests

| Job | Outcome | Successful / failed model calls | GPU allocation seconds |
|---|---|---:|---:|
|1200863|Coordinate preflight failure|0 / 0|4|
|1200864|Native Ca passes; checkpoint recomputation fails|1 / 1|22|
|1200865|Original core qualification passes|10 / 0|86|
|1200878|Derived-report preflight failure|0 / 0|6|
|1200884|Whole-GGR CUDA OOM|0 / 1|18|
|1200885|Revised core qualification passes|10 / 0|86|
|1200886|Whole-GGR qualification passes|6 / 0|188|

Total development:27 successful calls,2 failed calls,410 GPU allocation seconds,
6560 allocated core-seconds and470.555 reported CPU-seconds. Zero new DFT,
solvent-solver or training calls. The final whole gradient endpoints take
40.699184 and40.789196 model seconds; each signed energy-only displacement takes
about13 seconds. Six final calls total133.702123 model seconds. Preparation,
dry-runs, reports and tests have separate local resource receipts; interactive
edits and unprofiled reads are not included in these compute totals.

Seven final real-fixture regression tests pass (6 in12.883s; source mapping in
1.531s), none skipped. These include actual output serialization, force/gradient
sign and units, finite physical inputs, mandatory core qualification, rejection
of a corrupted criterion and complete whole-protein source/paired-gradient
mapping. Earlier native-regression and recovery tests are retained separately.
These parser/contract tests are distinct from the executed GPU qualifications.

## Meaning and next action

Numerical credibility: passed for the declared core and whole-GGR checks.
Affordability: full descriptor derivatives now fit a standard A5000 with useful
headroom and under a minute of inference per endpoint. Scientific usefulness:
source-mapped sensitivity is available, but no affinity improvement follows yet.
No relaxation, entropy, uncertainty covariance or hybrid gradient is supplied.
The preceding GGR structural and parvalbumin robustness failures remain failures.

Continue with [the declared response screen](MASKED_RESPONSE_SCREEN_PLAN.md):
compare40 masked-core evaluations with archived DFT displacement energies and
analytic gradients, without new DFT. This distinguishes usable direct response
from a need for DFT anchors or an unsuitable curvature model. Do not optimize
geometries or add numerical corrections simply because gradients are available.

Exact outputs and receipts: [compact result](MASKED_GRADIENT_RESULT.json).
Full artifacts live under `workspaces/mace_omol_20260917/`; the authoritative
core report is `masked_gradient_core_report_v5/result.json`, whole report
`masked_gradient_full_report_v1/result.json`. Earlier partial JSON files are
preserved and are not authoritative reports.

## Runnable receipt replay (no inference)

Run from the repository root. The report output directory must be new.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/masked_gradient_full_v3/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/masked_gradient_full_v3/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/masked_gradient_full_v3/implementation/mace_omol_gradient_run.py report --manifest workspaces/mace_omol_20260917/masked_gradient_full_v3/manifest.json --output workspaces/mace_omol_20260917/masked_gradient_review_report_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_omol_gradient_run.py -v
```

The completed full-stage `submission.json` contains exact allocation and execute
arguments. Its frozen implementation also supplies `collect` and `execute`
through `mace_hybrid.py`; successful receipts are reusable and failures remain
visible. Do not duplicate completed inference merely to review this result.
