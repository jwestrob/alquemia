# One failed strict native cell is recoverable from two compatible seeds

Both predefined native restarts converge and agree to **1.44×10⁻⁷ kcal/mol** at
the exact failed coordinates. The original failure is therefore avoidable for
this cell through its initial electronic state, without changing its Hamiltonian,
geometry, tolerance, temperature or cycle limit. This supports a separately named
recovery sensitivity. **The primary full100 benchmark still records the failure.**

The target is A0ACD6B9F2, La-conditioned sample0, La at adaptive_Ca, vacuum, in
`__envelope_2eda778295ffed58`. The original fresh calculation reached a500-cycle
oscillation/LEANSCF error with no normal termination; its process return code0
does not make its energy valid. No original input, output or matrix was changed.

## Exact two-start result

| Fixed seed order | Actual native guess | Cycles | Energy, Hartree | Final energy residual, Hartree |
|---|---|---:|---:|---:|
| Origin, designated recovery | XTBRESTART |67|−365.122955481166|+3.1719×10⁻¹¹|
| Own adaptive_La, agreement check | XTBRESTART |37|−365.122955480937|−4.7521×10⁻¹¹|

The signed check-minus-origin difference is+1.4371355921×10⁻⁷ kcal/mol, below the
predeclared0.1 kcal/mol gate. Both normal termination and SCF convergence are
confirmed. All three printed energy/max-density/RMS-density residual criteria
pass, and both calculations use the native mixer, TolE1e-10,300K,MaxIter500 and
one rank. The output Cartesian coordinates reproduce the actual target to the
printed six-decimal precision; no source-seed geometry replaced the target.

The only proposed energy is the first row, selected **before** execution by its
seed identity. Neither the lower-energy start nor a biological class was selected
afterward. This diagnostic does not itself compute a classifier or change a
reference. Root may join a separately labelled sensitivity using the unchanged
remaining cells. No third attempt or wider rescore is warranted by this result.

## Compatibility and preserved evidence

Both seeds are genuine paired GBW+xtbw files from successful same-context La
vacuum cells:204atoms in identical order, charge−3, singlet,636native electrons,
573basis functions and325shells. Each xtbw has17368bytes and each GBW3855516bytes.
Their native parameter exports match the failed target. Exact seed, original
input/output/receipt, context/preparation and target coordinates are pinned in
[SEEDS.json](SEEDS.json). Native binary contents are retained; the compatibility
check uses actual producer state/atom order/basis dimensions and observed seed
consumption, not a fabricated binary parser.

The existing strict-native recipe and executor are reused. Matching GBW+xtbw
are copied to the runtime basename to activate the documented native restart;
NoAutostart is the only removed input term. Pre/post seed hashes and actual
receipts are preserved. The earlier origin task's documented byte-identical
runner-path recovery is retained; inventory checks its actual runner hashes,
output, input and state rather than pretending a fresh origin calculation ran.

Result protocol:
`Nikasha_single_failed_native_cell_two_seed_recovery_sensitivity_v1`.
Status `qualified_recovery_sensitivity`; primary status
`unavailable_unchanged`. [Exact artifact pins](ARTIFACTS.json).

## Cost and checks

Job1211825 completed0:0 on node-224-2t-8gpu-1: **21allocation seconds,
42allocated CPU-seconds, zero GPUs**;2CPUs/8GiB requested. The nested executor
wall time was17.8597s. Actual job CPU was23.785s (distinct steps23.784s, rounding).
Scheduler MaxRSS0 is uninformative, not evidence of zero memory use. Exactly two
new nativeGFN calls ran; no new MACE,DFT,search,coordinates or protonation.
Original failure/seed creation costs belong to their historical runs and are not
added again. Local inventory/preparation/tests/reporting are additional unmetered
work. Full receipts and accounting are linked by ARTIFACTS.json.

Four actual-artifact tests pass in1.061s: source/seed state and order, exact
preflight/recipe, wrong-state/seed rejection, and final native restart/coordinate/
fixed-selection checks. Before execution the same collector was tested with
actual absent outputs and correctly returned unavailable. An earlier test-only
temporary manifest outside the task root hit the executor's path containment
check; its log is retained, and only the fixture location was corrected.
[Commands](COMMANDS.md); [predeclared scope](PLAN.md).
