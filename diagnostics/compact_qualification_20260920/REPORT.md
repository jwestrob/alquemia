# Compact solvent challenger: native precision supported

**All six tighter-native checks pass.** Tightening native GFN2's printed energy
criterion from1e-6 to1e-8Eh changes the solvent score correction by at most
**0.025545kcal/mol**, below the unchanged0.20 criterion. All24 tighter-native
endpoints converged, with identical nuclear/electronic states and parameter sets.
This supports the numerical precision of the working native recipe on these
consumed contexts. The independent ordinary-SCF comparison remains unavailable.

The separate [fresh scanner pilot](../compact_scanner_20260920/REPORT.md) measures
actual MACE-plus-solvent operation. Neither continuation changes the successful
formula, published bands, aquo reference or production default.

## Tighter native results

Values are new minus original; transfer and score differences are kcal/mol.

| Context | Ca solvent-transfer change | La solvent-transfer change | Ca−La correction change |
|---|---:|---:|---:|
| 1H4I | +0.000916 | -0.000150 | +0.001066 |
| Q9Z4J7 | +0.000680 | -0.003291 | +0.003971 |
| 4MAE | +0.000032 | +0.025577 | -0.025545 |
| A0ACD6B9F2 | -0.001009 | +0.002418 | -0.003427 |
| 1F6S | -0.000107 | -0.000886 | +0.000779 |
| 1GLG | -0.000348 | -0.000582 | +0.000234 |

Maximum endpoint-transfer change0.025577 is below0.10; maximum score change
0.025545 is below0.20. Actual outputs confirm the tighter TolE and continued use
of the native xTB mixer. Parameters, atom order, coordinates, charge, multiplicity,
electron count,300K smearing, water inventory and solvent selection were checked.
No improved score was selected as a replacement for a published result. Raw
energies, charge audits, receipts and exact differences are retained in
[RESULT.json](RESULT.json) and its pinned collections.

## What happened to the ordinary solver?

1. The earlier default ordinary-SCF attempt oscillated by large energies instead
   of approaching convergence. The native solver converged routinely.
2. Eight charge/multipole `.xtbw` restart attempts again used SAD, ignoring the
   proposed restart. All eight failed. Their first iteration reproduces the
   earlier unseeded result; these were not effective seeded comparisons.
3. Eight orbital-restart attempts using the simple MORead keyword printed an
   explicit warning that MOInp was ignored because Guess was not MORead. All
   eight failed. Native-xTB input setup overrode that simple-keyword request.
4. Explicit `%scf Guess MORead` and `MOInp` finally produced actual MOREAD starts
   in all eight outputs. **One of eight converged**:1H4I Ca in ALPB, at96 cycles.
   Its energy differs from the native primary by only-0.000505kcal/mol. The other
   seven remain nonconverged; there is no complete solvent-transfer or metal pair.

This distinguishes input-control failures from a correctly initialized solver
that still struggles. It supplies neither proof of a wrong native result nor a
passed independent-solver check. No extra solver search, parameter refit or
successful-backend substitution was launched. The actual failed attempts remain
visible, including their cost. The effective restart syntax is documented by
the [ORCA initial-guess manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#restarting-scf-calculations);
the observed override behavior comes from the saved installed6.1.1 outputs.

## Qualification and use

- **Native convergence precision:** supported on six contexts by the declared
  tolerances. This is a same-state numerical check, not fresh biological validation.
- **Fresh operational reproducibility:** measured separately by the four-site
  scanner pilot, with original model/recipe/bands and actual uncached evaluations.
- **Independent ordinary-SCF pair agreement:** unavailable. One agreeing endpoint
  is insufficient to qualify a transfer or Ca-minus-La score.
- **Physical/general biological validity:** remains limited by the approximate
  low-level solvent transfer, cluster boundary and consumed benchmark set.

Recommendation: use the completed **opt-in research scoring path** for controlled
comparisons alongside the released baseline. The current native precision result
is useful; resolving every alternative solver is not the sole condition for
progress. Preserve the original25/25PQQ plus3/3transfer and6/6structural-direction
evidence without claiming additional independent biological observations here.
There is no new model or calibration version and no automatic default promotion.

## Actual cost and implementation

| Job | Calls | Successful | Failed | Wall seconds | Allocated CPU-seconds |
|---|---:|---:|---:|---:|---:|
| 1203265 | 24 native-tight +8 xtbw-restart | 24 | 8 | 127 | 8,128 |
| 1203284 | 8 simple-keyword orbital restarts | 0 | 8 | 65 | 4,160 |
| 1203308 | 8 explicit-block orbital restarts | 1 | 7 | 60 | 3,840 |
| **Total** | **48** | **25** | **23** | **252 summed job-seconds** | **16,128** |

Zero new GPU time or DFT/MACE calls in this numerical track. All jobs requested
64CPUs/128GiB, eight concurrent eight-rank tasks. Scheduler memory/cost records
are in [COSTS.json](COSTS.json). The separate scanner pilot includes its own
failed scheduler-layout attempt and GPU allocation; neither cost is hidden here.
Original exploratory32-call plan and both technical fixes were recorded before
their corresponding submissions, with no changed chemistry or thresholds.

New `compact_solvation_qualification.py` reuses the existing runner, native input
parser and actual source receipts. It exposes prepare/validate/execute/collect,
preserves failed rows and verifies actual restart behavior. Five real-fixture
tests pass, including ignored versus effective restart output. These software
tests are distinct from the48 actual scientific attempts. No fake successful
output or numerical DFT gradient was generated. See [COMMANDS.md](COMMANDS.md).

The working primary scorer and immutable earlier panel remain in
[the original result](../compact_solvation_20260920/REPORT.md). Existing CC1202429
continues independently under its existing monitor; no new CC work was launched.
