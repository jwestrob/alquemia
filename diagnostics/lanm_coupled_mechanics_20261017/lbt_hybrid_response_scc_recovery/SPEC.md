# Four-call embedded SCC diagnostic at the failed response trial

## Question and preserved evidence

Does a declared change to initial Broyden damping recover the same embedded
GFN2 solution at the failed 7CCN/baseline/La response trial, while retaining the
accepted default solutions at two fixed controls?

Response job 13019 is unavailable after two accepted steps and a failed third
trial. Its point003 (`ddcaef5b`) has an unavailable embedded GFN2 component and
an accepted isolated component. The embedded SCC reached its unchanged 500
iteration limit; the recorded final RMSdq is 0.380 and dE is +1.94112 Eh. These
are unconverged iteration diagnostics, not accepted energies. The preceding
point002 (`2034ba36`) and original 7CCN/baseline/La origin are the controls.
All original results, failed requests and optimizer budgets remain unchanged.
This diagnostic neither resumes a search nor admits a replacement composite.

## Fixed matrix and method

| Case | Archived embedded input | Initial broydamp | Agreement reference |
| --- | --- | ---: | --- |
| failed_005 | Response point003 | 0.05 | failed_010 |
| failed_010 | The same point003 | 0.10 | failed_005 |
| preceding_005 | Response point002 | 0.05 | Its accepted default output |
| origin_005 | Original 7CCN/baseline/La origin | 0.05 | Its accepted default output |

**The primary candidate is failed_005**, declared before these outcomes. Exactly
four serial, single-thread native calls are permitted. There is no automatic
alternate damping, geometry, iteration limit, restart guess or subsequent run.

Reuse the published xTB 6.7.1 executable and full author parameter file from the
qualified hybrid evaluator: 172 QM atoms, 84 external point charges, QM charge
−2, 538 valence electrons, UHF 0/effective multiplicity 1, physical La
multiplicity 1, 300 K, accuracy 0.0001, 500 iterations, `--norestart`, populations
and analytic gradients. No solvent is requested. Exact archived core XYZ,
external PC coordinates/charges/hardness values and parameter bytes are copied.
The archived physical map and cap Jacobians are retained, not regenerated.

The existing embedding input is preserved and only the SCC block is appended:

```text
$embedding
 interface=orca
 input=external.pc
 gradient=pcgrad
$end
$scc
  broydamp=0.05
$end
```

The second target call changes only 0.05 to 0.10. The existing command already
has `--input control.inp`; it is retained exactly once. The installed binary,
package manifest, manual and archived parser/mixer source are pinned as in the
completed isolated diagnostic. The native output must print the requested
initial damping. Native adaptive damping behavior is unchanged.

## Unchanged admission and agreement gates

Reuse the existing native call, parser, GFN_ADMISSION and physical scatter
helpers. Verify all archived raw artifacts and re-admit both accepted default
embedded controls before calls. Bind the failed trial and preceding point to
the actual 13019 request list and terminal receipt, and the origin to the
complete original eight-origin result. The companion isolated result is
preserved and verified as accepted; it is not rerun.

All four new calls must meet unchanged convergence, state, population/charge,
energy/gradient consistency, analytic-gradient, field-inventory and native
termination gates. Compare the two failed-point damping results and each of
the two controls against its exact accepted default output. Each of the three
comparisons must meet both existing agreement thresholds:

- Absolute energy difference ≤ 1e−7 Eh.
- Maximum Cartesian gradient-component difference ≤ 0.005 kcal/mol/Å,
  independently for all 172 QM atoms, all 84 external PC sites and all 256
  physical atoms after cap and field scatter.

Native QM and PC gradients are dE/dR in Hartree/bohr, converted using the
unchanged admission constants and multiplication order. Physical scatter
includes both cap-anchor contributions and all actual PC gradients. The two
omitted cap anchors remain physical atoms although they are absent from the
84-charge field. Full charge vectors and their maximum difference are
descriptive; no new charge-agreement threshold is introduced.

Agreement is a finite solver-branch check. It does not establish a unique or
stable ground state, physical spin qualification, force accuracy, an acceptable
response endpoint or a general solver fallback. No energy sign is an admission
criterion. The spin-independent GFN2 model and its frozen-f limitation remain.

## Execution, failure and accounting

`SCOUT.py --stage preflight|run --config CONFIG.json --output FRESH`

Both stages require Slurm. Preflight performs only saved-artifact validation and
raw control replay. Run performs that same preflight and at most the four
declared native calls, serially, with one thread each. No HybridEvaluator,
optimizer, MACE, MM, ALPB or DFT is constructed or executed.

Retain all four slots. An SCC/admission failure leaves its slot unavailable;
other independent declared controls may still run. A native lifecycle failure
or mutated input/runtime stops further scheduling. Qualification requires all
four admissions and all three agreement checks. Save exact tasks, changed
control input, execution/runtime provenance, raw outputs, admitted native
payloads and actual per-call CPU/wall costs. The new NATIVE_RESULT receipts
explicitly identify their changed solver control and original component; old
failures are never overwritten. Root owns submission and any later scientific
decision. No response continuation is authorized by this script.
