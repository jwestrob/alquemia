# Whole-Mex force consistency: passed

19 October 2026. The whole2049 MACE vacuum plus fixed-charge ALPB model passes
all32 declared directional comparisons at the saved BEFORE geometry. Both
steps pass separately for MACE, solvent and their sum. This supplies the sampled
implementation check missing from the eight-cell whole-system panel.

It establishes numerical consistency of the tested forces and energies. Physical
accuracy, the correct aggregate electronic state, forces at other geometries,
converged accommodation and natural-protein metal preferences remain unresolved.

## Actual result

| Saved state input | Multiplicity | Direction/step checks | Outcome |
| --- | ---: | ---: | --- |
| La aligned | 10 | 8/8 | Pass |
| Dy aligned | 15 | 8/8 | Pass |
| La with one Nd reversed | 4 | 8/8 | Pass |
| Dy with one Nd reversed | 9 | 8/8 | Pass |

Four original directions cover a backbone carbon, the former partition-boundary
anchor, a complete water and a spectator Nd. Both original steps0.0002/0.0001Å
and the frozen componentwise gate were retained:

`abs(FD - analytic) <= 0.005 + 0.0001 * max(abs(FD), abs(analytic))` kcal/mol/Å.

For example, the Dy-aligned backbone comparison at0.0002Å has MACE residual
0.0000506667 against tolerance0.00856224; its total residual is0.0000506389
against0.00818884. All other declared comparisons also pass; no step or state
was selected after inspecting results. Exact per-component values are retained
in the [losslessly archived producer result](actual_run_v1/RESULT.json.gz).
The [result pin](actual_run_v1/RESULT_PIN.json) identifies the original readable
workspace JSON and its exact SHA256.

## Execution and preservation

Zero-call preflight14810 passed on1CPU in8s. It re-admitted all four MACE centers,
two ALPB centers and16La solvent displacements, including exact coordinates,
states, native receipts and gradient conversion. Result `9a66cf94`.

Run14811 completed64/64 new MACE evaluations and16/16 new Dy ALPB evaluations.
No central, DFT, GFN, MM or optimization calls occurred. One warm H200 model
loaded in6.75s; recorded peak allocation was81,684,463,616bytes. The job took81s
on12CPUs/oneH200. Including preflight, allocated cost is980CPU-seconds and81GPU-
seconds. [Exact accounting](actual_run_v1/COST.json).

Result SHA256: `c8c6166452996b862da4dc9fb35a758bc163ff44a60d60ceeba79051f27c4c07`.
The run and preflight share manifest `27c91c6c`; source `bf32bc23`, run config
`94718e39`, and all native requests/coordinates/forces/solvent outputs remain
bound at `workspaces/lanm_coupled_mechanics_20261017/mex_whole_force_run_v1/`.
The complete manifest stays in the workspace with an archived exact pin; no
generated coordinate copy is needed in source history. See the [specification](SPEC.md),
[source review](INDEPENDENT_REVIEW.md) and [actual result review](RESULT_REVIEW.md).

## Decision

Close this finite consistency check as passed. Preserve both aggregate-state
hypotheses and their previously observed energy sensitivity. This result does
not release an optimizer or natural-source scoring. It does not repair the old
partitioned model, qualify the closed embedded-DFT derivative recipe, or reopen
the failed electronic-base candidate. The broader mechanics goal remains unmet.

Independent7CCN whole-PBE0 reference14795 remains running under its existing
scout-admission and measured-cost gate. That reference addresses accuracy of
the fixed deformation on a second source; this check addresses implementation.
