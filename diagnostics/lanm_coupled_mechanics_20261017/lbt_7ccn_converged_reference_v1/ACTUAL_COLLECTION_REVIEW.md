# Independent completed comparison review

**Collector 15838 completes the declared two-endpoint comparison with zero new molecular calls. It supports the finite composite energy lowering, while retaining a material work discrepancy and a nonstationary reference force at the model-converged endpoint.**

## Binding and denominator

- Collector RESULT: `407bdf88e42c43cb2ae48ed1493cbbc008f5c0eb957628787a6aa1c1b795bb0d`.
- Executed source: `f251ca1e9514fb0c76659faf603b03ce2daf3684204fe4e98ef44e5d448bb240`.
- Collection config: `cd2189abe830ecfef2527718775bc72dc5db3033f0fe2104e7b117eccb125ef0`.
- Bound endpoint 15545 RESULT: `a53125b98e597ea0473e06d58316ac2b32c914e9647da595f163cf0aa76a8c19`.
- Old/new force-report hashes: `4f96db5dd86b2c11bfe2a4c647df04e59957d069d3eac1dfee10f4a1c5efcd68` / `e6c0d039500c9deb87a53f31f974dfdc741c1cd801e481b51dc3316d100b7a7c`.

The actual result is under `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_converged_reference_v1/collection/execution/`. I verified the selected source/config/parent/preflight/decision/direction/report pins. Both native rows, input evidence, solvent admissions, scientific binding and caveats match the native run exactly; both copied force-report hashes are unchanged. `full_two_available=true`, comparison status complete, no selected native execution tasks, and every molecular/geometry/optimizer count is zero. There is no unavailable survivor or substituted endpoint.

## Reported result

All values below are copied from the saved comparison, in kcal/mol, for **step 725 minus step 92**. This is not the shorter step 172→725 optimizer increment.

| Contribution | Work |
| --- | ---: |
| PBE0 vacuum | +15.133560 |
| MACE vacuum | +16.362172 |
| Shared fixed-FF ALPB | −45.442067 |
| PBE0 + ALPB | −30.308506 |
| MACE + ALPB | −29.079895 |
| MACE-minus-PBE0 work discrepancy | +1.228612 |

The reference supports the sign and substantial composite lowering. The discrepancy exceeds the original 1 kcal reporting resolution; that flag is not a failed native-admission gate, error bar or criterion to retune. Both vacuum contributions rise and both composites lower because the same solvent contribution outweighs them. This verifies neither solvent physics nor its missing density feedback.

At step 725, reference+ALPB maximum/RMS force is 7.868978/3.587757 kcal/mol/Å, whereas model+ALPB gives 0.075055/0.029611. Thus finite work support coexists with disagreement about stationarity. Full-force accuracy, electronic-state stability, a physical minimum and preference remain unqualified. Native details and retained failed residual flags are in ACTUAL_NATIVE_REVIEW.md.

## Bounded next physics question

The most direct next diagnostic is whether the **existing PBE0 analytic gradient is consistent with its own energy slope in the specific direction of this model/reference discrepancy**. A centered pair of vacuum singlepoints at the exact step 725 geometry, along one frozen normalized full 256-atom force-difference vector, is a coherent minimal test. Reuse 15545's central analytic gradient; do not repeat the EnGrad calculation or add solvent calls to test the vacuum derivative. The shared-solvent cancellation makes the same direction relevant to the model/reference force discrepancy.

Before execution, freeze the direction's exact all-atom construction/sign/normalization, displacement magnitude, original unit/sign conversion and numerical gate, with both displaced cells mandatory. Preserve the same molecular state, basis/ECP, SCF/COSX recipe, 12 MPI runtime and native admissions; record actual residuals and final integration diagnostics. Use native endpoint energy differences before conversion and account for printed-energy resolution. No geometry relaxation, state selection, solver sweep or automatic retry is justified.

One step and one direction test only that sampled derivative. Agreement would strengthen the interpretation that the admitted model/reference disagreement is in their energy response; it would not validate every force component or establish PBE0 as a physical oracle. Failure could reflect analytic consistency, incomplete SCF stationarity, root/path changes or finite-step effects; it would not uniquely prove a gradient bug. Do not change a gate after seeing the result. Any inconclusive step-stability check needs its own finite declaration.

This is more discriminating than another optimization or reference geometry because it tests the credibility of the force evidence already driving the decision. It is still expensive: 15545 reports 349.478 minutes for SCF work and 77.770 minutes for gradient evaluation. Two singlepoints avoid a new gradient but each still performs a full SCF. The existing documentation review in `mex_dft/gridx3_diagnosis_v1/DOCUMENTATION_REVIEW.json` records modern COSX derivative treatment and separates unknown external-point-charge guarantees; it does not establish a nuclear-gradient defect here.

This recommendation is source-only and conditional. Root owns any separate specification and execution decision. No new numerical values were calculated on login, and no new jobs, scientific imports, adapter changes or molecular calls were made by this review.
