# Native energy accounting passes; direct quantum source remains research

Job1201055 completed exactly12native direct-field queries and6constrained
no-response energies. Combined with the12previous mutual-polarization controls,
**all15accounting and rigid-transformation checks pass**. Production is unchanged.
These are ion-excluded protein/water frameworks, not new metal predictions.

The actual maintained Tinker26.2/GK calculation satisfies, with native units,

```
I_vac = −0.5*(electric/dielec)*sum(mu_vac,d · F_vac,p) = ep
I_GK  = −0.5*(electric/dielec)*sum(mu_solv,d · F_solv,p)
      = ep_mutual + es_mutual − es_no_response
```

Here fields are e/Å², induced dipoles eÅ and energies kcal/mol. The distinct
d/p field definitions remain explicit. Maximum accounting residual is
3.183231456205249e−11kcal/mol against the predeclared1e−7tolerance. Permanent
energy is unchanged. Maximum rigid-field component error is
3.986394547794703e−15e/Å² against1e−8. Actual no-response arrays are all zero;
no induced solve was performed in the new calls. Cavity and physical parameters
match before sharing static references across masks/convergence settings.

The d/p reciprocal contractions differ by up to0.003047kcal at the old standard
solver tolerance, decreasing to about2.46e−5kcal at the tight tolerance. These
are recorded diagnostics, not a new fit or acceptance criterion. Neither
field definition may be silently replaced with the other. The result verifies
native evaluated energy algebra; it does not by itself establish an arbitrary
external-field Hamiltonian or prove exact stationarity at finite SCF tolerance.

## Technical recovery and real tests

Six rotated outputs were initially rejected by a validation routine that
regenerated rotated coordinates. Each native parameter/coordinate record is
**exactly equal** to the pinned parent and actual serialized XYZ. The failure
occurred only in re-evaluation of the rotation on the compute node. Its exact
last-bit discrepancy was not logged; local regeneration agrees, so no numeric
magnitude or unique hardware cause is asserted.

Recovery verifies the actual frozen coordinates and full parameter record,
reparses the original native output, and preserves all failed receipts. No
scientific inputs/tolerances changed and no calculation reran. The recovery
implementation and output are frozen under `recovery_implementation_v1` and
`recovery_pinned_v1`; an earlier live-script reparse remains recorded separately.
Five real artifact/parser/algebra tests pass in19.316s, none skipped after
execution. Before execution, two integration checks were explicitly skipped.

## Cost and remaining work

The allocation was38wall-seconds on64CPUs:2432allocated core-seconds,
106.172actual jobCPU-seconds and zeroGPU. Native kernels totaled6.304154534s;
native process wall/CPU sums16.8593394943/89.561133s. Sampled batch peakRSS
349836KiB. Build0.706143wall/0.650211CPU-s; preparation1.607914wall/1.299933CPU-s.
Pinned read-only recovery/collection19.905204wall-s. Earlier unpinned reparsing
timing was not durably captured and is explicitly unavailable, not zero.

This establishes a cheap accounting prerequisite for a direct-field adapter.
No QM source has yet entered this native solver. Covalent-boundary charges,
damping, quantum reaction-field representation and matched full-cavity
subtraction remain unresolved; no full environmental correction is available.
Predictive usefulness is untested for this component. Retain the baseline.

Protocol `Tinker26_2_native_direct_field_accounting_v1`; workspace
`workspaces/mace_omol_20260917/native_field_accounting_v1`. The compact JSON
records exact manifest/result/cost hashes and all unrounded component values.
Next: qualify a per-site field-input adapter against these existing native
controls, then define a complete boundary/source model before metal scoring.
