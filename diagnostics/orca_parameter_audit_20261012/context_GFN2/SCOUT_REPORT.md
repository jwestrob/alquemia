# Linux GFN2 scout: primary reproduces archive; absolute independent gate failed

The four exact q9z4j7 original-context adaptive_La inputs reproduce the archived
corrected standalone Mac energies within1.1e-12Eh. Independent tblite0.7.0
electron populations agree within8.4e-8e, but its total energies differ by
6.47–6.73e-5Eh. **The original absolute-energy gate failed4/4**; its80-cell
collection retains76 unexecuted cells. No threshold or historical output changed.

| Metal / medium | Primary xTB energy / Eh | xTB minus tblite / Eh | Maximum charge error / e |
|---|---:|---:|---:|
| Ca / vacuum | -259.041699230004 | -0.000064681258209 | 2.6174e-8 |
| Ca / water | -259.680806949445 | -0.000064789247403 | 8.3142e-8 |
| La / vacuum | -259.743743660661 | -0.000067256124282 | 3.0585e-8 |
| La / water | -260.205179278232 | -0.000067348682137 | 7.9023e-8 |

Unrounded energies, atom charges, SCF traces, exact XYZ and software hashes are
in `workspaces/orca_parameter_audit_20261012/context_gfn2_v1/`.

## Cause and relevance to Nikasha

A four-call causal diagnostic changed only the non-self-consistent three-body
D4 factor s9 from5 to0 in separate Ca/La vacuum evaluations. Backend energy
residuals then fell to1.55246e-7 and1.48193e-7Eh. Thus the discrepancy lies in
that three-body contribution to the original gate's precision. It is
geometry-dependent, **not an element reference constant**. Exact-cutoff matching
was not executed; the inspected source implements40-Bohr versus25-Bohr cutoffs.
Both implementations construct that term from neutral-reference weights and
coordinates; ALPB and solved atomic charges do not enter it.

Nikasha uses GFN2(ALPB)−GFN2(vac) at the same geometry. The non-self-consistent
term therefore cancels algebraically within each backend. The actual scout
solvent-transfer discrepancies are1.07989e-7Eh for Ca and9.25579e-8Eh for La;
their paired contrast discrepancy is1.54313e-8Eh. These are direct differences
of computed totals, not estimated offsets. This diagnoses the interface
comparison; it does not demonstrate improved discrimination.

Root explicitly authorized a separately named component-qualified v2, with
the same1e-6Eh numerical tolerance applied to each metal's actual solvent
transfer, plus the original per-input population, convergence and state checks.
See [the frozen v2 plan](COMPONENT_V2_PLAN.md). It reuses the four scout inputs
and runs76 additional inputs per backend. All absolute residuals remain visible;
the failed v1 protocol is preserved. This changes the qualification target to
the component actually used, without changing either Hamiltonian or class bands.

## Executed work and tests

The scout used4 primary plus4 independent singlepoints: summed worker wall
73.057s, measured process CPU72.796s, maximum observed RSS82432KiB. The workers
ran at most two one-thread processes concurrently. Sum of worker times is not
batch latency or allocated Slurm time; this was specifically authorized login
execution. The dispersion diagnostic adds4 actual calls, with separate receipts.
No MACE, ORCA, optimization, DFT or Slurm job was run by this agent.

Four real-artifact parser/algebra tests pass: failed historical gate retention,
direct same-coordinate solvent subtraction, explicit missing components, and
the causal s9 diagnostic. These tests make no molecular calls. The actual scout
and dispersion diagnostic above are separately executed scientific evaluations.

```bash
python diagnostics/orca_parameter_audit_20261012/context_GFN2/TEST_REAL_SCOUT.py
```

Source accounting and standalone run commands are in PLAN.md,
DISPERSION_DIAGNOSTIC_PLAN.md and COMPONENT_V2_PLAN.md. No production promotion,
calibration or classification claim follows from this numerical qualification.
