# Density and permanent-field diagnostic — completed

**Atomic monopoles explain a substantial part of the interaction error; allowing
the quantum core to respond to the fixed protein field removes much less of the
remaining partition discrepancy. Neither result is a validated affinity score.**
Baseline/default and every prior experiment remain unchanged.

## What the real calculations show

Same four consumed1H4I qm33/qm36 La/Ca states, exact frozen coordinates and
permanent protein charges; native ORCA6.1.1 r2SCAN-3c/DefGrid3, no solvent.

| Contribution to qm36-minus-qm33 contrast | kcal/mol |
|---|---:|
| Vacuum quantum contrast | +61.738603357 |
| Exact-density permanent coupling | −51.091194742 |
| Additional self-consistent core response | −0.974010147 |
| Remaining embedded-core contrast | **+9.673398468** |

Replacing MBIS monopoles with actual saved-density coupling changes that
partition contrast by **−8.060657641 kcal/mol**, in a common atomic-unit
convention. The old TABI conversion convention is retained separately. This
shows that passing the exterior ESP criterion alone did not guarantee accurate
core/environment coupling. Residue-level contributions remain in the utility
collection; they are not an automatic rule for selecting quantum residues.

All four embedded states converged, and response energies satisfy the declared
variational accounting check:

| State | Embedded native energy, Hartree | Electronic response, kcal/mol |
|---|---:|---:|
| 1h4i_qm33_La | -1755.033401445475 | -32.511799630 |
| 1h4i_qm33_Ca | -2400.898454753969 | -30.881230469 |
| 1h4i_qm36_La | -1983.473868193360 | -31.219472620 |
| 1h4i_qm36_Ca | -2629.323505960147 | -30.562913606 |

The absolute response is about−31kcal/mol per endpoint, but its differential
partition effect is only−0.974kcal/mol. Thus core polarization alone did not
remove the remaining roughly10kcal/mol dependence. The remaining contributions
could involve the classical residue representation, response outside the core,
or short-range quantum interactions; these causes are not yet separated.

No old reaction-field energy is attached to the new density as a validated
global score. The previous v1 partition, surface and rotation failures remain
failed. No accuracy trial or changed threshold is hidden in this diagnostic.

## Cheaper charge representation

The separate [uniform CHELPG experiment](CHARGE_FIT_RESULT.md) used the same
four saved vacuum wavefunctions, without new SCF calculations. All four fits
improve exterior potential RMS and direct-coupling agreement versus MBIS on
these states. Standalone utility times are58.593–73.329seconds on one CPU each.
Paired coupling errors versus exact density are−0.120202/−0.527673kcal/mol.
This is a promising practical representation, not an automatically promoted
model or evidence of biological accuracy. No charge scheme was chosen per case.

## Scope, recovery, tests and cost

[Authorization and pre-execution scope](AGREEMENT.md), [energy expression](ACCOUNTING.md),
[new protocol and compact result](RESULT.json):
`native_r2scan3c_permanent_field_density_diagnostic_v1`.

- 1199979: four utility calls failed before computing potentials because copied
  files lacked ORCA's densitiesinfo index;8allocatedCPU-seconds.
- 1199983: four identical-density/probe utility retries completed;108CPU-seconds.
- 1199980: four embedded DFT/MBIS endpoints completed, no quantum retry;1259s ×
  64CPUs =80576CPU-seconds; peak batch RSS8010824KiB.
- 1199984: four CHELPG utility calls completed;74s ×4CPUs =296CPU-seconds.

**Total recorded allocation:80988CPU-seconds; zeroGPU.** All failed artifacts
remain intact. Preparation, downloads and local parser testing were not fully
timed; they are not treated as free or zero. The quantum node was an Intel
XeonE5-2683v4; comparison with other nodes is not a matched speedup benchmark.

Seven real-artifact/parser tests passed, including prepared geometry/field
identity, actual utility outputs, per-residue sum closure and corrupted native
control rejection. The completed quantum endpoint receipt/state/charge checks
also pass. The completed-output regression independently checks endpoint,
paired and partition component closure without double-counting coupling.
Scientific calculations and software checks remain distinct.

Workspace: `workspaces/density_embedding_20260916/`. Final collections:
`embedded_result_v1.json`, `potential_retry_v1/potentials_1199983.json`,
`chelpg_v1/comparison_v1.json`, `accounting_final_v1.json`. Implementation copies
preserve each executed revision. No job from this diagnostic remains running.

**Recommendation:** retain baseline. Prefer exact-density coupling for diagnosis;
pursue the cheaper fitted-charge representation only with further physical
validation. Investigate the remaining residue interaction before paying for a
larger global solvent sweep. Any such calculation gets a separately declared
scope under Jacob's later autonomous-research authorization.

## Read-only collection

Run from the repository root after the corresponding jobs finish:

```bash
DENSITY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
DENSITY_WORK=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/density_embedding_20260916
"$DENSITY_PY" scripts/density_embedding.py collect-potentials \
  --manifest "$DENSITY_WORK/potential_retry_v1/potential_manifest.json" \
  --output "$DENSITY_WORK/potentials_operator_v1.json"
"$DENSITY_PY" scripts/density_embedding.py collect-embedded \
  --manifest "$DENSITY_WORK/prepared_v1/manifest.json" \
  --potentials "$DENSITY_WORK/potentials_operator_v1.json" \
  --output "$DENSITY_WORK/embedded_operator_v1.json"
```

Writers refuse existing outputs. Missing/failed utility or quantum outputs
remain unavailable.

## Separate interaction follow-on

The [predeclared two-task interaction decomposition](INTERACTION_PLAN.md) tests
the remaining Asp303 interaction. Initial job 1199985 failed before any SCF:
ORCA requires the geometry declaration before the fragment block. The exact
same two inputs, reordered without changing scientific settings, run as
1199986 in `eda_order_retry_v1/`. Native ghost-basis references already fail
the frozen equivalence check; Asp303 fragment SCFs remain unstable at the
checkpoint. [Current status and automatic completion location](INTERACTION_STATUS.md)
retain that limitation and all attempts. Its results and cost are separate
from the completed density/field diagnostic above.
