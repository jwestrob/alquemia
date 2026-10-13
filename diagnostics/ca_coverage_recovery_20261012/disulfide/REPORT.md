# MMOL_2048 recovered: correct origin, inconclusive accommodation

**The preparation defect is repaired, but accommodation loses a correct La call
on this source.** Origin, complete shared pool and all six solvent components are
available. The origin is La under its frozen calibration; operational and
mathematical accommodation are both inconclusive under theirs. This is a
scientific limitation of the tested model/search, not another unavailable result.
No optimization, acceptance or calibration change follows this outcome.

The excluded source contains an intact disulfide. A malformed existing hydrogen
coordinate caused its preparation failure; this defect required no new protein
fold. The versioned one-H repair passes the unchanged199-atom envelope, physical
mapping and frozen chemical-state checks. Original primary failure is retained.

## Actual discrimination and accommodation

| Protocol output | R = E_Ca − E_La, kcal/mol | Frozen decision |
|---|---:|---|
| Repaired source origin | −405453.2995597891 | La, correct |
| Operational shared pool | −405454.5167513889 | Inconclusive |
| Mathematical row minima | −405454.5167513889 | Inconclusive |

The origin La boundary is−405455.45587457. The accommodated Ca/La boundaries are
−405457.564455343 and−405454.092284296. These are the existing canonical-only
calibrations, not refits to this source. The contrast is a protocol-specific
classifier quantity, not an affinity or binding free energy.

Both metals select their own proposal from the same three geometries. Composite
accommodation lowers Ca by5.331035580 and La by4.113843980kcal/mol. Thus the
Ca-minus-La contrast moves about1.217kcal/mol toward Ca and enters the inconclusive
region. Vacuum MACE work is−6.090204225 for Ca and−4.306729215 for La; the solvent
component opposes part of each decrease. The available calculations explain the
score change algebraically, without identifying a unique missing physical cause.

The shared selected torsions are Glu195 chi1/chi3/chi2 and Asn313 chi2. In stored
selector order `(195chi1,195chi3,313chi2,195chi2)`, final displacements in radians are:

- Ca: `[-0.046067433, -0.404428899, 0.155114241, 0.245482846]`.
- La: `[-0.143953987, -0.321569653, 0.082547963, 0.179535335]`.

Both original optimizers report success (15/14 iterations), with no active angular
or heavy-displacement boundary and no new severe overlaps. Maximum heavy motion
is0.736781Å for Ca and0.611757Å for La. All other physical atoms are fixed. These
are finite MACE-vacuum proposals evaluated with the composite; composite gradients
were not computed and unconstrained minima are not claimed. Residual omitted-mode
loads remain recorded; this result does not diagnose them as the cause.

## Solvent qualification

All six candidate/metal components pass the unchanged component-v2 gates:
absolute difference between xTB and tblite ALPB-minus-vacuum corrections≤1e−6Eh,
population difference≤1e−5e, converged/state-admitted endpoints and source charge
closure≤1e−5e. Paired-coordinate errors are exactly zero for all candidates.

| Component | Correction error, Eh | Maximum population error, e |
|---|---:|---:|
| Ca/origin | −1.019720912e−7 | 5.757314259e−8 |
| Ca/Ca proposal | −1.011989639e−7 | 6.724570284e−8 |
| Ca/La proposal | −1.014972213e−7 | 7.938494428e−8 |
| La/origin | −8.609185898e−8 | 5.719703139e−8 |
| La/Ca proposal | −8.566587439e−8 | 3.704880941e−8 |
| La/La proposal | −8.616348168e−8 | 4.386882596e−8 |

Absolute backend energies still differ by roughly1.60–1.67e−4Eh. The historical
absolute-v1 gate remains failed; component-v2 acceptance does not erase that
failure or establish physical accuracy from backend agreement.

The125-row annotated overlay, including the two actual MxaF technical recoveries,
is107correct/1wrong/2inconclusive/15unavailable. Original primary accommodation
remains105/1/1/18 and every original row remains present. This source adds coverage
and an explicit abstention, not a new correct accommodated prediction. No new
three-source aggregates were constructed.

## What actually failed

Slurm2864 inspected the rejected MMOL_2048 Ca seed1/sample0 source, accepted
Ca sample1 and the original La sample0 state anchor. The source topology contains
the real Cys127–Cys128 bond at2.117557 Å. Cys128 HB2 was present, but its CB–H
distance1.376989 Å exceeded the source graph's actual1.340000 Å parent-assignment
cutoff. The resulting missing H-parent entry triggered `unsupported disulfide
partner` before the S–S distance test. Both controls passed without alteration.

## What changed and passed

Slurm2870 restored only that existing HB2 along its original CB–H ray to the
pinned ff19SB CYX equilibrium1.09 Å. The explicit XML parameter is
`protein-H1`/`protein-2C`, length0.10900000000000001nm; topology SHA
`086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09`.
After PDB serialization, CB–H is1.089972477 Å and H displacement0.287017421 Å.
Every other coordinate and atom identity is unchanged. This is radial
preparation repair, not an angular optimization or a different protonation state.

The original source/protonation records remain immutable. New records point to
the repaired source, retaining the original core XYZ and fixed fragment union.
No builder cutoff or admission rule changed. Actual results:

- Target:199atoms, paired maps equal, exact original small-context replay.
- Frozen source-state key matches the existing canonical group:
  `f8eeff480ce4ee77f17fbccba0dd5f458922d734abd2320c94b3c997ee2cf2eb`.
- Two untouched controls reproduce their archived states and physical maps.
  Maximum coordinate-copy residual is2.22e−16 Å; the La control's XYZ hashes
  are identical, while Ca control serialization differs only at this rounding scale.
- Core coordinates, heavy geometry, proton/water inventory, charges, selection,
  canonical bands and original failed benchmark row remain unchanged.

Protocol: `Nikasha_disulfide_unmapped_source_CH_radial_repair_v1`.
Real receipt: `workspaces/ca_coverage_recovery_20261012/disulfide/repair_v1/RECEIPT.json`.
Result SHA: `147aedd094e2a43837c519d7a6a53fa7bbd3a7a7f5ac6260b43793b741d2b367`.

## Measured work and final disposition

Diagnosis:4.335s measured script wall time. Preparation:7.875s wall,
6.664processCPU-s, peak RSS177300KiB on one Slurm CPU. No energies, forces,
optimizations, new protons or waters were calculated.

Slurm2875 completed the one-source native stage on seven CPUs at
node-64-768g-5: both origins, both unchanged four-mode proposals and both cross
evaluations form a complete shared pool. The two existing runtime force fixtures
passed with exact energies and maximum force residual9.55e−14 eV/Å. This stage
used33 actual forwards,463.764s elapsed,2549.236 measured modelCPU-s and
3246.348 allocatedCPU-s. These are this CPU recovery's measured costs, not
GPU throughput or a full new-site latency estimate. All started model calls
have terminal receipts. The same-protein preparation controls were not rescored.

Corrected solvent Slurm2880 completed12 xTB and12 tblite calls, all receipted,
with784.661processCPU-s. The scientific stage including final collection took
80.224s; adapter preparation plus collection/overlay took83.558s and
1002.695allocatedCPU-s on12CPUs. These costs are separate from native evaluation
and source repair above. No model/solver call failed or remained partial.

`RESULT.json` copies the actual unrounded selections, bands, components, geometry
checks, source pins and receipts. Final solvent collection SHA:
`865224278a143023fceefc28dfe28c15d314185254a15fa174d2bd1ec9abecaa`.
Combined overlay SHA:
`d99d5c58c27b4a9260ad3dc8ba929658e2b72508d71e7715c367813649e7549d`.

Retain the preparation recovery and this negative accommodation evidence.
Do not tune the source, search or bands to make this consumed example La-like.
No additional computation, production promotion or independent push was done.

Original125-source results, frozen calibration and every original failure remain
visible. Later results must be appended as named preparation recovery.
