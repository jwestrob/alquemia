# LBT local electronic / MM evaluator

Module 3, 17 October 2026. Root owns execution. This is a distinct Hamiltonian
from whole-peptide MACE; no old whole-peptide energy is added to it.

## Frozen source and energy

Use the four source/region preparations from Slurm12953, with every physical
coordinate in the original 256-atom order. Source7CCO was La-conditioned and
source7CCN was Lu-conditioned. Neither is selected by its calculated score.

| Region | Complete residues | Real electronic atoms | Caps | External charges |
| --- | --- | ---: | ---: | ---: |
| baseline | 3–14 | 170 | 2 | 84 |
| expanded | 2–15 | 208 | 2 | 46 |

Each region has charge−2; the shifted external charge is zero within preparation
roundoff and whole-system charge is−2. The existing complete interaction ledger,
ff19SB topology and pure12-6 La/Dy parameters are immutable.

```
E(R) = MACE_vac(Q(R))
     + GFN2_emb(Q(R); shifted fixed MM charges at R) − GFN2_vac(Q(R))
     + MM_ledger(R) + ALPB(unshifted fixed FF charges; R)
```

The electronic region includes the target. Caps follow both real boundary
anchors. The full gradient includes electronic, cap, point-charge, every retained
MM-group and whole-system solvent derivatives. Fixed charges move with their
atoms; this is not a polarizable MM environment. Finite hardness-damped GFN2
embedding tends to zero at infinity, with no PC self-energy or duplicate MM–QM
Coulomb. Native cross-region electrostatic1–4 scaling is not represented in
the embedded field, matching the explicitly approximate Mex policy.

MACE uses physical multiplicity1 for La and6 for Dy in its graph spin input.
Published original GFN2 remains spin-independent, UHF0/effective multiplicity1;
valence electrons are calculated from the qualified full La+Dy profile. The
profile is metadata/parity evidence for standalone xTB, not silently imported
as an ORCA calculation. Standalone calls stage the exact qualified complete
parameter text and set XTBPATH to that staged directory; hashes and load-path
source evidence are retained. No ORCA defaults, parameter remapping or new ECP
is introduced. Physical-state accuracy remains unqualified.

## Interface and reuse

`HybridEvaluator(config, fresh_output).evaluate(source, region, metal, R)`
returns total kcal/mol energy, full256×3 gradient in kcal/mol/Å, component
energies/gradients, native units, state and immutable provenance. `difference`
subtracts native component energies before conversion. The complete cache key
includes all coordinates, partition, model/checkpoint, runtime/parameter profile,
force-field XML, charge/state metadata, solvent recipe and wrapper source.

Reuse the pinned preparation `update_geometry` and `scatter`; they preserve
gradient units. Reuse the warm MACE constructor and exact native forward call,
with explicit172/210-count and target0 gates. Generalize only Mex GFN output
admission dimensions and charges; preserve every numerical/convergence gate.
Read actual serialized OpenMM systems in the original physical order through
the pinned MM Python, querying their four declared groups. Reuse the existing
256-atom fixed-FF ALPB binary and its exact origin results only when coordinates,
elements, charges, state and binary lineage match. No whole-MACE or legacy48
electronic result is a substitute for a new local-region cell.

## Finite stages

1. `bind`: metadata/pin construction only, allowed on login.
2. `preflight`: Slurm, zero new molecular/MM/solvent calls. Verify all four
   preparations, parameter/runtime lineage and native MACE metadata. Replay all
   four actual Mex GFN energy/charge/analytic-gradient records with the generalized
   admission. Replay four existing LBT whole-system fixed-charge solvent records.
   Freeze actual direction IDs/vectors at the baseline7CCO/La origin.
3. `scout`: exactly nine composite evaluations: baseline7CCO/La central plus
   two directions × two step sizes × both signs. Its origin ALPB is exact reuse;
   eight displaced solvent calls are new. At most9 MACE,18 GFN and9 MM
   configurations (36 group queries). Cache-repeat tests add no chemistry.
4. `origins`: only after actual scout admission. Reuse its central point and
   evaluate the other seven declared cells:2sources×2regions×2metals=8 total.
   Retain every unavailable cell; never select regions/sources by the outcome.
   Compare metal-dependent **source responses within each region**; absolute
   energies from different partitions are not comparable.

The scout directions are selected before any new energy: (a) the first boundary's
omitted CA, perpendicular to its cut bond, using the Cartesian axis least parallel
to that bond; (b) the first lexicographic external backbone N/O with nonzero charge,
radially toward the target. Each direction has unit Cartesian norm. Freeze exact
physical indices, IDs and vectors in `FD_DIRECTIONS.json` during preflight.

Central differences use h=2e−4 and1e−4Å. For each component and total, require
absolute residual ≤0.005+0.0001×max(|FD|,|analytic|) kcal/mol/Å at both steps.
This resolves the planned0.10 max/0.03RMS optimization criteria; it is tighter
than the old coarse derivative check. Point-charge gradients retain the actual
native5e−9Eh/bohr per-component printing bound. Report all residuals and
step-halving differences, without adjusting tolerance after scoring.

New origin scores do not authorize optimization. Expanded-region and Dy
derivative qualification, saved-pool selection, and any optimization are separate
declared stages. This evaluator implements no optimizer or production promotion.
The first origin matrix answers representability/executability and source/region
sensitivity, not force-field accuracy, experimental transfer or binding affinity.

## Units, failures, resource accounting

MACE returns energy eV and force eV/Å; negate forces before projection and convert
using23.06054783061903 kcal/mol/eV. GFN2 QM/PC gradients and ALPB gradients are
Eh/bohr; convert with627.509474/0.529177210903. OpenMM returns kcal/mol energies
and physical gradients explicitly converted to kcal/mol/Å. Sum gradients in the
fixed order MACE, GFN field difference, MM, ALPB.

All numerical stages require Slurm. Root stages one allocated GPU and12CPUs;
no memory/time overrides. Native GFN calls use one thread each and run after
MACE, not concurrently with its CPU work. OpenMM uses its Reference platform in
the pinned environment. Preserve attempted calls, stdout/stderr, convergence
evidence, worker execution receipts, elapsed times and actual allocation costs.
An exact failed cache entry blocks automatic retry. Failures never become zero
energies or disappear from the denominator.
