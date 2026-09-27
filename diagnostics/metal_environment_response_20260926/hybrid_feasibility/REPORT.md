# 1H4I hybrid feasibility: reusable parameters exist; exact coupling remains unbuilt

**Pursue a bounded parameter-and-boundary preparation audit, not another broad
search or immediate optimization.** Protein mechanics and Ca/La Lennard–Jones
parameters are available locally. The missing pieces are narrower than “no
metal/PQQ parameters anywhere”: an exact-state PQQ cross-LJ typing, an explicit
choice of metal LJ family suitable for the embedding approximation, and a
verified additive boundary/exclusion ledger for this exact source. No complete
1H4I hybrid was located, built, or evaluated in this read-only audit.

## Reusable concrete assets

| Component | Actually available | Remaining issue |
|---|---|---|
| Exact-source protein parameters | `scripts/affordable_state.py:source_protein()` instantiates ff19SB with `NoCutoff`, no constraints and no rigid-water assumptions before extracting charges. The current successful environment therefore derives from an actual full standard-protein System. | Its caller returns charge/source records and discards the System's other parameters. Re-export bonded, LJ and exceptions for the same source, without changing protonation or terminal chemistry. |
| Existing serialized protein parent | `workspaces/scaffold_environment_20260922/inventory_v2/1H4I/{parent_system.xml,term_supports.json,parent_atoms.json}` | This is the other 1H4I preparation: 9,060 protein atoms, charge −9; current source is 9,113 atoms, charge −8. Do not transplant its indices or parameters blindly. |
| Protein boundary map | Current `scout_v3/core_mapping.json` and `boundary_mapping.json` identify Glu177, Asn261, Asp303 cut bonds, cap Jacobians, removed MM charges and redistribution recipients. | Force projection exists, but a complete chosen additive boundary Hamiltonian and corresponding MM exclusions have not been instantiated. |
| Metal LJ parameters | Installed Amber ion frcmod files explicitly include Ca2+/La3+. | They are parameter-family-specific empirical models, not automatically qualified cross-LJ terms for an embedded QM metal. |
| PQQ atom typer | Installed `antechamber`, `ATOMTYPE_GFF2.DEF`, and `gaff2.dat` version 2.2.20; `parmchk2`/`tleap` also present. | No exact-state PQQ atom-typed LJ artifact found in targeted legacy paths. A chemically explicit bond-order graph still needs pinning and inspection. |
| Legacy topology accounting | `benchmarks/hans_lanm_dy_qmmm_correction_v1/qmmm_calibration/low_level_subtraction.py` | Useful source/receipt/exclusion audit patterns, but hard-coded Hans Q155/periodic/capture contracts are not an executable 1H4I recipe. |

Exact file paths and SHA256 values are recorded in `ARTIFACTS.json`.

## What the metal files actually say

Under `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/dat/leap/parm/`:

- `frcmod.ionslm_1264_opc`: Ca2+ Rmin/2=1.602 Å,
  epsilon=0.08034231 kcal/mol; La3+ Rmin/2=1.719 Å,
  epsilon=0.15131351 kcal/mol. Comments cite Li et al., JCTC 2020
  16:4429 and JCTC 2021 17:2342. This is explicitly a **12-6-4 OPC family**.
- `frcmod.ions234lm_126_tip3p`: Ca2+ Rmin/2=1.649 Å,
  epsilon=0.10592870 kcal/mol (CM set, Li et al., JCTC 2013 9:2733);
  La3+ Rmin/2=1.718 Å, epsilon=0.15060822 kcal/mol (IOD set,
  Li et al., JPCB 2015 119:883). This is a **12-6 TIP3P family** and is
  concrete evidence that pure-LJ candidates also exist.

These numbers are read from installed files, not selected by any score. Neither
family has been chosen for the new model. Removing C4 from a 12-6-4 fit does not
turn the remaining A/B coefficients into a separately validated 12-6 model;
using the pure-12-6 family also does not establish transfer to QM/MM. Freeze and
justify the chosen physical approximation before evaluating scientific outcomes.

The legacy `make_zero_c4_topology()` explicitly zeros only LENNARD_JONES_CCOEF
and verifies that every other prmtop field stays unchanged. Its
`mm1_electrostatic_completion()` tracks topology identities, exclusions and
1–4 corrections when replacing QM and boundary charges. These are valuable
implementation precedents. They do not justify adding a C4 induction term to
a QM region already polarized by embedding, or importing its periodic/hard-coded
Hans correction into the finite 1H4I energy.

## A practical route to exact-state PQQ cross LJ

The current PQQ is wholly electronic, so a complete classical PQQ internal
bonded/charge model is not needed merely to define *nonbonded repulsion and
dispersion across QM/MM*. It still needs unambiguous LJ types for all 27 actual
PQQ atoms, including its three prepared hydrogens.

1. Pin the oxidized PQQ3− graph with the current 24 heavy atoms/three hydrogen
   identities and exact coordinates. Existing `pqq_microstates.py` supplies
   named connectivity and state, but its bond tuples do not encode bond orders.
   Obtain and check explicit bond orders/aromaticity/formal charges against the
   source chemistry; do not infer them solely from compressed metal geometry.
2. Use the installed standard GAFF2 atom-typing machinery on that fixed graph,
   without optimization or choosing atom types by preferred La/Ca outcomes.
   Map output atom types back one-to-one; require each LJ type to exist in the
   pinned GAFF2 table. This is a proposed preparation, **not an executed result**.
3. Extract only the established per-type LJ parameters needed for cross pairs.
   Keep QM partial charges absent from the classical Coulomb ledger; no dummy
   zero-charge/full-MM PQQ model is necessary or justified. If a converter
   demands additional MM records, resolve its contract explicitly rather than
   silently filling them with zeros.
4. Use the existing documented Amber/OpenMM combining convention and exact
   topology exceptions, not an invented mixing formula. Audit actual exported
   pair coefficients and exclusions before coupling to the electronic energy.

Targeted searches of `workspaces/mxaf_qm`, local PQQ workspaces, the sibling
`fep` and `tannase_fep` trees, legacy scripts and relevant reports found no
PQQ `.mol2`/`.frcmod`/`.prmtop` to reuse. This is a scoped negative search, not
proof no such file exists anywhere on the cluster. `mace_gb` system.xml matches
are not LJ parameterizations: that code explicitly uses zero epsilon for its
separate electrostatic/GB diagnostic.

## Boundary accounting is a model choice, not an impossible missing constant

The [ORCA additive QM/MM documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-molecules.html)
assigns MM internal and QM/MM bonded/LJ interactions to the force field and
requires the corresponding parameter file. Thus an existing additive scheme is
a coherent route once exact-source typing and exclusions are established.
Its actual link handling and input/output behavior still need verification.

The prior scaffold report's inability to subtract a matching *capped MACE
protein reference* should not automatically block ordinary additive QM/MM.
Those are different energy constructions. A new additive construction must
explicitly retain chosen crossing bonded terms, exclude wholly electronic MM
terms, avoid ordinary QM/MM Coulomb duplication, handle existing redistributed
boundary charges consistently, and account for cap forces through the source
Jacobian. Cap-dependent QM contributions remain an approximation to qualify,
not evidence that arbitrary FF penalties can be added freely.

Current r2SCAN-3c D4/gCP applies to the actual electronic region. Cross LJ must
not duplicate an already included cross-region dispersion term. The finite
point-charge scout contains electrostatic embedding only; its absence of
cross LJ is a limitation, not a reason to pretend those coefficients are zero.

## Decision and limits

The most useful next preparation would export the **current** ff19SB parent
parameters and term supports, audit a fixed-graph PQQ GAFF2 LJ assignment, and
write the chosen native additive boundary ledger. This is finite reusable work
with plausible payoff; no new optimization or blanket electronic campaign is
needed to establish it. Metal LJ selection remains an explicit approximation
to qualify, not a production default or a proven affinity model.

Read-only inspection and hashing only. No parameterization executable, new
System construction, energy/force call, dependency installation or job ran.
