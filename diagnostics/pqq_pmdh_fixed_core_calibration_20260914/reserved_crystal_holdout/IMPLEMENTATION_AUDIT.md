# Reserved crystal holdout preparation implementation

Date: 2026-09-15
Protocol: `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`
Status: pre-energy implementation amendment; no ORCA holdout energy has been
launched.

## Pre-energy release amendment

The first full preparation attempt revealed **130 missing standard-residue
heavy atoms in 6OC6**, while both primaries, 1H4I and 4MAE, had **zero**. The
fatal heavy-atom-repair gate is unchanged: 6OC6 is not repaired, rescued, or
silently accepted. Because 6OC6 is a nonindependent secondary geometry check
(its sequence is already represented by C5B120), its preparation defect must
not prevent the preregistered atomic primary pair from being tested. The
release may therefore contain exactly 1H4I+4MAE (four La/Ca legs), or those two
plus explicitly requested 6OC6 (six legs). The primary decision rule is
identical in both cases, and omission is recorded as `secondary-not-run`.

## Purpose and isolation

`prepare_holdouts.py` is a separate, holdout-only preparation path. It does not
modify or call the calibration-specific `prepare()` or `carve()` entry points,
which are intentionally hard-wired to the 25-member Protenix panel. It reuses
only hash-pinned protonation and fixed-core chemistry primitives. The original
calibration files and their hashes remain unchanged.

Preparation has no ORCA execution capability. It only emits coordinate and
input templates after every target-level check passes. The separate runner is
the only holdout component permitted to invoke ORCA. A fresh output tree
contains an `HOLDOUT_PREPARATION_INCOMPLETE` sentinel until the complete
requested set validates. Existing output paths are refused.

## Release gate

A calibration `result.json` is mandatory and has no default or bypass. Before
reading a holdout coordinate source or creating an output directory, the
preparer verifies its hash-bound calibration pins, preregistration, and
25-target preparation. It independently recomputes the 11/14 class coverage,
U, L, gap, midpoint, AUROC, and every leave-one-out midpoint classification
from the frozen R scores. All four preregistered gates, a gap of at least 5
kcal/mol, and internally identical released bands are required.

The primary operation is atomic: it always prepares both 1H4I and 4MAE or
leaves the whole tree incomplete. That four-leg primary preparation is
release-runnable. 6OC6 is added only with the explicit
`--include-secondary-6oc6` flag and remains a nonindependent secondary geometry
check that cannot replace either primary member.

## Raw-source selection and exclusions

Every row is read from the exact hashed `holdout_spec.tsv` source. Model `1`
and author chain `A` must resolve uniquely; there is no first-model,
first-chain, nearest-site, or cross-chain fallback. Before extraction the code
checks the source metal selector, element, occupancy, B factor, coordinate,
CCD PQQ graph, fixed roles, Asp/Arg contact, exact typed-donor set and
distances, water inventory, and complete source altloc inventory.

The protonation context retains the full standard-amino-acid chain A plus
exactly the selected PQQ and metal. Every source water, other chain, other
model, and noncore heterogen is excluded before protonation. For 4MAE this
explicitly removes every atom of `A:15P603`, including Ce-bound OXT at 2.747
A, and adds no replacement. Its eventual result must therefore be described
as a dry fixed-coordinate structural-transfer test with a ligand vacancy, not
as the intact crystallographic first shell.

Core atoms may not have alternate conformers. The remote 4MAE `A:ASP162`
alternate is resolved by the already pinned
`residue_consistent_highest_occupancy_blank_then_A_v1` policy (A wins its
occupancy tie); all discarded atoms and the choice are recorded.

## Metal normalization decision

The standard-only protonation wrapper requires one residue named `LA` and one
named `PQQ`. For raw Ca and Ce crystals, only the selected metal **residue
name** is therefore changed to `LA` before protonation. Its native atom name
and element remain Ca or Ce, and its coordinate is unchanged. This choice
minimizes alteration of the crystal during standard-residue H placement and
makes the authorized identity mapping explicit:

- 1H4I: `A:CA701/CA` -> `A:LA701/CA` (native Ca atom retained);
- 4MAE: `A:CE601/CE` -> `A:LA601/CE` (native Ce atom retained);
- 6OC6: `A:LA701/LA` -> `A:LA701/LA`.

Only after protonation and all coordinate checks does the carver replace that
single metal atom with the paired La and Ca QM arms. The arms use one common
protonated scaffold and therefore cannot differ in any nonmetal coordinate.
The source occupancy (notably 4MAE Ce occupancy 0.60) and B factor are asserted
and recorded as crystallographic provenance; neither is a QM-XYZ property.

## Chemistry and coordinate contract

PQQ must be the complete 24-heavy-atom `pdb_ccd_pqq_v1` graph with no source
hydrogen or altloc. The carver alone adds the three invariant hydrogens of
oxidized PQQ(3-), `pqq_ox_3minus_v1` (`C14H3N2O8`, singlet).

Standard amino acids use the unchanged deterministic calibration subprotocol:
pH 7.0, Python seed 20260914, OpenMM CPU with one thread, `forcefield=None`,
and no nonstandard CCD hydrogen definitions. Missing residues are not added;
any repaired/added heavy or terminal atom is fatal. PQQ and the metal must
receive zero hydrogens.

The fixed fragments and charges remain:

| Holdout | Included protein roles | La charge | Ca charge | Source CN |
|---|---|---:|---:|---:|
| 1H4I | Glu177, Asn261, Asp303, Arg331; Ala305 excluded | -1 | -2 | 6 |
| 4MAE | Glu172, Asn256, Asp299, Asp301, Arg326 | -2 | -3 | 9 |
| 6OC6 | Glu192, Asn276, Asp318, Asp320, Arg345 | -2 | -3 | 9 |

Fragment membership is role-defined, not distance-defined. The catalytic Asp
in 1H4I remains in the fixed core despite lying beyond 3.1 A. Arg is retained
as the fixed second-shell charge partner and never counted as a metal donor.

The preparer checks retained raw -> selected PDB -> protonated PDB heavy-atom
identity and coordinates at a maximum 0.001 A serialization tolerance. It
also compares every final fragment heavy atom directly back to the raw source.
Every generated PQQ H, protonation H, and C-beta link cap is recorded with its
origin, parent, and coordinate. The serialized La and Ca XYZ files must have
identical counts and byte-identical nonmetal lines; their metal coordinate is
identical and only metal identity, total charge, and electron count differ.

## Future invocation (only after review)

Run with the pinned Python and a fresh output directory:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/prepare_holdouts.py \
  --calibration-result /absolute/path/to/locked/result.json \
  --output /absolute/path/to/fresh/holdout_prepared
```

This prepares the primary release only. Add `--include-secondary-6oc6` solely
to request the optional secondary. The command prepares templates only;
submission and transfer scoring remain separate reviewed actions.

## Four- or six-leg execution and frozen scoring

`run_holdouts.sbatch` requests one high-memory node, invokes
`run_holdouts.py`, then invokes the frozen scorer against the receipt named with
that SLURM job ID even when one or more ORCA legs failed. The scorer
writes to a fresh result path (the default is
`reserved_crystal_holdout/result`); the batch script refuses an existing path
before ORCA begins. The runner accepts exactly the two-primary/four-leg release
or the same pair plus the optional secondary/six-leg release. It derives
simultaneous target count, leg count, MPI ranks per leg, assigned ranks, and
unused CPUs from that exact preparation and `$SLURM_CPUS_ON_NODE`. Each leg uses
at most 16 MPI ranks, matching the calibration cap. Per cluster policy, the
script specifies no memory, CPU count, wall time, or exclusivity. The runner
always writes its receipt before reporting per-target failure. The batch script
deliberately captures that status and invokes the scorer whenever the receipt
exists, so failed/unscorable primary arms still yield a durable FAIL
JSON/Markdown result. Catastrophic pre-receipt failures stop the job.

`score_holdouts.py` requires a receipt whose exact two- or three-target ledger
matches the preparation and independently validates every available task input, XYZ, ORCA output,
normal termination, SCF convergence, runner, renderer, executable, and hash
link. It derives the exact
frozen S bands from the locked calibration result and verifies their algebraic
identity to the frozen R bands and aquo gauge. It never fits a holdout cutoff.
The primary rule is 1H4I `S <= U_S` and 4MAE `S >= L_S`; a gap score, wrong-band
score, or unscorable arm fails and makes the scorer exit nonzero only after its
result files are written. If omitted, 6OC6 is explicitly reported as
`secondary-not-run`; if included but unscorable, it remains secondary-only and
does not make an otherwise passing primary verdict fail.

After reviewed preparation, the intended command is:

```bash
sbatch \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/run_holdouts.sbatch \
  /absolute/path/to/holdout_preparation.json
```

That one job performs execution followed by frozen scoring. It does not prepare
the holdouts and does not alter the calibration or its released bands.
