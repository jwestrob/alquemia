# Internal proton transfer: no competitive state on these frozen paths

**All 24 native r2SCAN-3c/CPCM analytic-gradient calculations completed normally.
Every sampled point is uphill from its neutral-water/carboxylate starting state.**
This experiment supplies no proton-transfer correction to promote. It does not
rule out a separately relaxed protonated site or establish pH populations.

| Source contact | Ca transfer cost | La transfer cost | Ca−La difference |
|---|---:|---:|---:|
| 1F6S water211 → Asp87 OD2 | +24.434540 | +20.796316 | +3.638224 |
| 6IP9 water310 → Asp82 OD1 | +62.860417 | +89.970689 | −27.110271 |
| 6IP9 water310 → Asp88 OD1 | +32.257900 | +37.962818 | −5.704918 |

Units are kcal/mol; these are electronic changes at the acid-H endpoint. All
intermediate values, gradients and components appear in [the tables](export_v1/TABLES.md)
and [unrounded record](export_v1/result.json). Positive Ca−La difference means
relative stabilization for La, but a high-cost chemical state cannot simply be
added to the discriminator because its relative shift looks useful.

## What the test establishes

We selected all three existing water/carboxylate-O contacts <=2.60 Å in the
full-water contextual preparations before computing proton-state energies.
Both structures are bovine alpha-lactalbumin: one biological group, not a new
alpha-versus-GGR affinity comparison. The four actual neutral endpoints were
reused. Four path fractions per contact and metal gave 24 new calculations.

One real H moves continuously along a donor-centered curved path to an acid O-H
distance of 0.98 Å. Every other atom remains fixed. The curved coordinate avoids
interpolating through water O when the initial H points away. The state change
`COO− + H2O → COOH + OH−` preserves atoms, charge and solvent inventory, so proton
and water reservoir terms cancel exactly. No H+ solvation energy, assumed G=E,
entropy offset or threshold fit was introduced.

The source water-H preparation differs between Ca and La; those spectator
coordinates remain explicit. Consequently these are comparisons between actual
prepared conditional paths, not pure electronic swaps at identical coordinates.
The CPCM dielectric component is retained inside SCF energy, not added twice.

All six source gradients point downhill for an infinitesimal displacement along
their paths (−15.37 to −31.47 kcal/mol per path fraction). All points at fraction
0.25 are already uphill. This suggests short-range H bond/angle response near
the neutral arrangement, distinct from the tested full proton transfer. It does
not quantify a relaxation correction or show that it improves discrimination.

## Important geometric limit

The 6IP9 Asp82 acid-H endpoint places H only **1.909 Å from the metal**; the other
two contacts place it at 2.908 and 2.289 Å. The conservative predeclared nonoverlap
gates passed, but do not establish a relaxed acid/OH geometry. The opposite metal
shifts therefore cannot be interpreted as universal intrinsic proton affinities.
Actual fixed orientation/coordination costs are part of these path energies.
[Geometry audit](GEOMETRY_DIAGNOSTIC.json) was descriptive, after outputs; no
selection, criterion, contact or denominator was changed to repair a result.

No path establishes a stable transferred-proton basin. No new protonation state
is selected for production. Populations, pKa, bound-state vibrational/quantum
corrections, nonpolar contributions and relaxed-state selectivity remain unavailable.
The failed wider-water entropy approach was not reused as a thermodynamic shortcut.

## Execution and cost

Job **1202444** used four concurrent 16-rank endpoints on node-64-768g-25:
**2418 wall seconds, 64 allocated CPUs, 154752 allocated CPU-seconds
(42.9867 CPU-hours), zero GPU-seconds**. Reported total CPU time was 146548 seconds;
batch MaxRSS was 13428488 KiB. Local preparation/tests/collection are additional
unmetered work. This is a development experiment, not a proposed routine scanner
cost. [Scheduler receipt](COSTS_sacct.tsv).

Slurm reports the batch FAILED because its original post-calculation collector
omitted a dependency. All 24 ORCA processes terminated normally with successful
receipts; no scientific retry occurred. Separate analyzers corrected that packaging
issue and reused the existing fixed numeric energy parser (the older copied parser
mistook diagnostic ellipses for numbers). Final analysis_implementation_v4 collected
all actual energies/gradients without changing or rerunning inputs. Earlier failed
and partial collections remain preserved.

Eight real-fixture tests pass with zero skips: physical coordinate mapping, exact
atom/charge/scaffold conservation, all 24 actual output/gradient receipts, partial
job accounting, and direct four-energy Ca−La algebra with one unit conversion.
Actual calculations and parser/geometry tests are distinct. Protocol:
`native_r2scan3c_fixed_scaffold_internal_proton_transfer_v1`.

**Recommendation:** keep this proton-transfer model as a negative, bounded research
result; do not add its shifts to the scanner or launch an automatic rescue. The
useful immediate work is the separately authorized independent parvalbumin test of
the already-demonstrated fixed-inventory water-H preparation. PQQ, original score
records and the protonation model remain unchanged. [Runnable commands](COMMANDS.md).
