# Compact frozen-core exchange: source reversal persists

**All six endpoints are available, but this static local model is not a robust
Hans/Mex discriminator.** Both Hans structures must remain visible. The changed
frozen-4f electronic treatment alone does not resolve the earlier source reversal.

D(H,M)=[E_H(Dy)-E_H(La)]-[E_M(Dy)-E_M(La)], kcal/mol:

| Hans source versus Mex8FNS | D |
|---|---:|
| Hans8DQ2 | +11.330886345 |
| Hans8FNR | -18.468182860 |

Positive means greater conditional relative La preference in Hans. No absolute
La/Dy threshold, affinity, or EF3 biological label is assigned. These are consumed
fixed-core development structures, not independent blind controls. Both metal
identities have exact paired coordinates within each source; different protein
compositions remain balanced in the exchange. No aquo offset is imported.

## What the forces add

The predeclared projection moves only the metal toward its nearest real oxygen.
Positive energy gradient means that approach is uphill locally; force has the
opposite sign. Gradients in kcal/mol/angstrom:

| Source | nearest O distance, angstrom | La | Dy |
|---|---:|---:|---:|
| Hans8DQ2 |2.45247|27.93011|19.00414|
| Hans8FNR |2.29269|64.32003|50.26367|
| Mex8FNS |2.35144|-56.66472|-56.63721|

Both Hans sources resist this local approach, more strongly for La, with stronger
loads in the more compact8FNR source. Mex has a different force direction. This
supports examining structural accommodation; it does not determine its integrated
energy, repair source dependence, or justify optimizing a single favorable direction.
Only the earlier Dy8DQ2 directional derivative has been numerically qualified.
All-source raw gradients, geometry, state and translation residuals are retained in
COMPACT_FORCE_RESULT.json; no physical protein/cap force mapping is claimed there.

## Verified preparation limitation and decision

COMPACT_HYDROGEN_REPORT.md finds source-derived H-angle defects and retained
~1.17–1.20angstrom H bonds in these archived compact inputs. Caps preserve rather
than create the directional defects. Their influence on exchange energies is
unmeasured. Do not attribute the entire reversal to hydrogen preparation or claim
that correcting H will produce the desired answer.

Proceed with the separately versioned, fixed-heavy-atom hydrogen repair diagnostic,
then compare matched repaired origins if admitted. Preserve this matrix unchanged.
The full-region embedded donor-response1220312 remains running on its original
normalized-bond-length inputs; it too contains neighboring H-direction strain.
No production change, broad classifier claim, or repeat of failed explicit-f jobs.

## Actual execution

Job1220308 completed five new analytic endpoints in1809seconds on40allocatedCPUs;
collector1220309 took1second on1CPU:72,361allocatedCPU-seconds, zeroGPU. One exact
Dy8DQ2 origin was reused, with its prior cost recorded in FROZEN_F_SCOUT_RESULT.md.
All electronic states, energies, coordinates and native gradients passed collection.
No failed endpoint was silently substituted. The timing includes energy and forces;
it is not a production screening benchmark. Slurm batch MaxRSS8,013,844KiB is a
reported step metric, not evidence of full-node peak memory or measured utilization.

Reproduce the force analysis with scripts/metal_environment_compact_forces.py and
COMPACT_FORCE_DESIGN.json against the actual FINAL_COLLECTION.json. The compact
manifest and collection are pinned in COMPACT_EXCHANGE_RESULT.json and workspace.
