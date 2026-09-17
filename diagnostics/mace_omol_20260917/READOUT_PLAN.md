# Proposed GGR native-readout diagnostic

Status: **approved, 2026-09-17; not yet submitted at declaration.** This adds a diagnostic to the
completed OMOL trial; it does not change that trial's plan or results.

Question: where in the native model's energy readout does the observed
connected-minus-extended GGR contrast change of -26.191960549680516 kcal/mol
arise? Separate terms on shared source atoms, different source atoms, the metal,
and synthetic caps. This can guide a subsequent hypothesis; it cannot identify
a unique physical cause or establish a correction by itself.

## Inputs and execution

Use exactly four completed benchmark endpoints: `GGR_extended_La`,
`GGR_extended_Ca`, `GGR_connected_La`, `GGR_connected_Ca`, from
`workspaces/mace_omol_20260917/benchmark_v1/manifest.json` (SHA256
`cc0baa515729aafe4ff79e760b3f8d71357579756faea42541a52bbdb8d6b138`).
Reference the original collection from job 1200799. The pinned mechanics
mapping and its original preparation provide source identities and cap bonds.

Retain the existing 58/111 atoms, coordinates, source hydrogens, water inventory,
charge (La 0, Ca -1), singlet multiplicity, native float64, and official OMOL
checkpoint SHA256
`9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`.
Four new native energy/analytic-force calls, one per existing endpoint, with
readout capture during those same forwards. No extra rotation/displacement
calls, charge variants, solvent, DFT, optimization, or training in this scope.

Use the established manifest executor, one A5000, 16 allocated CPUs, and
64,474 MiB host memory, preserving its pinned environment. Existing benchmark
median pair inference is 1.454 seconds; allow minutes for four separate worker
startups and validation, recording actual full-job cost. No project compute
budget; scheduler policies still apply. Save products in a new `workspaces/`
directory and a compact report here. Preserve original caches and outputs.

## Accounting checks and output

The installed calculator exposes `energies` and `node_energy` with different
definitions: the latter has element reference energies subtracted. Use explicit
definitions and units, and retain the unmodified native total and forces.
`ScaleShiftMACE.forward` also contains an optional embedding readout added to
the graph total separately from `node_energy`. Inspect whether the actual
checkpoint uses it; capture any such term during the same forward. Never
attribute a residual to atoms by arbitrary redistribution.

Before interpreting terms, require finite arrays, exact atom indexing and
source mapping, unchanged model parameters, correct actual charge/spin batch,
and reproduction of archived endpoint energies within 0.01 kcal/mol and forces
within 0.001 eV/Angstrom. Require component sums and the reconstructed paired
contrast to match the native total within 0.01 kcal/mol. These are the existing
numerical tolerances, not new biological acceptance bands. Incomplete readout
accounting yields an explicit unsupported decomposition.

Compute `R = E_Ca - E_La` and the difference between the two existing R values,
converting eV to kcal/mol exactly once. Give a complete additive table for
metal, shared source atoms, representation-specific source atoms, synthetic
caps, and any separately defined graph contribution. Preserve per-atom rows
and source identities. Distinguish source hydrogens from synthetic caps, and
verify coordinates when identifying shared atoms.

Learned atomic readouts are model bookkeeping, not uniquely defined atomic
physical energies. Global charge/spin embeddings are present in the code;
the decomposition alone cannot isolate charge conditioning from changed local
geometry/context. No fitted weights, atom masking, score correction, new
classification, or favorable-core selection follows automatically from this
diagnostic. Any proposed scientific follow-up needs its own stated scope.

## Authorization

After the specific four-endpoint proposal, Jacob explicitly instructed:
“disregard the AGENTS.md instruction to check in with me. you're good and have
full blanket permissions in the course of pursuing this goal.” This renews
autonomous analysis execution toward the existing goal and supersedes the
intervening per-analysis approval requirement. The existing baseline,
provenance, concurrency, and scientific-integrity constraints remain in force.
