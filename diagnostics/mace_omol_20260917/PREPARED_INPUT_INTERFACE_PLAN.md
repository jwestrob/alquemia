# Reusable interface for the qualified whole-chain descriptor

Declared 2026-09-17 while canonical job1200830 runs. Engineering within Jacob's
active goal; no new scientific model, label, geometry selection or energy run.
The baseline/default and running immutable manifests remain unchanged.

Add an opt-in operation to score a recorded whole-chain preparation through the
existing runner. Input is an explicit `mace_global_prepare` preparation.json with
source records, paired coordinates, atom mapping, protein/cofactor formal charge,
spin, water inventory, and provenance. Arbitrary XYZ files are insufficient.
The interface does not invent missing protonation, ligands, sequence or structure.
It initially supports the same single selected Ca/La site and chain-A policies;
multiple sites/assemblies remain unsupported, not silently reduced to one site.

Audit the input by replaying the existing protein and cofactor/water preparation
functions with the recorded sources, checking peptide connectivity and requiring
exact agreement with the recorded physical atoms/charges/paired coordinates.
Check every source pin and the qualified descriptor's checkpoint/adapters.
Create exactly four bound/detached endpoint tasks with the existing versioned
charge mask, geometry convention, finite manifest, immutable snapshot and cache.
Use the existing dry-run, allocated execute and receipt collection operations.

Report unrounded paired descriptor components in explicitly nonphysical model
units. Missing endpoints or references remain unavailable. A later compatible
passing canonical report can provide PQQ-only research bands; no generic affinity
zero or PQQ band for non-PQQ inputs. Default production scoring is unchanged.

Initial tests only replay already consumed real 1H4I and GGR inputs and reject
corrupted copies / the actual broken1KB0 input. Compare exact prepared tasks with
completed job1200828; do not launch duplicate forwards for interface validation.
No new model inference, DFT, solver, training, dynamics or structure prediction is
included in this engineering phase. Keep a usable raw-score interface even if the
canonical calibration fails, with explicit unvalidated/unavailable decisions.
