# Finite two-Ca occupancy refold

12 October 2026. Root-authorized recovery of the fifteen off-site Ca predictions.
Original sources, original125 outcomes and frozen corrected canonical bands stay
immutable. No model ranking selects or excludes a generated sample.

## Scientific distinction

The ExaA crystal structure 1FLG maps to Q9Z4J7 and contains two protein chains,
two PQQ molecules and four Ca ions. Its primary paper explicitly describes an
additional N-terminal Ca site contributing to stability, beyond the PQQ-site Ca.
This directly motivates testing two-Ca occupancy for ExaA. For PedE Q88JH5 and
ExaF C5AXV8, two-Ca occupancy remains a hypothesis based on related proteins and
the observed distant six-donor region; it is not experimentally established here.
Sources: https://www.rcsb.org/structure/1FLG and
https://doi.org/10.1006/jmbi.2000.3603.

## Frozen inputs and execution

Use the three exact original cached `*_data.json` files. Preserve every existing
field, including name, protein, modifications, paired/unpaired MSAs, templates,
PQQ SMILES, bondedAtomPairs, userCCD, dialect/version and modelSeeds=[1]. Append
exactly `{"ligand":{"id":"D","smiles":"[Ca+2]"}}`, copied from ligand B except
for ID. Do not change PQQ to a different charge representation at fold time.
New source/state IDs and output paths explicitly distinguish occupancy2Ca_v1.

Use the actual ADH9 AlphaFold3 3.0.4 container and weight SHA256 pins, five
diffusion samples and ten recycles, inference only, no data-pipeline/MSA search.
These are explicitly pinned new execution settings. The original source folds'
software version and recycle settings have not been verified, so comparison with
those predictions cannot isolate Ca count as the only execution difference.
First run Q9Z4J7's five samples as an executable scout. C5AXV8 and Q88JH5 require
the real completed five-model scout receipt and a separate root launch. No
automatic seed retries, extra samples, ranking-based selection or larger batch.
All fifteen samples and every failure remain in the final source ledger.
GPU scheduling/resource practicality, including the current12GB host-RAM cap,
belongs to root. No memory/time requests, policy bypass or submissions here.
The GPU wrapper requests twelve CPUs and one GPU, as scheduler-test accepted by
root; GPU numerical thread variables use those twelve allocated CPUs and process
affinity is capped to that allocation. CPU preparation stays at one thread.
The GPU receipt retains actual affinity and scheduler resource variables.
Large runtime hashes
are streamed in8MiB blocks instead of loading the container or weights into RAM.

## CPU preparation and reference mapping

On Slurm, verify pins, append D and prove that removing only D restores the
original input object exactly. Validate each input using the installed container's
`alphafold3.common.folding_input.Input.from_json`; no model is loaded. Record
input/MSA/template hashes and upstream schema acceptance.

Download the explicit 1FLG mmCIF once from RCSB, retain its raw bytes/hash, and
audit all four Ca ions with the existing3.1A typed donor and4.0A PQQ-association
definitions. Map each crystal protein chain to the exact cached Q9Z4J7 sequence
using Gemmi global alignment (match2, mismatch−1, gap-open−6, gap-extension−1,
no free gaps). Retain full alignment, residue mapping and mismatches. Report
the crystal identities corresponding to Q9 Asp45/Thr48/Asp51 and their actual
Ca contacts. Mapping is descriptive evidence, not an admission-score threshold.

## Later generated-source admission

Retain both Ca ions in each full source. Identify the scoring site only when
exactly one Ca associates with the unique schema-valid PQQ within4.0A; ambiguous
or absent sites remain unavailable. Report the other ion as the spectator with
its actual contacts; never silently delete or move it. Existing typed donor,
homologous-core, charge/microstate, dry-water, cap, frozen4.3A context and state
compatibility gates remain mandatory before any later score. The old protein
protonator accepts exactly one metal, so a separate explicit spectator-aware
adapter is still required; this refold wrapper does not claim prepared endpoints.
Use the existing corrected canonical score bands only, with no new fit.

All molecular parsing, input preparation, schema validation and mapping run on
Slurm. Login work is code/metadata/byte pinning only. Root owns execution, wakes,
admission and scoring. This package creates no production changes.

## Technical successor v3

CPU job2877 stopped before molecular imports because login and node5 have different
node-local Apptainer executable bytes. The frozen v2 attempt is retained. This
successor records the actual scheduler node launcher path, hash and version; it
does not equate that binary with the login launcher. Shared AF3 container, weights,
scientific inputs and model settings remain byte-pinned and unchanged. No molecular
work from2877 is repeated because none began.
