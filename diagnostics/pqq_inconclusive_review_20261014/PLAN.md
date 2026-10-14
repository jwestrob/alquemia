# Inconclusive benchmark cases: saved-result diagnosis

Jacob requested investigation of the inconclusive results after the completed
fixed-charge candidate. This is a development diagnosis, not a new calibration.
Production and all previously recorded predictions remain unchanged.

## Questions

1. Which inconclusives are near a correct-class boundary versus well inside the
   gap? How do all available same-protein structural peers behave?
2. Is ambiguity already present at the origin, or does accommodation move a
   source into the gap? Compare origin and selected scores using the same pooled
   bands as a diagnostic, while retaining separately calibrated origin decisions.
3. Which changes come from the MACE proposal energy and native solvent component?
   Reuse the completed model/selection/band decomposition; do not rerun chemistry.
4. Do actual donor geometry, cofactor geometry, local confidence or optimizer
   bounds distinguish the ambiguous cases from their same-state peers?

## Population and controls

Use the exact completed candidate result, SHA256
`82981de62dfa2f2ff5eaedd273e38d8b20bf36034dc39110c104908599470d25`, its immutable
inventory and native component collections. Include every inconclusive in the
expanded135 and current recovered125 cohorts. Deduplicate identical source and
selection identities while retaining every cohort alias. Retain the five
inconclusive existing triples and their exact members; do not create new votes.
All available same-protein peers from those cohorts are controls, including
different source-conditioning metals, reported separately by physical state.
Original unavailable inputs remain explicit and are not new inconclusives.

The expected union is nine source contexts across XoxF1 A0ACD6B9F2, XoxF4-2
MMOL_1770, ExaF C5AXV8, ExaA Q9Z4J7 and PedE Q88JH5. Recover exact identities
and state mappings from existing records. Source-conditioning metal is not the
experimental reference class. These are previously inspected development cases.

## Analysis and limits

Keep unrounded scores, bands, both selected geometries and individual metal
works. State distances to the correct and opposite class bands, without calling
them probabilities or error bars. Display all three candidate geometries where
helpful, but do not pick a favorable geometry or switch scoring policies by label.
Separate source-to-source composition differences from motion within a fixed
source. No across-source absolute endpoint minimum is allowed.

Reuse existing typed-contact/PQQ/source-mapping/confidence tools and archived
optimizer results. Report missing metrics as missing. Use their existing donor
cutoffs and geometric definitions; do not invent a score-dependent geometry
acceptance threshold. A high pLDDT/iPTM does not prove local ion geometry, and
neither disagreement nor a short contact alone identifies a faulty fold.

All numerical analysis runs on Slurm, one CPU per analysis, with explicit finite
input pins and completion wakes. No new MACE, xTB SCC, DFT, folding, optimization,
threshold fitting or production promotion in this diagnosis. Root owns execution
and score interpretation; separate agents handle geometry and source identities.
The deliverable is a case-by-case finding and prioritized physical next test,
not a demand to make every ambiguous case cross a band.
