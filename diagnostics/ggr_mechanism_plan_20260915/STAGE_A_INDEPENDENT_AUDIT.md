# Stage A final-input integrity audit

Read-only audit of the five pinned final preparations. No endpoint energies were inspected, no preparations or geometries changed, and no calculations launched. Exact preparation/XYZ hashes and parent comparisons are retained in `STAGE_A_INDEPENDENT_AUDIT.json`.

## Findings

- NMA extensions add only C/H atoms; no new typed donors.
- The connected GGR model adds six typed oxygen donors, all beyond the frozen 3.3 Å inclusion radius:

| Added source atom, chain A | Metal distance (Å) |
|---|---:|
| LYS137:O | 6.376739 |
| ASP138:O | 5.604942 |
| GLY139:O | 6.020902 |
| GLN140:OE1 | 7.815366 |
| ILE141:O | 5.899066 |
| GLN142:O | 6.332113 |

| Preparation | Shortest eligible cap contact (Å) | Parent (Å) |
|---|---:|---:|
| ggr_1glg_nma | 2.471388 | 2.396957 |
| aequorin_1sl8_EF3_nma | 2.486517 | 2.096599 |
| alacta_1f6s_nma | 2.238837 | 2.088290 |
| alacta_6ip9_nma | 2.307033 | 2.358765 |
| ggr_1glg_connected | 2.332224 | 2.396957 |

For cap contacts, the graph includes retained source heavy bonds, source hydrogen-parent bonds and each cap-parent bond. We exclude the cap itself, bonded 1–2 neighbors and angle 1–3 neighbors. We do not exclude 1–4 contacts or infer metal-coordination bonds. Distances include every remaining QM atom, including the metal.

No obvious cap overlap was identified. Source coordinates, source mappings, overlap merging, charge −3, waters and paired endpoint states remain consistent with the approved plan. Review of actual Stage A task packaging found no material state/provenance issue.

## Limits

These contact distances are an integrity check, not a force-field clash score or validation of cap forces. The audit cannot establish predictive correctness, protein response or physical convergence. The extended model changes local interactions and continuum cavity even without changing direct coordination.
