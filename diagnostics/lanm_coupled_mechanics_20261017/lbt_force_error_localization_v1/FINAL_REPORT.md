# Shared scaffold discrepancy and localized metal-response discrepancy

20 October 2026. Saved-only Slurm job **15797** completed all eight endpoints
and four matched La/Dy pairs. The broad force discrepancy shared by both metals
lies throughout the peptide's heavy atoms. The smaller metal-dependent
discrepancy is particularly concentrated in specific side-chain oxygens in 7CCN.
This separates two mechanics issues to investigate; it does not identify which
electronic model is wrong.

## What was measured

At identical coordinates, define the per-metal error as MACE minus PBE0.
The **shared error** is the mean of the La and Dy error vectors. The
**metal-dependent error** is the Dy error vector minus the La error vector.
The same saved ALPB term is used on both sides of each model comparison and
therefore cancels from these errors. Solvent accuracy is not tested here.

All 256 atoms remain included: one metal, 119 hydrogens, 69 named backbone heavy
atoms and 67 remaining protein heavy atoms. Existing physical IDs determine
the categories and all residue groups. No distance cutoff or selected region
was introduced. Forces are in kcal/mol/Å; RMS is over atom-vector magnitudes.

| Source / geometry | Shared error RMS | Metal-dependent error RMS | Largest metal-dependent atom error |
| --- | ---: | ---: | --- |
| 7CCO before | 3.511769 | 0.501500 | Metal, 2.704859 |
| 7CCO after | 3.516169 | 0.476265 | Asp7 OD1, 2.259919 |
| 7CCN before | 3.513145 | 0.726763 | Asp7 OD1, 6.242802 |
| 7CCN after | 3.581379 | 0.563471 | Asp3 OD1, 4.148825 |

These are the original 7CCO hybrid-La deformation and 7CCN whole-MACE-La
80-step deformation. They are two consumed examples of the same peptide.
Neither is the later 7CCN step725 endpoint.

## Where the discrepancies lie

The shared error is consistently broad: backbone heavy atoms carry
52.0–52.8% of its squared norm, and remaining protein heavy atoms carry
43.8–44.7%. Hydrogens carry 3.2–3.5%, and the selected metal less than 0.13%.
Trp9 has the largest residue contribution in every pair, 12.7–13.2%; its
CG atom has the largest shared atom error in all four. This distribution
does not support attributing the shared discrepancy only to the metal or
terminal atoms.

The distribution changes for the metal-dependent error:

| Source / geometry | Metal | Hydrogens | Backbone heavy | Remaining protein heavy |
| --- | ---: | ---: | ---: | ---: |
| 7CCO before | 11.36% | 4.34% | 37.91% | 46.38% |
| 7CCO after | 8.35% | 5.03% | 39.30% | 47.32% |
| 7CCN before | 1.02% | 1.46% | 11.52% | 86.01% |
| 7CCN after | 4.40% | 2.37% | 22.40% | 70.83% |

Percentages are fractions of the **squared force-error norm**, not energy
contributions, affinity errors, or fractions of a causal mechanism.

For 7CCN before, the largest atom discrepancies are Asp7 OD1 (6.242802),
Glu14 OE2 (4.788395), and Asn5 OD1 (4.684363). Their residue contributions
are respectively 31.57%, 22.27%, and 18.15%. After the original 80 steps,
Asp3 OD1 (4.148825) leads, followed by Asp7 OD1 (3.038746) and
Asn5 OD1 (2.586347). The largest discrepancy changes location; the recorded
motion has not eliminated metal-dependent force disagreement.

7CCO is less concentrated: backbone and remaining heavy atoms both contribute
substantially, and the metal itself accounts for 8.35–11.36%. The two structures
therefore do not support a universal single-atom correction.

## Decision and limits

The useful mechanics distinction is now clearer:

- Absolute scaffold forces have a broad discrepancy that largely cancels when
  comparing the metals at identical coordinates.
- The remaining metal-response discrepancy includes specific oxygen atoms and
  depends on the source and deformation. It warrants attention even when the
  total deformation energies agree closely.

This supports investigating electronic forces and their reference model,
while retaining the available structures. It does not justify fitting an
oxygen-specific correction, changing geometry to improve a score, or extending
minimization without checking what the reference says.

The already-running **15545** tests the whole-peptide La energy and forces at
the model-converged step725 endpoint. Collect that endpoint against the saved
step92 reference to assess finite work and residual reference forces separately.
Do not infer matched La/Dy convergence or a natural-protein preference from this
single La endpoint.

MACE's physical Dy multiplicity6 and PBE0's effective-singlet frozen-4f
representation differ. Native SCF/state and derivative limitations, conditional
termini, the vacuum-dianion assumption, and post-added solvent remain.
Neither model is established here as the physically correct force field.
No force-accuracy threshold, solvent validation, production promotion or
natural-system transfer is established.

The inherited 7CCN caveat text saying Dy MACE forces were unavailable describes
the older reference collection. Job15720 filled those missing arrays; it did
not remove the reference's physical limitations.

## Verification and preservation

All 47 input pins and eight saved admissions passed. Slurm performed
3,106 bookkeeping checks, including matched-vector replay within the
declared floating-point tolerances and the shared/difference squared-norm
identity; all passed. Root and independent review inspected the recorded checks
and immutable input bindings. There were zero molecular calls, model loads,
new geometry evaluations, optimizer steps or fitted parameters.

Scheduler cost: **1 allocated CPU-second, zero GPU-seconds**. Application elapsed
time was 0.286792s; measured TotalCPU was 00:00.355, a separate accounting field.

- Result SHA256:
  `c789abb0071a10b1b2dc06db2e24b1c087c44a15cd3fc2682add532cf24ca7d1`.
- Parent15720 result: `f42cffcf`; inventory: `e2c196d1`.
- Source/config/spec: `41502735` / `7297b790` / `588a423d`.
- [Independent actual review](ACTUAL_REVIEW.md), `17271e3f`.
- [Lossless execution archive](actual_result_v1/ARCHIVE.json), `aae39987`.
- [Supplemental root and independent reviews](actual_result_v1/SUPPLEMENTAL_REVIEWS.json),
  `47d3e0ea`; the base archive remains unchanged.
- Original execution:
  `workspaces/lanm_coupled_mechanics_20261017/lbt_force_error_localization_v1/run/`.

The archive preserves all per-atom vectors, all residue summaries, fixed top-ten
displays, checks, source bindings and scheduler evidence. Root's proactive
result review is separate from acknowledgment of a delivered notification.
