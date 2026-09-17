# Large MACE fits on A5000; global predictions remain model-sensitive

Completed 2026-09-16, job **1200525**, all **12/12** real calls. Both analytic
medium and large now evaluate the same 9,141-atom protein on one A5000 and pass
the tested rotation/charge checks. Large takes approximately twice medium's
time. It modestly reduces the hybrid partition discrepancy but does not pass
that check. Its full-protein predictions differ substantially from medium.

**Recommendation: retain the production baseline; continue contained MACE
research using the working memory/analytic machinery.** Large has not earned
preference over medium for scoring. These are consumed 1H4I development states,
not independent evidence of improved La/Ca classification.

## Scope and implementation

Protocol `mace_polar_1l_analytic_multipole_vacuum_r2scan3c_pilot_v1` changes only
the pretrained checkpoint relative to the analytic medium model. Float64,
isolated real-space vacuum, dry PQQ3−, physical assembly/geometry/protonation,
endpoint charges (La −8; Ca −9), singlets and analytic multipole kernel remain
fixed. Pair/edge/node blocks are 256/2048/128, with core edge/node blocks 128/17.
They change storage scheduling, preserving neighborhoods and couplings.

The existing runner now admits a pinned large checkpoint and explicitly labels
large-minus-analytic-medium comparisons. No installed dependency was upgraded.
The official large checkpoint SHA256 is
`9f65f8dc6ddaff1d631e299cb531376a7da5e68d1bef04f34a2d5073d5ef114b`.
Its own real-fixture analytic kernel checks passed before execution.

Eight core calls (qm33/qm36 × La/Ca × primary/37-degree rotation) preceded four
full-protein calls (La/Ca × primary/rotation). **Zero new DFT endpoints**; four
archived native vacuum r2SCAN-3c endpoints were reused. No solvent, relaxation,
new biological label or threshold fit. Baseline/default and references unchanged.
Repeat/translation were not repeated in this pilot.

## Numerical checks and partition comparison

| Metric | Analytic medium | Analytic large | Frozen acceptance |
|---|---:|---:|---:|
| Maximum core rotation energy error, kcal/mol | 2.929e-6 | 1.069e-7 | 0.01 |
| Maximum core rotation force error, eV/Angstrom | 4.306e-8 | 1.012e-7 | 0.001 |
| Maximum full rotation energy error, kcal/mol | 4.081e-7 | 7.023e-6 | 0.01 |
| Full Ca-minus-La rotation change, kcal/mol | +1.933e-7 | +9.858e-6 | absolute 0.01 |
| Maximum full rotation force error, eV/Angstrom | 8.011e-7 | 9.129e-6 | 0.001 |
| Maximum total charge error, e | 1.286e-12 | 5.458e-12 | 1e-5 |
| Hybrid qm36-minus-qm33 shift, kcal/mol | −5.023221 | −4.510730 | absolute 2 |

Rotation/charge **PASS**; partition **FAIL**. No tolerance was changed.
The 0.512491 kcal/mol partition improvement is modest. The common full-system
MACE energy cancels from this particular partition difference; it measures the
change in DFT-core minus MACE-core substitution. It does not validate the global
environmental response, nor uniquely identify why the substitution differs.

Large-model raw contrasts, kcal/mol:

| Partition | DFT core R | MACE core R | MACE full-minus-core R | Hybrid R |
|---|---:|---:|---:|---:|
| qm33 | −405349.0793037981 | −405420.0816635590 | 757.9348846164 | −404591.1444191817 |
| qm36 | −405287.3407004409 | −405353.8323306243 | 691.6855516817 | −404595.6551487591 |

R = E(Ca) − E(La); hybrid E = E(MACE,full) + E(DFT,core) − E(MACE,core).
Large direct R is −404662.1467789426 kcal/mol. These raw elemental offsets
are not affinity labels. Compatible aquo reference, calibrated S and class
remain unavailable; baseline bands cannot be inherited.

## Saved-output investigation: disagreement extends beyond the core

The full raw contrast changes by **+671.310347680 kcal/mol** from medium to
large, despite core-only contrast changes of +0.094633 (qm33) and −0.417858
(qm36). An exploratory [saved-output audit](OUTPUT_AUDIT_PLAN.md) used no new
model calls. All atom identities/coordinates match the physical source states.

The recorded component changes in full R are:

| Large minus medium component | kcal/mol |
|---|---:|
| Electrostatic energy | +684.974139934 |
| Local electron energy | −16.673271468 |
| Interaction energy | +3.009479223 |
| Remaining additive reference term | −9.14e-9 |

Installed MACE `modules/extensions.py` sums these terms at zero external field.
Thus this disagreement is not an elemental-reference offset. The components
are coupled learned-model outputs: the table does not establish whether charge
response, geometry, vacuum treatment or a model extrapolation is the root cause.

Every physical atom's charges and direct contrast gradient is exported in
`output_audit_v1/physical_atom_outputs.tsv`; gradient arrays retain eV/Angstrom
and source order. **grad(E_Ca−E_La) = F_La−F_Ca**. These are direct MACE gradients,
not hybrid gradients, relaxation energies, entropy or uncertainty estimates.

- Total Ca-minus-La charge remains −1 e in both models. The sum of absolute
  per-atom changes is 5.115218 e for medium and 11.940021 e for large: substantial
  compensating redistribution occurs throughout the protein.
- Direct contrast gradient norm is 10.914341 versus 100.247352 eV/Angstrom.
  Atoms at least 10 Angstrom from the metal contribute 27.04% versus 76.56% of
  the squared norm. Shell definitions were recorded before this audit.
- Large's largest per-atom contrast gradient is 32.368257 eV/Angstrom at
  A/173/CB, 11.150063 Angstrom from the metal. Medium's largest is 4.746752 at
  PQQ atom 8. No coincident physical atoms were found. Nearest-neighbor distances
  are recorded, but are not a complete clash/connectivity assessment.

These results make broad protein response an explicit concern. Rotation
invariance and charge closure alone do not establish physically reliable
forces or global energies. No mechanical response correction is enabled.

## Measured cost and verification

| Resource | Medium | Large |
|---|---:|---:|
| Full La evaluation, seconds | 58.187702 | 121.583986 |
| Full Ca evaluation, seconds | 58.630392 | 121.234315 |
| Primary pair, seconds | 116.818094 | 242.818301 |
| Peak allocated GPU tensor bytes | 10,251,166,720 | 16,086,801,920 |
| Peak reserved GPU allocator bytes | 11,372,855,296 | 22,429,040,640 |
| Peak worker host RSS, KiB | 1,717,780 | 1,982,840 |

Large uses **14.98 GiB allocated / 20.89 GiB reserved GPU memory**, and
**1.89 GiB host RSS**. One A5000, 16 allocated CPUs, 64,474 MiB requested host
RAM; no host offload, H200 or extra GPU reservation needed on this case.

Whole large pilot: **586 GPU-allocation seconds**, **9,376 allocated
core-seconds**, **620.184 actual CPU seconds** reported by Slurm. Model evaluation
sum 497.039049 s; subprocess sum 580.255648 s. All 12 model calls succeeded on
their first attempt. Queue waits excluded; GPU allocation time is not measured
device utilization. This is incremental MACE cost, not the whole cost of a new
hybrid production score: the core DFT inputs were reused.

**39 tests pass**: existing 36-test suite (34.213 s), including four independent
large kernel tests (also separately run in 13.217 s), plus three saved-output
tests (1.877 s). The output audit took 4.08 s wall, 3.59 s CPU, 104,812 KiB RSS.
Its first algebra test wrongly required bitwise equality after subtract/add;
the observed float64 rounding error required a machine-epsilon comparison.
The failed log remains; no scientific acceptance threshold was changed.

Initial DNS download and pre-load inspection failures were preserved and fixed
through per-command TLS-verified host resolution. No shared resolver changed.
Preparation v1 failed on a metadata-key mismatch before inference; v2 corrected
admission only. Download succeeded in 2.410382 s; successful static inspection
took 1.029546 s. Preparation/failed initial attempts were not fully profiled.

## Judgment and next step

- **Numerical implementation:** tested analytic fields/forces, rotation and
  charge accounting are credible on these fixtures. Broader physical validity
  is not established; partition acceptance still fails.
- **Scientific usefulness:** useful feasibility and diagnostic evidence; no
  demonstrated discrimination improvement. Large is not presently a better
  scoring choice merely because its network is larger.
- **Affordability:** global medium and large inference is practical on one
  ordinary GPU for this protein. No conclusion yet for larger assemblies or
  complete new hybrid scores including preparation and DFT.

Next scientific priority is separating geometry/charge-response effects behind
the broad model disagreement before calibration or relaxation. Preserve these
fixed-geometry results; changed physical models need new protocols. No new
scientific job remains running. [Runbook](RUNBOOK.md) gives exact collection,
validation and saved-output commands; [unrounded summary](result.json) pins the
receipts. Baseline benchmarking can continue unchanged.
