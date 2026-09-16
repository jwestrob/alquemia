# ADH9 protein-hydrogen recovery — 16 September 2026

## Completed result: metal assignment is sensitive to model/preparation — 16 September 2026

**Job1200182 completed/PASS in4m59s; the repaired preparation and both ORCA endpoints pass independent checks.** The same Hydrogenophilalia `4380_6` sequence scores **21.450083 kcal/mol with Protenix versus9.720464 with AF3 plus repaired protein hydrogens**, a **−11.729619 kcal/mol** difference. The latter lies in the existing Ca-supported band, but this is **not a robust calcium assignment**: the structural preparation changes the classification. Biological metal preference remains unresolved. Predictor geometry and hydrogen preparation have not been isolated as separate causes; no new controlled comparison was launched.

| Protein | Protenix S | Accepted AF3 S | Existing-band position of AF3 result | Biological interpretation |
|---|---:|---:|---|---|
| Rokubacteriales `32301_3` | 15.503190 | 19.984118 | Indeterminate | Unresolved |
| Hydrogenophilalia `4380_6` | 21.450083 | 9.720464, repaired H | Ca-supported band | Unresolved; model/preparation-sensitive |

The invalid earlier AF3 value17.134098 remains excluded. The repaired run preserves the same selected AF3 sample1 and all its heavy coordinates, PQQ hydrogens/caps, core membership and charges; only23 protein H coordinates changed relative to the failed AF3 preparation. **Protenix and AF3 have different predicted heavy coordinates**—the heavy-coordinate preservation statement refers to the AF3 hydrogen repair, not to equivalence between predictors. Independent review verified63 pinned files and both normal/converged ORCA outputs, plus core geometry, atom identities, charges and paired coordinates. No remaining gross overlap was found. The numerical hydrogen-method change and absence of independent recalibration are explicit.

Repaired energies: E(La)=−2426.60502454349 Hartree; E(Ca)=−3072.667079829895 Hartree; R=−646.0620552864052 Hartree. Post-completion review: PLM `energetics_queue/adh9/af3_comparison/hydrogen_repair/post_completion_review.json`. Current two-protein table: `energetics_queue/adh9/af3_comparison/reviewed_score_comparison.tsv`. Original execution results, initial failure review and preparation snapshots remain immutable. The watcher delivered completion; no further calculations are running for this repair.


The same selected Hydrogenophilalia ADH9 AF3 model is now prepared with valid
protein-hydrogen geometry. This repairs a preparation failure after folding.
It does not change the AF3 model, MSA, heavy coordinates, PQQ state, conserved
core, charges, electronic method or released score bands.

Jacob explicitly approved diagnosis, repair with AF3 heavy atoms fixed,
validation, and rerunning the affected La/Ca pair: “let’s do it”. The approval is
`EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/adh9/af3_comparison/hydrogen_repair/authorization.json`.
This preparation agent submitted no jobs. The parent agent owns the two-endpoint
execution and watcher.

**Endpoint rerun submitted:** the parent's independent preflight passed and
Slurm job **1200182** now owns exactly this repaired La/Ca pair. Watcher PID
**1901467** monitors it. The execution manifest is
`workspaces/plm_adh9_af3_20260916/hydrogen_repair/execution/manifest.json`,
SHA256 `28fc91ff52785d66f0e178a1347f0d9f0db74781f364f9df86b04845e624eda0`.
Preparation, code and pinned execution artifacts are frozen while it runs.
Job submission is not endpoint completion; consult the PLM repair
`job_status.json`/`completion.json` for the latest state.

## Selected input and failure

- Protein `PQQSEQ_13d74836d4b7a3e02140`, gene
  `PLM0_60_coex_redo_sep16_scaffold_4380_6`.
- AF3 seed101/sample1, exact source SHA256
  `f8517d8ed8042e941c5d09d1ca9b927382e0e8712f9ce4ff683a9777e4f1e23d`.
- Original fixed core: Glu201, Asn269, Asp320, Asp322, Lys344 and complete
  oxidized PQQ³⁻, with La/Ca charges −2/−3. Both endpoint cores have78 atoms.
- Original core contained13 H–H pairs below0.5 Å, closest Asp320 HB2/HB3
  0.058626 Å. All overlapping H coordinates were copied from the saved
  protonated PDB. PQQ hydrogens, link caps and AF3 heavy atoms were uninvolved.
- The original ORCA tasks converged and terminated normally, but the resulting
  raw S17.134098 kcal/mol is invalid because of the input geometry. Preserve
  it as failed-preparation evidence, not as an accepted score.

## What failed and what changed

`capture_original.py` observed the unchanged pinned protonator, including its
original CPU/one-thread `forcefield=None` OpenMM objective. It reproduced the
failed PDB byte-for-byte. Initial energy was8.29645×10¹⁹ kJ/mol, maximum H force
1.00156×10²⁶ kJ/mol/nm. The minimizer returned after one reported accepted
iteration, with enormous residual forces. This establishes unconverged early
termination; the exact internal line-search stopping reason is not exposed.

OpenMM places initial hydrogen coordinates in clusters around their parent
atoms, then minimizes them. The original wrapper checked atom counts and
protonation states but did not check minimizer convergence or H overlaps.
An exact pinned preparation retry had already reproduced the same defect.

Recovery uses the **exact serialized System captured from that original call**
(SHA256 `33eaa43870ad91e3b7ce3bb48d04b0168fe5cf076bd4ae1878517d33f1b75ce5`).
It restarts from the retained failed PDB coordinates and raises the numerical
iteration ceiling from50 to1000. The objective, tolerance1.0 kJ/mol/nm, CPU
platform and single thread remain the same. All heavy atoms have zero mass and
remain immobile; every original H identity and protonation state is preserved.
No new random draws occur. The source H state derives from the original
seed20260914. No waters, synthetic ligands or new hydrogenation states are added.

This is an explicitly changed numerical preparation. Restarting the optimizer
and its serialized coordinates is part of the repair; increasing the iteration
ceiling alone was not established as sufficient.

## Validation and deliverables

The exact-objective recovery converged after219 iterations in24.8 seconds:

| Check | Result |
|---|---:|
| Full-precision H force RMS | 0.650485 kJ/mol/nm |
| Heavy displacement | 0 Å |
| Retained core H–H minimum | 1.708046 Å |
| Retained protein H–parent range | 1.177657–1.198481 Å |
| Minimum distance between any nonmetal core atoms | 1.010000 Å |
| Changed core atoms | 23 existing protein hydrogens |
| PQQ hydrogens, caps, all heavy atoms, atom ordering | unchanged |
| Core formula excluding metal | C27H31N4O15 |

Convergence was measured on full-precision minimizer coordinates. Normal
three-decimal-Å PDB serialization produces a maximum H displacement0.000836 Å
and H force RMS2.982187 kJ/mol/nm. That measured rounding effect is recorded
separately; we do not claim the serialized file satisfies the unrounded force
criterion. Geometric checks are applied to the actual serialized/carved atoms.

All compute products are under:

`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/plm_adh9_af3_20260916/hydrogen_repair/`

- `original_capture/`: original System, observed failed minimization and failed
  protonated PDB; unchanged original helper replay.
- `diagnose_50/`: bounded50-step same-objective restart, improved but not
  converged; diagnostic only.
- `recover_1000/`: recovery with reconstructed mathematically identical System,
  verified against the captured System; diagnostic only.
- **`recovered_exact_objective/`**: authoritative recovery using the exact
  captured System; full-precision coordinates, PDB, minimization and
  post-serialization force records.
- **`prepared/prepared_pairs.json`**: one case, prepared for the existing
  manifested endpoint executor; no ORCA run by this preparation agent.
- **`prepared/geometry_validation.json`**: PASS, pinned by both prepared-pair
  inventory and carve manifest.
- `prepared/PQQSEQ_13d74836d4b7a3e02140_AF3_sample1_Hrepair/`: new immutable
  La/Ca XYZ/input pair, carve/protonation/heavy-coordinate records.

The new carve uses the unchanged pinned fixed-core fragment and PQQ helpers.
It records a new wrapper, hydrogen-recovery provenance, original failed carve,
new implementation pins and approval. It does not impersonate the original
protonation implementation.

Nine real-artifact regression checks pass in `test_repair.py`, including
rejection of the original overlapping coordinates; rejection of moved heavy,
PQQ-H or cap atoms; altered identity, charge or nonfinite coordinates; and
unconverged minimization. `software_validation.json` records the check.

## Interpretation boundary

The intended fixed-core chemistry and endpoint Hamiltonian are unchanged;
the numerical hydrogen-preparation repair is explicitly documented. No score
band was refitted and no separate calibration panel was rerun. Thus compare
the repaired score using the existing declared scale, retaining that preparation
deviation and the original calibration's limitations. The repair does not
establish biological metal use, substrate specificity or general transfer
accuracy. It does not authorize refolding or rescoring other ADH9/PLM proteins.
