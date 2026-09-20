# Structural causes of selected all250 preparation failures

2026-09-20. Read-only inspection requested by the parent during preparation
job1203744. No new preparation, energy, fold, replacement or repair was run.
The fixed denominator remains **250**. This is a targeted failure audit, not
a final count of successes in the still-running campaign.

## C5AXV8: five Ca-conditioned ions occupy another predicted site

All five failures say `exact selected PQQ required`. The actual selection
in `carve_with_pqq._prepare_nearby_pqq` requires a PQQ atom within **4.0 Å**
of the selected metal. The full PQQ is present, but the Ca is far away:

| AF3 sample | Ca-conditioned: nearest PQQ atom, Å | La-conditioned: nearest PQQ atom, Å |
|---|---:|---:|
| 0 | 25.0162 | 2.4469 |
| 1 | 25.1070 | 2.4587 |
| 2 | 25.0244 | 2.4937 |
| 3 | 25.0385 | 2.4640 |
| 4 | 25.2093 | 2.4716 |

In every Ca fold the three nearest protein O/N atoms are the backbone O atoms
of **Asp36, Thr39 and Asp42** (order varies). Asp36 O is 2.331–2.395 Å away.
The declared PQQ-site roles Glu198, Asn275, Asp317, Asp319 and Arg344 have
their nearest O/N atoms **27.9–32.7 Å** away. In the La folds, the selected
ion instead contacts the declared PQQ-site residues: nearest Glu198 O/N
2.455–2.658 Å; nearest Asn275 O/N 2.486–2.619 Å.

The raw metal coordinates are unchanged in the saved normalized/protonated
files. The Ca residue-name alias did not relocate the ion. Thus the rejection
identifies **off-site placement in these Ca-conditioned predictions**, not a
missing PQQ molecule or a scorer classification error. This geometry alone
does not establish an experimental Ca site or a biological specificity label.
No ion was moved to the PQQ site and no replacement sample was selected.

## MMOL1770 La sample4: malformed prepared hydrogens, before context expansion

The failure says `overlapping expanded atoms/caps`, but the existing
**80-atom core already contains 10 H/H pairs below 0.45 Å**. All 10 pair
members are recorded as `source_protonation_hydrogen`; **none is a link cap**.
Examples:

| Residue | Atom pair | Separation, Å |
|---|---|---:|
| Arg383 | HH11–HH12 | 0.104346 |
| Asn318 | HB2–HB3 | 0.111937 |
| Glu200 | HG2–HG3 | 0.139617 |
| Glu200 | HB2–HB3 | 0.145458 |

The full saved protonated protein has **1,759 atom pairs below 0.45 Å**;
the shortest is Pro492 HB2–HB3 at 0.004123 Å. These H atoms were added by
the pinned PDBFixer1.12/OpenMM8.5.1, pH7, RNG20260914, single-thread CPU
preparation; the original fold contains **zero H atoms**. All **4,895 raw
heavy atoms** are retained with **zero coordinate displacement**, with no
unmatched heavy atoms and no raw heavy-atom pairs below 0.45 Å.

Both saved Ca/La core XYZ files exactly reproduce all 79 nonmetal fragment
atom records. The context overlap check therefore catches a real existing
hydrogen-preparation defect; it does not demonstrate a defect in the source
heavy-atom fold or in synthetic cap placement. The precise internal cause
of the protonator's malformed H output is **not established** by this audit.
No retry, optimization or hydrogen repair was performed.

## Evidence and reproduction

- [Pinned all250 source configuration](PQQ_ALL250_SOURCES.json).
- [Read-only measurement script](audit_failed_fold_geometry.py).
- Full per-atom distances, overlaps, source/artifact hashes and protonation
  receipt: `workspaces/accommodation_controls_20260920/fold_failure_geometry_v1.json`.
- Original failures/artifacts:
  `workspaces/accommodation_goal_20260920/folds_v1/source_preparation/`.
- The `archive_coordinate_mismatch` last-bit issue was left to the parent,
  as requested. No score or biological label was derived from these failures.

From the repository root, use a new output path because the audit refuses
to overwrite an existing record:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/accommodation_controls_20260920/audit_failed_fold_geometry.py --sources diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json --prepared workspaces/accommodation_goal_20260920/folds_v1/source_preparation --output workspaces/accommodation_controls_20260920/fold_failure_geometry_recheck.json
```
