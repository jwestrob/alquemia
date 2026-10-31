# Declared implementation and metadata corrections

No feature additions, numerical acceptance changes or outcome-based descriptor tuning.

- Historical reproduction 31942: the original 101-record feature/prediction table
  replays exactly. Old outputs were never overwritten.
- First curation attempt made zero records because the native equilibrium text
  is `[ML]/[M][L]`, not `M.L`. Corrected native-string admission before any fit;
  raw constants were not reinterpreted or changed.
- Descriptor 31953 used a donor centroid for the center start. A one-donor ligand
  then places the center exactly on its donor, at the radial norm's singularity.
  Before inspecting predictive results, 31955 corrected the deterministic start
  to centroid+(0,0,1) A, identically for both radii. First numerical results are
  retained under `workspaces/ln_structure_increment_20261030/implementation_v1/`.
  Real acetate and EDTA analytic/finite-difference checks verify the correction.
- Initial grouped fits 31961 were inspected. A chemistry audit then identified a
  grouping-code error: the `N+` heuristic put nitrobenzoic acid (3198) in the
  pyridine family, apart from its benzoate relatives. The corrected rule explicitly
  recognizes pyridine N-oxide and gives phosphate priority over aromatic N (the
  ATP candidate was already excluded from primary measurements). Only ligand 3198
  changes primary group; no primary identities, labels, conditions or exclusions
  change. Original tables, freeze, fits and predictions remain in `grouping_v1/`.
  Current FREEZE is version 2 and flags this post-initial-results metadata repair.
  31968 repeats precisely the same models/grid/selection/uncertainty/criteria.
  Do not describe that revised grouping as untouched prospective registration.
- Check 31964 failed on RDKit rejecting a NumPy integer index, after saved energy
  replay; 31966 uses a Python integer. It repeats saved checks only, no new
  descriptor generation or fit. Both job logs/costs are retained.
- The optional 2025 peptide source is now verified: FIDTNNDGWYEGDELLA, free N
  terminus, amidated C terminus; N5D is FIDTDNDGWYEGDELLA. SI S3/S6–S7/Table S3
  reports pH6, 50mM MES, 100mM NaCl; competition-assay temperature is not stated.
  The fixed M4-S illustration therefore compares log preference, not an imputed
  experimental kcal/mol target. No peptide simulation or fitting is performed.
- Exposure audit 31958 was stopped during metadata normalization, before any
  inventory, application or primary inference result was produced. Its cache was
  keyed on the entire complex, repeatedly normalizing the same ligand for other
  metals. Recovery 31979 caches the
  ligand separately. Normalization rules and models are unchanged; only redundant
  work is removed. Original script/logs and scheduler costs remain preserved.
- Recovery 31979 was also stopped during the audit, before inventory/application/
  primary inference results. Full tautomer enumeration of irrelevant released
  ligands dominated runtime. Audit 32006 uses a heavy-element-composition filter:
  protonation and tautomerization cannot change that signature, so every possible
  exact evaluation-ligand match still receives full normalization. Other ligands
  remain charge-normalized and are all audited for families/close analogues;
  their row-level normalization status is explicit. The older LOGKPREDICT loop
  now uses the same ligand cache. No row is omitted from the released inventory,
  no model/input/score is selected, and no absence claim is based on this shortcut.
  Both interrupted metadata attempts' scripts/logs/accounting are preserved.
- 32006 completed the 121,078-row overlap inventory, then failed the modern
  LOGKPREDICT Python API example: its mask produces 21 features for a released
  checkpoint expecting 48. The original `LOGKPREDICT` executable in the same
  pinned release supplies 48 features. Recovery 32009 uses that executable
  unchanged and reuses the completed audit; it does not alter checkpoints or
  infer missing features. Pending metadata job 32007 was rebound to 32009.
- Original LOGKPREDICT executable 32009 runs but gives logK 5.71413952 versus
  bundled output 5.40875760 (difference 0.30538, failing the frozen 0.02 replay
  tolerance). This failure remains failed. Upstream history shows the example
  output changed in the 2025 refactor while the original executable did not.
  Do not retune versions/features or claim exact legacy reproduction. 32015
  reuses that completed example and audit, runs the separate 2025 M4-S uranyl
  application inputs (exact numerical ML outputs are not tabulated by the source),
  and collects the 2025 chemical-graph references. Their clean-comparison status
  remains unresolved. 32007 now depends on 32015.
- Saved summaries 32016/32021 wrote their aggregate tables, then failed a final
  grouping assertion because a local variable holding the primary index had been
  shadowed by the published-prediction table. 32027 gives that index a separate
  name and repeats saved aggregation only. No fits, descriptors or inference change.
