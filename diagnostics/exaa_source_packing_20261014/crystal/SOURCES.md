# Source evidence and limitations

- [Official1FLG deposition](https://www.rcsb.org/structure/1FLG): 2.60 Å,
  no mutations, two 582-residue chains; PQQ and Ca are the only nonwater ligand
  types. The deposited `_struct_conn` records contain metal coordination and
  the adjacent-Cys disulfide, with no PQQ external covalent connection listed.
  The Slurm audit checks the complete actual file and both sites independently.
- [Keitel et al.2000](https://pubmed.ncbi.nlm.nih.gov/10736230/),
  DOI10.1006/jmbi.2000.3603: the abstract identifies ATCC17933 enzyme and the
  extra N-terminal Ca. Current UniProt mapping is Q9Z4J7/PAO1; the saved exact
  sequence alignment supports identity over resolved residues, not identical
  strain/construct provenance.
- [Author-repository article](https://digitalcommons.unl.edu/chemistrystezowski/6/)
  was found; its linked PDF returned HTTP403 with both normal retrieval tools.
  No statement about the article's detailed redox assignment is inferred.
- Immutable coordinates and completed sequence/site audit:
  `workspaces/ca_coverage_recovery_20261012/refold/run_v3_h200/prepared/1FLG.cif`
  and `1FLG_AUDIT.json`.
- Completed cofactor metrics:
  `workspaces/pede_solvent_response_20261013/cofactor_review_v2/execution/RESULT.json`.
  Both1FLG instances have long carbonyl-labelled bonds/nonplanarity; the report
  expressly did not assign redox. These metrics are reused, not recomputed.
- Existing Q9 preparation manifest and completed output:
  `workspaces/ca_coverage_recovery_20261012/spectator_preparation_q9_v2/manifest.json`
  and `execution/FINAL_PREPARATION.json`.

The three already-scored crystals1H4I/4MAE/1KB0 are different proteins and are not
rerun.1FLG is a new same-sequence structural source, explicitly conditional on
the existing dry oxidized-PQQ model and the actual preparation gates.
