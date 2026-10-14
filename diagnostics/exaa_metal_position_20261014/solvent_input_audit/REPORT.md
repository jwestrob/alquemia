# Original solvent controls: all eight exact cells are reusable

The two ExaA sources, two metals and two saved geometries used by MACE4130 have
**8/8 available original frozen-charge ALPB receipts**. All eight also have
their original standalone xTB vacuum and water charge files, raw outputs,
execution receipts and accepted xTB/tblite component-v2 comparisons. No original
cell is missing. [EVIDENCE.json](EVIDENCE.json) records the complete joins and
pins; its eight rows are keyed by source, metal and candidate.

This audit checks saved metadata and byte hashes only. It makes no molecular
calls, recomputes no energies or charge sums, and changes no scientific admission
rule. The metadata audit verified 1,771 pinned artifacts, including the runtime
and native build provenance. Both source states retain 192 atoms, 588 GFN2
valence electrons, Ca charge −3 / La charge −2, and multiplicity 1.

## Exact input chain

- Original inventory: `workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json`,
  SHA `757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343`.
- Original SCC components and raw populations:
  `workspaces/ca_coverage_recovery_20261012/two_ca_solvent_q9_v1/FINAL_COLLECTION.json`
  and its `points/input_*/{RESULT.json,xtb/charges,xtb/stdout.txt,xtb/EXECUTION.json}`.
- Original native replay receipts:
  `workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/transfer_v1/tasks/native_*/RESULT.json`.
- Actual reused MACE/XYZ identities:
  `workspaces/exaa_metal_position_20261014/scoring_v1/execution/FINAL_COLLECTION.json`.

For every cell, the inventory MACE receipt/request and original XYZ exactly
match the original-cell record in MACE4130. The native input XYZ, executed SCC
XYZ and replay XYZ hashes match that same original geometry. State, source
preparation, maps and electron metadata also match. Each stored native water
diagonal has four passing informative component checks, each vacuum replay has
four informative precision bounds, and the water/vacuum replay cavity and Born
matrix files are byte-identical.

The original per-point `accepted` field retains the failed **absolute-energy**
xTB/tblite v1 gate. It is not the component-v2 acceptance flag and must not be
rewritten as passing. The separately stored component-v2 pairs are accepted;
both backend calculations are available with their original solver evidence.

## Qualified recipe and parameter identity

Standalone SCC uses xTB6.7.1, the original `--gfn 2` route, physical singlets
(`--uhf 0`), electronic temperature 300 K, accuracy 0.0001, maximum500 iterations
and `--norestart`. The exact parameter file is
`/home/jwestrob/miniconda3/envs/lanm_bench/share/xtb/param_gfn2-xtb.txt`, SHA
`f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545`.
The executable SHA is
`481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f4`.
The independent component check uses the separately installed tblite0.7.0;
the older libtblite0.4 linked into the xTB environment is not that independent
backend. Runtime and dependency pins are retained in the evidence.

The qualified native replay uses the **internal GFN2-xTB ALPB factory**, water,
P16, `gsolv`, 230 angular points per atom, factory temperature298.15 K and no
ionic screening. That factory temperature is distinct from the SCC electronic
temperature300 K. It includes all four native terms: Gelec, Ghb, Gsasa and Gshift.
The existing per-component rounded-charge bounds and 1e−6 Eh informative ceiling
remain unchanged. The exact qualified generic replay binary is
`frozen_distribution_candidate_v1/canonical_v1/replay`, SHA
`3d6a9d23aa1d1152c368b21f7f6d55a135279596b8d01a3a1b3ebc42617e75b6`;
its build, native static libraries and wrapper sources are pinned by the audit.

## What remains new

The eight original charge vectors and energies are exact comparator reuses.
**They are not reusable SCC solutions at the transplanted metal positions.**
The next eight geometries require their own fixed-geometry SCC/component
qualification and native frozen-charge replay under the root's finite solvent
plan. This audit supplies inputs; it does not evaluate those new cells or infer
a solvent-corrected response, classification or affinity.

## Existing coupled-search pointers, for later decisions

`scripts/adaptive_metal_proposals.py` already implements three Cartesian metal
translations plus four previously selected donor angles: `full_q`, `constraints`,
`physical_status`, `final_geometry`, `projected_gradient_eV`, `MixedProposal` and
`optimize`. It uses `scripts/mace_site_kinematics.py` and the physical cap Jacobian.
The actual frozen copy at
`workspaces/motion_envelope_20260923/pilot34_searches_v1/implementation/adaptive_metal_proposals.py`
has SHA `a89ac5fd974ceca2952d944264347c8fb41631eebab46979101b0f366192aac8`.

The existing checks retain the0.8 Å original-host heavy-atom displacement bound,
coordinate bounds, unchanged nonselected atoms and the source/cap/bond/overlap
checks. The present transplant preparation already reused that frozen physical
checker. The historical joint pilot is documented in
`diagnostics/adaptive_metal_20260922/FINAL_RESULT.md`: seven of eight candidates
were admitted, with force relief but no changed available class calls; the
oversized trial remained rejected. These are implementation/evidence pointers,
not a new search plan or permission to rerun that pilot.

Reproduce this metadata inventory with
`python diagnostics/exaa_metal_position_20261014/solvent_input_audit/AUDIT_INPUTS.py`.
