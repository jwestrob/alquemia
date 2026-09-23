# How much of the ten-fold pocket survives three-fold input?

Read-only inventory,23 September2026. Under Jacob's discretionary development
authorization, root compared every existing three-of-four noncanonical
La-conditioned source subset for all25 reference proteins. Membership comes
from the archived source fragment selections, without reading energies or labels.
No source was prepared anew and no molecular calculation ran.

Of100 declared triples,61 select exactly the same fragment list as the tested
ten-fold union;33 select a smaller list;6 contain an unavailable preparation.
The33 changed triples reduce to16 distinct protein/fragment selections and55
distinct source/selection pairs. These counts are an inventory, not a new score
comparison. Graph closure might make some differing lists physically identical;
that has not been tested here.

All four Q9Z4J7 and C5AXV8 triples have identical fragment selections. All four
A0A3F2YLY8 triples omit the Asn285 peptide and Cys129 sidechain present in the
ten-fold union. This matters because the successful full-fold method used a
broader physical context to repair A0A3's difficult Ca-conditioned structure.
Three La folds cannot simply be assumed to reproduce that representation.
No causal contribution is assigned to either omitted fragment separately.

The next practical question is whether a fixed, source-supported pocket can
retain the gain with fewer input folds. Preparation equivalence and exact cache
reuse can reduce the work, but this inventory does not authorize relabeling the
existing ten-fold result as a three-fold validation. Current source-group tooling
continues to expose the actually tested policy.

## Actual inputs and result

- Union preparation: `workspaces/consistent_context_20260922/prepared_v1/PREPARATION.json`.
- Source preparation: the hash-pinned `source_preparation` in that record.
- Each supported source's hash-pinned `context/preparation` supplies
  `added_fragments`; identities are kind,chain,residue number,insertion code and
  residue name, matching `consistent_context.fragment_id`.
- For each protein, the union of all supported source lists was checked against
  the saved ten-fold list. The four noncanonical La sources supply all four
  three-member combinations; missing preparation makes a triple unavailable.
- Full per-triple omissions, source IDs, counts and input hashes:
  `workspaces/union_subset_inventory_20260923/inventory_v1.json`.

All sources are previously consumed structural development examples. Neither
candidate counts nor fragment-list equality establish biological independence,
equilibrium populations or measured specificity.
