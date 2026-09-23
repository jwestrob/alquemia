# Existing disulfide/source topology policy

Before any contact distances are inspected, the graph implementation additionally
records use of the existing `environment_context_chemistry.cysteine_units` rule
for source cysteines: installed source-PDB topology, actual existing hydrogen
inventory and the already published1.8–2.3Å disulfide range. This reuses the
preparation's existing covalent-connectivity policy; it does not introduce a new
proximity-based bond graph or move/create any atom. Relevant missing template
neighbors, unsupported cysteine state or explicit unhandled source connections
produce an unsupported row. Internal PQQ connectivity is unnecessary for this
protein-versus-protein contact question because the entire PQQ heavy set is
mapped/fixed and excluded from both protein query and outside sets; verify that
condition explicitly and retain the cofactor handling in each result.
