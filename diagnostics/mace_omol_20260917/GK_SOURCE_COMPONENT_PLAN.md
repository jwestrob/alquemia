# Separate source self-solvation from protein cross-solvation

The saved GGR direct/induced localization passes closure and finds a distributed
protein contribution. All three metal effective Born radii are near 30 A. The
pinned Tinker `born.f:tanhrsc` explicitly bounds this model at 30 A; the OpenMM
reference implementation has the same bound and rescaling formula. This is an
approximation in the model, not yet evidence of a coding defect or a reason to
choose a favorable radius.

Before interpreting the remaining +12.6002 kcal/mol GK difference, compute its
source-source and source-environment pieces from the actual archived primary
outputs. Include all three GGR structures and both alpha structures. No score,
settings, state, density or cavity changes; no new scientific executable calls.

All frozen source moments are monopoles. With their actual projected charges q,
actual native effective radii a and fixed physical coordinates, the native
`esolv.f:egk0a` source-only term is

G_QQ = 0.5 * electric * (1-78.3)/78.3 * sum_ij q_i q_j /
       sqrt(r_ij^2 + a_i a_j exp[-r_ij^2/(2.455 a_i a_j)]).

Sum the diagonal and both off-diagonal directions exactly once. The environment
permanent term and the common nonpolar term have already been subtracted from
the archived GK transfer. Consequently G_QE = GK_transfer - G_QQ. This residual
contains all source/environment monopole and higher multipole cross terms; it
is not a residue attribution. Induced energy remains separate. Do not use a
bare Coulomb expression in place of the screened native kernel.

Pin native source/constants, physical identities, charges and Born radii; verify
source higher moments vanish and same-cavity endpoint radii match. Repeat using
the actual rigid outputs with <=1e-7 kcal/mol tolerance. Keep explicit status:
analytical decomposition of the native expression, not independently executed
source-only native energies. If this identifies a consequential term, validate
it with a separately declared native component test before changing a model.

Report every structure's G_QQ and G_QE Ca-minus-La contributions, their closure,
and the three GGR pairwise differences. No calibration, new classification,
post-hoc source selection or claim that the bound uniquely causes the failure.
