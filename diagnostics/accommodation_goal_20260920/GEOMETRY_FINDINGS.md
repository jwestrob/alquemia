# Geometry diagnosis and reconciliation

Read all25 canonical preparations and all176 completed PLM records:137 scored
coordinate maps reproduced,39 unavailable kept. No new molecular calculation or
score. Inputs and all162 available feature rows are in
`workspaces/accommodation_goal_20260920/geometry_v3/result.json` and `features.tsv`.
The initial two parser attempts failed on relative paths and different PQQ atom
names; their empty directories are retained. Corrected extraction maps the
documented Protenix/CCD PQQ atom conventions and checks actual endpoint coordinates.

Within the137 PLM predictions, the nearest anchor-Glu distance has descriptive
Spearman rho0.7299 with the existing DFT score. Extra-Asp nearest/second distance
medians are2.269/3.490A in the21 Ca-band predictions versus2.618/2.923A in the94
La-band predictions. All groups have median protein-ion iPTM0.98. These are
associations with a calculation using those same coordinates, not measured
metal-use labels or proof that a particular geometry is correct.

The two negative extremes have extra-Asp La–O contacts1.766 and1.847A, compared
with2.410–2.648A in the11 La calibration models. This is a useful mechanistic
lead, but **not a newly discovered defect**: the other PLM session already
completed a broader528-model geometry assessment. Its exact source is:

`/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/xoxf_all/geometry_quality/REPORT.md`

That report found26/137 warning-bearing selected structures, including15/21
Ca-band predictions. Both extremes remain compressed in every saved sample.
Ten selected models have warning-free saved alternatives; these alternatives
remain unscored. Reuse the existing all-sample descriptors instead of repeating
the review. Its separate RELAXATION_PROPOSAL.md proposed eight constrained DFT
optimizations but did not run them. No such campaign has been silently launched.

The next useful question is causal: does metal-dependent work along physically
allowed whole-carboxylate motion distinguish a strained representation from a
stable alternative organization? Generic vacuum-MACE relaxation and radial force
averages already failed to give robust gains. The parallel solvent-consistent
physical-mode gradient pilot tests a specific difference in Hamiltonian and
observable before the next PQQ intervention. No PLM prediction is relabeled,
no geometry is selected by desired score, and original results stay intact.
