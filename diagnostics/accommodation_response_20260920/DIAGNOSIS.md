# Why the earlier response experiments did not improve the scanner

The evidence supports several different failures, not a general conclusion that
metal-dependent organization carries no information.

1. **Optimizing the wrong target can lower energies without improving specificity.**
   All ten complete-context coordination searches optimized native vacuum OMOL,
   while their final readout was original-core CPCM DFT. Nine DFT endpoints went
   downhill, yet all three class/direction gaps narrowed. All ten proposals touched
   the imposed0.20 Å boundary; none establishes an unconstrained physical minimum.
   Those results reject that particular preparation policy, not the existence of
   solvent-aware donor response. PQQ and waters were fixed in that policy.
2. **A radial mean throws away physically linked organization.** The old native
   POLAR descriptor changes from+6.67365 for1GLG to−1.83307 for2FW0 and+0.85078 for
   2FVY kcal/mol/Å. Several residues contribute: Asn136 is+18.095/−2.923/+0.495;
   backbone Gln140 is+13.094/−1.093/+1.594. This is not just a single outlier.
   Alpha6IP9's Asp82 oxygens contribute+48.262 and−19.411, which a donor average
   obscures. These Cartesian radial components are not a covalently valid
   collective displacement of the protein. A source OD1/OD2 name switch alone
   is not a chemical-state change: equivalent carboxylate oxygens need consistent
   whole-residue mapping.
3. **The old forces do not validate one another.** At the exact old paired centers,
   mean radial POLAR-vacuum/DFT-CPCM differences are−2.990/+7.221 for1F6S and
   −1.414/+7.031 for6IP9. Their full donor-vector cosines are0.246 and0.900.
   Solvent, electronic model and sometimes representation differ; neither a
   matching order nor an opposite sign isolates the source of error. These values
   were independently replayed from12 actual POLAR force arrays and8 analytic DFT
   files in `workspaces/accommodation_response_20260920/archived_response_replay_v1.json`.
4. **The water result separates first-order usefulness from wider surface error.**
   Eleven anchored local steps lowered DFT energy, but one endpoint remained
   nonstationary and failed its soft-mode test. Wider-domain representative
   errors were0.331–0.685 kcal/mol weighted and up to2.156 kcal/mol at one point.
   A Cartesian DFT tangent fixes the force at one center; it does not fix solvent
   nonlinearity, curvature or another basin. More samples of the unchanged
   potential would not repair this. None of the data supplies entropy or occupancy.

The compact solvent-transfer composite is a useful new opportunity because its
*static* descriptor already improves the consumed core alpha/GGR ordering and
retains PQQ fidelity. Its gradient can be assembled consistently from the same
three energy terms, instead of using vacuum forces to guide a solvated score.
That is an explicit, testable change in Hamiltonian. Whether it yields useful
additional accommodation information is still a hypothesis.

The present test retains every signed derivative of actual donor torsions and
peptide modes. It does not fit a new force norm or radial coefficient to labels,
optimize until the desired class appears, or infer a relaxation energy from an
invented stiffness. A later fixed, whole-carboxylate displacement can then measure
actual metal-dependent work along a meaningful motion. Its useful outcome would
be new organization information; merely reducing an electronic energy is not
sufficient.

Primary project sources: `coordination_preparation_20260919/REPORT.md`,
`response_probe_20260919/{REPORT,CONTINUATION_REPORT}.md`,
`hydration_basin_20260919/REPORT.md`, `water_basins_20260919/REPORT.md`, and
`compact_solvation_20260920/REPORT.md`. The old force pairs are POLAR-medium;
the current compact composite uses OMOL0-100M. No Hamiltonian equivalence is implied.
