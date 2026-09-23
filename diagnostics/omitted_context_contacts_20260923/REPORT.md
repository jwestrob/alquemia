# No severe omitted-protein intrusion in the admitted proposals

The declared audit finds **no newly created protein-heavy contact below2.0Å**
in any of the110 threefold proposals or110 matched tenfold proposal records.
All55 context pairs are supported in both representations, with no pre-existing
<2.0Å contacts in these query/outside sets either. This does not support an
omitted-atom collision explanation for A8's regression. It does not establish
that omitted packing, electrostatics or solvent/cavity effects are negligible.

These55 pairs contain41 distinct folded sources from11 reference proteins;
tenfold counterparts repeat where one fold occurs in multiple threefold unions.
The220 proposal records are correlated comparisons, not independent biological
observations. The fixed set includes four successful A0A3 source examples.

## Moving donors remain separated from omitted protein atoms

Minimum observed distances across all55 paired records, after exclusions:

| Representation | Origin | Ca proposal | La proposal |
|---|---:|---:|---:|
| Threefold, moving-donor subset |2.7880Å|2.8423Å|2.9037Å|
| Tenfold, moving-donor subset |2.9896Å|2.9433Å|2.9448Å|

The minimum across **all** represented protein-heavy atoms is2.3278Å in both
representations and remains fixed across all three geometries. Thus the closest
contact is not created by these donor movements. Counts <3.5Å are descriptive:
threefold moving-donor totals are53/65/61 at origin/Ca/La; matched tenfold totals
are39/49/45. These sums include repeated sources and different context boundaries;
they are neither clash rejection rules nor comparable independent sample counts.
Every per-context minimum/count and atom-pair distance is retained in the outputs.

For the problematic168-atom A8 context, the nearest omitted-atom distance for its
moving donors **increases** during accommodation:

| A8 source | Origin | Ca proposal | La proposal |
|---|---:|---:|---:|
| La sample1 |3.5359Å|3.5764Å|3.6030Å|
| La sample3 |3.4494Å|3.4953Å|3.5174Å|

All seven A8 source/context pairs and both representations have zero severe
flags. This negative finding is consistent with keeping root's separate component
audit as the relevant current explanation: removal of neutral Ser352 mainly
changes the origin solvent contrast; differential accommodation slightly decreases.
The contact audit itself cannot assign an energy cause.

A0A3's four included sources also have zero severe contacts in either
representation. Its threefold sample2 approaches omitted atoms from3.5078Å at
origin to3.3667/3.4128Å, producing two descriptive <3.5Å contacts in each proposal;
this is not treated as a failure or a correction to its score.

## Actual source geometry and exclusions

Queries are real source-backed protein-heavy atoms in the scored context; the
moving-donor subset is the union of atoms displaced >1e-8Å in either proposal.
Outside atoms are raw-source protein-heavy atoms absent from that context,
frozen at their original coordinates. Exclusions cover represented atoms,
synthetic caps, H/D and covalent neighbors through three source-graph bonds.
The graph uses the pinned protein topology, observed peptide connectivity and
the existing validated source-disulfide policy, with1–3 disulfides per source.
No new bond policy, hydrogen reconstruction or coordinate preparation was used.

All source-heavy coordinates match the original structures exactly. No added
protein heavy atoms needed exclusion. Every saved proposal replays its actual
physical mapping/cross-metal coordinate pair. The complete24-heavy-atom PQQ set
is mapped and fixed, and lies outside this protein-versus-protein contact scope.
No unsupported cofactor/graph case occurred. A separate check of all330 context/
geometry records confirms zero motion of unscored physical heavy atoms, so the
frozen-outside assumption agrees with the actual saved kinematics.

The local audit took113.704s wall time; zero Slurm, GPU, molecular or protonation
calls. Four final real-artifact tests pass, covering actual coordinate distances,
source/graph coverage, three-bond exclusions, thresholds and outside stationarity.
The first test log retains an incorrect hand-written PQQ count assertion22;
it was corrected to the actual pinned24 without changing any analysis output.
Existing source structures, energies, proposals, selections, references and
classifications are untouched. No energy correction or rejection filter follows.

Recommendation: close the severe-intrusion hypothesis for this fixed population
under the declared2.0Å criterion. Keep the existing environmental/membership
investigation separate. [Actual results/pins](RESULT.json), [commands](COMMANDS.md)
and the complete TSV/JSON preserve all55 rows and both representations.
