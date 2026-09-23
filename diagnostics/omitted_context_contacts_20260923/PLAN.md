# Frozen geometric audit of omitted protein contacts

Root explicitly assigned this bounded read-only diagnostic under Jacob's
ongoing discretionary research scope on2026-09-23. Question: do the admitted
proposals move real donor atoms into protein atoms omitted from their scored
context? Audit all55 newly scored source/context pairs in the completed
threefold transfer, and each exact source's tenfold-precision counterpart.
Population is fixed by `pool_reused=false` in the existing comparison; no extra
sources, score-based selection, molecular calls or Slurm submission.

For each context replay the actual saved origin, Ca proposal and La proposal.
Use actual source-heavy coordinates and existing source/cap Kinematics mapping.
The query set is source-backed protein heavy atoms included in the scored
context; also report the subset that moves in either proposal (>1e-8Å, only to
exclude floating-point identity). The outside set is original source protein
heavy atoms absent from that context. It is fixed at source coordinates.
Synthetic caps, hydrogen/deuterium, cofactor atoms and metal are not protein
contact queries. Mapped context atoms cannot also occur in the outside set.
Any source-preparation-only added heavy atom is explicitly inventoried rather
than silently treated as an observed raw-source atom.

Exclude actual source-graph pairs separated by one, two or three covalent bonds,
using the same pinned protein topology and observed peptide connectivity.
Unsupported relevant graphs/cofactors, unmatched source atoms, nonidentity source
heavy coordinates or inconsistent saved candidate mappings remain explicit;
no inferred new bonds, atom reconstruction or geometry rescue. PQQ must remain
complete, fully source-mapped and fixed. Other unexpected cofactors or explicit
nonstandard covalent connections are reported as unsupported rather than guessed.

Predeclare severe geometric flags as newly created distances <2.0Å whose exact
origin pair was >=2.0Å. Preserve pre-existing <2.0Å contacts separately. Report
minimum distances and counts <3.5Å as descriptions, not admission or rejection
rules. Compare all context protein-heavy queries and the moving-donor subset,
with exact pair identities for severe/new contacts. These metrics do not change
any score, selected geometry, reference, threshold or classification.

The A8 regression is not presumed a steric failure: root's existing component
audit attributes most of its origin contrast shift to solvent after removal of
neutral Ser352, with slightly reduced differential accommodation. This contact
audit can test one physical-plausibility mechanism, not assign a unique energy
cause. Include A0A3 wherever it occurs in these same55; do not add label-selected
controls. All55 paired rows/unsupported counts and repeated-source/group
identities stay explicit; they are correlated structural repetitions.
