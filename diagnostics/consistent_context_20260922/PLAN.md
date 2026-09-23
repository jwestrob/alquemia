# Consistent protein context across saved folds

## Agreed question and fixed scope

Jacob approved the three proposed next directions on 22 September. This branch
asks whether changing context membership between folds contributes to unreliable
PQQ calls. It holds surrounding fragment membership fixed within each protein.
Root assigned this branch independently; production/defaults and immutable source
preparations remain unchanged. No new folds, DFT, water, protonation, optimization,
PLM rescore or automatic promotion.

Use all 250 saved reference sources (25 proteins, five Ca- and five La-conditioned
folds each) from the existing PQQ_ALL250_SOURCES record. The 25 exact canonical
matches remain calibration; the 225 other sources remain structural transfers.
Keep all 17 earlier preparation failures in the denominator. Geometry selection
does not read score/label outcomes. The four pilot groups were explicitly fixed
by root: A0A3F2YLY8, A0ACD6B9F2, Q9Z4J7 and Q88JH5.

## Membership and chemistry, frozen before energies

For each protein take the union of the existing selected, source-backed fragment
identities across all previously supported saved folds. A fragment identity is
kind (complete sidechain or peptide), chain, residue number/insertion and actual
residue name. Reuse the original fixed core, complete amide graph closure,
proline and source-disulfide handling, cap construction and formal-charge rules.
Do not apply the geometric neighbor cutoff again to remove union members.

Validate actual protein sequence/residue numbering across the group, rather than
assuming that an accession supplies an atom map. The designated canonical source
defines the common selected-state signature after union preparation: exact source
atom identities/elements, H identities and covalent parents, retained bonds,
cut-bond/cap identities, fragment formal charges, cofactor/proton/water inventory,
paired electronic charges and multiplicities. Every admitted fold must match.
No renamed/guessed H, proton transfer, missing-atom construction or microstate
selection is allowed to make a fold agree. A source or canonical signature failure
is explicitly unavailable; no substitute canonical or reduced pool is chosen.

All source heavy atoms and H stay at their archived coordinates. Existing core
atoms and surviving caps are exact copies; graph closure may replace a previous
cap with its actual source atom. Each paired Ca/La preparation has identical
nonmetal coordinates and differs only in metal and its charge. All retained/cut
bonds and source-to-QM mappings are saved. Reproducing the old per-source fragment
list must reproduce that source's archived context within the existing 1e-12 Å
roundoff allowance before treating the new union adapter as compatible.

This is a new opt-in protocol, not a silent geometry/threshold migration. The
union uses all available unlabeled fold geometries, so this is a development
cohort test of fixed membership, not prospective blind validation or evidence
that a future three-fold scanner already has the same context.

## Preparation, finite scoring pilot and calibration

Prepare all 250 denominator rows, but initially score only the 25 designated
canonical sources plus all nine noncanonical sources in each of the four pilot
groups: 61 declared unique sources. The prior inputs support at most 50 of those
sources (11 prior failures), giving at most 100 native MACE and 200 native GFN2
singlepoints before new state failures or exact compatible reuse. No energies
start until the actual prepared inventory and finite manifest are sent to root
for engineering coordination. Preparation alone does not establish improvement.

Use the unchanged native OMOL checkpoint and

`E_M = E_OMOL,vac,M + E_GFN2,ALPBwater,M - E_GFN2,vac,M`;

`R = E_Ca - E_La`, in explicit once-converted units. No aquo offset, fitted mixture
or new physical term. Use the already qualified native GFN2 MaxIter500 ceiling,
unchanged Hamiltonian/mixer/smearing/convergence and explicit qualification pin.
Failures remain failures. Existing warm-MACE and ORCA runners handle execution.

Old static bands are a separately identified transfer check only. If all 25
canonical union preparations/endpoints succeed, the same previously used
class-extrema rule can define this new protocol's reference: Ca upper band is
maximum canonical Ca R, La lower band minimum canonical La R, requiring a gap
greater than the existing 0.02 model-kcal/mol qualification. No transfer fold may
set or adjust bands. Any missing member makes the new reference unavailable.
Report raw values/components, within-group spreads, errors/inconclusives, exact
coverage and static-versus-union matched changes. Do not interpret a larger
calibration gap alone as improvement or use biological labels to repair geometry.

Resources: one warm GPU worker with the existing 32-CPU/200000-MiB convention;
native solvent tasks use existing 8-rank runners and finite disjoint manifests.
No new project compute/time budget. Record actual preparation, allocation and
failed-attempt receipts. Root owns broader scoring/compare coordination; this
branch owns its new adapter, tests, diagnostic records and vault note only.
