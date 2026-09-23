# A8R3S4: loss of starting solvent margin, not increased relaxation

The two new inconclusive triples arise mainly from a **context-dependent shift
in the starting solvent correction**. Removing the Ser352 sidechain fragment
raises the Ca−La contrast before accommodation. The subsequent differential
accommodation is slightly smaller than with the ten-fold context, yet crosses
the frozen Ca edge because less initial margin remains. This refines the earlier
transfer report's shorthand about response in a smaller context; it does not
support a claim that the new optimizer produces excessive additional motion.

This is a read-only audit of all four declared A8 triples, seven distinct
source/membership pairs, four old ten-fold pools, and the changed canonical pool.
No new energies, forces, restarts, labels, thresholds or parameters were used.
A8R3S4 is a Ca-associated reference. All structures and labels were already
consumed during development.

## Exact physical membership

The ten-fold union contains174 atoms and12 caps. The168-atom triple union has11
caps and omits the neutral **chain A Ser352 sidechain fragment**. The149-atom
union has10 caps and additionally omits the **chain A Trp282 sidechain fragment**.
These are capped fragments, not deletion of whole source residues. The source
role map identifies Ser352 as the neutral homolog of the extra acidic ligand;
this audit does not reassign its chemistry. La charge stays−2, Ca charge−3;
water inventory, core membership and retained protonation rules are unchanged.

| Triple | La samples | Atoms | Ten-fold median R | Three-fold median R | New call |
|---|---|---:|---:|---:|---|
|008|1,2,3|168|−405467.346652|−405460.145287|Inconclusive|
|009|1,2,4|149|−405473.639675|−405470.244625|Ca|
|010|1,3,4|168|−405467.346652|−405460.145287|Inconclusive|
|011|2,3,4|168|−405473.639675|−405468.176499|Ca|

All four ten-fold calls are Ca. The two new inconclusive medians lie3.5670131596
model kcal/mol above the frozen Ca maximum,−405463.71230032056. The still smaller
149-atom triple remains correct, so context size alone is not a monotonic remedy.
No triple/member was chosen or discarded after scoring.

## Starting offset and accommodation are distinct

Define R=E(Ca)−E(La), with each endpoint's native OMOL plus native GFN2
ALPB−vacuum correction. Selected work W is candidate energy minus its own origin;
its effect on R is W(Ca)−W(La). Positive changes move toward the La region.
All table values below are model kcal/mol.

|168-atom source|Origin ΔR, new−old|Origin native ΔR|Origin solvent ΔR|Old accommodation ΔR|New accommodation ΔR|Final ΔR, new−old|
|---|---:|---:|---:|---:|---:|---:|
|sample1|+8.690253|+0.524511|+8.165742|+9.290776|+9.022723|+8.422199|
|sample3|+7.703223|+0.721941|+6.981281|+8.553178|+8.051320|+7.201364|

Both new origins still fall in the Ca region. Their existing La-favoring
relaxation then uses the lost margin. For sample1 the selected native work is
Ca−7.074644/La−16.050028 versus old−7.141000/−16.076759; sample3 is
Ca−8.830204/La−20.483650 versus old−8.727977/−20.344289. The corresponding
native differential work changes only+0.039626/+0.037134. The solvent correction
opposes the sample3 relaxation more strongly after contraction, not less.

The final pooled shift has the following component accounting across **all seven
pairs**, including the correct outcomes:

|Sample|Atoms|Native ΔR|Solvent ΔR|Total ΔR|New member call|
|---|---:|---:|---:|---:|---|
|1|168|+0.564137|+7.858063|+8.422199|Inconclusive|
|2|168|+1.464570|+3.998606|+5.463176|Ca|
|3|168|+0.759075|+6.442290|+7.201364|Inconclusive|
|4|168|+4.331474|+2.795871|+7.127345|Ca|
|1|149|+3.588289|+3.955372|+7.543661|Inconclusive|
|2|149|+7.396094|−4.001045|+3.395049|Ca|
|4|149|+7.214869|−4.558567|+2.656301|Ca|

For failing samples1/3, the ALPB metal contrast changes−25.456296/−29.299062;
the vacuum contrast changes−33.314359/−35.741351. Their difference gives the
+7.858063/+6.442290 solvent shift. These compare changed compositions/cavities
and their selected geometries; they are not an isolated interaction energy for
Ser352. A specific hydrogen bond, dielectric effect or electronic-solution
artifact is not established by subtraction alone.

## Selected motions and candidates

For both failing168-atom members the selected physical subspace is unchanged:

- sample1: Glu213 χ1/χ3/χ2 and catalytic Asp350 χ1. Both metals select the
  Ca-proposed geometry before and after contraction.
- sample3: Glu213 χ1/χ3/χ2 and Asn300 χ2. Ca selects the Ca proposal and La the
  La proposal before and after contraction.

The largest change in the corresponding selected torsion coordinates is
0.003487rad for sample1's Ca proposal and0.003817rad for sample3's La proposal;
sample3's Ca proposal changes by at most0.011639rad. Full native proposal angles
and physical checks are pinned in the audit. Both cases retain their original
Ca-proposal nonboundary/La-proposal boundary statuses. No exact minimum is claimed.

Among the other members, sample2 changes its fourth mode from Asp350 χ1 to χ2
in both smaller contexts; its168-atom selection changes both metals from the La
proposal to the Ca proposal. Sample4 changes Asp350 χ1 to χ2 only in the149-atom
context; both metals still select its La proposal. These mode switches do not
explain the two failing sample1/3 outcomes.

## Native scalar convergence and continuation identity

All144 distinct scalar outputs considered here terminated normally using the
native energy convergence check. The actual printed density criteria were
usually unmet, in old and new contexts alike:

|Inputs|Cells|SCF cycle range|Printed energy criterion passes|MAX density criterion exceeded|RMS density criterion exceeded|Retained GBW+xtbw pairs|
|---|---:|---:|---:|---:|---:|---:|
|Old four ten-fold pools|48|17–108|48|48|48|48|
|New seven source/context pools|84|19–190|84|84|82|84|
|Changed three-fold canonical|12|22–71|12|12|12|12|

Raw unrounded energies, actual input/runtime controls, charge/electron counts,
occupations, last five printed iterations, output/receipt pins and both saved
seed hashes are retained. Normal energy-only termination is not evidence of
well-converged densities or stable derivatives. Its presence in both methods
also does not demonstrate that the A8 regression is numerical. No entropy term
or favorable electronic branch is inferred.

The separate32-source continuation uses **old174-atom A8 canonical sample0** from
the precision34 archive. It contains neither the changed168-atom canonical pool
nor noncanonical samples1/3. None of its A8 XYZ hashes equals a changed pool's
input. Thus that continuation cannot directly qualify this representation change.
The actual saved seeds make a separately declared matched-state continuation
possible; none was executed in this audit, and no new follow-on is assumed.

## Reproducibility and interpretation

Run the exact read-only command in [A8_AUDIT_COMMANDS.md](A8_AUDIT_COMMANDS.md).
[Compact pins](A8_AUDIT_PINS.json) identify the complete workspace audit, source
comparison, implementation and continuation inventory. The original v1 audit is
preserved; v2 adds actual proposal-angle receipts, without any new molecular work.

The supported diagnosis is loss of the starting solvent margin upon fragment
contraction, followed by nearly the same local accommodation. Restore neither a
hand-picked residue nor a threshold on this evidence alone. Keep the negative
three-fold fidelity result visible while separately resolving whether the solvent
shift reflects reproducible changed-context physics or native solver sensitivity.

## Source-distance check for a general motion envelope

A separate read-only geometry pass uses the **actual frozen donor-anchor sets**
from each source's original preparation, including completed functional-group
partners. Ser352 OG's nearest anchor is Asn300 OD1 for every saved origin and
admitted proposal considered. The omitted OG stays at its exact source coordinate:
none of the actual selected motions moves it. Distances below are Å; Ca/La columns
refer to the respective proposed geometries, not newly relaxed structures.

|Source|Context|Origin|Ca proposal|La proposal|
|---|---|---:|---:|---:|
|La sample1|ten-fold /168 /149|3.514751|3.514751|3.514751|
|La sample2|ten-fold|3.546596|3.562384|3.570734|
|La sample2|168|3.546596|3.562866|3.574018|
|La sample2|149|3.546596|3.560620|3.573312|
|La sample3|ten-fold|3.512943|3.515343|3.522099|
|La sample3|168|3.512943|3.516482|3.522716|
|La sample4|ten-fold /168 /149|3.502589|3.502589|3.502589|

Every origin lies just beyond the existing3.5 Å polar-neighbor cutoff. None of
the actual admitted endpoints brings this pair inside3.5 Å. These data therefore
show a sharp membership cutoff near a persistent polar contact, not evidence of
a newly formed contact during accommodation.

A **general proposed question** is whether freezing the complete fragment union
from a motion envelope reduces this membership sensitivity across all proteins.
For original anchors that may move by at most0.8 Å relative to an omitted fixed
neighbor, the triangle inequality gives3.5+0.8=4.3 Å as a conservative initial
search envelope. All four observed Ser contacts lie inside it. If both members
can move independently by0.8 Å, the analogous bound is5.1 Å; the4.3 Å argument
must not be advertised for all moving pairs. It also does not include newly
eligible donor-anchor identities or establish continuum-solvent convergence.

This is a physical, label-independent preparation proposal, not a radius fitted
to A8 outcomes. The next preparation must apply one frozen rule to the whole
existing panel, report changed charges/graph closure/unsupported inputs, and
create its own canonical reference before interpreting new transfer scores.
The geometry audit adds zero molecular calls and changes no current preparation.
Four focused actual-artifact tests cover complete scope, energy decomposition,
retained seed identities, continuation mismatch and exact geometry replay.
