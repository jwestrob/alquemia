# LanM La/Dy conditional selectivity pilot — preparation, not launched

Jacob requested within-lanthanide affinity testing. Root is discussing La/Dy
solution affinity versus La/Lu resin selectivity with him; these are different
experimental labels. This document prepares **La/Dy only**. Molecular execution
awaits that choice and the explicit electronic-state policy below. No production
default, historical input or La/Ca reference changes.

## Question and declared sources

Does the current native OMOL plus native GFN2 ALPB-minus-vacuum descriptor capture
Hans LanM's stronger light-over-heavy preference relative to Mex LanM? Use the
three primary sites in sequence order, without favorable-site selection or
weights: Hans 8DQ2 author chain A, La201–203; Mex 8FNS author chain A, Nd201–203.
EF4 is not one of these three high-affinity sites. Crystal conditioning remains
explicit: Hans is La-bound, Mex Nd-bound. Both test metals use exactly the same
nuclear coordinates, ligand protonation, caps and source-water inventory within
each site. This is a vertical comparison, not a sampled affinity calculation.

Reuse the actual Hans pH5 peptide-amide-v3 preparations from September15. Prepare
Mex with the same pH5 standard-residue protonation, deposited first-shell-water
policy and 3.3 Å generic selection followed by the same amide repair. No missing
residues are synthesized. Source heavy atoms must agree, supported standard
hydrogens/caps are explicit, and missing nearby heavy atoms or unsupported states
remain unavailable. Nd source selection is a private adapter extension of the
existing geometric selector; it is not a change to production metal support or
a computed Nd endpoint. All newly created intermediates are preparation-only.

## Two distinct electronic-state representations

Native float64 OMOL uses its released checkpoint, unmodified charge and spin
embeddings, physical La(III) singlet and Dy(III) sextet (4f9, S=5/2). Formal total
charge is the same for La and Dy in each site. Check element coverage, exact
input batch multiplicity and all-electron parity; element coverage alone is not
accuracy validation.

Installed native GFN2 exports three valence electrons with 5d/6s/6p shells for
both elements; it has no explicit 4f shell and its original energy is spin
independent. Use an explicitly named **effective valence singlet** for both in
this f-in-core solvent-transfer approximation. Record physical multiplicity and
effective solver multiplicity separately. Never pass sextet as if its five
unpaired f electrons belonged to this valence model. Check actual output electron
count, parameters, charge closure and state. This approximation does not model
Dy spin-dependent ligand interactions or spin–orbit thermodynamics.

Each GFN2 cell uses native mixer, ALPB(Water) or vacuum, 300 K electronic
smearing, MaxIter500 and ordinary unchanged tolerances. First run NoAutostart;
then exactly one same-geometry/same-medium continuation using that cell's matching
GBW and xtbw, with AutoStart positively confirmed by `INITIAL GUESS: XTBRESTART`.
The continuation is the declared reported value, not the lowest observed value.
Keep both results and their difference. Failed or unconfirmed continuations stay
unavailable; no escalation or silent primary-value fallback. Continuation alone
does not prove complete SCF branch convergence.

## Balanced observable and interpretation

For protein P/site i/metal M, in consistent energy units:

`E(P,i,M) = E_OMOL,vac(P,i,M) + E_GFN2,ALPB(P,i,M) − E_GFN2,vac(P,i,M)`.

Report each ordered site's balanced double difference:

`D_i = [E(Hans,i,Dy) − E(Hans,i,La)] − [E(Mex,i,Dy) − E(Mex,i,La)]`.

Positive D means Hans is more La-selective relative to Mex in this electronic
descriptor. Both proteins occur once on both sides of a metal-exchange cycle;
aqueous ion references and element-dependent atomic offsets cancel. Different
protein compositions/water counts do not require guessed water chemical
potentials because each protein's composition is unchanged by its metal swap.
No absolute aquo reference, absolute affinity, Kd, calibrated band or population
is claimed. Keep EF1/2/3 as an ordered vector. The experimental comparison is
protein-level apparent affinity coupled to folding/dimerization; these are not
three independent site-affinity labels. Crystal conditioning and incomplete
ensemble/entropy remain limitations.

## Finite calls and capability gate

Six sites ×two metals: **12 native MACE evaluations**. Two media per endpoint:
**24 GFN2 cells ×initial plus one continuation =48 GFN2 calls**. No DFT, dynamics,
optimization, extra electronic states, extra folds or fitted site weights.
First actual capability checks are the declared Hans EF1 Dy MACE endpoint and
its vacuum GFN2 cell, reused in those totals. If physical batch spin or effective
valence accounting fails, stop with unsupported status rather than substitute a
different Hamiltonian. No extra capability calculation is hidden in preparation.

Warm MACE uses one H200,32CPUs/200000MiB; native GFN2 uses existing8-rank workers
up to8concurrent,64CPUs/128GiB, through the existing runner. Actual recent native
restart job1210185 used27s on64CPUs for8calls on a larger128-atom PQQ context;
that gives an indicative six-wave48-call scale of minutes, not a guaranteed
LanM duration. New51–60-ish atom counts must be measured from preparation before
final resource prediction. Preserve all failures and allocated costs; scheduler
waiting is separate from molecular runtime. No submission is part of this file.

## Primary evidence

[2023 Hans/Mex primary study](https://www.nature.com/articles/s41586-023-05945-5):
Hans apparent La affinity68pM versus Dy main response2.6nM atpH5, with partial
Dy-induced folding and heterogeneous sites. Mex has weaker variation across the
series. Exact same-assay Mex La/Dy numbers/uncertainties are being recovered;
do not fit to approximated values. Retain site and construct qualifications.

[ORCA6.1 native semiempirical documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html)
and the installed parameter exports establish the actual implementation;
[GFN2 primary method](https://doi.org/10.1021/acs.jctc.8b01176) describes its
f-in-core lanthanides. The existing La/Ca protocols remain unchanged.
