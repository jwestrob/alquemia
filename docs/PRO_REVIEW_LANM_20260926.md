# Nikasha / LanM: model-selection review for GPT-6 Pro

Prepared 26 September 2026. **Start here.** This supersedes the older root
`GPT6_PRO_REVIEW_BRIEF.md` for this session's current question. That older brief
remains historical. Local reading anchor before this document: `6d80022`.

## What Jacob wants from this review

Help us choose a scientifically credible model **before another substantial
engineering campaign**. The current discussion concerns within-lanthanide
specificity in lanmodulin, beginning with Hans-LanM and a Mex-LanM comparator.
Jacob suspects coupled changes outside the direct coordination sites are central
and wants whole-protein structural response treated seriously. Treat that as a
mechanistic hypothesis to evaluate, not a result our calculations have proven.

The broader Nikasha La/Ca classifier is already useful for PQQ proteins and is
being applied to the PLM manuscript in **another session**. Preserve that work.
Accuracy matters more than throughput in this LanM branch, although failed,
hours-long calculations are not progress by themselves. A pocket-only method
that agrees on one structure but reverses on another is insufficient.

**No new LanM model has been chosen after the failures below.** Our owned
calculations are terminal; no automatic continuation or new search is queued by
this session. We are discussing models, not requesting another undirected survey
or automatic relaunch. Existing broad research authorization does not change
that immediate request for scientific discussion.

## Critical correction: the current model already includes the whole protein

It is inaccurate to describe this experiment as a carve with a frozen exterior.
The evaluated source is the complete Hans8DQ2 chain-A monomer, two ions in EF1/EF2,
all74 retained source waters and hydrogens:1889atoms. Every atom was allowed to
move. There are no synthetic caps or link atoms in this whole-chain representation.

### Current preparation and planned comparisons

Prepared, not all calculated: Hans8DQ2 (La-conditioned), Hans8FNR
(Dy-conditioned), Mex8FNS (Nd-conditioned), each at EF12, EF23 and EF1234 occupancy:
nine conditional systems, two metals per system. **Only Hans8DQ2 EF12 reached
whole-protein molecular evaluation.** Do not infer results for four-ion occupancy,
Mex or the other source from preparation manifests.

Source waters, proton inventories, actual constructs and coordinates are fixed
within each paired comparison. Standard source pH5 protonation was reused.
8DQ2's EF4 Na position is an explicit hypothetical fourth-Ln starting position;
it is not an observed four-Ln state. Complete monomers omit Hans dimerization,
unfolded/partially folded populations and occupancy equilibrium.

La inputs use physical singlet. Dy uses the declared maximum-spin hypothesis,
multiplicity1+5n for n Dy ions (11 or21); this is not validation of the ground
state, exchange coupling or spin-orbit treatment. Native GFN2 uses its distinct
f-in-core effective singlet representation. Physical/effective states are logged
separately; never silently replace one by the other.

### The actual proposal energy and motion

MACE-OMOL supplies vacuum energy and analytic forces on the entire system,
using the existing pinned checkpoint in float64 with physical charge/spin inputs.
Exact model/software hashes are in the compact evidence export linked below.
We independently minimize each metal's **vacuum MACE energy**, using L-BFGS-B:
60iterations maximum,180function evaluations maximum, per-Cartesian-coordinate
box ±0.45Å from origin, ftol1e-10, projected-gradient tolerance1e-3eV/Å.
Maximum possible displacement per atom is0.779423Å.

Thus protein mechanics are supplied by the learned potential, including its
predicted coupling between atom motions. There is no separate classical scaffold
force field and no harmonic spring added to the energy. However, there are also
no explicit bond-preserving constraints during search. Final admission checks
require source covalent lengths within0.8–1.2 of their starting values, preserved
C-alpha chirality/volume, and no new heavy-atom clash below1.2Å.

These are bounded local proposals throughout the molecule, not sampling of large
folding transitions, conformational populations or relaxation in solvent.
Elemental coverage and successful force evaluation do not validate metal
selectivity or long-range/cooperative response.

### The intended scoring energy

For metal M at geometry q:

```text
E_M(q) = E_MACE,vac,M(q) + E_GFN2,ALPB-water,M(q) - E_GFN2,vac,M(q)
Q = {origin, admitted La proposal, admitted Dy proposal}
R_P,n = selected E_Dy(P,n) - selected E_La(P,n)
D_n = R_Hans,n - R_Mex,n
```

Both metals must access the same compatible geometry pool. Retain mathematical
minima and the operational policy (retain origin for improvements below0.1model
kcal/mol). Larger D means relatively more La-selective Hans versus Mex. Equal
ion counts in this balanced exchange cancel common aqueous/element-reference
terms; raw cross-element R is not an affinity. Different occupancies require
separate state/reservoir accounting and cannot be compared as raw populations.

**Solvent affects final selection but supplies no forces to the proposal search.**
The subtraction requires a converged whole-protein vacuum GFN2 endpoint. Native
ORCA6.1.1 GFN2 uses the native mixer,300K smearing,TolE1e-10Hartree,MaxIter500.
This is a perturbative composite; it is not self-consistent solvated MACE or a
binding-free-energy calculation. A failed component makes the composite missing.

## What actually happened

### Earlier pocket result failed structural transfer

The original same-pH, protein-level Hans/Mex comparison gave the ordered site
vector D=(+4.866368,−0.264440,+12.499116)model kcal/mol using La-conditioned Hans.
With the actual Dy-conditioned Hans structure, it became
(−9.568461,−28.187049,−20.673113), with matched composition and the same Mex
reference. The reversal already appears in MACE; solvent does not remove it.
This establishes source sensitivity, not that all pocket models are impossible
or that a particular missing global mechanism is uniquely identified.

Read [the transfer report](../diagnostics/lanm_series_followup_20260923/DY_TRANSFER_REPORT.md)
before interpreting the earlier, apparently encouraging
[initial report](../diagnostics/lanm_series_followup_20260923/REPORT.md).
No independent site-affinity labels were invented or selected.

### Whole-protein MACE ran, but both relaxed geometries were invalid

129 actual MACE evaluations finished in451.420worker seconds. Exact native versus
checkpointed-adapter comparison gave identical energy and maximum force discrepancy
6.84e-14eV/Å. Peak framework GPU allocation was73,678,407,680bytes, including the
native full-system comparison. This establishes executable feasibility only.

Both60-iteration searches elongated covalent carboxylate C–C bonds from~1.52Å to
2.6–2.9Å (Asp85/Glu91/Asp107). Both were rejected. Hundreds of Cartesian
components reached movement bounds. The large energy decreases are not physical
accommodation evidence. With neither proposal admitted, only origin energies
remain; there is no successful relaxed whole-protein discriminator.

We have not established whether the cause is primarily model applicability,
starting chemistry/geometry, vacuum conditions, optimization coordinates, or a
combination. Merely constraining the bad bonds would prevent that symptom but
would not itself validate the remaining energy surface.

### Solvent GFN2 converged; vacuum GFN2 failed twice

| Attempt | Actual outcome | Allocated resources/time |
|---|---|---|
|1213018|MACE completed; all native calls failed MPI startup|458s,32CPUs,1H200|
|1216461|MPI fixed; ORCA MaxCore2000MB below~5.7GB requirement|308s,64CPUs|
|1216547|Optional Slurm memory envvar absent; no molecular work|0s recorded|
|1216564|Both ALPB endpoints converged; both vacuum endpoints failed500cycles|21141s,112CPUs|
|1217219|Same-metal ALPB electronic seeds really activated; both vacuum retries failed500cycles|39377s,112CPUs|

Valid origin ALPB energies: La−3040.435409827431Eh; Dy−3039.193803581625Eh.
State/charge/strict convergence checks pass. These values alone are not the
missing composite or an affinity comparison.

The retries preserved geometry, Hamiltonian, spin/charge, temperature, tolerances
and iteration count. Matching native GBW+xtbw files from each metal's ALPB result
seeded its vacuum run; both printed `INITIAL GUESS: XTBRESTART`. They still failed.
The final seeded energy residuals were approximately−0.174Eh (La) and−0.495Eh
(Dy), with nontrivial density residuals. These are not nearly converged energies
that we can accept by rounding or relaxed thresholds.

Higher process count did not demonstrate a speedup: first112CPU/four-cell run
used28ranks per endpoint; the two-cell retry used56ranks. The systems took many
iterations and the electronic paths differ. Do not infer scaling from job CPU
counts. Allocation costs, including the listed failures, total6,812,384CPU-seconds
and458GPU-allocation seconds; these are not measured CPU utilization or a
production benchmark. No new DFT was executed in this whole-protein phase.

## What is established, and what is not

- Whole-molecule energies/forces can run at this size. The current search did
  not yield chemically acceptable response, and the current composite cannot
  be evaluated because of its vacuum term.
- ALPB convergence does not validate its lanthanide energetics. Vacuum failure
  does not prove an incorrect structure or a unique physical explanation.
- Changing the initial guess alone did not fix the vacuum branch. We have not
  tested every possible solver, solvent Hamiltonian or protein mechanics model.
- No four-ion outcome, whole-protein Hans/Mex selectivity, Kd, occupancy
  distribution, binding entropy or cooperative free energy has been obtained.
- Spicy-Lams inventory found616assay sequences and complete mature-chain models
  for16reserved-panel proteins. Numerical reserved outcomes remain unopened.
  See [inventory](../diagnostics/spicy_lams_inventory_20260923/REPORT.md) and
  [evidence mapping](../diagnostics/spicy_lams_inventory_20260923/EVIDENCE.md).

## Requested scientific judgment

Please recommend a primary model and a small discriminating feasibility test,
with explicit reasons to reject alternatives. Address these questions:

1. What protein response/assembly/ensemble must be represented to compare the
   actual experimental observables? How far can fixed-occupancy monomer minima
   take us before missing folding or dimerization dominates interpretation?
2. Should we use a parameterized protein/solvent mechanics model with reliable
   local metal interactions, a polarizable whole-system potential, or a coherent
   ML/MM hybrid? What La/Dy parameters/training coverage actually support it?
   AMOEBA, MACE/POLAR and hybrid approaches are candidates, not validated choices.
3. How should solvent and polarization enter the **same energy/forces used for
   relaxation**, without requiring this failing whole-protein vacuum subtraction
   or double-counting interactions already represented by another component?
4. What local chemical checks would distinguish faulty metal energetics from
   inadequate scaffold freedom, without fitting springs to the desired labels?
   Covalent-preserving motion is useful but is not an energetic validation.
5. Which paired structures/occupancies would provide the cheapest informative
   rejection test before a larger series? Use the consumed Hans/Mex sources
   first; preserve source-conditioned results and reserve library validation.

A useful answer can reject our current model. Please avoid simply prescribing
more iterations, bigger CPU allocations, arbitrary springs, favorable-source
selection, or an unsupported model because it nominally includes lanthanides.
MACE need not supply every term, but an AI component must demonstrate useful
prediction or structural-response capability to merit inclusion.

## PQQ context and communication lesson (not this session's execution task)

The PQQ improvement is real development evidence: the strongest consumed
225-structure comparison improved same-solver static199correct/1wrong/8inconclusive/
17unavailable to207/0/1/17. The practical three-source4.3Å protocol separately
achieved91/100correct, or94/100 with named numerical recovery; correlated
structures are not independent biological observations.

The earlier~42seconds/source was an H200+32CPU static source-to-score benchmark,
excluding folding. The later practical accommodation pilot scored six sources
in251allocation seconds, but that small pilot is not cohort throughput.
The separate PLM session observed~8.9minutes/prepared three-source protein in an
older A5000/one-CPU job; scalar work was serialized. Its split CPU/GPU execution
is separate ongoing work, not a speedup established by this LanM session.
Do not use these different denominators/hardware as an exact DFT speedup claim.

## Reviewer evidence and reproducibility

Read these next, only as needed:

- [Compact whole-chain evidence export](../diagnostics/lanm_global_occupancy_20260923/PRO_REVIEW_EVIDENCE.json):
  actual hashes, preparation/model settings and collection summaries.
- [Original scope](../diagnostics/lanm_global_occupancy_20260923/PLAN.md),
  [MACE physical failures](../diagnostics/lanm_global_occupancy_20260923/FIRST_RUN.md),
  [first native result](../diagnostics/lanm_global_occupancy_20260923/NATIVE_RESULT_1216564.md),
  [seeded result](../diagnostics/lanm_global_occupancy_20260923/NATIVE_RESULT_1217219.md).
- [Whole-chain implementation](../scripts/lanm_global_occupancy.py),
  [seed recovery](../diagnostics/lanm_global_occupancy_20260923/seed_vacuum_recovery.py),
  [current status](../diagnostics/lanm_global_occupancy_20260923/CURRENT.md).
- [PQQ release](../diagnostics/pqq_fast_release_20260920/REPORT.md),
  [strongest structural transfer](../diagnostics/strict_native_transfer_20260923/REPORT.md),
  [practical SOP](PLM_PQQ_SOP.md).

Large workspaces/wavefunctions and parts of the legacy live environment are not
in GitHub. Compact exports preserve provenance; they do not pretend the public
checkout contains every executable/artifact. Actual pinned implementations under
`workspaces/lanm_global_occupancy_20260923/scoring_v2/implementation` control the
reported executions. Do not reconstruct missing outputs from prose or silently
attribute a pinned run to a different committed carver. Older checkpoint headers
may describe then-live jobs; terminal reports above supersede them.
