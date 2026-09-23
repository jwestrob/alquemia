# Small whole-protein LanM occupancy pilot

Jacob: “Try a few with two-ion occupancy and a few with four. Tinker around a
bit with a small number of structures—maybe start with Hans-LanM—and let's see
if we can get specificity within the lanthanide series before we blast.”
Accuracy takes priority over metagenomic throughput in this branch. This is
authorization for the contained pilot below; not a whole-library rescore.

## Question and fixed starting material

Does bounded motion of the complete LanM chain improve source-structure
consistency and the relative La/Dy selectivity signal that failed with separate
44–50-atom EF-hand cores? Retain the existing PQQ pipeline unchanged.

Use Hans chain A from both 8DQ2 (La-conditioned) and 8FNR (Dy-conditioned), plus
the existing Mex 8FNS chain A (Nd-conditioned) as a same-occupancy exchange
comparison. Actual mature constructs are retained, not sequence-normalized.
These are three already consumed real structures, not reserved SpyCI-LAMBS inputs.

For each source prepare exactly three conditional states:

- Two ions in EF1 and EF2.
- Two ions in EF2 and EF3.
- Four ions in EF1–EF4.

Both La and Dy evaluate each state: nine source/occupancy systems, 18 endpoints.
The EF4 Na position in Hans 8DQ2 is an explicitly proposed fourth-lanthanide
starting site, not an observed four-Ln state; Hans 8FNR actually deposits Dy at
all four sites, including EF4. Source identities remain explicit. All chain-associated source waters
are retained within a source across occupancies; no score-dependent water choice.
Exact missing-atom/construct/terminal conditions and any unsupported source are
reported before scoring. Existing standard pH 5 protonation is reused where
compatible. No synthetic carve/link atoms and no fixed pocket boundary.

This first representation is the complete **monomer**. It tests coupled motion
within the chain, not Hans dimerization or a folding/unfolding equilibrium. The
conditional monomer limitation remains relevant to protein-level affinity labels.

## Method and response

Use the existing native MACE-OMOL checkpoint in float64 with physical total
charge. La is singlet. Dy uses the declared maximum-spin coupling hypothesis,
multiplicity `1 + 5*n` for n Dy(III) ions; other exchange-coupling and spin–orbit
states are not surveyed. Native GFN2 retains its documented f-in-core singlet
representation. Check physical and effective electron parity separately.

Each metal proposes lower vacuum-MACE energy independently from the same
coordinates. All protein atoms, metals and retained waters can move. Initial
finite search: L-BFGS-B, at most 60 iterations/180 evaluations, per-coordinate
box ±0.45 Å (maximum atom displacement 0.779423 Å), objective relative to q0,
`ftol=1e-10`, projected-gradient tolerance `1e-3 eV/Å`. No artificial spring
energy is added. Covalent connectivity, source C-alpha stereochemistry and
new severe nonbonded clashes are checked before admitting final candidates.
These are bounded proposals, not certified composite minima or thermal samples.

Pre-execution geometry admission gates: each source covalent bond length must
remain within 0.8–1.2 of its prepared length; each non-glycine C-alpha oriented
N/C/CB volume must retain its sign and at least 0.2 of its prepared magnitude;
no new nonbonded heavy-atom pair may approach below 1.2 Å; and the declared
coordinate box must be respected. These reject changed connectivity/stereochemistry
and severe clashes without fitting bond lengths to a selectivity outcome. Identity
checks pass on all nine actual prepared states, with no initial severe clashes.
Exact checkpointed/native MACE qualification tolerances are 1e-6 eV in energy
and 1e-6 eV/Å in maximum force component on the first full-chain source.

Both metals score the SAME finite pool {origin, La proposal, Dy proposal}:

```
E_M(q) = E_native-MACE,vac,M(q)
       + E_native-GFN2,ALPB(water),M(q) - E_native-GFN2,vac,M(q)
R_P,n = selected E_Dy(P,n) - selected E_La(P,n)
D_n = R_Hans,n - R_Mex,n
```

Positive D means stronger *relative* La preference for Hans than Mex under the
same declared occupancy pattern. Equal ion counts cancel element-reference and
common aqueous-reference terms. Do not interpret raw R as absolute affinity or
compare two- and four-ion raw totals as populations. Keep each Hans starting
structure separate; do not pool different source-water inventories.

Retain mathematical minima and operational selections (origin retained for an
improvement below 0.1 model kcal/mol), static and accommodated components, chosen
geometries, boundary hits, per-site donor/water distances and full-protein motion.
A required failed cell leaves the result unavailable. No label-specific repair,
threshold fitting, or substitution of a pocket result.

## Execution and numerical feasibility

Reuse the working warm MACE and native ORCA runner/memory/thread infrastructure.
Begin with one actual prepared Hans two-ion system to measure full-chain gradient
memory/runtime and native GFN2 feasibility. Numerical equivalence of the existing
exact checkpointed MACE adapter is checked against a native call on this same
full source if it fits. Preserve failed attempts and exact receipts.

Native GFN2: ORCA 6.1.1, fresh NoAutostart, ALPB water/vacuum pair, native mixer,
TolE 1e-10 Hartree, MaxIter 500, electronic smearing 300 K. Start with explicit
resource layout after source atom count is known; no external xTB substitution.
Initial numerical checks retain 0.1 kcal/mol component and 0.2 kcal/mol paired
contrast tolerances. Memory-only implementation changes preserve the Hamiltonian.

The full declared matrix has at most 54 composite endpoints / 108 scalar calls,
with exact same-run reuse/deduplication. MACE optimization calls are additional
and counted. Prepare the full finite manifest, but execute the first feasibility
cells before committing the remaining global solvent work. This is staging based
on measured feasibility, not a project-wide time or compute cap. No new DFT,
full quantum Hessians, new folds, FEP, or long trajectories in this pilot.

## Decision

First judge physical admissibility and source consistency, then the relative
Hans/Mex direction under each fixed occupancy. A favorable single state does not
establish within-series prediction, and an unfavorable state is not discarded.
Retain all occupancies, missing cases and costs. An interpretable result is needed
before extending to other ions or experimental orthologs. The reserved numerical
SpyCI-LAMBS profiles remain unopened. Root owns scoring/integration; Khoury owns
source-backed whole-protein preparation only.
