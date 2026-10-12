# Matched context and accommodation: result

**Larger contexts change the score offset, but do not suppress the problematic
La-favoring donor response. Proposal choice also mattered in the earlier apparent
repair.** With every saved candidate available to both metals, the expanded
context keeps these two reference anchors correctly ordered by1.265264kcal/mol,
less than the3.472171 margin from its own three-candidate pool. This is a two-case
mechanistic result, not a new canonical calibration or accuracy claim.

All40 required composite endpoints are available. No optimization, new structure,
proton/water inventory, threshold or MACE weights were introduced. The original
strict208 correction remains complete and its accommodated calibration remains
unsupported. [Independent numerical review](../context_comparison_review/REPORT.md).

## What the matched experiment separates

Both source proteins use one pinned parent structure. q9z4j7 is a Ca reference
with145 versus192 atoms; c5axv8 is a La reference with172 versus192. Source heavy
atoms and hydrogens map exactly. Missing original entries are boundary caps,
whose replacements are explicit. Added fragments, caps and solvent cavity change
together, so a context effect is not uniquely electrostatic.

Original/expanded searches also selected a different fourth donor torsion:
q9 A350chi1 versuschi2; c5 A317chi2 versusA198chi1. The other three selected
modes are retained. We replayed every actual saved physical torsion vector in
each context, then scored the same common geometry pool for Ca and La. No cap
was independently moved; all context-exclusive atoms stayed fixed. Preparation
and exact maps: [context/REPORT.md](../context/REPORT.md).

The five physical candidates are origin, original-context Ca/La proposals and
expanded-context Ca/La proposals. Within each context, the score is

    E_M(q) = E_OMOL,vac,M(q) + E_GFN2,ALPB,M(q) − E_GFN2,vac,M(q)
    R_pool = min_q E_Ca(q) − min_q E_La(q).

The operational policy retains the origin when improvement is≤0.1kcal/mol.
Mathematical minima are reported separately; both policies give the same selected
candidates here. Ties follow the predeclared common geometry order. Larger R is
more La-like on this model's raw scale. No aquo offset or universal zero is
introduced, and no candidate count is interpreted as a population.

## Context × candidate-pool result

Values below are **R(c5 La reference)−R(q9 Ca reference)** in model kcal/mol.
Positive values mean the expected ordering of these two anchors only. The
origin/original-three values reproduce the corresponding prior boundary values;
other canonical members were not recomputed in this experiment.

| Common candidate pool within each source | Original context | Expanded context |
|---|---:|---:|
| Origin only |+8.867782|+12.175484|
| Origin + original-context Ca/La proposals |−4.556602|−1.212842|
| Origin + expanded-context Ca/La proposals |+0.437421|+3.472171|
| Origin + all four saved proposals |−3.714195|+1.265264|

Expanding the context while retaining the original proposals improves their
ordering by3.343760kcal/mol but leaves it reversed. Changing only the candidate
set can change the sign too. The common-five result preserves these findings;
we do not select the favorable reduced pool and claim search independence.

## The physical response that still needs explaining

At q9, the original La proposal moves R toward La by12.086147kcal/mol in the
original context and12.111762 in the expanded context: an interaction of only
+0.025615kcal/mol. Its MACE response changes+8.517684→+8.831655; the solvent
response changes in the opposite direction. Adding these fragments does not
provide appreciably greater resistance to that differential motion.

For the same q9 proposal in the expanded context, Ca lowers its energy by1.979582
and La by14.091344kcal/mol. La still selects that original proposal when it is
admitted to the common pool. Ca selects the expanded La proposal, whose Ca
energy is lower. Consequently the shared pool reduces the apparent expanded
reference margin by2.206907kcal/mol. Both metals have equal candidate access;
we do not force them to select the same geometry.

For c5, both metals select the expanded La proposal in the common pool. Its
differential motion response is−0.584073 in the original versus−0.514215 in the
expanded context. Across all eight nonorigin geometry/context comparisons,
changes in differential motion work are0.0256–0.3901kcal/mol. In contrast, the
origin-only context shift in R is−2.976718 for q9 and+0.330984 for c5.

These observations separate a context-dependent score offset from a correction
to the suspect relaxation response. They do not determine whether the remaining
problem is local electronic energetics, missing scaffold resistance, solvent,
state preparation, or a mixture. Geometry and existing compatible reference
calculations should guide the next discriminating test.

## Numerical qualification and actual work

MACE:24 verified archive reuses;16 new cross-context calls; four additional CPU
origin replays match archived GPU energies exactly. Same checkpoint/head/float64
and physical total-charge/multiplicity inputs; no neutral-charge adapter.

GFN2:80 actual primary xTB6.7.1 calls plus80 independent tblite0.7.0 calls,
including the scout; a separate four-call causal diagnostic isolates a nonsc
three-body dispersion difference. All80 absolute-energy comparisons fail the
original1e-6Eh parity gate. That failure remains visible. Because this term is
density independent and identical between same-geometry vacuum/water endpoints,
it cancels from each backend's solvent transfer. A separately frozen, inspected
development protocol qualifies the component actually used by this composite.

All40 transfers pass that1e-6Eh component gate; maximum residual1.24519e-7Eh
(0.00007814kcal/mol). All-atom population differences≤9.194e-8e pass1e-5e.
All16 available exact archived primary energies reproduce. No fitted correction
or estimated dispersion offset is used. See [solver report](../context_GFN2/REPORT.md).

The collector initially rejected the two producers' different XYZ file paths.
Direct checking showed only comment/decimal serialization differences, at most
5.5511e-17angstrom; the corrected join verifies elements/order and coordinates
under the already declared1e-12angstrom reuse gate, preserves both hashes and
still rejects state changes. Three real-artifact regression tests pass, including
an explicitly corrupted missing-cell copy that keeps the common-pool result
unavailable. Four MACE reuse tests and four solver parser/algebra tests also
pass. These file/algebra checks are separate from molecular execution.

| Component | Actual work | Measured process CPU | Other measured cost |
|---|---|---:|---|
| MACE |20 new calls,24 reuses|590.696s|113.884s elapsed;4,451,228KiB peak RSS|
| Primary xTB |80 calls|486.154s|Included in solver batch|
| Independent tblite |80 calls|1960.283s|Qualification overhead|
| Causal dispersion diagnostic |4 calls|43.859s|Separate failed-gate diagnosis|

The76-input post-scout solver batch elapsed1195.658s on two one-thread workers.
Total solver CPU2490.296s; peak single-process RSS111104KiB. ZeroGPU/Slurm work.
Thread limits were8 for root MACE and2 for this solver branch, within the user's
aggregate16 allowance alongside the separate native audit. These are development
measurements, not a new whole-protocol throughput benchmark.

## Decision and next work

Do not promote the expanded model from this two-anchor result or simply widen
the calibration band. Keep the unchanged geometry comparator and separate
corrected origin-only evidence visible. Inventory the existing expanded
candidate's known-reference calibration/transfer pools for corrected evaluation
without new optimization; retain the shared-pool sensitivity reported here.
In parallel, inspect the actual q9 donor/second-shell response and any compatible
archived DFT evidence to select a physical test. Both tasks are contained follow-
ups under Jacob's existing authorization. Larger candidate searches or tuning
against these two labels are not a demonstrated remedy.

Machine-readable unrounded data: [RESULT.json](RESULT.json),
[contrasts.tsv](contrasts.tsv), [effects.tsv](effects.tsv),
[anchor_ordering.tsv](anchor_ordering.tsv). For reproduction, invoke
`scripts/compare_context_factorial.py --help` with the exact MACE and solver
paths recorded in RESULT.json and a fresh output directory. All consumed
structures remain development evidence.
