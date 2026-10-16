# LanM reciprocal-source force diagnosis — 16 October 2026

## Question and decision

Do the saved La/Dy forces identify a coherent structural response along the
observed Hans8DQ2 versus Hans8FNR protein differences that merits a finite test?
Jacob authorized return to La/Dy profiling after preserving the selected PQQ
policy. This first step reuses completed calculations and makes **zero new
molecular calls**. It is a diagnostic for choosing a useful physical experiment,
not another attempt to solve the already-converged electronic states.

The completed static-source comparison gives opposite conditional ordering:
Hans8DQ2 versus Mex +11.689857 kcal/mol; Hans8FNR versus Mex −18.269452.
Mex cancels in the difference between these two comparisons. A larger source gap
is not, by itself, evidence for either a particular geometry error or missing
protein relaxation. The completed Glu91 finite test changed the pool contrast
by only −0.173940 kcal/mol and is closed.

## Frozen inputs and state

Use only the exact six-origin provenance in
`diagnostics/claude_lanm_20261012/03_comparison_root_recovery/FINAL_PROVENANCE.json`
and its pinned final result; profile the two Hans sources and both metals.
Their current repaired source inputs are
`workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/`
`{Hans8DQ2,Hans8FNR}/INPUTS.json`. Reuse the admitted nuclear and point-charge
gradients, cap Jacobians and matching classical component gradients through
`scripts/metal_environment_force_assembly.py`. Synthetic caps are not independent
physical atoms. Verify native and source coordinate mappings before assembly.

Both Hans structures have the same 1,665 protein atom IDs/charge map (A24–133),
with the existing pH5 proton inventory. EF3 is exchanged; all source-specific
waters and EF1/EF2/EF4 spectators remain fixed. DQ2 has 2/74 QM/total waters and
La/La/Na spectators; FNR has 3/179 and Dy/Dy/Dy. There is no cross-source water
correspondence or matched-hydration claim. Electronic spectators use formal
monopoles; the classical ledger also retains ion-specific pure 12–6 terms.

The existing frozen-f PBE0-D4 electronic approximation and state limitations
remain: La ECP46; Dy ECP55 effective singlet valence, physical 4f9/sextet recorded
separately. No new electronic method, protonation, metal occupancy, cavity,
relaxation or solvent model is introduced. Full moving-scaffold hybrid forces
are not yet qualified as an optimization engine.

## Declared analysis

1. Assemble both metals' mapped physical gradients for each Hans source using
   the existing implementation. Reuse exact accepted results; fail explicitly on
   missing/incompatible input or component evidence.
2. Align FNR to DQ2 by equal-weight least-squares fitting of all exactly matched
   protein heavy atoms. Save identities and the rigid transform. Transform
   gradient vectors consistently into the same frame.
3. Define reciprocal Cartesian source-difference directions on the common
   protein atoms, including attached hydrogens. These are local diagnostic
   vectors, **not covalent-preserving finite moves**. Normalize by their physical
   Euclidean displacement norm before comparing loads across groups.
4. Report individual La and Dy gradient loads and Dy-minus-La loads, both whole
   protein and per residue. Keep electronic and classical contributions visible
   where available. Also report individual and differential force norms for
   protein, waters, target and spectator atoms without inventing water alignment.
   The EF3 sequence segment A83–94 is fixed before analysis; report remaining
   protein separately. Do not select residues by desired preference sign.
5. Any inventory of existing donor/hinge directions must use real source axes,
   atom maps and native covalent connectivity. Unsupported directions are absent,
   not filled with guessed Jacobians. No mode search or Hessian is authorized
   by this diagnostic.

For a unit direction u in physical Cartesian coordinates, report
`g_M dot u` and `(g_Dy − g_La) dot u` in kcal/mol/angstrom. Force has the opposite
sign. A negative directional gradient means local downhill energy under that
source's fixed state. Report both source-specific directions and signs clearly.
Do not multiply these local slopes by the full inter-crystal displacement and
call the result an accommodation energy. Do not sum different source states
into a population or infer an affinity from the source gap.

## Interpretation and next gate

This can localize where source geometry and differential loads meet. Agreement
of useful directions across sources would motivate one separately declared,
small covalent-preserving displacement per source and paired evaluation of
both metals. Divergent or diffuse loads may instead motivate a controlled
hydration or spectator-field question. No four-cell quantum batch is prequeued
from this plan. Neither load magnitude nor a favorable predicted sign establishes
an improvement in La/Dy discrimination. Preserve common strain and electronic/
state uncertainty; do not tune the displacement to an expected Hans label.

## Execution and preservation

Root owns submission. One standard Slurm CPU; no memory/time overrides, no GPU,
no molecular executable. Claude/Codex orchestration stays on the login host.
Pin the plan, adapter and source receipts in a fresh workspace. Emit JSON plus
a compact per-residue TSV and an actual execution receipt. Record scheduler
cost separately from reused molecular costs; do not report CPU utilization from
allocated slots. Arm the established completion wake. No reserved SpyCI-LAMBS
numerical profiles are read. No PQQ production or historical output is modified.
