# Parallel discriminator development, 2026-09-19

Jacob's current instruction: “Second shell modeling is crucial. This should be a
parallel track. ... let's get this and the water occupancy/entropy process each
going in parallel.” He also nominated the proposed structural-response track for
“additional parallel investigation.” The preceding cheap contextual preparation
by scorer comparison was agreed as first on the roadmap. Existing discretionary
execution authorization applies; preserve production and immutable experiments.

## Independent work streams

1. **Immediate scanner utility (root):** reuse both alpha-lactalbumin structures,
   all three GGR references, archived original/normalized/context-prepared DFT,
   and the frozen masked whole-protein MACE scorer. Transfer only actual prepared
   water H coordinates by source identity into the existing whole proteins.
   Eight new MACE energy calls: bound/detached Ca/La for each alpha structure.
   No new DFT, water inventory, protein geometry or threshold. Dry PQQ and GGR
   inputs are identity operations and retain archived scores. Report normalized
   water geometry separately: original DFT had stretched H whereas whole-chain
   MACE already normalized them. This is one consumed biological comparison.

   The four-state MACE interaction descriptor is
   (E_Ca,bound-E_Ca,detached)-(E_La,bound-E_La,detached). Metal-specific water
   orientations invalidate the common-environment cancellation used by the
   two-call shortcut; calculate all four terms. Also retain the separately named
   atom-reference-subtracted total bound contrast for component interpretation;
   never select a definition by outcome. Expected allocation: one GPU/16CPUs,
   minutes from prior small whole-protein calls, no arbitrary compute cutoff.

2. **Second shell (second_shell agent):** direct matched native r2SCAN-3c/CPCM
   and native OMOL on graph-completed source hydrogen-bond neighborhoods.
   Freeze exact cases, neighbor rule and endpoints in its own plan before runs;
   reuse previously computed alpha expansions. Fixed source coordinates/water
   inventory isolate context from geometry preparation. Changed core composition
   does not inherit an absolute reference or classification bands.

3. **Water basins (water_basins agent):** four existing 1F6S11/6IP9110 Ca/La
   development states, archived DFT-gradient-anchored OMOL potential, 298.15K.
   Two Sobol scrambles of1024 configurations per center,8192 cheap evaluations;
   up to16 deterministic representative native DFT EnGrad validations. Shared
   physical domains per Ca/La pair, nested COM radii0.30/0.45/0.60Å and rotational
   radii0.50/0.80/1.10rad. Exact proposal density and SO(3) measure, fixed scaffold,
   rigid water geometry/inventory. Necessary accuracy/convergence gates and
   representatives are frozen in the agent's plan. Wider conditional basin
   integrals are the question; missing internal-water/nonpolar terms remain
   unavailable. No unrestricted entropy or occupancy claim from local sampling.

4. **Response discriminator (khoury_benchmark agent):** reuse56 PQQ and8 direct
   archived MACE-POLAR-medium analytic core endpoints; existing compatible DFT
   direct gradients. Signed mean radial donor load and transverse-load fraction
   of grad(E_Ca-E_La), with physical source/cap mapping. Existing grouped PQQ
   folds; fixed ridge lambda1, each observable alone and each plus DFT contrast,
   versus DFT alone. No feature search or direct-site training with two groups.
   Zero new molecular energy calculations; saved compute serves other tracks.

New methods remain opt-in development. Accuracy and PQQ preservation are primary;
runtime is reported separately. No promises of improvement before results.
Independent agents own separate scripts/diagnostics/workspaces and vault notes.
No jobs killed, defaults changed, broad rescore, push or new package installation.
