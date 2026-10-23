# La and Dy selectivity through validated free energy calculations

21 October 2026. **Implementation plan and goal handoff.** Read [CURRENT.md](CURRENT.md) for authorization history, ownership and actual evidence. This revision answers the request to write down the complete strategy; it does not itself submit calculations or change goal state.

## At a glance

**Question:** What determines La versus Dy binding when the protein, water and electronic interactions can all respond?

**Calculation:** Replace La with Dy in the protein and in water, then subtract the two free-energy changes. This requires ensembles of configurations, a consistent interaction model and a defined experimental comparison.

1. **Specify the target.** Fix the construct, protonation, solution conditions and measured binding quantity. Keep unresolved assay assumptions explicit.
2. **Establish trustworthy interactions.** Compare both metals on the same water and donor configurations. Test energies, forces and reference-method uncertainty before choosing or repairing the sampling model.
3. **Build the complete simulated system.** Use the intact peptide or protein, explicit water and ions with consistent polarization and protein mechanics. Validate the numerical implementation and metal-transformation path.
4. **Sample and calculate exchange.** Run independent replicas and intermediate metal states in both environments. Resolve hydration and conformational transitions, demonstrate statistical overlap, and add quantum corrections only with adequate sampling.
5. **Explain and validate.** Identify which physical changes repair the previous failure, freeze the method, and test independent experimental cases before extending to multi-site LanM.

**Done means:** a reproducible exchange prediction with converged sampling, quantified model limitations, an experimental mapping and independent validation. A matching sign or a completed calculation alone does not satisfy the goal.

**Short goal to type:** Follow the saved La/Dy exchange plan through a validated free-energy prediction and explain the previous model's failure.

## Objective

Develop a physically consistent prediction of La/Dy preference by calculating metal exchange between a protein and water. Establish credible interaction energies and forces, sample the relevant states, and connect their free energies to an experimentally measured quantity.

The first milestone is one defensible prediction for an experimentally characterized monomer, preferably with a single binding site. Independent validation follows before extending the method to multi-site natural LanM and assembly-dependent assays. Success requires useful uncertainty and an explanation of the mechanisms; matching one known label is insufficient.

## The approach in plain language

Calculate the cost of replacing La with Dy **inside the protein**, then subtract the cost of making the same replacement **in water**. Both calculations must allow the surrounding molecules to respond and sample their relevant states. The difference gives the metal preference under the specified chemical conditions.

Our current static energies and limited relaxations do not establish that quantity. The new work must determine whether the main failures come from local metal interactions, solvent and protein response, missing state populations, or the comparison with the assay. Changing MACE checkpoints or repeatedly optimizing a starting structure does not answer those questions by itself.

The practical route is a qualified polarizable model for extensive explicit-water sampling, with quantum reference calculations testing its interactions and supplying statistically valid corrections when needed. DFT is a reference candidate with its own uncertainty, not an automatic ground truth. The affordable model is selected by these tests; no replacement MACE model is a prerequisite.

## Proposed goal

Develop and validate a coherent explicit-solvent free-energy method for La/Dy selectivity, starting with an experimentally characterized monomer. Qualify the metal, ligand, protein and water interactions against independent reference evidence; establish converged bound-versus-aqueous metal-exchange free energies; identify why the previous approach failed; and test the frozen method on independent experimental cases. Preserve the existing baseline and report both sampling and model uncertainty before extending to multi-site or assembly-dependent LanM predictions.

## Scientific starting point

The existing calculations establish conditional energies and forces for selected chemical states. They do not yet establish the complete equilibrium preference measured experimentally.

On the saved late La relaxation, both MACE and PBE0 vacuum energies increase while the shared fixed-charge ALPB term makes the total decrease. Reusing that solvent term in both comparisons does not independently validate it. This motivates testing solvent response, but does not establish that solvent is the sole error. See the [completed reference comparison](../lanm_coupled_mechanics_20261017/lbt_7ccn_converged_reference_v1/FINAL_REPORT.md).

The [experimental target record](../lanm_coupled_response_20261016/EXPERIMENTAL_TARGET.md) distinguishes apparent Hans/Mex affinities from microscopic site affinities. Folding, occupancy and assembly can affect the former. The [new primary-source review](experimental_target_v1/REPORT.md) identifies LBT3 as the development target, with explicit qualifications:

- Primary construct: **NH3+–FIDTNNDGWIEGDELLA–CONH2**, peptide charge −4 and metal-complex charge −1 for the declared protonation state. The 9CEQ solution structure supports this amidated construct.
- Historical free-carboxylate LBT3 remains a separate sensitivity calculation, not a silently corrected baseline.
- The 2020 ITC data provide matched La/Dy free energies at 283.15 K, pH 6, but terminal chemistry is unresolved. The 2025 assay specifies amidation but leaves temperature and numerical uncertainty incomplete. Keep the assays separate.
- These are consumed development data. Independent validation remains to be registered; another structure of LBT3 is not an independent test.

## Stages and progression criteria

| Stage | Deliverable | Criterion for progressing |
| --- | --- | --- |
| 1. Define the observable | Experimental target and state specification | The calculated quantity can be mapped to the measurement |
| 2. Qualify the chemistry | Matched La/Dy reference set with uncertainty | Interaction errors are acceptable for the intended selectivity resolution |
| 3. Implement the coupled model | Consistent protein, metal and explicit-water energies and forces | Numerical correctness and physical reference tests pass |
| 4. Establish sampling | Reproducible conformational and hydration populations | Independent simulations agree and relevant states communicate |
| 5. Calculate exchange | Bound-minus-aqueous exchange free energy | Statistical overlap, convergence and correction checks pass |
| 6. Diagnose the mechanism | Controlled comparisons explaining failures and improvements | Conclusions survive uncertainty and alternative state assumptions |
| 7. Validate and transfer | Predictions on independent cases | Useful accuracy without fitting to those outcomes |

These are scientific progression criteria. During authorized goal execution, passing a stage allows autonomous continuation within its scope; it does not introduce a new per-analysis permission requirement. A failed criterion directs recovery or method development and remains recorded.

## 1. Define the experimental observable

Select a monomer with verified La and Dy equilibrium measurements, preferably a single binding site. Record sequence and construct, temperature, pH, ionic strength, buffer and chelator conditions, free-metal activity conventions, stoichiometry, assay interpretation and experimental uncertainty.

Determine whether the measurement represents a microscopic binding constant, an apparent conformational affinity or a competitive capture response. If a simple target is unavailable, explicitly construct the required occupancy or assay model before claiming quantitative agreement. Do not substitute an isolated-site score for a protein-level observable.

Separate development cases from genuinely uninspected validation cases. Historical LBT and Hans/Mex outcomes are consumed evidence, and structures of the same protein are not independent validation proteins. Set the desired resolution and acceptance criteria before evaluating validation outcomes.

Prepare an explicit observable contract: free aqua-ion activity versus analytical metal concentration, buffer and anion complexation, protonation assumptions, and which conformational states are included. A box containing water and NaCl alone does not establish equivalence to a MES-buffered competition or ITC assay. Resolve or bound those differences before claiming quantitative experimental validation.

Include possible buffer or anion binding to the metal–peptide complex as well as to unbound metal. A correction for free-metal speciation alone does not account for a populated ternary complex; test or bound its contribution under the assay conditions.

The initial production design target is a **95% sampling-confidence half-width no larger than 0.5 kcal/mol** for the final exchange estimate, with agreement across independent starts. This is a planned precision requirement, not an achieved uncertainty or a tolerance for every individual interaction. Reference and model-error limits must be justified against the experimental separation and recorded before production or held-out evaluation. Report systematic sensitivities separately; do not shrink them into the sampling error bar.

## 2. Qualify the metal interaction physics

Build a finite, matched reference set for both metals covering hydrated ions, carboxylate and amide donors, ligand competition with water, coordination changes, and representative protein configurations in their surrounding electric fields. Include configurations away from optimized minima.

Use balanced ligand-exchange reactions and forces to test metal-dependent response. Preserve chemical composition or account explicitly for exchanged water and protons. Absolute electronic energies of different elements or differently capped clusters are not binding selectivities.

Qualify the electronic treatment: relativistic basis and effective-core conventions, Dy electronic state, SCF convergence, numerical grids and boundary sensitivity. Compare appropriate small cases with independent electronic calculations and experimental hydration or coordination evidence. Treat PBE0 and other DFT methods as candidates with uncertainty. Use existing admissible results where compatible.

Reference admission must also include basis and basis-set-superposition sensitivity, force checks against finite energy differences along several relevant motions, and electronic-state sensitivity. Establish whether frozen-f versus explicit-f treatment and spin–orbit effects change the exchange contrast at the required resolution. A comparison that also changes basis and relativistic representation cannot isolate the f-electron contribution by itself.

If multiple Dy electronic states contribute at the experimental temperature, include their energies and degeneracies in the bound and aqueous state populations, or justify neglect or cancellation at the required resolution. Comparing selected electronic states alone does not establish this thermal contribution.

Start with the [prepared matched reference panel](reference_chemistry_v1/READY.md): water, carboxylate and amide dimers, radial distortions, and hydrated ions. The first eight native cells qualify execution and provide initial comparisons; they do not qualify protein selectivity. Extend to hydrated mixed-donor configurations, competing coordination numbers, external electric fields and actual sampled protein configurations. Keep reference configurations used for fitting separate from those used for validation.

Prioritize errors in the La/Dy difference and in physically relevant motions. A large shared force error need not imply an equally large selectivity error. Freeze a reference set for validation separately from any data used to parameterize the affordable model.

## 3. Implement one consistent explicit-solvent model

The proposed architecture uses a polarizable molecular model for extensive sampling, with embedded quantum calculations validating and, where needed, correcting binding-site physics. Polarization allows the charge distribution to respond to the environment.

Published AMOEBA La parameters provide a concrete starting candidate for water, acetate and acetamide interactions. They do not establish a qualified paired La/Dy protein model. Check parameter compatibility and both metals against Stage 2 before selecting the production research Hamiltonian. [Primary La parameter study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9957963/)

Implementation must preserve the actual parameter functional form, including the published pair-specific polarization damping. The current native Tinker route can represent it. Available Dy parameters came with a different water model; transferring them onto the La water/protein model is a hypothesis to test. Do not copy La donor-specific parameters onto Dy or combine parameter sets merely because they execute. Reproducing a published La-water dimer is an implementation control, not qualification of the paired model.

Represent the full chosen construct, explicit water and ions in one consistent energy model. Validate energy derivatives, electrostatics, protein deformation costs and any quantum/classical coupling. Account for induced polarization and avoid missing or duplicated interactions across boundaries. A quantum region must include the chemistry needed for the question; test its convergence with region size and treatment of exchanging waters. Do not silently change water identities across a discontinuous boundary.

Explicit solvent replaces the implicit-solvent contribution for the represented environment. Do not add the existing ALPB term to the same explicit-water system. If the sampling model fails the chemistry tests, revise its parameters or coupling and repeat the affected qualification tests before scaling up.

## 4. Sample the relevant states

Use existing structures as starting configurations for independent La and Dy simulations. Include multiple initial solvent arrangements and independent replicas. Sample protein conformations, donor orientations, coordination changes and waters entering or leaving the binding region.

Use enhanced sampling when ordinary trajectories remain trapped. Choose the biased coordinates from the physical transitions being tested, remove the bias correctly, and verify that important hidden barriers remain sampled. Simulation duration alone is not an acceptance criterion.

Handle plausible protonation and occupancy alternatives through explicit state free energies and the appropriate proton or metal chemical potentials. Raw energies of different chemical inventories cannot provide their populations. Apply the same care to alternative folding or assembly states when they are part of the observable.

Require agreement between independent replicas and starting states, repeated relevant transitions or otherwise validated state connections, and stable population and free-energy estimates as sampling increases. Declare convergence criteria before examining the resulting preference.

Begin production design with at least three independent replicas per thermodynamic leg, with independent solvent configurations and seeds; pilot throughput determines the initial trajectory lengths and intermediate-state count. Inspect hydration number, donor identity, ion association and protein conformational populations as well as the scalar free energy. Increase sampling where these diagnostics expose missing transitions. Account for time correlation when estimating uncertainty.

## 5. Calculate balanced metal exchange

The single-site reaction is:

**protein–La + Dy(aqueous) ⇌ protein–Dy + La(aqueous)**

Define the exchange free energy as:

**ΔG_exchange = ΔG_bound(La→Dy) − ΔG_aqueous(La→Dy).**

For this convention, a negative value favors Dy binding relative to La. For matched microscopic dissociation constants, ΔG_exchange = RT ln[Kd(Dy)/Kd(La)]. Apparent multi-site or conformational affinities require the additional mapping in Stage 7.

Use intermediate interaction states to connect La and Dy in the bound and bulk-water environments. Both ions retain charge +3. Use consistent physical conventions, replicas, overlap diagnostics, forward/reverse consistency and statistical analysis in both legs. Account for restraints, standard states and finite-size effects where applicable; equal charge alone does not establish that every correction cancels.

Specify how every metal parameter and induced polarization changes along the intermediate states. Test polarization convergence and finite energies/forces throughout the path before production. Use an estimator such as BAR or MBAR with adequate adjacent-state overlap; require stability under added intermediate states and longer trajectories. Test periodic box size and counterion association, and account for any restraint used to define the bound state. Diagnose a failed path rather than accepting a numerically returned free energy.

If a quantum correction is required, evaluate sampled configurations in both bound and aqueous environments and calculate the corresponding free-energy corrections. Demonstrate that the affordable and quantum models sample sufficiently overlapping distributions. A mean singlepoint energy difference is not generally the free-energy correction. Failed overlap calls for additional intermediate models, improved sampling or an improved sampling potential. [Primary QM/MM correction methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC6834343/)

Define each correction as the free energy of changing the sampling Hamiltonian into the specified target Hamiltonian for that environment and metal. The complete cycle is:

**ΔG_target = ΔG_model + δG_bound,Dy − δG_bound,La − δG_water,Dy + δG_water,La.**

All four terms need compatible definitions. Check effective sample size, concentration of reweighting weights, independent replicas and boundary/region convergence. If a few frames dominate, bridge the models or resample; do not present the correction as reliable. The final target must remain a consistently defined solvated system, not an unexplained mixture of vacuum DFT, classical energy and the old ALPB correction.

Report sampling uncertainty separately from model, electronic-state and protonation sensitivity. For relative selectivity between proteins, the common aqueous leg may cancel algebraically under matching conventions, but protein-specific hydration and state populations still require treatment.

## 6. Identify what failed in the previous approach

Compare controlled model changes on matched chemical states and, where possible, common sampled configurations. Distinguish incorrect local interactions, insufficient sampling, solvent-treatment sensitivity and incorrect mapping to the measured observable.

| Evidence | Interpretation and response |
| --- | --- |
| Reference ligand or water exchange is wrong | Repair local interaction physics before longer protein simulations |
| Independent simulations disagree or overlap fails | Improve sampling or intermediate states before reporting preference |
| Explicit hydration reproducibly changes selectivity | Test which hydration populations or interactions explain the change |
| Conditional site exchange agrees but the assay does not | Investigate occupancy, protonation, folding or assembly mapping |

Compare against the preserved baseline without fitting to a desired sign. Component interactions prevent unique additive attribution from an energy decomposition alone. Use controlled perturbations and independent structural or thermodynamic observations to strengthen causal conclusions.

## 7. Validate and extend the model

Freeze the qualified protocol and predict independent proteins or mutations. Require useful quantitative agreement, robust preference within reported uncertainty and consistency across reasonable starting states. Preserve failures rather than retuning the same validation cases.

For multi-site LanM, calculate the relevant coupled occupancy free energies and their equilibrium populations. Predict the experimental binding response under its actual solution conditions. Include protonation, folding and oligomerization when required; summing isolated-site energies is insufficient. Extend to Hans assembly only when that is the chosen experimental target.

The final delivery should include the experimental mapping, executable protocol and exact model versions, reference evidence, convergence diagnostics, uncertainty, independent predictions, explanation of improvements, remaining failures and measured computational cost. A failed scientific gate is informative, but is not completion of a validated predictor.

## Execution after goal activation

Begin parallel work on target verification, small-system reference chemistry and model implementation. Keep the coupling and reference gates ahead of expensive production sampling. Then distribute work across metals, environments, intermediate states and independent replicas. Measure pilot throughput before setting a larger resource allocation or promising a completion date.

Declare each question, inputs, settings, finite run inventory and interpretation before execution. Keep results immutable, preserve failed attempts and reuse compatible completed calculations. Use `standard_expand` for production and GPU work under verified scheduler limits. The recorded login allowance is 16 aggregate numerical CPU threads across owned processes, with explicit agent shares; record local and scheduler costs separately. Maintain completion notifications while keeping routine scheduler messages out of scientific discussion.

Use the established root-owned submission model with disjoint implementation and independent review tasks. Preserve PQQ production and the existing La/Dy baseline. Do not promote defaults or treat checkpoint changes as physical validation.

### Concrete work sequence

1. **Parallel foundation:** finalize the experimental/state contract; execute the small paired reference panel; reproduce native polarizable-model controls and qualify both metals on the same fixtures. Existing prepared inputs are linked above. Root owns submissions, and independent agents own chemistry, implementation and review.
2. **Choose and freeze the sampling model:** evaluate local exchange energies and forces, hydration and donor competition. Repair demonstrated failures using development data, then pass reserved reference tests. Do not begin long protein production with an unqualified mixed parameter set.
3. **Build both thermodynamic legs:** intact amidated peptide with explicit water/ions, and the corresponding aqueous-ion system. Prepare the terminal-chemistry sensitivity separately. Test numerical derivatives, periodic electrostatics, polarization and the alchemical path.
4. **Run parallel pilot ensembles:** metals, environments, intermediate states and independent replicas provide useful parallel work. Measure throughput, state transitions and overlap; expand the productive calculations rather than waiting on a single diagnostic job.
5. **Complete exchange and any QM corrections:** converge both legs, reconcile assay conditions, and report the precision and model-sensitivity ranges. An unconverged or model-sensitive result remains unresolved even if its sign matches experiment.
6. **Explain and independently test:** compare controlled physical-model changes with the preserved baseline, freeze the protocol, and predict registered independent cases. Only then extend to the occupancy and assembly model needed for natural LanM.

### Definition of done

Deliver an executable, versioned protocol with: a chemically specified target; qualified paired interactions; converged bound and aqueous ensembles; the balanced exchange result and uncertainty; explicit assay mapping; controlled evidence explaining the earlier failures; and independent quantitative validation without tuning to its outcomes. Preserve failures and actual costs. If an essential model or experimental ambiguity cannot be resolved, identify the missing evidence and report the resulting limitation rather than calling the full predictor validated.

## Authorization and earlier diagnostic branches

This document records the strategy, not a launch receipt. Follow the latest user instruction and actual goal state, with authorization history retained in CURRENT.md. This documentation update does not create, pause, resume or complete a goal, and does not authorize new work beyond the user's existing scope.

The [earlier CPCM and reciprocal-Dy plan](../lanm_mechanics_parallel_20261021/PLAN.md) is paused. The CPCM branch has an unreviewed source draft and ten-cell inventory, without a finalized executable configuration. The Dy branch has a recovered checkpoint inventory but no new executable bundle. Neither branch was submitted. They may inform this program if later justified, but are not automatic prerequisites or queued steps.

Existing calculation 16680 was left unchanged under its prior authorization. Its outcome is not a dependency for planning or beginning the new program once activated. Historical receipts and checkpoints remain the authority for its actual status.
