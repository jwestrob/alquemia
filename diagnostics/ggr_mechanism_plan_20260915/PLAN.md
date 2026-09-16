# GGR mechanism investigation: proposed autonomous scope

**Status: PROPOSED — awaiting Jacob's agreement. No new preparations or scientific runs have started.**

Requested by Jacob: “Let's make a detailed plan … when we agree i'll lower the effort meter a bit so you can proceed with autonomy.” Planning used existing artifacts, source-atom inventories and primary documentation. Repository reference: `c94b5f0`; existing working-tree changes belong to earlier work and remain untouched.

## 1. Recommendation and the actual question

Run one staged investigation of **representation, source geometry and local structural sensitivity**, keeping the electronic Hamiltonian fixed. The intended outcome is an explanation of which assumptions matter, plus a defensible next discriminator experiment. Obtaining a favorable GGR sign is not the success criterion.

There are two distinct observations:

1. Original GGR v2 passed its frozen sign test at S=−1.183146. Correcting a formaldehyde-like backbone fragment to an amide moved S to +3.057473, a +4.240618 kcal/mol shift. The repaired protocol has no calibrated zero; this is a sensitivity finding, not an independently calibrated failed classification.
2. Under the same repaired protocol, both alpha-lactalbumin structures rank below Ca-favoring GGR despite condition-qualified La-favoring evidence. This ordering conflict survives any common reference shift and is the stronger generalization challenge.

Use GGR to investigate the first issue while checking whether changes help, harm or leave the second issue unchanged. All these cases have already been inspected and are development cases. None becomes a blind test.

The plan contains **38 named high-level evaluations**, plus a specific half-step consistency branch of up to 16 evaluations. These counts specify the experiment; they are **not CPU-time, wall-time or spending limits**. Technical retries are recorded separately. No production rescore or default change is included.

## 2. What remains fixed

- Native ORCA 6.1.1 r2SCAN-3c, native basis/ECP, CPCM(Water), DefGrid3. No extra composite correction or alternative functional.
- Existing protonation states and exact protonated source coordinates for existing structures. La/Ca pairs use identical coordinates, atom ordering, waters and ligand composition; charges differ only by the metal's +1 charge difference.
- No heavy-atom relaxation in stages A/B. No water addition/removal, microstate search, donor cutoff tuning or experimental relabeling.
- Existing donor qualification at 3.1 Å and inclusion at 3.3 Å, including the established residue-consistent alternate-conformer policy.
- Existing aquo reference remains untouched. Report raw R=E_Ca−E_La and differences in R; S is an explicitly uncalibrated reporting gauge for new models. Convert Hartree once with the released 627.509474 factor.
- Current baseline/default and PQQ bands remain unchanged. New generic models inherit no PQQ threshold.

The stronger evidence is an effect that persists under physically comparable preparations. Comparing absolute energies of different atom inventories is not meaningful; comparing their Ca-minus-La contrasts is.

## 3. Stage A — Does the nearby fragment boundary drive the result?

### A1. Extend the amide boundary consistently: eight endpoints

Replace each repaired formamide-like backbone unit with an **N-methylacetamide-like unit** using actual source atoms. This retains the two Cα atoms on either side of the peptide bond. Preserve the source carbonyl C/O, amide N/H and Cα hydrogens. Cap the four outward Cα bonds to omitted N/C/Cβ atoms with the existing 1.09 Å C–H convention. Bonded neighbors must come from connectivity, not residue-number arithmetic.

| Existing preparation | Amide donor unit(s) | Atoms, current → proposed | New evaluations |
|---|---|---:|---:|
| GGR 1GLG | Gln140–Ile141 | 52 → 58 | Ca + La |
| Aequorin 1SL8 EF3 | Ala125–Ile126 | 43 → 49 | Ca + La |
| Alpha-lactalbumin 1F6S | Lys79–Phe80; Asp84–Leu85 | 40 → 52 | Ca + La |
| Alpha-lactalbumin 6IP9 | Same two peptide units | 43 → 55 | Ca + La |

Each retains ligand charge −3, La core charge 0 and Ca core charge −1, singlets. Alpha's two source states retain their respective two/three waters. They remain one biological observation.

**Why these controls:** aequorin EF3 had almost exactly GGR's repair shift, so it tests whether the boundary effect is shared. It supplies no site-specific affinity truth. Both alpha structures test the existing relative-ordering conflict without choosing a favorable geometry. No new biological panel is needed.

Implement this as a reusable, versioned selection policy on the existing source graph. Overlaps must merge before capping; preserve source-to-QM maps, retained/cut bonds, charges and electron parity. Unsupported chemistry must remain explicit. Proposed ID: `generic_peptide_alpha_caps_native_r2scan3c_dev_v1`.

### A2. Connect GGR's native local peptide: two endpoints

Prepare one additional GGR representation, estimated **111 atoms** from the existing source inventory:

- Keep all source atoms in Asp138–Gly139–Gln140–Ile141–Gln142.
- Retain Lys137 Cα/C/O and its Cα H as the left acetyl-like terminus; cap the omitted N and Cβ directions. Do not import Lys137's charged side chain.
- Retain Phe143 N/Cα, N–H/Cα–H as the right methylamide-like terminus; cap omitted C and Cβ directions.
- Merge existing Asp138 and Gln142 fragments into that peptide; preserve the existing Asp134, Asn136 and Glu205 fragments.
- Charge remains −3: added complete residues are neutral, and Asp138 was already present and charged.

This is a topology-defined local segment between the nearest selected side-chain donors bracketing Gln140's carbonyl, with chemically complete terminal peptide bonds. Verify the actual source connectivity and atom count before scoring. Proposed ID: `ggr_connected_segment_native_r2scan3c_dev_v1`.

Do **not** automatically apply this expansion to the other controls: the analogous aequorin segment imports Asp129, while alpha expansion imports excluded charged side chains. Those would change the scientific comparison. Do not delete or neutralize those groups to force matching.

### A readouts and limits

Report the GGR ladder: historical defective v2; repaired formamide; extended amide; connected peptide. Only the last two require new endpoints. Report each change in R, all four A1 comparisons, and alpha-minus-GGR ordering under the common extended-amide policy. Reuse existing component extraction, with the existing warning that SCF minus CPCM is not a vacuum energy.

These models add native electronic interactions and alter the continuum cavity. **This is local representation sensitivity, not a matched-system QM/MM partition test and not an isolated measurement of protein electrostatics.** The frozen connected peptide includes no relaxation free energy. Two small successive changes would support stability over the tested representations, not prove convergence.

## 4. Stage B — Does GGR's source conformation matter? Eight endpoints

Preselect two structures from one primary study:

| Structure | State | Resolution |
|---|---|---:|
| [2FW0](https://www.rcsb.org/structure/2FW0) | E. coli WT, sugar-free/open, Ca retained | 1.55 Å |
| [2FVY](https://www.rcsb.org/structure/2FVY) | E. coli WT, glucose-bound/closed, Ca retained | 0.92 Å |

The structures share a WT construct, but crystallization conditions differ. 2FW0 contains sodium/citrate in the sugar cleft; 2FVY has reported radiation damage at Glu149. Thus this tests **structural robustness**, not a causal sugar-only effect. Both structures lack terminal residues 1 and 307–309. [Primary study](https://kiesslinglab.com/sites/default/files/labs/kiessling/pdfs/2007BorrokProtSci.pdf).

After approval, acquire/pin the source CIFs and inspect the actual Ca-site mapping, alternate conformers, missing atoms, covalent modifications, waters and nearby nonstandard species before energies. Use the sequence-matched chain A monomer and the existing altloc policy; do not select another chain or structure based on a score. Missing terminal residues outside the representation are recorded, not reconstructed.

Build **both current repaired v3 and the A1 extended-amide representation** for each source: two structures × two representations × two metals = eight evaluations. Use one frozen protonation preparation per source at the established pH 7; reuse it across representations. Map residues by sequence/connectivity before selection, rather than assuming 1GLG's author numbering. The frozen selector and qualification gates must recover the same six donor residues and zero-water GGR composition; do not manually force that membership. Different denticity among atoms of the retained residues can be recorded as geometry variation if the existing qualification gates still pass. An added/lost donor residue, water, missing local atom or incompatible microstate makes the matched comparison unsupported. Continue other independent work.

Compare these with 1GLG under both representations. Report the range of R across the three sources and the representation-by-source differences. All remain **one GGR observation**. The sugar is outside these cores, so this measures geometry-associated changes, not its direct electrostatic contribution or an allosteric free energy.

## 5. Stage C — What real motions change the electronic preference? Twenty endpoints

### C1. Analytic gradients: four center evaluations

On GGR 1GLG's **extended-amide and connected-peptide models**, obtain Ca and La analytic gradients at the frozen centers. Use the same native Hamiltonian/CPCM/grid, explicitly adding **EnGrad and TightSCF** for derivative accuracy. Preserve the four stage-A normal-SCF center energies separately; this is a declared numerical recipe, not a silent baseline update.

Proposed sensitivity protocol: `ggr_local_sensitivity_tightscf_dev_v1`, with the underlying representation recorded separately. Cache identity must include source geometry, representation, numerical recipe, displacement coordinate/amplitude and task type. A single-point receipt cannot stand in for a gradient receipt.

Record the TightSCF-to-normal-SCF change in each endpoint and in R. Require each endpoint energy change and the change in R to be ≤0.05 kcal/mol before treating the normal/tight results as interchangeable at the interpreted scale. If this bridge fails, report the numerical discrepancy; do not pool them or reinterpret it as a boundary effect. C can still report internally consistent tight-SCF results if its own checks pass.

Existing `affordable_response.py` already supports analytic extraction, link-atom chain rules and paired gradients. The exact integration is unvalidated. Required implementation checks include real `.engrad` row/scalar formatting, units, atom order, coordinate identity, La ECP conventions, explicit charges/multiplicities, normal termination and complete native composite/solvent derivatives. The current La endpoint uses a 46-electron Def2-ECP. Do not infer coordinate-file atomic numbers from that electron count.

ORCA documents analytic DFT gradients and relevant CPCM support, but this does not substitute for an executed check of our combination. [Gradient documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html), [solvation documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html).

Do not invoke the larger `site_mechanics.py` scan workflow. Adapt the existing runner and gradient extractor only for the named tasks. Any numerical-gradient fallback makes the gradient task unsupported; it does not authorize a numerical DFT gradient campaign.

### C2. Validate two physical directions: sixteen single points

Use both metals and both representations, at each positive/negative displacement. **All C2/C3 single points use TightSCF**, matching the C1 centers and the same native Hamiltonian, solvent and grid:

1. **Metal motion:** ±0.02 Å along the frozen metal-to-Gln140-carbonyl-O direction. All protein atoms remain fixed.
2. **Peptide crankshaft:** ±1° rotation of the real Gln140 C/O–Ile141 N/H peptide unit about the line connecting its two fixed source Cα atoms. This preserves the peptide unit and its anchor bond lengths while changing local angles. Every other source atom remains fixed.

Directions and amplitudes are fixed from geometry before gradients/scores are inspected. Generate both endpoint geometries from the same displaced physical source. Recompute link caps from their source bond endpoints, never move caps as independent coordinates. Require maximum displaced source-heavy-atom motion ≤0.05 Å and unchanged coordination membership, protonation and water inventory. A violated condition makes that direction unsupported, not an invitation to pick a favorable alternative.

Project each analytic gradient through the physical coordinate Jacobian, including cap chain rules where applicable, and report gradient rather than force. Units are kcal/mol/Å for metal motion and kcal/mol/radian for the peptide rotation. Export the paired derivative of R. These two particular motions leave all link-cap anchors fixed in the extended/connected models, so they validate the selected physical derivatives, **not moving-cap force transfer**. Retain the existing geometric Jacobian check on a real GGR bond as a separate software test; do not claim an end-to-end displaced-cap energy validation.

For coordinate q and amplitude h, compare the odd energy change

`[E(+h) − E(−h)] / 2`

with `h × dE/dq` from the analytic gradient. Check each metal and the paired R. At the nominal amplitude, the proposed tolerance is **max(0.02 kcal/mol, 5% of the magnitude of the predicted odd energy change)**. This tests energy/gradient consistency at a scale well below the observed several-kcal/mol representation shifts; it is not a biological classification threshold or an uncertainty estimate. Retain the even energy change separately as curvature-sensitive response: a perfectly quadratic energy surface already produces it. Do not equate it with anharmonic breakdown or fit a relaxation correction.

### C3. Explicit conditional half-step check

If that consistency test fails for a validly prepared representation/direction, repeat **both metals at ±half amplitude** for that block: ±0.01 Å or ±0.5°. Four evaluations per triggered block, at most 16 across the two representations/two directions. At half amplitude, halve the absolute energy-tolerance floor to **0.01 kcal/mol**, retaining the 5% relative term; this preserves the absolute derivative-consistency requirement. Report the residual and residual/h at both amplitudes. This branch probes step dependence; a half-step pass alone does not identify the error's cause. It cannot be triggered by an inconvenient score. Do not keep shrinking steps or changing grids until the result agrees. Persistent failure stays a failed numerical check.

### C interpretation boundary

This is sensitivity of the isolated CPCM electronic model. Even the connected fragment omits much of the protein's mechanical resistance. Large gradients are not relaxation energies; small gradients are not proof of correct affinity. No Hessian, fitted spring, free optimization, entropy term, thermal covariance or environmental correction is included. `response_model_not_validated` remains true and numerical relaxation corrections remain unavailable.

## 6. How we will distinguish outcomes

| Observation | Supported conclusion | What it does not establish |
|---|---|---|
| Large formamide→extended shift across controls | Nearby donor boundary materially changes the descriptor | Which representation predicts affinity correctly |
| Small extended→connected change for GGR | Limited sensitivity over those tested local models | Whole-protein convergence or global-field irrelevance |
| Large extended→connected change | More local chemistry/cavity context matters | A uniquely electrostatic or mechanical cause |
| Large spread across the three source structures | A single crystal geometry is a fragile input for GGR | A causal sugar effect or three independent failures |
| Validated large paired derivative for physical motion | R is locally sensitive to that physical coordinate | An affordable, accurate relaxation free energy |
| Similar paired derivatives across representations | Sensitivity is less likely to be solely a nearby cap artifact | Adequate scaffold curvature or thermodynamics |
| Alpha–GGR ordering remains reversed under A1 | The local amide extension alone does not resolve the known challenge | Failure of every possible environmental model |
| All tested sensitivities are small but disagreement remains | Prioritize evidence-state alignment and omitted affinity physics | Permission to tune a reference or relabel GGR |

Report continuous effects and all unsupported cases. Do not select a winning structure, threshold or representation because it gives the desired sign. Proposed scientific triage scale: differences around 1 kcal/mol already matter for GGR's narrow original margin; shifts of several kcal/mol rival the known repair effect. These are interpretation scales, not fitted decision bands or pass/fail accuracy claims.

## 7. Execution, accounting and autonomous boundaries after agreement

| Task group | Purpose | Named evaluations |
|---|---|---:|
| A1 | Four extended-amide Ca/La pairs | 8 |
| A2 | One connected-GGR pair | 2 |
| B | Two source structures × two representations × two metals | 8 |
| C1 | Two representations × two analytic-gradient centers | 4 |
| C2 | Two representations × two directions × two signs × two metals | 16 |
| **Main sequence** | | **38** |
| C3, conditional | Four per failed representation/direction check | **0–16** |

A precedes B/C preparation; B and C can then proceed independently. Complete the named comparisons regardless of whether earlier scores look favorable, subject to the stated chemistry/numerical validity checks. Existing matching outputs can replace a task only with full provenance and scientific-setting agreement. Technical recovery may fix parsing, quoting, memory sizing, MPI layout, missing receipts or interrupted execution without changing the scientific task. Preserve all attempts. A different Hamiltonian, state, structure, label or additional experiment requires a revised agreed scope.

Reuse `SourceGraph`, the established manifests/cache checks, `run_orca_task_manifest.py`, the baseline collector and `affordable_response.py`. Candidate inputs, coordinates, outputs and receipts belong under `workspaces/ggr_mechanism_20260915/`; compact plans/results under this diagnostic directory. Implementation snapshots accompany prepared manifests. No watcher, queue priority, shared environment or other user's job is changed.

Use the established small-job allocation shape, initially four concurrent 16-rank workers on 64 CPUs when appropriate. Account for the larger connected model's memory and distribute actual endpoint work across workers; do not reserve a large node for a serial task. No GPU work is needed.

Measured reference cost: the 52-atom GGR repaired pair consumed **264.288805 summed endpoint wall-seconds and 4,228.620880 assigned rank-seconds**, using 16 ranks per endpoint. These are endpoint sums, not elapsed batch time. The 111-atom model and analytic gradients have no measured matched runtime yet. Their first receipts will establish the cost projection for the rest. No CPU-time, wall-time or spending stop is imposed. Research diagnostic cost and proposed routine per-site cost must be reported separately; no automatic production promotion follows.

After approval, autonomy includes implementation, tests against real artifacts, source acquisition, the exact preparations, submission of the named tasks/conditional branch, technical recovery, collection, comparison, vault notes and scoped commits. Do not pause for individual commands. If one scientific check is unsupported, retain that status and finish independent approved work.

## 8. Deliverables and stopping point

1. Versioned experimental preparation policies and graph/coordinate/charge regression checks; unaffected PQQ fixtures stay identical.
2. All actual unrounded energies, costs, failed attempts and physical mappings in the benchmark ledger, marked development and grouped by protein.
3. One compact report with a GGR representation ladder, representation-by-source table, matched alpha–GGR differences, and physical-gradient/energy-check plots. No pooled accuracy or calibration fit.
4. A validated analytic-gradient integration at the level actually achieved, retaining disabled relaxation/entropy corrections.
5. A handoff with exact resumable commands, completed/pending/unsupported status, preserved implementation hashes and one recommended next scientific step.

The investigation is complete when these defined comparisons and checks are reported, including honest failures. It need not produce a replacement scorer.

## 9. Follow-on ideas: documented, outside this execution scope

- **GGR gateway mutants:** published Q142 substitutions are promising same-protein contrasts. Exact same-assay Ca/La tables and matching structures were not verified, so no specific mutant is ready to score. A neutral substitution would be particularly useful for avoiding gross-charge confounding. Obtain primary evidence before selecting a test. [Primary publication](https://doi.org/10.1021/bi952430l).
- **Hydration/protonation:** a justified matched state/reference cycle could test these, but adding a water or proton to rescue GGR would not diagnose them. The current plan preserves both states. The source-sugar mismatch remains a qualification, not a changed label.
- **Mechanical correction:** if stage C shows reproducible physical sensitivity, the next question is an independently defensible scaffold stiffness model. Neither an arbitrary spring nor the previously tried 12-6-4 curvature is automatically acceptable.
- **Environmental/global representation:** the failed APBS and whole-protein xTB routes remain archived. Local representation/gradient evidence can inform a later alternative, but this plan does not repeat those experiments or infer that they must succeed.
- **Electronic-method sensitivity and cheap comparators:** inspect the existing historical comparisons before proposing more calculations. They remain useful for the broader paper benchmark, but are not prerequisites for this specific investigation.

Approval of this document would cover sections 2–8, including the explicit conditional branch, and not these follow-on experiments. Record Jacob's actual agreeing message before launching scientific work; the existence of this agent-written plan is not approval.
