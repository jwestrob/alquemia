# Nonlinear donor accommodation: substantial response, incomplete qualification

The two compressed PLM sites have large, metal-dependent accommodation energies
under the declared model. Terminal carboxylate rotations relieve their very short
metal–oxygen contacts and shift Ca−La contrasts by **+60.93 and +37.54 kcal/mol**.
The two known PQQ controls change by less than0.56kcal/mol and retain their expected
regions in the explicitly developmental transfer of the old bands.

All eight optimizers produce candidates; seven pass the separate interior/raw
gradient/geometry gate, but only **three individual endpoints and zero complete
metal pairs** pass every curvature check. Relaxed classifier accuracy and binding
free energies remain unestablished. The baseline and released fast mode are unchanged.

## Fixed experiment and actual outcome

The predeclared [plan](PLAN.md) uses the original four source contexts, one start
per metal, anchor-Glu chi3 plus actual extra-Asp chi2 where present, bounds±0.8rad.
All other physical coordinates, atom identities, charges, protonation, explicit
waters and caps stay fixed. There are no label-selected starts or rescues.

The evaluated energy is native vacuum OMOL plus native GFN2 ALPB(water)−vacuum,
using Tight analytic gradients and the current physical Jacobian. The optimizer
receives component-wise energies relative to its own q0; eV and hartree each convert
to kcal/mol once. No entropy, fitted spring, density-functional gradient or aquo
reference is added. The exact settings and hashes are in the manifest and
[implementation note](IMPLEMENTATION.md).

| Case | Ca work | La work | Change in R=Ca−La | Old-band transfer, q0→candidate |
|---|---:|---:|---:|---|
| 1H4I, Ca-associated control | −0.13104 | −0.05016 | −0.08088 | Ca→Ca |
| 4MAE, La-associated control | −0.62350 | −0.06522 | −0.55828 | La→La |
| PLM8344, unknown label | −35.76648 | −96.69794 | +60.93146 | Ca→Ca |
| PLM07ab, unknown label | −18.02357 | −55.56410 | +37.54053 | Ca→indeterminate |

Energies are kcal/mol on this model's scale. Full identifiers, unrounded energies,
native/solvent components, raw contrasts, coordinate pins and donor distances are
in [RESULT.json](RESULT.json). Changed geometries and Tight calculations do not
inherit the released calibration. The table is a developmental transfer check;
the PLM sites have no experimental labels. Both controls are consumed development
cases. Their raw gap contracts64.95440→64.47701kcal/mol. The full25-reference panel
was not rerun by this pilot, and no accuracy improvement is established here.

For PLM8344, the compressed extra-Asp contact changes1.766→2.250Å for Ca and
→2.423Å for La. For PLM07ab, it changes1.847→2.234/2.400Å. The metal is fixed;
these are complete terminal-group rotations. Maximum candidate heavy-atom
displacement is0.76653Å. All eight candidate geometries pass the declared physical
checks, and none lies on the angular boundary.

## Qualification and component audit

| Endpoint | Raw stationarity | Curvature | Final status |
|---|---|---|---|
| 1H4I Ca | pass | pass | qualified constrained minimum |
| 1H4I La | pass | refinement fails | stationary candidate |
| 4MAE Ca | pass | pass | qualified constrained minimum |
| 4MAE La | pass | ALPB probe unavailable | stationary candidate |
| PLM8344 Ca | pass | symmetry/refinement fail | stationary candidate |
| PLM8344 La | fails | not evaluated | nonstationary candidate |
| PLM07ab Ca | pass | pass | qualified constrained minimum |
| PLM07ab La | pass | refinement fails | stationary candidate |

PLM8344 La satisfies scipy's projected stopping test but its absolute derivative
is0.41159kcal/mol/rad, above the independent0.2 gate. That gate prevents the
optimizer's success flag from being mistaken for physical stationarity.

The three available failed Hessian checks are dominated by the solvent derivative
component. Native versus solvent refinement errors are0.00724/0.31081 for1H4I La,
0.00335/2.64649 for8344 Ca, and0.00260/2.45058 for07ab La, in
kcal/mol/rad². The combined errors exceed their frozen tolerances; matrices are
not repaired or eigenvalues clipped. This audit localizes the observed discrepancy,
but does not distinguish convergence precision, higher-order behavior or an
approximation in the solvent derivative. The4MAE La ALPB curvature probe actually
fails SCF convergence after125cycles; its matching vacuum calculation succeeds.
No failed point was retried. Candidate energies themselves are all available.

Primary-to-Tight q0 changes are separate from accommodation. Native MACE repeats
exactly. The largest solvent endpoint shift is0.08078068kcal/mol, and the largest
paired contrast shift is0.08112124kcal/mol, both for07ab. These numerical-policy
differences do not explain the tens-of-kcal accommodation effects.

The separate native-candidate execution gate requires at least one completely
qualified pair. **It is not met; no native DFT candidate calculations were launched
under this pilot.** A genuine prepare/dry-run adapter can preserve the candidates
and failed statuses without overriding that decision. No harmonic entropy or
relaxation free-energy scalar is available.

## Execution, checks and practical cost

Job1204162 completed all eight starts. It made89 native MACE energy/force calls
and178 native GFN2 analytic endpoint attempts;177 GFN2 endpoints completed, with
the single failed curvature probe retained. There were88 complete composite
evaluations and one failed evaluation. No numerical DFT gradients or new DFT calls
ran. Seven focused real-fixture tests pass, including actual component algebra,
mapping, the observed curvature failure and projected-versus-raw stopping behavior
([TESTS_v6.txt](TESTS_v6.txt)). Parser/algebra replay and executed scientific calls
are separately counted. Original failed test-harness logs remain preserved.

Measured operation time is1723.192s. The actual allocation is**1725s ×32CPUs =
55,200core-seconds**, plus**1725GPU-seconds** on one H200. All failed attempts are
included. MACE requests total14.272s after a1.201s model load; most elapsed time is
the repeated analytic solvent operations and their execution overhead. GPU tensor
peak is6,055,856,640bytes; worker host RSS is1,848,868KiB. A whole-job aggregate
host-memory peak is unavailable. These are prepared-context research timings,
excluding folding and prior preparation. Nested ORCA allocation records are not
added again. Exact receipts are pinned in [COST.json](COST.json).

The full composite optimizer has not earned routine scanner use. The result
supports testing cheaper physical geometry proposals while measuring the complete
scoring energy separately, as the independent proposal experiment is doing.
Keep the released classifier in place until a common preparation rule demonstrates
reference accuracy and structural robustness at useful cost. No new scientific
run, wider bound or threshold change is proposed by this report.

See [COMMANDS.md](COMMANDS.md) for read-only validation/collection. Final full
collection: `workspaces/accommodation_nonlinear_20260920/collection_final_v1.json`.
The immutable automatic collection and every attempted calculation remain under
`pilot_v2/`. Live reporting adds donor distances and native/solvent curvature
components without changing the executed optimizer or its results.
