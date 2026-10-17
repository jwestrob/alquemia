# Mex full-coordinate energy/force qualification

17 October 2026 UTC. **The existing Mex composite passes the declared checks for newly movable atoms.** This permits testing broader relaxation under the same energy. It establishes numerical coupling consistency, not accurate La/Dy preference or scaffold stiffness.

## Actual result

Slurm12952 completed all34 composite evaluations: two fresh central replays and32 central-difference displacements. Both metals pass all four directions at both step sizes; every individual energy component and their sum pass the unchanged tolerances in SPEC.md. Fresh central component energies reproduce exactly. Largest componentwise difference between fresh and cached total gradients is1.098e−12kcal/mol/Å.

| Direction | La total residual at0.0001Å | Dy total residual at0.0001Å |
| --- | ---: | ---: |
| Backbone carbon | +1.49569e−5 | +8.07273e−6 |
| Boundary anchor | +2.11832e−6 | −1.39106e−7 |
| External water | −1.94245e−6 | +2.39616e−5 |
| Nd spectator | +2.43038e−6 | −3.17887e−6 |

Residuals are finite difference minus analytic projection, in kcal/mol/Å. Both0.0002Å and0.0001Å steps pass separately; this table does not select a favorable step. Directions and physical identities were frozen before any new evaluation. Complete component residuals and explicit ten-check denominator are in the immutable RESULT below.

The full2,049-atom coordinate update propagates through the209-atom electronic region,1,840 external charges, two-anchor link atoms, native MM ledger and full fixed-charge ALPB. The same complete energy supplies the gradient. Environment charges remain fixed in magnitude. Physical Dy multiplicity6 and spin-independent GFN effective multiplicity1/UHF0 remain distinct recorded assumptions.

## Cost and reproducibility

- One warm MACE load;34 MACE calls,68 single-thread GFN calls,34 MM configurations/136 component queries,34 ALPB calls. No DFT or optimization.
- Application elapsed305.998s; scheduler307s on12 allocated CPUs and oneGPU:3,684 allocated CPU-seconds and307 allocated GPU-seconds. These are allocations, not measured utilization.
- Native GFN child CPU time311.784s; scheduler batch MaxRSS2,294,556KiB. GPU memory was not separately measured.
- Zero-call interface preflight12951 completed9s/1CPU; search guard preflight12955 completed7s/1CPU with all protein/water/target/spectator identities admitted.
- Frozen evaluator SHA256 `5153440c685dc2c1c335c35bdae42b19ce93f15de38d4706bd5b616365a32a8c`; search driver `35d12483cc72105b1562296892416b9e8d80ebf2f93fb37e05cef82ddc2a7104`.

Authoritative artifacts relative to repository root:

- `workspaces/lanm_coupled_mechanics_20261017/mex_preflight_v1/CONFIG.json`
- `workspaces/lanm_coupled_mechanics_20261017/mex_qualification_v1/execution/RESULT.json`
- `workspaces/lanm_coupled_mechanics_20261017/mex_qualification_v1/SUBMISSION.json`
- `workspaces/lanm_coupled_mechanics_20261017/mex_search_preflight_v1/execution/RESULT.json`

SUBMISSION records the exact command, allocation and completion watcher. Do not repeat these completed evaluations. The main search stage still follows the plan's LBT scout gate. Independent DFT response and scaffold-stiffness checks remain outstanding.
