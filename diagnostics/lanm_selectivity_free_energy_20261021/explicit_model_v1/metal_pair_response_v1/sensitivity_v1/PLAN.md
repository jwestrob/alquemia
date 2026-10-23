# Finite consumed D/M response sensitivity

## Question and fixed profiles

Do independent metal–water oxygen direct and mutual cubic-Thole controls change the known collective energy/dipole errors in useful independent directions while retaining the original isolated metal–water force controls? This is a fixed sensitivity experiment, not a fit, candidate selection, reference promotion or validation on held-out data.

Freeze the 17676 six water values, intrinsic alphas, all other parameters, graphs and coordinates. Use the qualified source_v5/prototype binaries and exact same pair selectors, kernels, native settings and numerical admission. Four profiles, applied identically to both allowed metal/O pairs:

| ID | delta_direct | delta_mutual |
| --- | ---: | ---: |
| D_plus | +0.05 | 0 |
| D_minus | −0.05 | 0 |
| M_plus | 0 | +0.05 |
| M_minus | 0 | −0.05 |

The zero baseline is the completed selected 17676 and admitted moment/fragment results. Absent/zero/pure-water numerical identity was independently established by 18383. No new zero-profile molecular calls are needed. Every profile and every declared outcome remains available even if its physical errors increase. No sign-dependent continuation or selected profile is defined.

## Exact finite cells and reuse

Each profile has 30 cells:

- Ten original consumed water cells: water_minus/central/plus × La/Dy, and aquo8/aquo9 central × La/Dy. This retains all ten energy and fifteen rigid-force controls with the unchanged 413-entry objective.
- All eight consumed M+2-water trimers: W8/W9 × La/Dy × near/far, at the original geometry-only selected water indices.
- Twelve distinct literal M+water dimers required by those eight seven-term cycles: W8 water indices 1/5/6 and W9 indices 1/3/5, each metal. These are exact extracted coordinates, not the geometrically different published dimers.

Thus there are 120 profile/cell outcomes. Each requests analyze `system.xyz EM` with the existing reporting-only `polar-print`/`save-uinduce`, plus analytic `testgrad system.xyz Y N`. The D_plus and M_minus near-W9 trimer gradients for both metals already exist in 18383. Reuse those four complete analytic arrays only after exact XYZ/key/parameter/binary equality and raw admission; their original FD provenance is retained without another FD request. There are **120 new analyze calls +116 new analytic testgrad calls =236 new native programs**, four single-thread workers, four CPUs. No new quantum calls or Cartesian/directional FD energies are authorized.

Unchanged ion, literal water, two-water and shell moment/energy terms are reused from the admitted model78/trimer results. This is justified by the exact unmatched-pair fallback and qualified pure-water/ion identity, not by setting missing values to zero. Full historical receipts remain attached. The native final induction residual must be finite, nonnegative and ≤1e−8; complete energy/components/moments/forces and coefficient readbacks are required. The first native or parser failure stops new starts, retaining all 120 slots. Scientific errors are reported and do not stop the finite diagnostic.

Baseline forces are available for the ten original training fixtures and the six W9 qualification fragments (two literal dimers and four trimers). The other fourteen fragment baseline force arrays are explicitly unavailable. They are not inferred from energy/moment-only runs or rerun to fill a reporting column. Reference forces exist for the original ten fixtures; no matched trimer/literal-fragment quantum force dataset is asserted.

## Observables and unchanged references

Recompute each of the eight exact seven-term E3/mu3 cycles using its original manifest and unchanged full-parent-basis PBE0 electronic reference. The permanent multipoles/local frames and atom masks are unchanged; induced sums retain the original six-decimal Debye bound (14e−6/component per trimer cycle). Preserve the lower-precision COM-corrected M check and original reference printing bounds. These bounds exclude physical-reference and SCF uncertainty.

Report the four center energies, analytic forces and C interaction dipole vectors, using the same bare-ion and shell subtraction and existing full-precision quantum C evidence. Report all original 25 development checks separately and retain all four profiles. Do not claim a new full-shell P or B: completing those would require another 22 literal dimers per profile, outside this inventory. The eight sampled E3 values are not a reconstruction of full B.

The isolated-pair force guard is the unchanged minus/central/plus published-dimer comparison against the original consumed PBE0-D4 E/F, with existing 2/5 rigid-force and 1 kcal/mol energy criteria. Keep initial failures and newly lost/passed checks explicit. Trimer E3 uses the previous electronic reference without adding D4; center development uses original targets. No target replacement, reweighting or extra fitting criterion is introduced.

Raw results support later central sensitivity estimates between ±0.05 at fixed kernel roles. A useful response is not proof of the mechanism or a production repair. Quantum E3/mu3 do not uniquely identify polarization, penetration or charge transfer. Separate correlated-reference uncertainty remains recorded.

## Execution boundary

`RUN.py prepare` performs metadata staging and saved raw reuse admission only. Runtime requires the actual passing 18383 independent review, exact executable/source pins, frozen preparation, and root's finite release. Root owns submission. Sources, inputs and immutable result paths are versioned; production, prior fits and held-out outcomes remain untouched.
