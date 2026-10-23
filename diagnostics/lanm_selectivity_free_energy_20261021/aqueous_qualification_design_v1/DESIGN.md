# Next aqueous qualification steps

Source-only design, 21 October 2026. The frozen17605 water candidate passed the
four reserved orientation cells:11/11 checks, with largest energy error0.8244
kcal/mol, rigid RMS1.6406 and body-component error3.6302 kcal/mol/Å. That tests
small water configurations. It does not yet establish a credible bulk hydration
ensemble. The full donor-orientation denominator remains incomplete.

## Smallest useful sequence

| Step | Existing inputs and finite work | What a pass establishes |
| --- | --- | --- |
| 1. Resolve the water mapping | `aqueous_water_mapping_v1/prepared_v1/MANIFEST.json` (`0018cc60`): exact two BIO aquo8 controls and two neutral boxes. Root's two-cell numerical comparison reuses SI results and accounts explicitly for the changed Urey–Bradley term. | Correct namespace transfer and accounted SI→BIO energy/force difference. BIO H–H equilibrium distance1.5537Å stays unchanged; SI used1.5337Å. Atomic forces need not be identical. |
| 2. Test coordination transfer | [Frozen coordination screen](coordination_screen_v1/SPEC.md): two existing W9 geometries, four native programs; reuse four ECP reference results and W8 model anchors. One balanced W8→W9 paired-energy gate and three rigid-force gates. | Consumed development transfer to these unoptimized coordination/geometry changes, against provisional DFT. It is not a ninth-water insertion free energy or a coordination population. |
| 3. Qualify the periodic implementation | [Periodic numerical outline](../aqueous_periodic_numerical_v1/OUTLINE.md):24 programs, exact La/Dy endpoints and intermediate states, CPU/GPU parity, positive polarization residual admission and two directional finite differences. | Numerical behavior for the mapped BIO model and this static periodic path. It does not establish salt or hydration accuracy. |
| 4. Check finite dynamics | After steps1–3 pass, adapt the already declared `paired_legs_v1/REPORT_v2.md` aqueous smoke pilot: both metals, capped minimization, then4000 steps at0.25fs with fully flexible BIO water. Preserve its finite limits; freeze actual integrator/thermostat/readback before release. | Execution stability and measured throughput over1ps, not equilibration or a hydration-accuracy pass. |
| 5. Establish hydration and salt populations | Use independent solvent starts and at least three replicas per metal, with explicit hydration/ion-association diagnostics and a declared finite sampling budget informed by step4. | Only agreement across starts, demonstrated state communication and appropriate external/reference evidence can support an aqueous ensemble for exchange. No production duration is chosen from the current static results. |

The coordination screen deliberately retains SI water. BIO qualification must
retain its different internal water term; small-cell mapping arithmetic does not
prove that this difference has no effect on flexible bulk-water populations.

## Concrete missing physical gates

- **Salt cross interactions:** the mapping lists metal–water O/H, Na and Cl
  pairs explicitly. Water O uses the fitted override; Na/Cl retain native mixing.
  Their coverage prevents omitted or misassigned pairs but supplies no validation
  of contact or solvent-separated ion association. Record metal–Cl and metal–Na
  distance distributions and residence/transition evidence in later trajectories.
  Persistent association must be tested against matched reference/experimental
  evidence or reported as an unresolved sensitivity before aqueous free energy
  is promoted. No new salt reference campaign is released by this design.
- **Hydration:** record the full metal–O radial distribution, first-shell
  population distribution, water identities and exchange events. Freeze the
  coordination counting convention before comparing model outcomes. Do not force
  La/Dy into integer9/8 states or call lack of exchange in a short trajectory
  convergence. Distinct initial shell populations must communicate or be joined
  by a validated enhanced-sampling calculation at fixed global composition.
- **Physical reference and conditions:** the saved PBE0-D4/frozen-f comparison
  is provisional. The W9 screen has no direct native force finite-difference
  qualification of its quantum references. Independent hydration structure,
  electronic-state/method sensitivity and any quantum correction still require
  explicit evidence. The primary La study reports a nine-coordinate interpretation
  in aqueous perchlorate; the later series study examines gradual hydration
  changes. These are leads, not matched acceptance limits for0.1M NaCl at283.15K:
  solution concentrations, exact Dy estimates and uncertainty conventions need
  recovery from the original methods/data before numerical thresholds are set.
  [Näslund et al.,2000](https://doi.org/10.1021/ic991208s),
  [Persson et al.,2008](https://doi.org/10.1002/chem.200701281).
- **Finite-size and assay mapping:** both prepared boxes remain neutral with
  a fixed+3 metal. No charge-changing correction is applied. Residual polarization,
  salt atmosphere, cutoff and box-size sensitivities remain. The aqueous box has
 3980waters,Na7/Cl10,49.324Å edges; seven background salt pairs represent about
 0.0969M at that starting volume. This is not the bound box's total ionic strength
  or a proof of equal salt chemical potential. Neither box includes MES, and the
  2020/2025 assay conventions remain separate.

Reuse the pinned Tinker CPU/double-precision CUDA12 engines, existing native
parsers, `alchemical_path_v1` interpolation equations and the exact mapped BIO
parameter bytes. Do not rerun the passed small-system suites or infer bulk
physics from their success. Any failure identifies the affected numerical,
sampling or physical question; the completed fit and reserved results stay
immutable. Only root binds a finite execution declaration and submits jobs.
