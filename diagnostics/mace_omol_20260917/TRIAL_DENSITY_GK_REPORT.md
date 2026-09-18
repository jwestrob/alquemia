# Responsive-density / AMOEBA-GK / MACE development result

**Complete, 2026-09-18 UTC. All predeclared development gates pass.** The
responsive-density candidate orders both alpha-lactalbumin structures above
both GGR representations and reduces the GGR partition discrepancy to
0.541233 kcal/mol. All 32 numerical checks pass. These are **two consumed
biological groups**, not four independent or blind observations. The production
baseline remains unchanged; the active discriminator goal is not complete.

Protocol: `saved_responsive_trial_density_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1`.
The [plan](TRIAL_DENSITY_GK_PLAN.md) was recorded before these calculations.
No labels, geometry, waters, thresholds, cavity parameters or MACE readouts
were changed after inspecting the new results.

## Matched results

Here R = A_Ca − A_La; larger R is more La-like. Differences between sites
cancel the common metal offset. No compatible aquo reference or absolute
classification bands have been established for this protocol.

| Check, kcal/mol | Frozen density | Responsive trial density | Unchanged criterion |
|---|---:|---:|---|
| Alpha 1F6S − GGR extended | −0.437800 | 10.127062 | >0.02 |
| Alpha 1F6S − GGR connected | 3.345908 | 10.668296 | >0.02 |
| Alpha 6IP9 − GGR extended | −3.546267 | 9.252283 | >0.02 |
| Alpha 6IP9 − GGR connected | 0.237440 | 9.793516 | >0.02 |
| GGR connected − extended | −3.783707 | −0.541233 | absolute value ≤2 |

GGR carries a direct same-assay Ca-favoring direction. Alpha carries a qualified
cross-study strong-site La-favoring direction. Their ordering is a development
challenge, not a quantitative common-condition affinity comparison. Both alpha
structures and both GGR partitions must remain in the denominator.

## What changed, and what it teaches

Every density-dependent term was recomputed together using eight saved
protein-field-responsive quantum densities. Exact matched comparisons verify
unchanged physical coordinates, keys, masks, MACE terms and evidence records.

| Change in R from frozen density, kcal/mol | Intrinsic core | Direct field | GK transfer | Induction | Total |
|---|---:|---:|---:|---:|---:|
| Alpha 1F6S | 0.031510 | 0.077382 | 4.228387 | −2.262909 | 2.074370 |
| Alpha 6IP9 | 0.041909 | 0.114625 | 6.928529 | −2.777005 | 4.308058 |
| GGR connected | 0.094416 | −0.450311 | −5.747166 | 0.855044 | −5.248018 |
| GGR extended | −0.139595 | −0.070867 | −9.657419 | 1.377390 | −8.490492 |

Small changes in the *differential intrinsic core energy* did not imply that
density response was negligible: the consistently updated reaction-field and
induction terms materially changed the ordering and partition discrepancy.
This is useful evidence about this model's sensitivity, not a unique biological
explanation or proof of an exact solvation treatment.

The unchanged MACE short component contributes +5.115575 kcal/mol to the GGR
partition contrast. Algebraically removing it leaves −5.656808 kcal/mol,
outside the partition criterion. This supports its usefulness for representation
robustness in this development example. It does not establish that MACE alone
caused the improved biological ordering. This omission calculation is an
algebraic audit of inspected outputs, not a separately validated alternative.

## Energy accounting and remaining approximation

```
E_core[rho*] = E_DFT,old_embedded[rho*] − sum(q_old * phi_rho*(old_positions))
A_M[rho*] = E_core,M[rho*] + V_current_AMOEBA[rho*]
           + es_static(Qproxy,E;B) − es_static(0,E;B)
           + I_current_density+proxy_RF − I_environment
           + T_short(full,M) − T_short(core,M)
I = −0.5*(electric/dielec)*sum(mu_solv,d * F_solv,p)
```

The original embedded calculation used a permanent ff19SB field and
`DoEQ false`; subtracting its actual density interaction once removes that
generating field without adding external-charge self energy. The resulting
trial density is **not self-consistent with the current AMOEBA/GK functional**.
Every endpoint uses its preassigned trial density; no favorable density was
selected per protein. Intrinsic polarization costs are positive in all eight
same-geometry comparisons with the archived vacuum calculations.

Actual quantum potential, electric field and spatial potential Hessian supply
direct multipole coupling and induction. Projected endpoint CHELPG charges
supply the approximate GK reaction field. The native AMOEBA2018 environment,
common Ca2018 metal cavity, source exclusions and full physical boundary are
unchanged. This common cavity is not a validated La force-field parameter.
Environment-only and nonpolar terms cancel within that common boundary.
MACE's jointly trained short readout is not uniquely nonelectrostatic.

No CPCM energy, second complete solvation term, geometry relaxation, entropy,
or unsupported combined gradient is added. No absolute affinity is reported.

## Executed physical and numerical checks

- All 32 full-model checks pass. Maximum change on response refinement from
  1e−7 to 1e−9 Debye: 2.419978e−7 kcal/mol; maximum rigid component difference:
  2.046363e−12 kcal/mol. Both full vacuum identity corrections are exactly zero.
- The prescribed ±5% common-metal-radius perturbation changes any tested
  relative contrast by at most 0.000331971 kcal/mol. This tests that specific
  radius perturbation, not all dielectric or cavity uncertainty.
- All source induced dipoles remain zero; maximum environment atomic induced
  dipole is 0.336748 eÅ, below the declared 1 eÅ flag. Actual unrounded solver
  residuals pass. No omitted or zero-filled failed contribution exists.
- Eight real saved-density queries pass their independent refinements: maximum
  electric-field change 1.434678e−8 au; quadrupole contraction change
  0.000227741 kcal/mol. These are spatial potential derivatives, not nuclear
  DFT Hessians. Sixteen real boundary initializations pass all 12 paired checks.

## Execution, tests and measured cost

Jobs 1201135 and 1201137 completed without scientific failure or retry.

| Job | Work | Wall s | CPUs | Allocated core-s | Actual CPU-s |
|---|---|---:|---:|---:|---:|
| 1201135 | Eight saved-density utilities | 642 | 8 | 5,136 | 2,215.715588 |
| 1201137 | Full hybrid plus collection | 580 | 64 | 37,120 | 2,630.625999 |
| Total | Incremental development | 1,222 | — | 42,256 | 4,846.341587 |

Zero new DFT, MACE or charge-fit calls; zero GPU time. The 111 logical native
tasks contain 37 qualified environment-only reuses. New executions comprise
34 static energies, 32 field queries and 40 response solves. Native kernels
sum 82.792887 wall seconds; complete child processes sum 140.591004 wall and
2,216.401163 CPU seconds. Sampled batch peak RSS: 1,598,144 and 3,183,280 KiB.
Reuse checks actual inputs, source implementation, field bytes and receipts;
reused costs are excluded from these new-call totals.

Local field preparation: 13.711318 wall /12.687261 CPU seconds. Boundary
preparation: 65.835948 /47.364624; its native initializations: 16.733242 /
16.553698. Hybrid preparation inside job 1201137: 48.770911 /40.601441,
already included in that job's cost. Prior source audit: 15.358869 /14.253791.
Reused quantum, CHELPG and MACE costs are additional. Hardware differs between
the two field pilots; no matched speedup is claimed. **Fresh end-to-end
production-pair cost is still unmeasured.**

Twelve distinct new parser/source/cache/algebra/actual-integration tests and
13 regression tests pass across the recorded focused runs. Final integration
tests use actual scientific outputs; none remain skipped. An initial expected
reuse count of 29 was corrected to the actual 37 before new native energies;
the scientific inventory and criteria did not change. Corrupted-input tests
are explicitly corrupted copies of real fixtures, not scientific evidence.

## Decision and next work

Numerical credibility: **passes the executed controls**. Scientific usefulness:
**promising on the two development groups; wider validation required**.
Affordability: **native component practical; complete fresh cost unmeasured**.

Freeze this candidate and test wider real inputs without retuning these cases.
Retain production baseline. No new absolute threshold may be fitted on the
validation set; functional PQQ labels, affinity directions and protein-level
multisite evidence must remain separate. Unsupported preparation is explicit.

Unrounded components, provenance, receipts and tests:
[TRIAL_DENSITY_GK_RESULT.json](TRIAL_DENSITY_GK_RESULT.json).
Runnable audit/collection/comparison commands:
[TRIAL_DENSITY_GK_COMMANDS.md](TRIAL_DENSITY_GK_COMMANDS.md).
Immutable products: `workspaces/mace_omol_20260917/trial_density_gk_*`.
