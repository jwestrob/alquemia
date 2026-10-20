# Independent terminal review — 15458

## Actual outcome

The declared attempt completed and stopped at its **accepted-step limit**. It did not converge and did not report execution or physical/state failure. All80 new requests completed; all80 steps accepted their first trial with multiplier1.0 and passing physical/state guards. Cumulative counts are172 accepted steps/172 evaluations, within172/292 limits. New calls are80 MACE and80 ALPB, with no DFT/GFN/MM, other-state search or cache hits.

| Evidence | SHA-256 |
| --- | --- |
| Actual `RESULT.json` | `9159ea43cf76692fb9ad40acbf199c1a6b146acdf8dc21ff68a045350a764a10` |
| Actual `SEARCH.json` | `4d59031becf63f38381a33733bbda62fef4ebe48e8ba0cceef8266bd627304ca` |
| Final checkpoint | `14cd220d924040caf046060c0b6fe3869941bff4bfeb9b210cc1a1e74a5eeacc` |
| Final POINT | `67690a33f7233a3e1b17b1abb687d41717fefd5903e3c04199a5060b1d6b7d81` |
| Run configuration | `8a78864e5baa4836c1fef00240023cdc941493f41a9dad2e2a47aef31ea9a997` |

Source remains `551d964b`, optimizer `c5092f24`. Actual preflight/decision bindings, original start/source, settings, initial resumed history/trust, historical request prefix and all80 new point/request/trace joins match. Eighty new POINT recipes match the admitted recipe. This review verified87 selected pins and literal saved fields; it performed no numerical replay, scientific imports, model calls or submissions.

The stored step92-to172 work is MACE **+5.308810226**, ALPB **−20.894473621**, composite **−15.585663395 kcal/mol**. Final maximum/RMS forces are **4.824511183/1.045572389 kcal/mol/Å**, above the unchanged0.10/0.03 criteria. Force changes are nonmonotone: stored RMS is0.629776 at step151,2.121072 at167 and1.045572 at172. Both final physical admission and bound-state labeling pass; model convergence, force accuracy and electronic-state qualification do not.

## What the saved evidence already distinguishes

- Energy lowering continues through the final accepted step: its stored composite work is **−0.247798946 kcal/mol**. The attempt is not shown to be a repeated-cache, zero-step or rejection-loop stall. This establishes decreasing recorded model energies, not reference accuracy or a measured numerical-noise bound.
- Stored trust radii are0.05/0.1Å, ending at0.1Å; the trust floor was not reached. A trust radius is a ceiling, not the actual step length, so these fields alone do not prove clipping or efficient displacement.
- Final raw and rigid-removed RMS are1.045572389 and1.045572171. Removing global translation/rotation does not explain the current large residual. The source correctly checks convergence on the full force.
- A favorable solvent contribution opposing a positive vacuum contribution is a permitted tradeoff in the stated composite objective. These signs alone neither identify a force inconsistency nor show that ALPB or MACE is physically wrong. The new terminal has no independent reference.

## Static optimizer interpretation

The unchanged search uses eight-pair L-BFGS, rigid projection, a descent fallback, a maximum-atom step cap and Armijo sufficient energy decrease. **It does not require decreasing force norms or a Wolfe curvature condition.** A full accepted trial therefore need not reduce the residual. Coupled soft motion and stiff local coordinates, step clipping, inaccurate limited-history curvature, skipped curvature updates or fallback resets could sustain oscillating forces while energy decreases. These are hypotheses for the saved audit, not diagnosed causes.

The helper scales its initial inverse-Hessian estimate with the latest accepted history pair. The optimizer only stores a new pair when its curvature checks pass; it clears history when the proposed direction fails its descent check. Neither final history length8 nor all-first-trial acceptance reveals the frequency of these events. Exact checkpoint replay can resolve them without proposing a new trajectory.

## Minimum useful saved-only audit

Use all81 saved states and all80 transitions, preserving the original denominator and any failed replay. Run numerical analysis on Slurm only.

1. **Replay the actual step.** Bind each archived checkpoint's history/trust and the original helper. Reconstruct projected gradient, raw L-BFGS direction, fallback, second projection and trust scaling. Compare the predicted accepted coordinates and next checkpoint/history with the saved ones. Report actual step norms, clipping ratio, descent slope, actual/predicted decrease and trust updates. Distinguish floating-point replay tolerances from unchanged scientific force gates.
2. **Inspect curvature and motion.** Report history resets/skipped updates, curvature quantities and helper scale; compare successive directions/gradients and actual displacement. This can distinguish small or oscillatory progress from a continuing directed path. Such path-derived curvature is not a physical Hessian or proof of stiffness.
3. **Replay raw components and locate residuals.** Preserve the original MACE-force sign and ALPB conversion order, verify saved composite gradient/energy closure, then report component norms, relative alignment/cancellation and persistent residual atoms/residues. Project both component gradients onto the actual saved displacement alongside each native component energy change. Finite-step mismatch can reflect curvature; it is not a finite-difference qualification or an exact integrated-work identity.
4. **Keep numerical and physical conclusions separate.** Original qualification12958 passed both declared directions at both steps for each of four starts, including initial7CCN–La POINT `50a597c8`; it supports the earlier local force interface. It does not certify every later geometry or every Cartesian direction. A zero-call audit can find replay, bookkeeping or optimizer-efficiency problems, but cannot independently bound terminal force error or establish physical accuracy.

The finite attempt remains closed. Budget exhaustion alone proves neither a bad model nor numerical stagnation. No additional optimization, modified solver/settings, tolerance change or DFT call follows from this review; the declared saved-only audit is the appropriate next discriminator.
