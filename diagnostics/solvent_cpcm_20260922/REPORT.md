# CPCM pilot: numerical route failed, discrimination remains untested

**Retain the released ALPB scorer.** Installed ORCA6.1.1 genuinely applies CPCM
to native GFN2, but automatically switches to ordinary ORCA SCF. That solver did
not produce any complete Ca/La contrast in this fixed four-context pilot.
There is no evidence here that CPCM improves or worsens biological discrimination.

All 16 declared calculations ran; **5 converged and 11 failed SCF convergence**.
No failure was caused by scheduler allocation exhaustion. Every failed output
explicitly reports `SCF NOT CONVERGED AFTER 499 CYCLES` under the requested
MaxIter500 and then error termination in LEANSCF. Outer ORCA return codes were
zero even for these failures; the existing runner correctly rejected them using
the scientific termination/convergence checks. No last-iteration energy was scored.

## Actual coverage and what it establishes

| Consumed context | CPCM Ca / La | Matched vacuum Ca / La | Ca/La contrast |
|---|---|---|---|
| 1H4I | converged / failed | failed / failed | unavailable |
| 4MAE | converged / failed | converged / failed | unavailable |
| Q9Z4J7 canonical | converged / failed | failed / failed | unavailable |
| A0A3F2YLY8 canonical | converged / failed | failed / failed | unavailable |

The four successful CPCM endpoints confirm Gaussian-vdW cavity construction,
water epsilon 80.151, surface density 5/Å², unchanged native parameters and a
nonzero printed reaction-field contribution. The frozen installed radii include
Ca 2.772 Å and La 2.400 Å; no cavity settings were tuned. Thus this is a real backend,
not an ignored keyword. Its numerical route is presently unreliable on these
actual contexts. The concentration of CPCM failures in La does **not** establish
an affinity direction or explain a biological failure.

The eight matching vacuum controls used genuine explicit MOREAD from their
archived native vacuum orbitals, ordinary SCF, the same parameters, 300 K state,
coordinates, charge/spin and default convergence tolerance. All eight restart
requests were actually honored; seven still failed. The sole successful
4MAE-Ca vacuum differs from its native-mixer result by **+0.0007018304 kcal/mol**.
That supports numerical agreement for one endpoint, not the unavailable La states.

Only 4MAE-Ca has a complete CPCM-minus-vacuum transfer:
**−351.2564804293 kcal/mol**. This is one endpoint contribution, not a metal contrast
or affinity. The intended score remains
`R = [MACEvac + GFN2CPCM − GFN2vac]Ca − [MACEvac + GFN2CPCM − GFN2vac]La`.
All four R values, new classifications and calibration remain unavailable. Archived
ALPB scores remain separate in the collection. No successful CPCM cell was rerun.

## Measured cost

| Job | Calls | Scheduler wall | Scheduler allocated core-seconds |
|---|---:|---:|---:|
| 1209970 CPCM | 8 | 4237 s | 271168 |
| 1209980 matched vacuum | 8 | 295 s | 18880 |
| **Total** | **16** | overlapping jobs | **290048** |

The inner executor records **289550.752817 core-seconds**; the scheduler total also
includes launch and collection overhead. Sum of individual endpoint wall×ranks is
117529.026472 rank-seconds; allocation totals include idle reserved ranks while
slower endpoints finish. **GPU use: 0 seconds.** Batch-step MaxRSS was 18478780 K and
4698204 K respectively; these are scheduler-reported batch measurements, not an
independently measured whole-node peak. No matched-hardware speed claim is made.

The observed cost and zero completed metal pairs do not support an affordable
production replacement. This rejects the tested ordinary-SCF route for immediate
benchmark expansion; it does not reject continuum solvent corrections generally.

## Disposition and reproducibility

No retries, full30 calibration, 225-fold scoring or production/default change.
Further CPCM work would require a separately justified numerical strategy before
spending on another biological panel. The independent observer completed once,
retained all outputs/receipts and wrote the vault note. Its immutable collection
was inspected rather than recollected.

Final existing-fixture checks: 7 passed; 1 phase-specific unrun-fixture test skipped
because execution has completed. That missing-result test passed before launch.
These include actual CPCM cavity/solver output, real ordinary-vacuum restart,
unchanged source algebra and explicit failure/missing behavior.

Artifact pins: [FINAL_ARTIFACTS.json](FINAL_ARTIFACTS.json). Actual products under
`workspaces/solvent_cpcm_20260922/matched_vacuum_v1/`:

- `final_collection.json`: authoritative matched results, SHA
  `60a80f536d13a0745fee1e0cd535e3f2b1e1422631140e465b5e4dbaedc952b2`.
- `terminal_audit.json`: 16 actual attempts, exact failure messages and receipt pins.
- `final_costs.json` and `scheduler_accounting_final.txt`: endpoint and allocation costs.

Runnable inspection remains in [COMMANDS.md](COMMANDS.md); no new calculation is
needed to reproduce this report from the preserved artifacts.
