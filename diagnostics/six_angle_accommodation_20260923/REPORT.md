# Six-angle extension: more relaxation, no useful discrimination gain

Completed 2026-09-23. **Close this branch and retain the four-angle candidate.**
All 18 searches, 18 cross-metal MACE evaluations and 72 strict-native scalar
cells succeeded. All nine five-geometry pools are available. The intended
C5AX sample-3 improvement did not occur: added Ca and La composite stabilization
almost cancels. No more searches or canonical expansion are recommended here.

## Decisions and raw ordering

Under the **unchanged four-angle strict32 reference**, calls change from
8 correct/1 inconclusive to 6 correct/3 inconclusive, with zero wrong calls.
This is old-band transfer, not evidence that a separately calibrated six-angle
classifier would intrinsically perform worse. No nine-source calibration was
performed. All sources were already inspected development cases; the five C5AX
replicas are one protein. Ca/La in source names means fold conditioning, not
the biological class: A0A3/A0AC are La-class and A8 is Ca-class.

| Source | Four-angle R | Expanded-pool R | ΔR | Old-band call, before → after |
|---|---:|---:|---:|---|
| C5AX La0 | −405452.236296 | −405451.240411 | +0.995885 | La → La |
| C5AX La1, canonical | −405456.379940 | −405457.048000 | −0.668060 | La → Inconclusive |
| C5AX La2 | −405445.107240 | −405443.124752 | +1.982488 | La → La |
| C5AX La3 | −405457.848335 | −405457.860024 | −0.011689 | Inconclusive → Inconclusive |
| C5AX La4 | −405435.358279 | −405435.236918 | +0.121361 | La → La |
| A0A3 Ca1 | −405450.716619 | −405449.152245 | +1.564375 | La → La |
| A0AC Ca4 | −405454.843450 | −405454.741857 | +0.101593 | La → La |
| A8 La1 | −405466.928339 | −405466.759638 | +0.168701 | Ca → Ca |
| A8 La3 | −405467.349478 | −405462.461602 | +4.887876 | Ca → Inconclusive |

R is the model Ca−La electronic contrast in kcal/mol; larger values are more
La-like on this model scale. Even without any calibration, the result gives
little reason to expand: the smallest La-class contrast minus largest Ca-class
contrast **shrinks from 9.080004 to 4.601578 kcal/mol**. Raw ordering remains
correct across these probes, but the margin contracts. C5AX's five-source range
**increases from 22.490056 to 22.623106 kcal/mol**. These are descriptive values
on a consumed, correlated panel, not held-out accuracy statistics.

## Added physical work

Each row below compares the new metal-specific proposal to that same metal's
own old four-angle proposal, with unchanged chemistry and source-q0 bounds.
Values are kcal/mol. The final common-pool score can select a cross-metal
proposal, so these own-proposal works need not equal the selected row-minimum
change. Both mathematical minima and the existing operational selection are
retained in the complete result.

| Source | Ca native | Ca solvent | Ca composite | La native | La solvent | La composite |
|---|---:|---:|---:|---:|---:|---:|
| C5AX La0 | −3.876951 | +3.213260 | −0.663691 | −4.772587 | +1.628137 | −3.144450 |
| C5AX La1 | −7.538008 | +4.438907 | −3.099100 | −6.230590 | +3.845447 | −2.385144 |
| C5AX La2 | −1.638710 | −0.089519 | −1.728229 | −3.032920 | −0.089338 | −3.122257 |
| C5AX La3 | −1.899818 | +1.309999 | −0.589819 | −2.446260 | +1.868130 | −0.578130 |
| C5AX La4 | −5.190603 | −0.733923 | −5.924526 | −6.886464 | +0.145869 | −6.740595 |
| A0A3 Ca1 | −1.263962 | +0.972422 | −0.291539 | −1.778374 | −0.077540 | −1.855914 |
| A0AC Ca4 | −2.063224 | +0.541113 | −1.522111 | −2.924969 | +1.301265 | −1.623704 |
| A8 La1 | −0.035589 | +0.305748 | +0.270159 | −0.850237 | −0.386555 | −1.236793 |
| A8 La3 | −1.583017 | +0.543502 | −1.039515 | −3.814423 | −2.112968 | −5.927391 |

For C5AX La3, the unchanged force-ranking rule adds Asp319 χ2 and Asn275 χ1.
The native differential benefit is +0.546442, canceled by −0.558131 of solvent
differential work. Its final Ca/La composite works are −0.589819/−0.578130.
The native gradient therefore identified real missing relaxation, but **not a
useful omitted discriminatory signal** in this test. No donor was forced into
the selector and no favorable endpoint was selected after seeing the outcome.

For A8 La3, extra La stabilization is much larger than Ca stabilization:
−5.927391 versus −1.039515 composite work. Both native and solvent contributions
favor that change, reducing the raw Ca-class margin. This is not repaired by
keeping the old candidates: the common-pool minima select both new proposals.

Both metals select the new La proposal on C5AX La0/La1/La4. They select their
respective new proposals on C5AX La2/La3, A0A3, A0AC and A8 La3. A8 La1 keeps
the old Ca proposal and selects the new La proposal. Old candidates were never
discarded; all mathematical endpoint row minima are nonincreasing.

## Physical and numerical checks

- All 18 native energies decrease from their warm starts, by 0.035589–7.538008
  kcal/mol. Every search converges under the original precision policy.
- All source/cap/bond/overlap and fixed-atom checks pass. Nine endpoints touch
  the unchanged physical boundary. Maximum final heavy extent is
  0.800000005702 Å, within the existing 1e-7 Å numerical allowance.
- There are 64 recorded infeasible trial requests; the largest intermediate
  extent is 1.089307 Å. The frozen policy constrains final candidates, not every
  SLSQP trial. No trial was clipped or assigned artificial energy.
- C5AX La3 final selected normalized native loads are 0.176700 Ca (boundary)
  and 0.005417 La (interior) kcal/mol/Å; omitted angular maxima remain 4.338860
  and 2.488373. These are vacuum-MACE gradients, not composite forces, entropy
  or proof of a fully relaxed physical site. Full vectors are retained.
- Existing omitted-protein-contact policy was applied read-only to all nine
  sources and all five geometries: complete PQQ mapping/fixed coordinates,
  raw source-heavy mapping and exclusion of atom pairs separated by ≤3 bonds.
  **Zero of the 18 new proposals creates a new <2.0 Å contact** with omitted
  protein heavy atoms. This is a physical flag only; it did not filter or alter
  scores and does not validate the missing scaffold energetics.
- All 72 scalar cells pass actual fresh-SAD initialization, native state,
  parameter, charge, TolE1e-10 and rank1 checks. Reused old108cells remain
  distinct. Scalar qualification is not a solvent-force qualification.

Seven real-artifact tests pass in 6.590 s, no skips. The six-column coupled
constraint Jacobian agrees with coordinate differences to 3.31×10⁻⁹ Å²/rad;
Cartesian-force work independently checks the new gradient projection. Final
tests replay exact old cells, pool signs/unit conversion, required-cell failure
handling and unchanged geometry bounds. No scientific results were fabricated.

## Cost, failures and delivery

Actual new work: **256 search MACE evaluations +18 cross-MACE =274 MACE calls,
72 native GFN2 calls, zero DFT**. q0 and warm-start calls were reused. Summed
search worker time75.979 s; warm GPU worker79.564 s, peak CUDA allocation
5.582 GB. Summed native ORCA elapsed time1386.854 s; scalar executor52.282 s.

Whole allocations, including setup/collection and failed postcompute staging:
GPU1211537:101 s×32CPU; CPU1211603:59 s×32CPU. **Total5120allocated core-seconds,
101requested GPU-seconds.** The canceled pending1211538 used zero runtime.
Local preparation/reporting time is additional and not scheduler-metered.

Preserved technical failures: preparedrun_v1 caught the legacy q0 receipt
schema before molecular execution; run_v2 completed all GPU chemistry, then its
scalar staging required complete archived task metadata (XYZ/input/medium).
Both empty staging directories remain. Recovered `scalar_v3` uses unchanged
actual GPU results and the same72declared cells; no molecular calculation was
retried. GPU job's final FAILED status therefore describes postcompute staging,
not failed native energies. Successful and failed receipts are both retained.

Main artifacts:

- `workspaces/six_angle_accommodation_20260923/run_v2/GPU_COLLECTION.json`
- `workspaces/six_angle_accommodation_20260923/run_v2/scalar_v3/COLLECTION.json`
- `workspaces/six_angle_accommodation_20260923/run_v2/COSTS.json`
- `workspaces/six_angle_accommodation_20260923/CONTACTS_v1.json`
- Compact component/geometry ledger: `RESULT.json`; replay: `COMMANDS.md`.

Protocol `common_six_angular_warm_native_OMOL_strict_GFN2_v1` remains an
experimental completed result. Four-angle candidates, strict32 reference,
historical failures and production defaults remain unchanged. This branch
does not warrant full canonical qualification or a new classifier claim.
