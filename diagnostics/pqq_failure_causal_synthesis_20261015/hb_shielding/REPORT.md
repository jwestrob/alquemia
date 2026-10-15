# Hydrogen-bond shielding supplies no useful ExaA remedy

**Close this finite component question.** Saved-output job 5147 retained all twelve endpoints and reproduced the previous electrostatic-only comparison. Hydrogen-bond burial raises the selected endpoint energies by roughly 12.5–12.8 kcal/mol, largely in common across metals and sources. Its accommodation effect is small, and no selected geometry changes under either mathematical or operational policy. No higher-precision replay or new cavity calculation is justified by this result.

Both sources retain `adaptive_Ca` for Ca and `adaptive_La` for La across all three component sensitivities. The HB-only paired accommodation changes are −0.033159 kcal/mol for sample0 and −0.057864 for sample3. These individual source comparisons pass their recorded numerical gates.

## Raw source-dependent effects

These are **sample3 minus sample0** changes, in kcal/mol. The full paired-pool change includes the origin component shift; the accommodation column subtracts it. Every aggregate value below remains explicitly **unqualified** under the unchanged charge-precision ceiling.

| Component transfer | Paired-pool change | Accommodation change |
|---|---:|---:|
| Gelec | −0.171134 | −0.076458 |
| Ghb | −0.098818 | −0.024705 |
| Gelec + Ghb | −0.269952 | −0.101164 |

Mathematical and operational results agree. The small raw HB source difference provides no supported discrimination gain. Large individual burial shifts alone do not imply a useful metal/source response.

## Qualification and limits

The saved HB translation, rotation and permutation controls pass; their largest absolute residual is 3.25261e−19 Eh against 1e−9 Eh. These are shared controls from sample0/Ca origin, not twelve independently rotated endpoints. Saved native repeat checks pass. HB per-source bounds are 6.333419e−7 and 6.351275e−7 Eh, below the existing 1e−6 ceiling; their aggregate **1.268469400e−6 Eh exceeds it**. The aggregate HB comparison therefore remains raw/unqualified. Gelec and combined aggregate bounds also fail, at 3.738758144e−6 and 3.748323601e−6 Eh. No tolerance was changed.

Union cap-zero and cap-removed controls remain recorded. Local cap-zero values exist only for sample0/Ca origin; **11/12 are unavailable**, and separate cap-zero repeat checks are unavailable. These omissions are not passes. The controls do not establish a physical source-only cavity or harmless charged caps.

This calculation changes only the declared native components at saved vacuum charges. It includes no full-protein SASA, permanent protein electrostatics, electronic response, new geometry, complete solvation model, calibrated class decision or production change.

## Evidence and cost

Authoritative [RESULT.json](../../../workspaces/pqq_failure_causal_synthesis_20261015/hb_shielding_v1/execution/RESULT.json), SHA256 `62b00e1e885f2aa42e0591cc828ef7e9f85fd0ffdec42353b3412b457c071a2a`, retains all twelve rows, raw values, exact source/state pins, both policies, controls and explicit rejection reasons. [EXECUTION.json](../../../workspaces/pqq_failure_causal_synthesis_20261015/hb_shielding_v1/execution/EXECUTION.json) pins config `6905ada8…` and the result.

Job 5147 completed in 3 scheduler seconds on one CPU: 3 allocated CPU-seconds. The script recorded 3.346539 s elapsed and 2.001231 s process CPU over its own measurement interval. **Zero new MACE, SCC, DFT, cavity or optimization calls; zero GPU use.** The report reads saved results only. Previous outputs and the canonical branch remain unchanged.
