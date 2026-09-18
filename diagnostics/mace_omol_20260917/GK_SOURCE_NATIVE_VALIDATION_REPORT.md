# Native source-solvation accounting confirmed

All15 prescribed static calculations completed in job1201172: five GGR/alpha
physical cavities, each with Ca source only, La source only and no permanent
charges. No DFT, MACE, density query, induced solve or field calculation.
The analytical source-self expression agrees with every actual source-only
minus empty-cavity energy and Ca-minus-La contrast: maximum error
6.32383034826e-13 kcal/mol, below the frozen1e-7 limit.

Every actual physical coordinate, base cavity parameter and Born radius agrees
with its parent. Source moments and zero exterior moments match exactly;
all induced arrays remain zero. Empty-cavity subtraction removes the native
nonpolar term. The isolated frontend differs only in explicit zeroing modes,
with the same pinned native library and original static energy routine.

This confirms the **+18.129021 source-self / -5.528812 cross** split of the
+12.600209kcal GK contribution to2FW0-minus2FVY. It validates component
accounting, not the physical accuracy of GK or its projected charges, and it
does not prove that the30A Born bound causes the failure. No score, calibration,
state, radius or baseline changed. The broader hybrid remains rejected for
promotion after its parvalbumin and GGR transfer failures.

Three real-input/output tests pass in35.054s, including actual energy parsing
and an explicitly corrupted completion marker. No final integration skip.
v1was never executed: review caught the moment parser's requirement for the
pre-energy prefix. The corrected immutablev2 was used for all15calls; no
scientific retry. Initial preflight2tests passed before the scientific run.

Cost:49wall seconds,64CPUs,3136allocatedcore-seconds,141.559actualCPU seconds,
zeroGPU. Native processes23.502209wall/118.974290CPU; kernels12.004682wall.
Batch peakRSS108472KiB. Build/preparation/analysis costs are separate receipts.
All15calls completed on their first attempt. Full outputs/receipts:
`workspaces/mace_omol_20260917/gk_source_native_v2/collection_job_1201172.json`.
Compact checks and cost: [RESULT](GK_SOURCE_NATIVE_VALIDATION_RESULT.json).

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/gk_source_native_v2/implementation/mace_gk_source_native.py dry-run --manifest workspaces/mace_omol_20260917/gk_source_native_v2/manifest.json
```

Next: the separately declared GK_SOURCE_FACTORIZATION_PLAN uses saved arrays
to distinguish projected-charge changes from distance/effective-radius changes.
No new score or scientific executor is required by that diagnostic.
