# PQQ-MDH fixed-core v3 preparation audit

**Audited:** 2026-09-15T04:09:59Z  
**Protocol:** `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`  
**Pinned implementation commit:** `31e9b99`  
**Preparation manifest:** `prepared/preparation.json`  
**Preparation SHA-256:** `42f5b9e6326d176047929e6c899fbc8bf86526349ffbf8f50a528c957a4e0781`

## Result

The fresh preparation is accepted for the preregistered calibration run:

- 25/25 frozen targets prepared: 11 experimentally Ln-class and 14 Ca-class;
- exactly 50 single-point tasks: one La and one Ca vertical-swap leg per target;
- all source heavy atoms retained, none added, and coordinates preserved within the declared PDB serialization tolerance;
- every La/Ca pair has byte-identical nonmetal coordinates;
- all task input and XYZ hashes match their target manifests;
- no source or synthetic water, point charges, geometry relaxation, or reserved holdout structure was used;
- no untyped nearby O/N atoms or unexpected direct donors occur in any target;
- typed coordination numbers are 7–8 for the Ln class and 6–7 for the Ca class;
- Ln-class cores contain six fragments and 80 atoms per metal leg; Ca-class cores contain five fragments and 71–73 atoms per leg;
- charges are exactly Ln-class La/Ca = −2/−3 and Ca-class La/Ca = −1/−2;
- the closest atom pair belonging to different retained fragments is 1.893 Å; no cap or fragment collision was detected.

The D+2 acidic Asp is present in all 11 Ln-class cores and absent in all 14 Ca-class cores, exactly as frozen. Consequently, a successful energetic separation calibrates the operational energy scale but does not establish energetic information independent of the already class-informative motif.

No ORCA calculation had been launched when this audit was written.
