# Complete the missing Dy forces on the two existing whole-peptide reference pairs

## Scientific question and fixed scope

Do the previously reported Dy deformation-energy agreements also hold for forces? Are the MACE–PBE0 force discrepancies shared by both metals, or do they change the force difference between La and Dy?

Use the already consumed eight cells: 7CCO and 7CCN, each with its original before/after coordinate pair and both metals. The pairs have different provenance: 7CCO follows the saved hybrid La search; 7CCN follows the saved whole-MACE La search through step80. They are two deformations of the same peptide, not independent protein validation. Neither pair is replaced by the later step725 endpoint.

Supply only the three missing whole256 MACE force arrays: **7CCO/before_Dy, 7CCN/before_Dy and 7CCN/after_Dy**. Their archived MACE requests were energy-only. Reuse the other five MACE force arrays, all eight admitted whole256 PBE0 gradients and all eight saved fixed-charge ALPB gradients. Preserve the original energy-only records; new force records are separate supplements.

The source-only inventory must pin the exact complete collectors, configurations, force tables, native endpoint identities, POINT/task/coordinate/physical-ID records, MACE model and solvent receipts before execution. If this inventory finds an exact compatible saved force array, reuse it and reduce the call count; do not repeat available calculations.

## Execution and admission

One finite Slurm GPU job, at most **three** serial native MACE energy-and-force requests using one warm model. Run all saved input/state/geometry/provenance checks before model loading. Use the existing `LBT3_FD.py::NativeMACE` implementation `0e48b8c6c2f271512f2c445f826036d8a68f7dca0aabc5923697ab3ecb7d9896`, original whole-force configuration `d3c28300c8c7ef88ef747b733603a344790b932db56aa53b0ccd4e999a7a4023`, original model checkpoint, float64 and omol head.

Use the original MACE environment: torch2.8.0, mace-torch0.3.16, ase3.29.0. Retain charge−2, physical Dy multiplicity6, target index255 and the exact256 physical atoms. Preserve all whole-peptide state checks, force sign convention, model/parameter checks and coordinate units. The new energy must agree with its original saved MACE energy within the existing **1e−7 eV** force-replay tolerance. Record each difference. This is a reproducibility gate, not a physical-accuracy tolerance.

Use one H200 on the qualified host with12 allocated CPUs, `standard_expand`, existing scheduler defaults and original thread settings. No memory/time overrides. All numerical work, including saved-vector comparisons, runs on Slurm. No scientific import or numerical replay on login.

There are zero new DFT, ALPB, GFN, MM, geometry, optimization or other-state calls. No new native reference admission rule: consume the pinned complete collectors and their admitted gradient tables, retain their actual native receipts, states, residual flags and qualification limitations. Do not modify old adapters/configurations to bypass their environment or provenance checks.

Require finite full256×3 forces, exact physical IDs/order, unchanged coordinates/state and original saved-energy replay. Preserve partial/failed results and stop further molecular requests on an admission failure. No automatic retry or adjustment of gates.

## Output and interpretation

Report the complete eight-cell denominator, marking unavailable cells explicitly. Retain full per-atom force vectors, maximum and atom-vector RMS residual forces, and MACE-minus-PBE0 force-error distributions for vacuum and for the shared-solvent composite. Keep each source and geometry separate.

Also compare the **Dy-minus-La force vector** at each of the four matched source/geometries under the two models. Report its full vectors and maximum/atom-vector RMS discrepancy. These differences test whether errors cancel between the metals; they are not an affinity or preference. Shared ALPB cancels from the model-versus-reference force error, including its matched metal difference, so it receives no independent validation.

Reuse only the original frozen directions for supplementary projections. Do not select directions, atoms, tolerances or favorable survivors to improve agreement. Retain the original saved energy/work comparisons unchanged. There is no new accuracy-pass threshold and no force-qualification promotion.

Physical Dy multiplicity6 in MACE and the native frozen-f effective-singlet Dy reference remain distinct representations. Retain conditional termini, vacuum dianion, post-added unpolarized solvent, ECP/no-SOC assumptions, nonstationary endpoints and all native SCF/derivative limitations. Completing these comparisons does not establish a physical minimum, validated optimizer or natural-LanM transfer.

This independent fixed-force completion does not depend on the running step725 La reference15545. It does not restart the challenged hybrid or closed Mex derivative branches. Root owns execution and subsequent interpretation; the broader saved mechanics plan remains active.
