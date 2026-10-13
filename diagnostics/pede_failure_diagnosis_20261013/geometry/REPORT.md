# PedE: saved source geometry does not isolate the two errors

**The measured geometry and confidence do not support rejecting only wrong samples1/2 as obviously bad folds.** All five recovered two-Ca PedE models preserve the same PQQ binding motif, similar local confidence and intact broad PQQ connectivity. They differ at fine donor distances. Correct samples share the conspicuous contact and accommodation changes seen in wrong samples. This does not prove that any predicted structure is experimentally accurate, or uniquely distinguish structure uncertainty from electronic/solvent-model sensitivity.

## Actual scope and execution

Slurm3340 measured all five new sources and all five original La-conditioned controls, including the excluded La sample1. It completed raw geometry/confidence and canonical alignment, then retained nine local-reader failures: the scored XYZ has192 atoms while its physical map has189 atoms and separate cap links. Slurm3344 completed only the nineteen unavailable saved local geometries using the existing Kinematics map and archived full-q values. Raw results were copied unchanged. Every saved proposal replays exactly; origin residuals are at most4.45e-16 Å. All ten source rows remain present; La1 remains raw-only with its original exclusion.

Both jobs used one real Slurm CPU on node-64-768g-5, with zero energy, force, optimization or folding calls. Actual process wall/CPU times:3340 3.652/3.895s;3344 0.932/1.091s. Peak RSS136872/78832KiB. These are saved-artifact measurements, not new scientific endpoint evaluations.

## Starting donor geometry

Distances in Å. CN uses the existing typed3.1Å definition. PQQ RMSD is relative to the consumed canonical La sample2 after the predeclared twenty backbone anchors, allowing only schema-equivalent carboxyl-oxygen swaps. It is a comparison with another prediction, not an experimental RMSD.

| New sample / final call | CN | Glu221 OE1 / OE2 | Asn308 OD1 | Asp358 OD1 / OD2 | PQQ symmetry RMSD | PQQ ring-plane RMS |
|---|---:|---|---:|---|---:|---:|
|0 / correct Ca|7|2.271 /3.001|2.371|3.051 /3.905|0.162|0.0120|
|1 / wrong La|7|2.295 /2.959|2.383|3.053 /3.860|0.147|0.0139|
|2 / wrong La|6|2.331 /2.938|2.376|3.107 /4.044|0.210|0.0130|
|3 / correct Ca|6|2.239 /2.989|2.355|3.179 /4.068|0.259|0.0126|
|4 / correct Ca|6|2.368 /2.928|2.383|3.108 /4.012|0.189|0.0129|

All five retain two typed Glu oxygen contacts, one Asn contact and three PQQ contacts. PQQ O1/N2/near carboxyl-O distances remain approximately2.405–2.445/2.409–2.457/2.531–2.580Å. The weak Asp contact crosses3.1Å in wrong2 **and correct3/4**; wrong1 closely resembles correct0. Correct3 has the largest PQQ displacement from the canonical control. Backbone-anchor RMSD is0.022–0.059Å across all five. None of these observations gives a wrong-only failure criterion.

All PQQ schema bonds pass the existing0.90–1.85Å connectivity gate. Ring planarity overlaps the saved La controls (0.0133–0.0153Å for the four admitted controls). Fine covalent distances remain available individually in the result; for example C5–O1 is1.337–1.343Å in new sources versus1.344Å in canonical La2. This shows similar predicted geometry, not independent validation of quinone bond lengths. Named-atom PQQ RMSDs of0.64–0.93Å largely reflect equivalent oxygen label swaps; symmetry-aware values in the table prevent treating those swaps as a different cofactor pose.

## What the saved proposals actually move

Mapped source heavy atoms reproduce each prepared origin exactly. Metal and PQQ remain fixed through every proposal; no PQQ covalent distortion is introduced. Mapped bond-length changes are below3e-15Å. Glu221 and Asp358 dominate the motion, with Asn308 moving when selected by the unchanged mode selector.

| Sample | La-proposal Glu OE1 / OE2 | La-proposal Asp OD1 | Glu chi1 / chi3 change (°) | Asp chi1 / chi2 change (°) | Selected La geometry max heavy motion (Å) |
|---|---|---:|---|---|---:|
|0|2.654 /2.718|2.995|−8.26 /+4.47|+7.17 /0|0.509|
|1|2.680 /2.711|2.886|−9.62 /−5.29|+10.14 /+10.70|0.610|
|2|2.834 /2.575|3.018|−9.76 /+20.92|+11.75 /0|0.784|
|3|2.803 /2.561|3.113|−8.44 /+20.66|+10.08 /0|0.720|
|4|2.597 /2.764|3.035|−5.75 /−1.17|+10.29 /0|0.433|

All La scores select the La proposal. Both metal scores select that same interior proposal for samples2 and3. Those two have similar Glu reorientation and final dentate geometry despite opposite decisions; their unselected Ca proposals both touch the existing0.8Å boundary. Sample1's additional Asp chi2 motion is real, but it cannot explain sample2 as a shared defect. Sample0 also makes Glu distances more equal and shortens the Asp contact while remaining correct. The complete record retains both proposals for every sample, including all Ca-side motions and unselected boundary candidates.

## Confidence and controls

| New sample | Active-Ca pLDDT | PQQ mean pLDDT | Typed donor atoms mean pLDDT | Mean PAE metal→PQQ / PQQ→metal (Å) |
|---|---:|---:|---:|---|
|0|91.24|97.84|97.26|2.71 /2.15|
|1|91.30|97.86|97.29|2.77 /2.23|
|2|91.34|97.90|97.29|2.75 /2.17|
|3|91.46|97.85|97.25|2.67 /2.23|
|4|91.74|97.94|97.36|2.65 /2.08|

No wrong-only local confidence deficit appears. Exact confidence arrays match the model's atom-chain ordering; site-specific residue, spectator and contact-probability statistics are retained. Confidence belongs to the raw fold and is not assigned to optimized structures.

The four admitted original La controls have CN7 and shorter Asp OD1 distances2.818–2.871Å, compared with3.051–3.179Å in the new sources. Their PQQ positions are closer to canonical La2 (other admitted controls0.045–0.078Å). This is a source-family difference. Occupancy, conditioning metal and incompletely known historical folding settings prevent attributing it solely to adding the second Ca. The excluded original La1 has its ion at the remote Asp53/Thr56/Asp59 site, roughly29Å from the metal-facing PQQ atoms; its original exclusion remains justified and unchanged. It also has high local pLDDT, illustrating why confidence alone is insufficient for occupancy admission.

## What admission establishes—and what remains unresolved

Existing preparation checks identity, broad connectivity, typed CN≥6/≤2N, included donor membership, homologous roles, fixed state/maps/caps/waters, no source heavy-atom repair and the remote-spectator boundary. Existing motion checks preserve bonds/fixed atoms and reject new severe overlaps and oversized trials. They do not test fine metal coordination or PQQ bond geometry against experiment.

The measured fine source differences can change an electronic/solvent calculation, but this geometry analysis supplies no causal perturbation or independent electronic reference. Root's separate saved-component result locates the large sample3→1 origin-score contrast primarily in the solvent component; the present measurements do not prove whether that sensitivity reflects physical geometry effects, electronic-state response or model error. Preserve both structure and model uncertainty. No new quality threshold, exclusion, favorable fold selection, fitted boundary or production promotion follows.

## Evidence

- [Frozen measurements](PLAN.md), [local-reader recovery](LOCAL_RECOVERY.md).
- [Compact copied-field evidence](EVIDENCE.json) pins manifests, actual receipts and results and contains all ten source rows' directly measured fields.
- Complete result: `workspaces/pede_failure_diagnosis_20261013/geometry_local_v2/execution/RESULT.json`, SHA256 `b4d1d99ca1e7b7471323312f63aa5ab370df38c3c72e957f4c0cc17d398c0995`.
- Original partial result: `workspaces/pede_failure_diagnosis_20261013/geometry_v1/execution/RESULT.json`, SHA256 `9045340de5214591230352c90159a2d0faa87bea8caea9ee8175807fd85a5bf5`.
