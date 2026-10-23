# Exact matched water removal in the current AMOEBA candidate

Declared before preparation or native execution. Root owns submission.

## Question and fixed six-cell denominator

Evaluate the already prepared water-deletion controls for monodentate acetate, bidentate acetate and acetamide, each with La and Dy. The source panel is `matched_water_removal_v1/prepared_v1/PANEL.json`, SHA256 `0352a54293309d02e84be9a99fd34835f5f0ba2641eef52be078b5f7b295d004`. It deletes only water atoms from the six completed hydrated central geometries. Preserve every surviving coordinate, symbol, atom order, covalent bond and state. Acetate complexes have eight atoms/charge+2; acetamide complexes have ten atoms/charge+3. Physical La1/Dy6 metadata is retained, although AMOEBA has no electronic spin variable.

Bind the completed hydrated AMOEBA result `2d55465328abcd60dcf0ca5b14761d3f4c96314310c7ec0381a2696d8c3f8674`. Copy each matched parent's parameter and key files byte-for-byte. Retain types500/501 for La/Dy and all ligand types; no new POLPAIR, VDWPR, charge, multipole, damping or polarizability. The asymmetric Dy common-water/donor transfer remains unqualified. No ALPB, implicit solvent, optimization, geometry search, fitting or protein scoring.

## Preparation and native execution

Preparation runs serial saved-text/metadata checks only. Verify exact prepared XYZ against the surviving parent lines/coordinates, graph restriction, deleted-water counts, paired coordinates, charges and parameter hashes. Write six Tinker inputs and a pinned manifest. There are zero native energy or force calls in preparation.

After actual preparation review, propose one root-owned Slurm job: standard partition, standard_expand QoS, one CPU/thread, scheduler-default memory/walltime, no GPU and no requeue. Reuse exact admitted Tinker26.2 executables and `RUN_PANEL.py` native parser/invocation helpers. Perform six `analyze E` plus six `testgrad Y N` invocations, at most12 native program invocations, with unchanged300-second per-invocation timeout. Each invocation evaluates the stated potential; invocation counts are not generalized counts of internal energy operations. First native/admission/nonfinite failure stops new scheduling, preserving all six rows and attempted receipts.

The key retains `polar-eps 1e-8`, `digits 10` and `openmp-threads 1`. Full finite analytic gradients, component and total energies, raw native outputs, and actual CPU/wall costs are retained. Admission retains exact agreement of energies printed by the two native programs, full gradient shape, finite values and absence of undefined-parameter/nonconvergence/error messages. No new derivative or physical-accuracy gate is added. Printed polarization iteration residuals remain unavailable.

## Interpretation and collection

Use the separately declared saved comparator to combine these six dry results with existing hydrated and aquo8 endpoints. The metal-dependent water contribution is `[E(Dy,wet)-E(Dy,dry)]-[E(La,wet)-E(La,dry)]` at unchanged surviving geometry. All terms and both metals are required. This total many-body response can contain water–water, ligand–water, metal–water, dispersion and polarization effects; it is not a unique pair energy or a hydration/binding free energy. Shared intraligand bonded terms cancel in the matched wet-minus-dry comparison. No automatic parameter fit or sampling release follows.

No new aquo, hydrated or original dry dimer call is authorized here. Independent protein donor qualification remains a separate next decision.
