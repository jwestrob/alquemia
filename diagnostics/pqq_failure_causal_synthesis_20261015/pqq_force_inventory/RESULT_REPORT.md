# PQQ forces: shared shortening load, with a separate pose signal

Job4598 completed all eight saved-array analyses and passed the frozen frame, normalization, rigid-transform, paired-coordinate and fixed-PQQ checks. It made zero molecular calls. [Exact copied values](RESULT_SUMMARY.json) retain the result/execution pins; elapsed analysis time was0.2430502469651401s and process CPU time0.09686571199999999s. No new numerical calculation was performed while writing this interpretation.

## Main finding

Both quinone oxygens experience shortening-directed MACE force in both metals, both sources and both donor geometries. Positive projection meant elongation, so the negative sign is unambiguous. For C5–O1, the origin projections are:

| Source | Ca, eV/Å | La, eV/Å |
|---|---:|---:|
| sample0 |−4.738040593289152|−2.9817968237938057|
| sample3 |−4.827402276222341|−3.335718514441759|

The corresponding fixed C5–O1 lengths are1.3347501638883574Å and1.3365287127480663Å. This establishes the proposal Hamiltonian's local preference at those coordinates; it does not establish an experimental bond length, wrong redox state or defective AF3 structure.

The classifier-relevant distinction is much less specific than the large individual loads suggest. The saved source3−source0 **Ca−La** C5–O1 projection is+0.26456000771476385eV/Å at origin and+0.04603095557557113eV/Å after donor accommodation. For C6–O2 the corresponding values are+0.06037166923593462 and+0.046594860504517666eV/Å. Thus the strongest signal primarily indicates a shared constrained cofactor load. It does not single out the ambiguous source as uniquely needing carbonyl shortening.

## Pose deserves separate consideration

Rigid PQQ rotation about chemical-frame e3 has a persistent source-dependent metal response:

| Geometry | sample0 Ca−La | sample3 Ca−La |
|---|---:|---:|
| origin |+0.018575221546511672|−0.21733540506047455|
| adaptive_La |+0.052202049477703114|−0.16809750339072496|

These are normalized generalized forces in eV/Å of collective heavy-atom displacement. They are **not** the raw torque, which is separately retained in eV/rad. All27 physical PQQ atoms, including three hydrogens, contribute; normalization uses24 heavy atoms. The result supplies an exploratory pose-coupling hypothesis that is more source-specific than the common quinone shortening load. It does not establish the energetic size, curvature, favorable finite direction, solvent response or predictive benefit of that motion.

Donor accommodation changes C5–O1 force more for target La than control La, but substantial shortening force persists everywhere. This is evidence of coupling within MACE, not proof that fixed cofactor strain caused the donor-work contrast. Matched DFT nuclear forces remain unavailable8/8; the earlier DFT energy works do not validate these PQQ projections.

## Smallest useful follow-up, if pursued

Do not promote a carbonyl relaxation rule or add free PQQ optimization on these data. A narrowly informative exploratory test would use **one symmetric, predeclared rigid e3 pose perturbation** of intact PQQ, on both saved origin/adaptive_La arrangements in both sources and both metals. Two signs give16 new fixed MACE points, reusing the eight centers; all other atoms and chemistry stay fixed. Freeze a common small angular magnitude and existing geometry checks before execution. Both signs and both sources must remain in the denominator.

That finite test can measure whether the omitted pose coordinate actually changes matched Ca−La donor work, rather than infer energy from force. It is not yet a classifier test: any scoring interpretation would still require the same matched solvent treatment and independent reference support. No displacement magnitude, new coordinates, executable molecular manifest or submission has been prepared here. The broad discrimination task may reasonably prioritize other directly supported error mechanisms first.

The alternative shared-carbonyl hypothesis would first need independent force/finite-response support under the same physical state. This result alone cannot distinguish a common source-preparation limitation from a model-specific preference, and it supplies no license to change protonation, oxidation state or label-dependent geometry.
