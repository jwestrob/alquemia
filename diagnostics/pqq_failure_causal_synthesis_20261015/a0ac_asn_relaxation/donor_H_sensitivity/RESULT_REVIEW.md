# Two-H sensitivity: the matched vacuum response survives

Job5039 completed all12 projected-H MACE endpoints and reused all12 original endpoints. Shortening only Asn290 HD21/HD22 to the explicitly mapped ff19SB equilibrium of1.01Å lowers every endpoint by about20–21kcal/mol. Much of that correction is shared across metals and geometries. The alternative-amide path response remains strong in both sample4 sources; inherited N–H stretching is not its sole cause in this model.

## Actual matched path responses

All numbers are saved5039 values in kcal/mol. `D` is Ca work minus La work, with each work measured from the matching original or projected-H origin. “Relaxed-La” is the saved La-generated geometry evaluated with both metals, not a newly relaxed projected-H structure.

|Source / destination|Original D|Projected-H D|Change in D|At least declared1kcal resolution?|
|---|---:|---:|---:|---|
|Ca4 / fixed flip|31.825581321|31.058360456|−0.767220865|No|
|Ca4 / relaxed-La|32.343120936|31.042278153|−1.300842783|Yes|
|La4 / fixed flip|31.715107138|31.209853390|−0.505253747|No|
|La4 / relaxed-La|33.560259733|32.614291817|−0.945967916|No|

The Ca4 relaxed path is meaningfully sensitive at the predeclared interpretive resolution; it should not be described as invariant. Nevertheless, its large positive differential persists. Projected-H fixed flips still raise Ca energy by4.764532825/4.864528240 while lowering La by26.293827632/26.345325151 for Ca4/La4 respectively. Projected-H relaxed-La geometries lower both metals: Ca by12.011967730/13.991429439 and La by43.054245883/46.605721256. These preserve the original qualitative vacuum response.

## Per-metal absolute H corrections

Correction is projected-H minus original-H energy, in kcal/mol.

|Source / geometry|Ca correction|La correction|
|---|---:|---:|
|Ca4 origin|−20.382731244|−21.263341170|
|Ca4 fixed flip|−20.303712929|−20.417101990|
|Ca4 relaxed-La|−20.784354261|−20.364121405|
|La4 origin|−20.504935870|−21.254647017|
|La4 fixed flip|−20.078425462|−20.322882860|
|La4 relaxed-La|−20.432274573|−20.236017803|

This supports retaining the missed-local-basin explanation **as a finite MACE response**, while keeping hydrogen preparation as a quantitative limitation. The response occurs in both source conditionings; it is not unique to the earlier inconclusive fold.

## Scope and admission

Actual5029 preparation preserved all heavy atoms, caps, cofactor atoms, proton inventory, charges and N–H directions, changing only the two physical H distances. Both neutral ASN template bonds and exact source H parents were verified. All six physical frames/twelve exchanged-metal cells passed the unchanged full-source relative severe-contact gate; no new severe contact was recorded. The contact union retains the corresponding saved204-atom frame and all original exterior source atoms. This is not a clash-free or equilibrium certificate.

No optimization, SCC, solvent or DFT followed H projection. These are fixed same-path comparisons, not a replacement common pool. They do not establish projected-H optimizer minima, a composite response, classifier accuracy, or that1.01Å is the exact equilibrium under MACE. The original `forcefield=None` preparation and matched original-H results remain preserved; no class/band calculation was performed here.

## Evidence and actual cost

- [Preparation5029](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_sensitivity_v1/execution/PREPARATION.json), SHA256 `1320a71b0420542029723b64906b57a415807052cdce01f3d8f5dbd29f461f76`.
- [5039 configuration](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/CONFIG.json), SHA256 `f5532b200c53ed92b0fd43eda8330cfaffd7f794af344931046ccbf2b7708909`.
- [5039 collection](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/execution/FINAL_COLLECTION.json), SHA256 `6d2d49b868ebaeb5bfe2de8e2358c690296d498bf221308ed815cfee1cf0be67`.
- [Preflight](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/execution/PREFLIGHT.json), SHA256 `e2989740bcff20c32a385f141ec6a7f67935a8c238f5ca6c0cab477a8d6228ad`; [execution](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/execution/EXECUTION.json), SHA256 `dd33aafd26ed91ebbaa0813cd775e0e37f081f878e8ff232b80225d87e7664ea`.

Actual scheduler accounting: COMPLETED0:0,211s×7CPUs=1477 allocated CPU-s. Execution records12 actual model forwards,211.629774s internal wall,2.696698s parent CPU; zero optimizer/SCC/DFT calls. This review copied saved values only; no new numerical calculation or molecular evaluation ran on login.
