# Solvent preserves metal-specific selection; source sensitivity remains

All eight whole-peptide standalone GFN2 energy/analytic-gradient cells passed.
The same four native MACE energies were reproduced exactly while obtaining
analytic forces. Both metals were evaluated on both intact 256-atom sources.

The composite still selects La-derived 7CCO for La and Lu-derived 7CCN for Dy.
Solvent reduces Dy's energetic advantage for 7CCN substantially; selection alone
should not conceal that sensitivity. The source difference in paired Dy−La
contrast changes from −72.049259 to −72.966739 kcal/mol. Thus adding this solvent
model does not repair source dependence. The common two-source minimum contrast
changes by +19.237990 kcal/mol; this is not a measured affinity or population.

The force is the derivative of the declared MACE-vacuum plus GFN2-water minus
GFN2-vacuum energy. Four complete composite force arrays are available. Their
independent directional derivative check has not yet run; convergence and parser
admission alone do not establish force correctness or biological transfer.

Geometry-only job 12509 produced 64 mapped, independent motions per source:
peptide phi/psi, whole sidechain torsions and target translations. There are no
synthetic caps or added mechanical energies. The frozen LBT3_RESPONSE_DESIGN.md
selects a small common subspace from individual and differential forces, then
checks one real displacement numerically before optimization.

## Actual incremental work

- Native MACE force job 12500: four calls, 12 s on 12 CPUs + one H200;
  energies identical to prior singlepoints, peak CUDA allocation 8,512,327,168 bytes.
- Geometry-only job 12509: 2 s on one CPU, zero energy calls.
- Solvent job 12510: eight analytic energy/gradient calls, 58 s on four CPUs;
  measured total CPU 181.875 s, scheduler peak host RSS 598,248 KiB.

These three allocations total 378 CPU-seconds and 12 GPU-seconds. Including the
previously reported work, this branch has used 786 allocated CPU-seconds and
38 allocated GPU-seconds through job 12510. They are allocations, not utilization.
No DFT, optimization, Mex solvent calculation or production change occurred.

Charged termini, absence of explicit waters, NMR coordination restraints,
MACE physical spin versus original GFN2 f-in-core approximation and unmatched
assay conditions remain explicit limitations. Mex comparison remains vacuum-only.
