# Same-pool electronic-model sensitivity

Declared17October2026UTC before reading/comparing outcomes. Protocol
`lanm_monomer_fixed_pool_electronic_substitution_v1`. Root owns execution.

Question: is the weak/source-sensitive Mex-La/LBT-Dy contrast specific to MACE's
local electronic surface, or does it persist with the already-executed GFN2
surface under identical geometry pools, mechanics and fixed-charge solvent?
This is a model-sensitivity diagnostic, not a competing validated classifier.

Use exactly all36 LBT cells from12720 (18geometries) and10 Mex cells from12731
(fivegeometries). Keep all15 existing LBT source/menu combinations. Existing
MACE components, GFN vacuum/embedded energies, MM ledger and fixed-FF ALPB
components are archived and admitted. Pin their state, coordinates, charge,
environment and source records; missing evidence leaves its menu unavailable.
No new molecular calls, geometry, optimization, preparation or reference offset.

Replace electronics only:

- LBT: E_GFN = E_GFN2,vacuum(whole256,q=-2) + E_fixedFF_ALPB.
- Mex: E_GFN = E_GFN2,embedded(core209,q=0;1840PC) + E_MM + E_fixedFF_ALPB.

The Mex embedded-minus-isolated correction cancels when MACE vacuum is replaced
by isolated GFN2. Do not add it again. Do not add GFN water solvent or a peptide
MM term. All charges, state inventories, caps, cavity and other interactions stay
unchanged. Use the same archived standalone xTB6.7.1/author parameters, effective
singlet/UHF0, accuracy0.0001 and300K. This spin-independent/f-in-core Hamiltonian
does not validate physical Dy sextet or spin-orbit behavior. No whole-protein SCC.

Report R=(Dy-La)_Mex-(Dy-La)_LBT in kcal/mol. For each existing menu compare
original MACE-composite R, GFN substitution at the SAME MACE-selected geometries,
and GFN reselection within the SAME compatible pool. Keep positive/negative and
within1kcal/mol resolution visible, plus source-local spread and selected geometry
identities. Do not hide source sensitivity behind the global pool. Use exact native
component differences before conversion/summation to limit cancellation. Replay
old selections/components with the existing arithmetic policy; report residuals.

Use one Slurm CPU and actual execution receipt; no time/memory override. Exclude
future matched-objective search geometries lacking archived GFN energies. A sign
or magnitude difference establishes sensitivity to the electronic model; agreement
does not identify which model is biologically right. No new thresholds, fitting,
model selection by favorable sign, production change or automatic new calculation.
