# Fixed-geometry CPCM-water solvent response

## Question and frozen selection

Does self-consistent PBE0 CPCM-water support the solvent lowering and metal-dependent deformation response seen with the existing fixed-FF ALPB correction? Exactly ten whole256 states are declared before execution: the original 7CCO and 7CCN before/after geometries, each with La and Dy, plus the exact 7CCN-La step92 and step725 endpoints. The late pair directly tests the previously reported approximately −45.44 kcal/mol ALPB lowering. `INVENTORY.json` binds all ten admitted vacuum references, physical coordinates, electronic states and cached ALPB receipts; selection is independent of new CPCM values.

## Scientific and native settings

Copy each admitted vacuum XYZ byte-for-byte. The complete scientific input must equal its original input with the unique ` EnGrad` token replaced by ` CPCM(water)`. Thus the only scientific changes are removing analytic gradients and adding CPCM water. Preserve PBE0-D4, def2-TZVP/def2-J for light atoms, the exact custom lanthanide basis/auxiliary/ECP assets, RIJCOSX, NoAutostart, DefGrid3, VeryTightSCF, PModel, charge −2 and effective RHF singlet. La uses ECP46; Dy uses frozen-f ECP55 with physical multiplicity6 retained separately. No reused orbital seed, optimization, changed geometry, parameter fitting, new vacuum calculation, ALPB, MACE, GFN or MM request.

Use the pinned ORCA 6.1.1 executable, original native runner and allocation-derived renderer, twelve actual MPI slots on one node, one numerical thread/rank, default memory and time. The original allocation helper observes the scheduler memory grant; it requests no memory override. The intended default is 48,000 MiB and approximately 3,100 MB MaxCore per rank. Actual grants and costs are recorded. Preparation and collection use one CPU; there is no GPU requirement.

## Stages and finite calls

`prepare` runs under Slurm with zero molecular calls. It verifies every declared source/hash, original collector admission, state, all256 coordinates/IDs, exact native XYZ, basis/ECP recipe, and cached ALPB authority before writing ten decks and a manifest. An actual successful preparation and a pinned root release are required before any native call.

`run-one --task CELL` makes one fresh native runner invocation for the named frozen cell, with `task_type=single_point` and `engrad_path=null`. Root may submit all ten independent singleton jobs after the shared preparation passes. The adapter schedules no jobs and performs no retries. Each failure retains its native artifacts and all ten result slots. The total ceiling is ten new CPCM whole-molecule endpoints; no central or old vacuum repeat is permitted. Internal native initialization work is part of those requested endpoints and its cost is retained.

`collect` is saved-only. It admits available singleton receipts independently, retaining unavailable/failed cells and every declared comparison slot. A complete interpretation requires all ten. There is no waiting collector or automatic continuation.

## Admission and comparison

Require native normal termination, SCF convergence, unique finite final energy, the exact ORCA executable/version, RHF charge/electron/ECP state, D4/PModel evidence, original convergence semantics, and CPCM `WATER`, epsilon80.1510, refractive index1.3328, solvent radius1.3000 and GAUSSIAN VDW output. Retain the actual element radii, cavity and solvent-energy lines, all SCF residual flags and warnings. Do not require every printed SCF residual flag to pass. Missing/contradictory native evidence is unavailable, with no rerun for parser formatting. No gradient artifact is expected or admitted.

For each cell report `S_CPCM = E_PBE0,CPCM − E_PBE0,vac`, the cached `S_ALPB`, and their difference. Native printed energy differences use Decimal before the frozen 627.509474 kcal/mol/Hartree conversion. Report before→after works for each source/metal, Dy−La contrasts at each matched geometry, and Dy−La deformation-work differences under vacuum, vacuum+ALPB and CPCM. The mandatory late-La step92→725 table reports the same three works and both solvent increments. Keep all ten cell values and all comparison slots, even on failure. Existing 1 kcal/mol reporting resolution is descriptive, never an admission tolerance or uncertainty estimate.

This compares solvent models at fixed geometries. CPCM's self-consistent density/cavity treatment and fixed-FF ALPB's charges, cavity and non-electrostatic terms differ, so disagreement does not uniquely isolate polarization or prove either model correct. Absolute metal energy contrasts are conditional bookkeeping, not affinities or binding preferences. Conditional termini/protonation, dianion/electronic-state limitations and frozen-f/no-SOC caveats remain. Energy-only results provide no force qualification, equilibrium or production-promotion claim.

## Review and preservation

Root owns releases/submissions and final reports; the adapter performs no scheduler submissions. Source, configuration, specification, inventory and parent evidence are hash-bound. All numerical imports/checks and calculations run only on Slurm. Preserve original references, prior failures and all new native outputs unchanged. A source-only AST check is permitted; no mock or numerical login test is part of this declaration.
