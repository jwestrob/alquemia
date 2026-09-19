# Electronic audit, 2026-09-19

The working baseline is native ORCA6.1.1 r2SCAN-3c/CPCM(Water)/DefGrid3,
with its own def2-mTZVPP basis, D4 and gCP terms. Archived real La output
confirms def2-ECP replacing46 electrons; Ca is all-electron. The new test keeps
that relativistic core representation but changes basis and electronic method.
It therefore tests the combined electronic recipe, not exchange-correlation
error in isolation. A changed score cannot by itself identify a unique cause.

Current native evidence was inspected in
`workspaces/ggr_mechanism_20260915/stage_c_tasks_v1/extended_center_La/endpoint.out`.
The actual first-panel baselines and source inputs are pinned in the new
`direct_v1/preparation.json`; the six copied XYZ files are byte-identical to
their real source files. Existing amide repair remains in force.

## What was already tested

Actual legacy B97 outputs terminate normally and reproduce published old scores:
XoxF21.460266, MxaF−7.417821, tannase22.717914, calexcitin11.321841kcal/mol.
See [machine-readable extraction](LEGACY_B97_AUDIT.json). Inputs explicitly
override La with def2-TZVP/def2-ECP and use an8-coordinate aquo reference.
These are historical sensitivity results, not independent validation of the
current native fixed-core method or current PQQ bands. No completed CC
benchmark of today's selected cores was found in the relevant local records.

## Basis and correlation accounting

The ORCA6.1 manual's default Ca frozen core is10 electrons, not18; thus its
default already retains Ca3s/3p correlation. The proposed input makes this
explicit and supplies Ca cc-pwCVTZ functions. La's46 frozen-core count includes
its46 ECP electrons, leaving all explicit La electrons correlated. No global
NoFrozenCore shortcut is used with ligand valence-only bases.

ORCA's documented CC/CPCM PTES includes correlation-dependent solvent response;
it is not simply gas-phase correlation added to an HF-CPCM number. Actual
parsed support with iterative triples will be checked in execution. The
finite-basis/PNO/CC response approximations remain, as does the possibility of
inadequate single-reference treatment. This check is an independent diagnostic,
not automatic promotion of a more expensive method.

Source guidance and frozen comparisons: [PLAN.md](PLAN.md).
