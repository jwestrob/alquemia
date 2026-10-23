# Reference evidence and boundaries

## Published La ligand model

[Wait et al., J. Phys. Chem. B 127, 1367–1375 (2023)](https://doi.org/10.1021/acs.jpcb.2c07237)
provides the exact four starting dimer structures and AMOEBA09/22 parameters.
The authors used relativistic MP2 geometries and counterpoise-corrected
DLPNO-CCSD(T) interaction curves. Supplement s002 states the SARC-ZORA-TZVP
metal basis alongside ZORA-def2-TZVP for light atoms. Its Table S3 reports a
49.09 kcal/mol AMOEBA22 underbinding error for monodentate acetate, despite much
closer bidentate agreement. The supplied mono structure came from7CCO. These
are consumed development fixtures, not independent validation data.

Author-supplied assets, unchanged:

- [Coordinates and parameters, s001](https://acs.figshare.com/articles/dataset/22008454), file39062827.
- [Supporting discussion and tables, s002](https://acs.figshare.com/articles/journal_contribution/22008457), file39062830.
- [Public manuscript](https://www.osti.gov/servlets/purl/2311656).

Original download responses and text extractions are retained in `assets/`.
The Figshare article JSON verifies DOI, filename and download URL. An initial
unverified download URL returned an empty response; it was not retained as an
input. Only the DOI-resolved Figshare assets enter the panel.

## Relativistic recipes

[ORCA6.1 basis documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html)
supports all-electron SARC-ZORA-TZVP with ZORA, SARC/J and explicit removal of
ECPs. `NoECP` and `DelECP` prevent an implicit default ECP from surviving a basis
replacement. The matched frozen-f basis/ECP coefficients are copied byte-for-byte
from the accepted local La/Dy reference lineage and pinned in PANEL.json.

The frozen-f model has effective singlets for both ions. The all-electron scalar
model has physical La singlet/Dy sextet. Neither includes spin-orbit coupling;
the latter still needs determinant/state qualification. This initial comparison
also changes basis and scalar-relativistic representation, so it cannot uniquely
attribute a discrepancy to 4f freezing.

## What is not yet supplied

There are no equilibrium hydration populations, protonation corrections,
counterpoise corrections, correlated references, protein embedding results or
experimental affinity predictions in this initial panel. No result has been
fitted to an expected selectivity. Later parameter fitting must reserve distinct
configurations for validation rather than rebranding these published fixtures.
