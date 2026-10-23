# First explicit polarizable model implementation

21 October 2026. Stage 3 implementation pilot, under the newly activated free-energy goal. This is a candidate interaction model, not a qualified selectivity predictor.

## Question and exact starting model

Can an isolated, native Tinker implementation evaluate real La/Dy aquo and donor systems with finite, internally consistent induced-polarization energies and forces, retaining the published La donor corrections?

Use the complete AMOEBA09 parameter file and AMOEBA22 key block supplied in Wait et al., DOI [10.1021/acs.jpcb.2c07237](https://doi.org/10.1021/acs.jpcb.2c07237), SI `jp2c07237_si_001.txt`. Preserve its atom types, multipoles, ligand polarizabilities, buffered 14-7 vdW, and La POLPAIR/VDWPR interactions. Do not substitute the current repository's differently numbered `amoeba09.prm`.

For Dy use the exact ion parameters in Cisneros's author-contributed `amoebail22.prm` in Tinker commit `87050685eff8840d312e2a332cc82c33f63c7c3d`, remapping only its ion type/class to an unused type. The Dy ion has charge +3, polarizability 0.7280 Å³, Thole parameter 0.1390, buffered-14-7 diameter 3.2300 Å, and well depth 15.0000 kcal/mol. Its permanent dipole/quadrupole are zero. It receives the native combining rules and no invented La-derived donor override.

**Transfer hypothesis:** the repository's Dy asset accompanies GEM-derived water, whose multipoles differ from the La SI water. Applying Dy to the common La-SI water and AMOEBA22 donors is a deliberately unqualified transfer candidate. A matched Hamiltonian is executable but published paired protein accuracy is not established. Do not label this hybrid parameter provenance as a previously validated published La/Dy model. Original Dy/GEM-water controls may be prepared separately, clearly outside the common-water candidate.

## Fixtures and bounded execution

Reuse the chemistry agent's exact `reference_chemistry_v1/prepared/PANEL.json`: four published ligand dimers (water, acetamide, mono- and bidentate acetate), three fixed radial placements, both metals; plus paired explicit 8- and 9-water clusters. Preserve all 28 geometry slots. Reuse the corresponding four fixed free-ligand geometries to define interaction energies; the isolated spherical ion has zero self energy in this model. There is no ALPB, implicit solvent, MM/MACE sum, external field, geometry optimization or fitted parameter.

Local allowance: at most four CPU build processes, followed by serial one-thread tiny smoke/derivative checks on the central water and bidentate-acetate pairs only. No login GPU. Record every command, output, failure and elapsed cost. Root owns the full Slurm panel and any subsequent dynamics. Compile numerical sources unchanged with portable `-O2`; preserve the original failed CMake build and the build-list-only repair adding omitted `uatom`/`uniquetyp` source files.

The root-ready pilot evaluates each selected structure's total/component energy and full analytical Cartesian gradient. On the four central water/acetate-bidentate La/Dy structures also run native central finite differences at 1e-5 Å. Use native mutual induced dipoles, `polar-eps 1e-8`, and one numerical thread. A fixed numerical check accepts maximum absolute analytic-minus-numeric gradient residual at most 1e-3 kcal/mol/Å. This is an implementation gate only, frozen before the first energy call. Native failure/nonfinite output/missing parameters fails the affected slot; retain the full denominator and do not silently substitute values or retry with changed physics.

## Interpretation and next gate

Report full potential terms, induced polarization, force convention, net charge, pairing/provenance, and measured cost. Shared water/ligand intramolecular energies cancel in fixed-geometry interaction differences, but polarization is many-body and does not admit a unique additive physical attribution. The published monodentate-acetate underbinding is already known and must remain visible.

No protein simulation, hydration free energy, preferred coordination number, experimental affinity or production promotion follows from finite small-cluster energies. The reference-chemistry comparisons must decide whether a common La/Dy donor model needs further parameter development. Protein terminal chemistry and atom typing must be resolved separately before assembling the experimental construct.
