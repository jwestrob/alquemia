# Possible next prerequisite: import the actual saved density

Read-only source inspection,2026-09-18. No dependency install, orbital export,
new density evaluation or new scientific pilot has been performed for this note.

The current ddPCM pilot consistently uses projected source monopoles for both
surface potential and the density-integral vector. It can test a resolved
dielectric boundary, but cannot establish whether the monopole representation
causes the GGR source-self discrepancy. Replacing only its surface potential
with exact ORCA values would mix representations and is not a valid next step.

The pinned Psi4/ddX interface computes the electronic integral by atom-centered
quadrature of the density against `scaled_ylm`; its nuclear contribution uses
the same cavity centers. It separately computes the actual source potential.
A prospective ORCA bridge needs all these pieces, including real cap positions
that are not necessarily physical cavity centers. Moving cap nuclei to a parent
atom is not exact density coupling. Partition/integration weights and support
must be derived consistently before such a bridge is declared runnable.
[Upstream interface](https://raw.githubusercontent.com/psi4/psi4/master/psi4/driver/procrouting/solvent/ddx.py).

ORCA6.1 documents `orca_2mkl` orbital export and a MOKIT path to PySCF, plus
`orca_2aim` WFN/WFX export. PySCF's Molden loader handles angular ordering and
normalization, but Molden does not carry the full ECP definition. Its CORE
section supplies a core-electron count when present. The MOKIT documentation
also warns that MKL lacks ECP information. None of these formats alone proves
correct La effective-nuclear charge or an unchanged density.
[ORCA utilities](https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/utilities.html),
[PySCF loader](https://pyscf.org/_modules/pyscf/tools/molden.html),
[MOKIT utilities](https://jeanwsr.gitlab.io/mokit-doc-mdbook/chap4-5.html).

IOData exposes effective core charges and a WFX parser; this is another import
capability, not a validated choice. PySCF and Psi4 are absent from the current
base analysis interpreter. Any dependency must go into a separate pinned
environment. Avoid implementing several competing import backends.
[IOData effective core charges](https://iodata.readthedocs.io/en/latest/pyapi/iodata.iodata.html),
[WFX parser](https://iodata.readthedocs.io/en/stable/_modules/iodata/formats/wfx.html).

Before an imported density is used scientifically, declare a small real paired
Ca/La qualification with immutable saved orbitals, exact basis/coordinate/ECP
accounting, electron-number and orbital-overlap checks, and comparison to actual
native ORCA potentials at archived points. A file loading successfully is not
a density-identity test. No new SCF, physical-state change, fitted rescue,
full-protein calculation or numerical score is implied by this assessment.
