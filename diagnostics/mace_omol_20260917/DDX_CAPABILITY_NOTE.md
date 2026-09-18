# ddX coupling assessment: source inspection, no solvent calculation yet

The native GK decomposition localized an18.129kcal GGR source-self shift;
18.453kcal is attributed algebraically to changed projected charges. The
same-density CHELPG sampling experiment is running. A separate question is
whether a resolved protein dielectric boundary behaves similarly to the GK
approximation. This would help distinguish charge representation from solvent
model effects, without selecting a model by its biological classification.

ddX0.9.0 provides domain-decomposition PCM/COSMO/linear-PB solvers and supports
multipolar sources and quantum-density host coupling. It exposes iterative,
matrix-free and fast-multipole options. Those features make a whole-protein
capability test worth considering; no local runtime or accuracy is established.
[Project documentation](https://ddsolvation.github.io/ddX/).

The interface needs **both** surface potential phi and a density integral psi.
Its energy is0.5*dot(psi,x), with x the solved reaction-potential expansion.
Passing exact ORCA potentials together with unrelated fitted-charge psi would
mix source representations. That is not an exact-density implementation.
The upstream Psi4 interface builds psi by atom-centered Becke quadrature and
separately evaluates surface potentials. An ORCA bridge therefore needs a
validated density import/integration route, including La ECP conventions.
[Theory](https://ddsolvation.github.io/ddX/md_docs_theory.html),
[Psi4 coupling source](https://raw.githubusercontent.com/psi4/psi4/master/psi4/driver/procrouting/solvent/ddx.py).

The rendered Python documentation has stale constructor argument order compared
with0.9.0 source/examples. Follow the pinned actual source: Model(model,centres,
radii,...), State(model,psi,phi,...). The model's multipole helpers require
sources centered on its spheres. No cap nucleus should be silently relocated
to make this convenient. A projected-monopole solvent diagnostic can use the
already validated physical source map, but remains a projected-charge model.
[Python interface](https://ddsolvation.github.io/ddX/label_python_interface.html),
[upstream example](https://raw.githubusercontent.com/ddsolvation/ddX/master/examples/run_ddx.py).

Source archive and three primary references are pinned under
`workspaces/mace_omol_20260917/ddx_capability_source_v1/`. PyPI release0.9.0 sdist
SHA2561f8faf41d61483f063d29934ea969747ea5e060442eed188adbf629e2a945561.
GitHub clone failed DNS resolution; the verified PyPI source download succeeded.
No existing environment modified; no solver installed or scientific solve run.
No model replacement, new score, calibration or baseline change follows.

Next engineering prerequisite: isolated reproducible build/import. A useful
subsequent pilot would hold the actual full GGR cavity, default source charges
and nuclear geometry fixed, comparing source-only PCM reaction energies against
the verified GK source-self values. It must declare radii and convergence/
rigid-transformation checks before solves. Different cavity constructions and
approximations mean agreement with GK is not itself an acceptance criterion.
Exact-density coupling and full permanent/polarizable environment accounting
remain separate required work before any new hybrid prediction.
