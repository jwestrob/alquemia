# Native dynamics implementation check

Declared21October2026 before execution. Question: does the successful pinned
double-precision GPU build pass its existing integration and periodic pressure
control regression tests before we use its dynamics machinery? This complements
the completed energy/force and POLPAIR checks; it is not a new physical model.

Exactly one native invocation selects four existing upstream cases:
NVE-Verlet-ArBox, NPT-Monte-Iso, NPT-Bussi-Iso, and NPT-VRESPA. Retain every
native reference, tolerance, source input and random seed. No source changes,
parameter fitting, protein trajectory or exchange calculation. Some upstream
fixtures use different water parameters; passing is numerical implementation
evidence only, not qualification of our chosen protein/water Hamiltonian.
NPT-VRESPA checks initialization/settings only; the other three exercise native
trajectory steps. Preserve that distinction when reporting a pass. Capture
successful assertions as XML and require their actual source-declared counts:
40 each for Verlet, Monte and Bussi, and7 for VRESPA initialization.

Use the already qualified double/arch90 CUDA12 build and staged user-space
runtime on one H200/four allocated CPUs, native CPU threads1, standard_expand
and default memory/time. Run in a fresh owned directory. Require process exit0,
all four exact testcase names, nonzero assertions, zero failures and actual
native XML completion. Pin executable, test source and result artifacts.
Preserve failure; do not weaken upstream tolerances or silently omit a test.
