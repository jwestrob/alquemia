# Execution repair before full-system solvent outputs

Job 1200695 completed four real Reference-platform core checks. Native versus
custom OBC-II differences: La 0.000039819, Ca 0.000046033 kcal/mol; maximum
force difference 4.34e-7 eV/A. Both pass the frozen numerical gates.

The first full GPU context failed before an energy evaluation with
CUDA_ERROR_UNSUPPORTED_PTX_VERSION. The installed Conda OpenMM 8.5.1 plugin
links libnvrtc.so.13 (13.2.78); the allocated A5000 node has driver 570.195.03.
No full GB energies or classifications were available when this repair was
chosen. Job cost: 8 allocated GPU seconds, 128 allocated core-seconds,
12.959 actual CPU seconds. No scientific input was changed.

Keep the failed v1 immutable. Prepare v2 with the same OpenMM 8.5.1 release
from its official PyPI wheels, selecting the CUDA12 extra and pinning CUDA
12.8 components in a new workspace venv. Inherit the existing environment's
Python/NumPy without modifying it. Wheel files, install logs, package inventory,
Python, solver source and native libraries are pinned. CUDA components:
runtime/cupti 12.8.90, nvcc/nvrtc 12.8.93, cufft 11.3.3.83.

Rerun all 19 declared solver calls because the native build changed, including
the four small core validations before full execution. Thus total planned
successful solver calls become 23, with one recorded failed context attempt.
No additional scientific variant, MACE or DFT calculation is introduced. Keep
the exact geometry, charge distribution, energy expression and tolerances.
