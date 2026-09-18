# Build-tool portability correction

1201202 failed its source/executable hash preflight before creating a venv or
running a compiler. `/usr/bin/cmake` differs between the login and allocated
node. The GNU compiler checks preceding it passed. Preserve this failed
manifest and submit script unchanged. Actual cost1wall second,64allocated
core-seconds,0.474CPU seconds; no scientific calls or dependency installation.

Use a new `ddx_software_v2` manifest and `run_ddx_build_v2.sbatch`. Add the
pinned official CMake3.31.6 binary wheel to the isolated venv, and use that
venv's command through an explicit process PATH. All source/other dependency
versions, GNU compiler pins,64CPU parallel build,32GiB request, import-only
validation and no-scientific-calculation scope remain unchanged. Record CMake's
actual version and hash through the preserved wheel and installed artifacts.
This is a technical build recovery, not a changed scientific model.
