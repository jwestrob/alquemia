# Pair-scoped direct/mutual CPU prototype

Prepared source only. Root owns building and all numerical execution. The original source tree and binaries are unchanged.

## Exact contract

```
POLAR-PAIR-RESPONSE 500 36 delta_direct delta_mutual
POLAR-PAIR-RESPONSE 501 36 delta_direct delta_mutual
```

Only these unordered atom-type pairs are supported. Each coefficient is the original resolved cubic-Thole POLPAIR coefficient multiplied by `1+delta_role`. The native exponent remains three. D acts on both permanent driving-field channels and their derivatives; M acts on both mutual-solve channels, the preconditioner and the mixed induced-dipole derivative. Neither the intrinsic alpha nor any other force-field parameter changes. No coefficient is fitted here.

Absent records, all-zero records and every unmatched pair retain the original kernel arithmetic. All-zero declarations bypass the new execution-mode guard. Exact metadata parsing rejects null, missing, duplicate, nonfinite, nonpositive-result and unsupported pair values. Missing atom types are allowed so the same key can be used on bare-ion and pure-water fragments.

The five changed files are `polar.f`, `kpolar.f`, `damping.f`, `induce.f`, `epolar1.f`. Energy and analysis sources remain byte-identical and receive the new direct field/kernel through existing calls. The exact stationary functional and source-line map are in [DUAL_FIELD_IMPLEMENTATION_MAP.md](../collective_mechanics_options_v1/DUAL_FIELD_IMPLEMENTATION_MAP.md). The first prototype is confined to the original-Thole MUTUAL, vacuum, nonperiodic double-loop path; unsupported routes fail before induction. It does not enable global alternate Thole.

## Public native audit interface

- `pairresponse(i,k,mode,gamma,selected)`: atom indices; mode 1=D, 0=M. A selected pair returns its actual coefficient. An unmatched/zero pair returns false and gamma=0 (a sentinel, not an undamped coefficient); the original branch follows.
- `dampresponse(i,k,mode,rorder,r,dm,selected)`: selected-pair factors; an unmatched call leaves dm untouched. The existing `damptholed` and `dampthole` expose complete direct/mutual factors with original fallback.
- `derivresponse(i,k,mode,r,xr,yr,zr,rr3,rr5,rr7,sc3,sc5,sc7,rc3,rc5,rc7)`: selected-pair derivative factors, with `rcN=grad(scN)*rrN`.
- `kpolar` prints `PAIR RESPONSE TYPES/GAMMA0/D/M:` followed by the two input types and three native coefficients, or `PAIR RESPONSE ABSENT TYPES:`. The accepted qualification plan combines this actual coefficient readback with the shared-selector source audit and numerical identity/derivative controls; it adds no separate probe executable.

Planck owns `QUALIFICATION_PLAN.md` and its finite numerical harness, including symmetry/unmatched readback. No numerical pass is implied by source preparation. Positive gamma is not proof that the response matrix is positive definite.

## Separate build command, not executed here

From the repository root, after source review:

```bash
cmake -S workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1/source_v3/cmake \
  -B workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1/build_v1 \
  -DTINKER_DIR="$PWD/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1/source_v3" \
  -DCMAKE_BUILD_TYPE=RELWITHDEBINFO \
  -DCMAKE_Fortran_COMPILER=/usr/bin/f95 \
  -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO='-O2 -g -fbacktrace' \
  -DFFTW_LIB=/usr/lib/x86_64-linux-gnu/libfftw3.a \
  -DFFTW_THREADED_LIB=/usr/lib/x86_64-linux-gnu/libfftw3_threads.a
cmake --build workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1/build_v1 \
  --target analyze testgrad --parallel 1
```

This is the existing baseline CMake/compiler/flag/library recipe in a separate build directory. Do not use the upstream RELEASE/Ofast default. There is no install or production-path change. Compiler diagnostics, binary pins and subsequent numerical outcomes must be retained before claiming executable qualification.

The first source-preparation draft (`source_v1`) stopped on one newly added fixed-form line longer than 72 columns, before compilation. `source_v2` fixes that line; `source_v3` additionally initializes all keyword values to invalid sentinels before list-directed reads. Those earlier source copies are preserved, and no native calculation occurred in any preparation.
