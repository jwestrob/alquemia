# Stage1 source review: complete periodic operator parity

No blocking source issue was found in the frozen Stage1 package. This review supports one separately released compile/link and one separately released explicit-operator program at the old e41f La state. It executes neither action and grants no numerical or physical admission.

## Lineage and isolation

The prepared files join the exact e41f selected artifact, BIO composition d8622d40, La geometry b59954d7, job20053 input4b166082 and initialized PME export6689dd39. The original composed key is retained byte-for-byte with only the already admitted PME settings appended. Both original response records remain present. The literal expectations reproduce all11958 coordinates/types/intrinsic alphas and the saved normalized metal stiffness; native h/H are calculated by the original geometry routine, not copied into the operator from the certificate.

The setup clone is exactly the original checked routine except its separate symbol and replacement of the final backend-guard call with comments. The immutable library61e2fab0, normal setup and both fatal guards remain unchanged. The executable calls initial/getxyz/mechanic, frame initialization, the separate checked setup, one dfield0c initialization and the bounded mutual/reciprocal paths. Inspection of dfield0c/ufield0c/umutual1 confirms no induced solve, energy or gradient call is hidden in those field paths. No original energy/force entry is overridden or enabled.

## Operator correspondence

The full action is correctly z−sqrt(alpha)*field(sqrt(alpha)z) for both unequal channels. Native umutual1 subtracts reciprocal dipfield, so the independently replayed reciprocal contribution enters the normalized matrix with a plus sign. Common ufield0 subtracts beta times each induced channel exactly once. The separate primitive chain uses the same Cartesian/fractional transforms, grid_uind, forward FFT, unsymmetrized qfac, unnormalized backward FFT and fphi_uind. The Python sparse scatter and complex FFT replay retain both packing channels and odd-q leakage; the symmetric lower CSC action subtracts its diagonal once.

Ten probes in three packed modes produce30 complete calls and30 separate reciprocal-reference calls:60 mutual reciprocal primitive applications, plus one direct/permanent initialization. There are zero induced solves, E/G calls, sparse assemblies or factors. Per-channel comparisons retain the inherited scaled1e−11 bound. Existing norm admission is reused only after exact snapshot-byte equality; large norm sums are not silently recalculated or relaxed.

The native snapshot/readback checks exact map/alpha/xyz/coefficient identities, PME state and raw finite outputs; qfac/splines/grid map/alpha/h/H remain checked throughout. Explicit u/v and FFT scratch are restored between applications. The compiled field path has no coordinate/cell mutation; changed inputs cannot reuse this one-state lineage. Saved replay requires complete ordered rows and terminal count/scope tokens and preserves failure. Runtime input/executable/release hash validation and receipt admission remain the root launcher's responsibility, explicitly stated in the source handoff.

## Coverage limits and build

Stage1 deliberately initializes raw Ed/Ep with dfield0c and reports h separately. It does not numerically exercise the common dfield0 addition of h during a solve. The common wrapper's equal addition is unchanged source; later Stage2 owns that solve-level check. This narrower coverage is explicit in SPEC and is not a blocker for the operator-only question.

The recipe is one GNU Fortran13 external compile/link against the original module/archive/FFTW dependencies, with unchanged generic flags and no archive rebuild. Source line lengths fit their fixed/free forms. Prescribed native/runtime threads are one. Allocated-node loader compatibility, executable pin, actual limits and releases are still unbound and must be recorded before application. The32 existing source/manufactured controls were inspected without rerun.

This is a qualified-source handoff only. It does not establish native parity, a uniform application-error bound, energy/force correctness, another candidate or Dy state, trajectory stability, MD throughput or physical accuracy. Planck authored the sparse assembly/static scheduling adapters and independently reviewed Singer's new native source; Arendt owns the separate numerical contract.
