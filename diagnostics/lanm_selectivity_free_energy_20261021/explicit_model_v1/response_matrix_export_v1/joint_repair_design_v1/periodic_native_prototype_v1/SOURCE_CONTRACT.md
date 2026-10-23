# Source handoff and execution boundary

## Implemented in the isolated tree

`source_v1/MANIFEST.json` lists every source output and the exact parent, compact-kernel and review authorities. Only `polar.f`, `induce.f`, `epolar1.f` and `epolar3.f` change. The patch includes complete source; no external dynamic library, dependency installation, candidate file or new build-system target is required to read it.

The BIO namespace dispatcher accepts the existing composition adapter's physical `a,b,c,d` record. It validates roles by type **and** element, graph reciprocity, literal water neutrality, unbonded charged ions, population, active positive-alpha sites, coefficient domain, and native mode. This first scope is one metal, at most 3980 waters, seven Na and ten Cl; it does not generalize to arbitrary electrolyte, zero-alpha sites or donor topology. The response arrays are atom-sized for h/H, and water-sized for pair geometry and both contraction passes. They rebuild at every induction setup; there is no stale geometry cache or image counter.

The compact helper returns the admitted native F/F' arithmetic for r≤20 Å and the product-rule continuation through 21 Å. Each oxygen's single nearest image supplies its environment weight and h. All primary waters remain in the selected-water table; outside support their F, F', phi contractions and explicit derivatives are zero. Other-water effects enter through positive prefix/suffix sums. The stiffness sits once at the primary metal. Native salt polarization remains in the full operator; salts receive no added h/H.

The original `dfield0` and `ufield0` wrappers are byte-identical. The common gradient hook runs after native real/reciprocal/frame derivatives; the periodic helper uses distinct u and v, and adds oxygen derivatives, metal reactions and `sym(r_image ⊗ G)` virial. The common pairwise analysis hook adds only the existing `−K/2 u·h` term. The alternative nonpairwise energy path receives no duplicate addition. The original vacuum geometry, spectrum, report, native field kernels and permanent-frame/moment routines are unchanged; the original vacuum gradient body follows a separate periodic dispatch.

## Explicit disabled execution contract

1. Active BIO setup computes/validates the current compact geometry, then calls `periodicresponsebackendguard`. It always reports `backend binding UNIMPLEMENTED` and calls native `fatal`. There is no key, environment variable, runtime flag or successful-PCG escape.
2. A second fatal `periodicresponsestability` hook is inserted after `dfield0` has prepared the current native direct/PME state, before native direct guesses and the mutual iterations. It cannot be reached until a later reviewed source change supplies the first backend admission. Static uniform-field addition follows that point; it does not alter the operator being certified.
3. The original dense vacuum spectrum/export routines are unchanged and never selected by the BIO branch. Active BIO rejects `RESPONSE-MATRIX-EXPORT`. No full-box dense matrix is allocated.

A future implementation must bind the exact composed candidate, intrinsic/damping assets and native periodic settings; certify the **current** normalized full operator under the existing relative spectral requirement; establish reciprocal/self contribution signs, PME application/adjoint conditions and finite arithmetic error; preserve both RHS and trial state while certifying; and invalidate admission on any coordinate/cell/grid/damping/coefficient change. Installing CHOLMOD or observing PCG convergence supplies none of these joins by itself. The backend cannot return a fabricated certificate/status. Unsupported or inconclusive states must remain fatal.

The compact scalar approximation bounds are separately reviewed mathematical bounds. They are not a full self-consistent force/stress error guarantee. The frozen h/H amplification, induced-response bounds, numerical derivative/virial qualification, approximation budget and actual PME/stability qualification are still prerequisites to execution/sampling. No physical-accuracy or transfer claim follows from this source.

## Controls and costs

The saved 59 source/artificial controls verify exact parent/output pins, all unaffected induction routines, original vacuum setup/gradient bodies after their narrow dispatch additions, native field/frame source identity, unique common hooks, fixed-form statement widths, topology/domain/fatal guards, dynamic storage, scalar native-expression agreement with the compact prototype, unequal-dual two-pass contraction, equal/opposite reactions and the affine virial identity. They use three artificial displacement vectors, with one in the switch region; they do not run Fortran or evaluate a physical box. The exact result records the one-thread wall/CPU cost. There was no compiler, native call, candidate read or fit inspection. Compilation and native qualification remain unperformed.
