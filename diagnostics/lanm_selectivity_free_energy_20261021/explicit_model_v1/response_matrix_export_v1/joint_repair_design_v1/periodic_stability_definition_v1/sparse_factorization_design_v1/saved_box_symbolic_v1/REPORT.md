# Saved BIO graph: symbolic resource result

**The single symbolic check completed; its declared numerical-memory forecast rejects proceeding.** No numerical factorization or solve occurred. This is a resource-feasibility finding, not instability or a model-accuracy result.

| Quantity | Actual saved result |
|---|---:|
| Active sites / coordinates |11958 /35874 |
| Unique real-space pairs |1815015 |
| Lower-CSC entries |16371009 |
| Input CSC array storage |196595608 bytes (0.183GiB) |
| Symbolic L entries |399831681 |
| Structural fill/input entry ratio |24.4232 |
| Fraction of fully dense lower triangle |0.6213 |
| Packed factor-array payload |4798841168 bytes (4.4693GiB) |
| Numerical-stage planning forecast |10264187872 bytes (9.5593GiB) |
| Declared planning budget |4294967296 bytes (4GiB) |
| Index-width check |int32 capacity passes |
| Backend numeric state |false; rcond=-1 |

The packed factor-array payload alone exceeds the4GiB planning budget; the broader forecast also rejects. The forecast includes explicit extra-factor/input/workspace/reserve allowances and remains an engineering estimate, not a proven peak bound. It does not prove that an8GiB process limit would necessarily fail a different implementation, and it does not release a larger budget or numerical attempt. The sum-of-squared-column-count work proxy is5942759498717; this is not measured CHOLMOD Common.fl or a numerical timing prediction.

## Scope and provenance

Exactly one La lower-CSC assembly and one AMD/simplicial symbolic analysis were executed. The frozen La/Dy coordinates, topology and nonmetal types are identical; all9 entries per pair are retained structurally, so the same graph covers both metal identities without a second assembly/analysis. This does not make their numerical matrices identical.

Literal saved BIO alpha/thole/POLPAIR parameters were parsed; La/Dy metal–O gamma remains0.28782327521600637/0.28375642691227121. The earlier prospective gamma family was not substituted. H=0, shift=0; Ewald0.4/A and9A cutoff are the declared diagnostic assembly settings. All active alpha values are positive, all mutual scaling factors are1. No candidate or periodic custom-response Hamiltonian was admitted. The unshifted site bound remains negative; neither that bound nor this symbolic result determines full-PME stability.

## Actual resource limits and cost

The worker applied hard+soft8GiB address-space and180s CPU limits before numerical imports/assembly. The parent applied180s wall timeout; no timeout occurred. The full pair query was included inside these external limits despite preceding the assembly estimator. Actual outer wall: 14.163174s; child CPU: 13.625333s; peakRSS579696KiB (566.109MiB). One numerical thread. The assembly storage forecast was817427872 bytes, within4GiB. Recorded virtual mappings before and after symbolic analysis were644739072 bytes; these snapshots are not peak mappings.

Empty stdout/stderr, one symbolic attempt/completion, zero numeric factors/solves/native/PME calls, no retries. Full input/source/config/release pins, array hashes, symbolic column counts/permutation, partial-state protocol and actual resource readbacks are retained. Raw C common memory/minor/status, numeric factor performance, complete native operator parity, physical stability and MD remain unqualified.
