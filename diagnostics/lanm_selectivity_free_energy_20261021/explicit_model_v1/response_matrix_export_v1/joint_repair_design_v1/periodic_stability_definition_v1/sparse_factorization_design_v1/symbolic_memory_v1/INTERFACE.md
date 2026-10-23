# Public symbolic fill and memory preflight

**The installed public API suffices; no C shim or installation is needed.** All34 declared tiny controls pass. Six dimension-six symbolic analyses produced structural L counts6,11,11,15,21,15, matching independent graph elimination and public nnz. The numerically indefinite cycle has the same symbolic counts as its positive-valued counterpart. No numeric factor or solve was performed.

## API and storage meaning

`SYMBOLIC.analyze_only(lower, max_numeric_estimate_bytes=..., reserve_bytes=...)` directly accepts the real float64 canonical lower CSC returned by sparse_assembly_v1. It uses lower=True, sym_kind=sym, order=amd, supernodal_mode=simplicial. Public colcount, permutation, nnz, integer width, and is_numeric are available before numeric allocation; the tested integer width is32bits. The installed constructor documentation and cholmod.h agree: simplicial symbolic factors retain Perm and ColCount, without numeric factor arrays.

For n coordinates and m=sum(colcount), packed simplicial real-double factor arrays require `(8+index_bytes)*m + index_bytes*(6*n+5)` bytes, excluding struct headers, alignment, workspace and optional inverse permutation. The latter adds index_bytes*n. This is array payload accounting, not measured peak allocation. The report also exposes sum(colcount²) as a work proxy, expressly not the unavailable Common.fl.

The explicit forecast is two packed-factor payloads (including optional inverse permutation), two additional input CSC copies,128 workspace bytes per coordinate, and a caller-provided reserve. It is an engineering forecast, not a proven upper bound. If it exceeds the caller's budget, the helper returns None for the factor and a rejection report; it never invokes factorize. Before a large numerical release, also reject structural fill beyond the reported index type's capacity. Symbolic analysis itself has allocation costs and is not made safe by a later forecast.

## Actual enforcement and future worker ordering

Use a disposable worker and call `set_worker_address_space_limit(explicit_bytes)` before pair enumeration/assembly. This refuses to relax an inherited soft or hard limit. The tested child bound current mapped virtual space plus8MiB; a1MiB allocation succeeded and a16MiB allocation raised MemoryError. Parent limits were unchanged. RLIMIT_AS bounds virtual address space, not RSS or other processes; it is distinct from a scheduler/cgroup memory grant. This small allocation test is not a qualified CHOLMOD out-of-memory recovery path.

Future order: explicit root memory/time budget → isolated worker hard ceiling → pair query/assembly → symbolic analysis → inspect fill/index capacity and live remaining address-space headroom → forecast rejection/acceptance → separately released numerical factor. Preserve a rejected/failed worker and never reuse a failed partial factor. The pair query precedes sparse_assembly_v1's own entry-count estimate; the worker limit must cover it too.

## Preserved limitations

Common lnz/fl/memory_usage/memory_inuse and raw minor/status are not exposed. No private pointer access, C shim, model mutation, eigenvalue repair, or physical-threshold change occurred. Native C status/minor integration, real-box symbolic fill, peak memory, numerical factorization, PME parity and periodic stability remain unqualified. A symbolic fill pass is neither SPD nor engine qualification.

Cost: controls body0.593511s wall/0.017790s own CPU, with one disposable child0.563924s wall/0.444785s CPU. Full process tree1.570838s wall/1.231904s child CPU, maxRSS89,660KiB. One numerical thread; six symbolic analyses; zero numeric factors, solves, saved molecular matrices or native calls.
