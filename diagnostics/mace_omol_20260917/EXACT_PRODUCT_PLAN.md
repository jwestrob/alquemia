# Exact atom batching for larger intact models

Declared 2026-09-17 before implementing or running this adapter. The completed
9088-atom PQQ inference peaks at 21.72 GB even after edge batching and allocator
recovery. The canonical inventory includes sources above 10,000 atoms. The
actual earlier OOM traceback points to symmetric contraction inside the native
EquivariantProductBasisBlock, after full neighbor aggregation.

Batch that block over atoms while invoking its original forward on every slice.
The installed MACE 0.3.16 product block contracts and transforms each atom's
features separately; it has no reduction between atoms. Concatenate outputs in
original atom order before the next full interaction. Preserve the qualified
1024-edge adapter, graph, weights, buffers, float64 precision, global charge/spin,
and all physical inputs. No force or gradient support is introduced.

New execution adapter: `omol_exact_edge_and_product_batches_v1`. Scientific
energy and intact descriptor protocols remain unchanged. This is a technical
qualification, independent of classifications or score direction.

Two finite stages use the existing native reference manifests and executor:

1. Eight calls: the four real 1H4I/4MAE core endpoints, with product batches of
   32 and 1024 atoms. The smaller batch exercises multiple slices on real cores.
   Require every endpoint and Ca-minus-La contrast to match the archived native
   energy within 0.01 kcal/mol, with unchanged state and complete atom/edge counts.
2. Only after stage 1 passes, fourteen ALPHA_1F6S bound/detached, repeat, rotation
   and farther-detachment calls with product batches of 1024 atoms. Require the
   same 31 numerical checks and all 24 direct comparisons against the completed
   native CPU reference, each within the existing 0.01 kcal/mol tolerance.

Twenty-two declared new MACE energy-only calls; zero DFT, optimization, training,
solvent or force calls. Reuse all native reference calculations. Record actual
memory, timing and failed attempts. One A5000, 16 CPUs, 64474 MiB host memory,
with the now working expandable-segments allocator. Keep finite manifests and
scheduler QOS; there is no project compute-time budget.

Do not use this adapter for further scientific scoring unless its actual
qualification passes. Do not change existing benchmark manifests or replace
their successful outputs. Training-domain/charge limitations of the model are
separate from this memory implementation and remain explicit.
