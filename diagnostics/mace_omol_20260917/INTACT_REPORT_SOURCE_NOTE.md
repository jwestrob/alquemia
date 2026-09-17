# Report-source provenance clarification

The comparison ran from a live script. While it was running, the reporting code
received an operation-local checksum cache; the reporter copied its source at
output time. Consequently the `implementation` record identifies that writer-time
source snapshot, not an exact snapshot taken when the process began. Preserve
the original result and this distinction; no scientific endpoint was recomputed
or changed by the reporting edit.

An independent invocation of `verified()` using the updated source reproduced
the complete saved comparison, including all24 passing checks and zero energy
discrepancy. Its measured receipt is `../edge_equivalence_cached.resources.txt`;
the original invocation is `../edge_equivalence.resources.txt`. The subsequent
intact benchmark snapshots its entire verifier and independently recomputes the
comparison while validating its preparation. Both underlying native and batched
inference campaigns ran from immutable implementation snapshots throughout.

Future report executions should keep their source unchanged while running, or
run directly from a frozen implementation directory. This note does not claim
that writer-time source metadata is an exact historical runtime pin.
