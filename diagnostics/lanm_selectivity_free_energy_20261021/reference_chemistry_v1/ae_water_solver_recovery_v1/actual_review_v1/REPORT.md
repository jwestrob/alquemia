# Actual Dy–water solver recovery: job17648

**The one declared recovery converged and supplies a native-admitted energy and
gradient.** Eighty saved checks reproduce both the original native collection and
the additional recovery admission. This is a proactive terminal review, not an
acknowledgement of a delivered watcher event.

The actual deck differs from failed17260 only by nested
`QuadRegionStart 1e-6`. Coordinates, charge+3, unrestricted sextet,73electrons,
164basis functions, scalar-ZORA/PBE0-D4 recipe, MaxIter300 and TolG2e−6 remain
unchanged. The4MPI/16,000MiB grant gives the original MaxCore3000 per rank.
All raw input/output/gradient, wrapper, native helper, release and preparation
pins replay. Original17260 remains failed and has no admitted energy or force.

The native solver reports534cycles including microiterations and40macro rows
(0–39). Its final **macro** orbital-gradient norm is4.543915e−7, below the
unchanged active2e−6 threshold. The final row uses NR. Native normal termination,
finite4-atom gradient and printed-final-energy/EnGrad closure all pass. Density,
DIIS and orbital-rotation figures are retained in the receipt; they are not
substituted for the active TRAH/NR macro criterion. Success with the changed
solver setting does not prove it was the unique cause or that the electronic
root is globally preferred.

The raw output retains39alpha/34beta electrons, `<S²>=8.753692` versus the
sextet ideal8.75, and Mulliken spin4.986453 onDy with total5. These describe this
determinant. They do not qualify its state, spin–orbit physics or force accuracy.
Native final COSX integration changes exchange energy by+0.002117190Eh
(about1.329kcal/mol). Energy–gradient file closure is not an independent
displacement derivative test; that qualification remains absent.

Cost:1131scheduler seconds×4CPUs = **4524allocated CPU-seconds**. The wrapper
measured1130.799263wall seconds; process-CPU utilization is unavailable in the
bound receipts. Saved review used one thread,0.465106wall/0.417205CPU seconds
and zero new native calls. A prior saved-review regex mistakenly matched a
progress ellipsis as dispersion energy; its source/failure are preserved under
`attempt_v1`, and no molecular calculation was involved.

Exact independent receipt: matching workspace `actual_review_v1/REVIEW.json`,
SHA`56d26833e4e7946de61df8865e28bc46707c1b062760eaf12aecd19f6b283c7c`.
Raw native result SHA`d7e070cfd82f41af5a4820d77b1a4677c1b5a05443e89b254f3f15a1ef5d3f99`;
native output SHA`8356acac8a4444446bfda3b021d54a1ada6382df3728afe564e8ae508b8f7808`.
See the adjacent `paired_comparison_v1/REPORT.md` for the matched force comparison.
