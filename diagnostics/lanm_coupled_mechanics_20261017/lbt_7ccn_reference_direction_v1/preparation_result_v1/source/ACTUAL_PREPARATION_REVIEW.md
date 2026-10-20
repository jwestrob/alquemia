# Independent actual preparation review

**Preparation 15845 passes. No blocker was found to root releasing the two declared singlepoints with unchanged scientific binding.** This is an admission of prepared inputs, not a derivative result.

| Artifact | SHA256 |
| --- | --- |
| Actual RESULT | `779412e93cd5945f44bb5564bb6bf8af305021bb0671d235ddc9d5bdf5d8e8d9` |
| MANIFEST | `585fef986ee9e127169b304a926784caa1e7dbfcbef86d7cac4a2aa23d296e9a` |
| DIRECTION | `f20be9142797fbcbe6f6a468380dc617880f534de1bbf2d64644fc6d9a7165c4` |
| Executed source | `cd4396c02cf50905307acbf90772d062ae4a6a12aaf9914e01562167dca1ab32` |
| Configuration | `cec2f6e393f191c656f7dd6fd57441895c14f2bb5fd6a75014608f2b5428f43f` |
| Shared SP deck | `e99306dc36ed29ae86742b380a03442c352642d96b6cef6fd864a8eef5b6e8ed` |
| Minus XYZ | `8e66d66cf3c57e9e8a87e74ee9b7a3fc025fa26e6f2fef9103677b57bfb28a7a` |
| Plus XYZ | `875b2c9a612e647209c30802ebc1b81566c35ff7f6fcac62032024c7c1e54852` |

The result is `prepared`, on node-64-768g-1, with both minus/plus tasks selected and no requested singleton. All molecular, runner, central, other-metal and optimizer counters are zero. Both native rows remain explicitly unavailable/not attempted, and the comparison remains unavailable. Reported application time is 5.491746 seconds; root records the authoritative scheduler cost.

Sixteen selected artifact pins were independently hash-verified. RESULT, manifest and CONFIG scientific bindings agree. Each task contains the same 256 physical IDs and chemical metadata as the central endpoint, with fresh displaced physical/position records that agree literally. Whole charge −2, physical/effective singlet, 1030 explicit electrons, La ECP46 and the target index are unchanged. Tasks explicitly declare `single_point`, `engrad_path=null` and ±0.01 Å direction signs. Both decks exactly match the original central deck after removal of its sole EnGrad token. The manifest retains one worker and 12 MPI ranks for subsequent native execution.

The saved direction binds the exact central native result, collection, force report, gradient and state. It reports discrepancy norm 57.4261573845 and PBE0/MACE force projections −66.4531999104/−9.02704252585 kcal/mol/Å. These are copied Slurm outputs, not new calculations by this review. The declared preparation executed the finite-vector, norm/sign and displacement checks; they were not replayed on login.

Central SCF evidence is retained exactly: native convergence after 18 cycles, passing orbital-gradient/rotation flags, failing printed energy/density/DIIS flags, and `all_printed_residuals_pass=false`. Physical-state and force-accuracy qualification remain false. Only the complete admitted native pair can determine the frozen one-direction, one-step consistency verdict; no full-gradient, physical-Hamiltonian, solvent or preference claim follows from preparation success.

Review used stdlib JSON, literal metadata equality, source text and hashes only. No scientific imports, vector arithmetic, jobs or source/config changes were made.
