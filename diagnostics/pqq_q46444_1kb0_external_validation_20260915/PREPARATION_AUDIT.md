# Q46444 / 1KB0 preparation audit

Prepared: 2026-09-15  
Status: **ready for ORCA; not submitted**

The pinned fixed_core_v3 preparer completed without repair or fallback. The
release record is:

`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_q46444_1kb0_external_validation_20260915/prepared/external_preparation.json`

Observed release invariants:

- One 1KB0 target and exactly two tasks, La and Ca.
- Typed CN = 7 (6 O, 1 N); fixed-core total charges La -1 and Ca -2.
- All 5,138 retained source heavy atoms survived normalization and protonation
  with zero displacement at serialized precision.
- PDBFixer added no missing heavy atom and no missing terminal heavy atom.
- All 1,016 source waters and six noncore heterogens were removed. The
  normalization ledger explicitly contains excluded `A:TFB1810`.
- No source or synthetic water, replacement ligand, point charge, or geometry
  relaxation enters the QM pair.
- La/Ca nonmetal coordinates are byte-identical.
- Both inputs contain exactly
  `! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3`; runtime `%pal` is supplied
  only by the pinned manifested ORCA runner.
- No ORCA output exists in the preparation tree.

Key artifact SHA-256 values at release:

- external preparation:
  `3c8a944f4462f63670bb2c9069df3cc45a9bcf00c29e73e8b85c743c159155fc`
- carve manifest:
  `9ecc4229ffa077f8756a0a9f5a99cc3329a2d47f0416e2beb5e4989d28460114`
- La input / XYZ:
  `308a89f3a37d94c0262ae7871ffe6a44fb7b1f38df5e83b59dd8aea038b93411` /
  `df0b8083d182f3c71eedf0ad4a518846daaee0e4c83f85d93169413be1303750`
- Ca input / XYZ:
  `7bca3336d3aa5fd0ee2935e7c69a512848683207ee3d8d5154bee99e771ed4c1` /
  `c80259cbaf814e3262610eeb71363fdae7eebd2481eb1a3eb7e4cd9567453252`

This audit establishes input readiness only. The preregistered biological
verdict remains unopened because neither ORCA arm has been run.

