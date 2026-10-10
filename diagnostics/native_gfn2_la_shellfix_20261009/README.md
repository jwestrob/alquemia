# Native GFN2 La shell assignment correction — 9 October 2026

Jacob reports bidirectional causal verification: reproducing the erroneous La
shell assignment in standalone xTB recovers tested native ORCA, while fixing the
assignment in ORCA recovers published-parameter xTB. The motivating Ca-minus-La
solvent correction changes by +22.60 kcal/mol. That is one example, not a
universal offset, recalibration rule, or resolution of accommodation errors.

This session verified the two full JSON files and copied the corrected file
byte-for-byte. Exactly `/element/La/{kcn,shpoly,lgam}` differ semantically; shells
remain `[5d,6s,6p]` and every other JSON field is unchanged. No new molecular
calculation ran here. Primary file/hash pins and exact changes are in
[PARAMETER_RECORD.json](PARAMETER_RECORD.json).

## Pinned use in new research calculations

Protocol component: `orca611-native-gfn2-la-shellfix-v1`.
Corrected full-file SHA256:
`a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c`.
Staged artifact:
`workspaces/native_gfn2_la_shellfix_20261009/parameter_source/control.xtb.json`.

For each new native-GFN2 La calculation directory, copy the entire pinned JSON
into that directory as `control.xtb.json`. Do not reconstruct a partial table or
reorder any other fields. Use this exact input block:

```text
%method
  XTBParamFile "control.xtb.json"
  ReadXTBParam true
  WriteXTBParam false
end
```

Apply the corrected parameter identity consistently to both vacuum and solvent
legs of any new La solvent contrast. Include the full parameter hash, executable
hash and correction protocol in input/manifests/cache identities. Confirm actual
parameter import in output. An old cached native-La endpoint cannot satisfy a
corrected task. Keep original and corrected results separately named; do not
silently apply historical calibration bands to the changed scorer. No MACE
weights or declared electronic states change with this repair.

Tested executable, from the supplied campaign's manifest and receipts:
`/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca`,
SHA256 `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83`.
The archived corrected output declares ORCA 6.1.1 RELEASE. Other executables,
versions and metals remain unverified; successful table loading is not their
qualification.

## Impact on current work

The already prepared two-cell Hans8FNR continuation is PBE0-D4 frozen-f DFT,
with explicit La/Dy basis/ECP definitions, not native GFN2. Preserve that manifest
and its four accepted reference endpoints. This correction does not explain its
SCF history or establish LanM discrimination.

Historical native-GFN2-based comparisons may be confounded by this parameter
assignment. Their archived energies remain records of the old Hamiltonian.
Reassess conclusions using matched corrected calculations when that branch is
resumed; do not label every prior physical hypothesis rejected solely from the
uncorrected native-La implementation. The independent accommodation/calibration
problem remains open. PQQ production belongs to another session; this update
changes no production runner/default, starts no rescore and changes no archive.

## Evidence locations

Supplied campaign root:
`/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-10-09_Nikasha_native_parity/`.
Original and corrected tables: `inputs/La_vacuum_{original,corrected}/control.xtb.json`.
Executable pins/import receipts: `bundle_manifest.json`, `results/preflight.json`,
`results/La_vacuum_corrected/receipt.json` and its `endpoint.out`.
