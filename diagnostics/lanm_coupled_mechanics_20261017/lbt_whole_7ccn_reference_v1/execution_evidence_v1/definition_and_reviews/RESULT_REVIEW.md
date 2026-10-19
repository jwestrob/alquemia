# Actual7CCN preflight14793 result review

**PASS — saved-input preparation only.** The actual one-CPU Slurm job completed with exit0:0 in3 scheduler seconds (application3.5622883720789105s), with all molecular and native-runner counters zero. No native output or gradient file was produced. This review used source/text/JSON metadata and hashes only; it performed no numerical replay, imports or jobs.

## Reviewed pins

- RESULT: `64db262eadb942246cff51a8bccca5b5b52520bf802cf93ff6f241610691be39`
- MANIFEST: `b74e9aed49eb9078274cbfd9418c8f684c083bdcd05008679584f7ad4e8a71a0`
- DIRECTIONS: `5c773e988ca0cb0750b652739f9fbdf6135bf468190f9d726b02ec7c416cca20`
- Executed source: `c0cca6dd4f44619fc2e3201f024b7913f6b7e9d8e8c4fca694fe17bd942d0ecf`
- Staged config: `4768bc5825b8012ae2eae9985b9db9d73d859c63ce0f686f0d55df42132303ec`
- SPEC: `a8cb688c8d5f8f1f109a6c1436da24e68410cb92528f503d123cc83650307447`
- Frozen inventory: `4203d1aabfc4f312b8627469e0b94450927e06b73d38294b3db6a5a6739cb0e2`

428 unique result/manifest/config-linked file pins verified. All four task IDs and exact physical/NPY/ALPB parents match the declared geometry_20/cartesian_03 inventory. The source selection remains exact12968 7CCN-La start→80-step terminal, bound/admitted but nonconverged; Lu conditioning and terminal-chemistry caveats are retained.

## Actual admitted preparation

All four native tasks contain256 atoms, target255, charge−2, no caps/external charges, and analytic-gradient output paths. La physical multiplicity1 / ECP46 gives1076 all-electron and1030 explicit electrons; Dy physical multiplicity6 / ECP55 gives1085 all-electron and1030 explicit electrons. All four effective DFT states are singlets,515 alpha/515 beta. These values are copied from actual preflight state records, not recalculated here.

All four decks were checked byte-for-byte against the unchanged pinned PBE0-D4/DefGrid3/VeryTightSCF/EnGrad/PModel whole-vacuum template and corresponding basis/auxiliary bytes. No native SCF run has occurred; the prepared12MPI resource declaration is for later singleton executions, not a claim that this one-CPU metadata job used12 ranks.

The executed frozen input path passed all four raw ALPB scalar/gradient/charge/coordinate/state replays and exact whole-MACE model/task/cache/native-energy checks. The two available La raw/composite gradient equality gates precede preflight_passed and passed in the actual Slurm environment. Both Dy MACE requests remain energy-only: raw forces and composite gradients are absent, not filled. Saved source, physical ID, paired-coordinate and topology gates remain unchanged. One full-deformation chord is frozen from the prior coordinate pair, not a new reference outcome.

## Consequence and limits

No blocker remains for root's declared single before-La scout at12 actual MPI slots under unchanged native settings. This does not admit any new DFT energy/force, qualify physical spin/electron binding or force accuracy, or clear the remaining-three continuation gate. The latter still requires an actually admitted scout and separately pinned measured-cost decision. No source/config changes or new submissions were made by this review.
