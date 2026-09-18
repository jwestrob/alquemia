# Matched PQQ utility benchmark — running

Approved goal: GOAL.md. No new model, calibration, physical preparation or
unlabelled application test. Baseline/default unchanged.

- The two parser/algebra regressions pass on actual DFT and masked-MACE results.
- Prepared-input interface exposed a chain-selector bug: canonical AF inputs
  have protein A, La B, PQQ C. The old single-case guard counted metals only in
  A. SOURCE_SELECTOR_FIX.patch now checks the explicitly selected source metal
  and still rejects multiple metals in the selected/protein chains. All 25
  exact source/coordinate/charge audits pass; no physical atom changed.
- Isolated interface_source_v2 copies the prior frozen implementation, changing
  only that guard. Scientific adapters, model, states and calibration unchanged.
- Preflight_v2 prepared a valid two-task manifest, zero inference, and audited
  all 25 cases in 126.1672248840332 s. This is development overhead.
- Regression testing also found the older crystal records use
  `normalized_selector` instead of separate chain/residue fields. Source V3 and
  the current prepared-input script accept both recorded schemas explicitly.
  The running V2 snapshot is preserved; all of its 25 canonical cases use the
  first schema and already passed. No scientific model change or new inference.
  All three selector regressions pass in 12.595 s, including the real crystal
  and an explicitly corrupted copy of the canonical selected-residue record.
  Final V4 additionally preserves the generic `insertion_code` field alias;
  four selector regressions now pass in14.902s. V3/crystal-only test history is
  retained. TimingV2 remains immutable and unaffected by this compatibility fix.
  The four existing prepared-interface regressions also pass in28.420s, none
  skipped; these reuse actual PQQ/GGR outputs and check invalid1KB0 rejection.
- DFT attempt 1201561 failed before SCF: `ntasks=1,cpus-per-task=32` exposed
  one MPI slot. Four failed launches/two cases are retained in timing_v1.
  Recovery uses 32 scheduler tasks, one CPU each, two 16-rank endpoints.
- **DFT 1201562** and **MACE 1201566** are submitted against timing_v2. Inspect
  live scheduler state; never resubmit a running campaign or overwrite results.

Manifest timing_v2 SHA:
7a57fea328c1e11668ec23dab4587acd8bf38c42069c206e693d1126c908df22.
50 fresh endpoint calls per method; no detached MACE tasks. Separate allocations
on node-128-512g-8gpu-1: DFT 32 CPU/no GPU; MACE 16 CPU/one A5000, 64474 MiB each.
No project runtime cutoff. Complete paired timing/classification result pending.

**Additional technical check, same frozen scientific model:** parser-cache
job1201609 runs50freshMACE endpoints, no newDFT, against timing_cached_v1 SHA
794e45b9475d87d37441cead94c26e335e60e31ebd3dde219f5554326c128953.
One additionalA5000/16CPU/64474MiB on the same host. SourceV5 alters only parsed
XYZ reuse inside existing verification operations; all54real endpoint coordinate
arrays are bitwise equal and four cache tests pass18.552s. Matched report-only
job1201589 measured62.06s→44.67s with identical scientific output,107core-s,noGPU.
Full optimized workflow speed remains unproven until1201609 completes. Original
timings remain authoritative and separately reported. See REPORT_CACHE.md and
CACHED_RUNTIME_CHECK.md; these are engineering checks, no new biological/model
comparison, scientific protocol, calibration or default change.

First completed pair: DFT239.74962226301432s, MACE161.11215551942587s. The MACE
score exactly matches51.455254788976355. DFT differs by−0.00024073708349kcal/mol,
within the declared0.01reproduction tolerance but across its exact calibration
extremum, giving a literal inconclusive call. Keep that numerical boundary
crossing visible; no adjusted threshold or biological regression claim.

The benchmark driver now returns nonzero for an incomplete execution. The
already pinned timing_v2 driver is unchanged and its per-case success/failure
receipts, rather than Slurm COMPLETED alone, determine scientific completion.

Original accuracy evidence remains 25/25 canonical and 2/2 supported crystal
transfers for masked MACE; 1KB0 unsupported. DFT has 25/25 plus 3/3. These are
consumed reference tests, not fresh independent biological validation.
