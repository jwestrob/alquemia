# Partition preparation 12950: inactive MM parameter comparison

Job12950 constructed all four source/region fixtures and failed their final
La/Dy common-MM equality check. It ran zero energy or molecular force calls.
The original source, config, result and all XML remain immutable; exact pins are
in `RECOVERY_V2_EVIDENCE.json`.

Text/XML inspection found the same single difference in every fixture: group1
`MM_LJ`, physical particle255 (`metal`), has La sigma
`.30611279955302056` nm versus Dy `.2866912074975612` nm. Both particle charge
and epsilon serialize as exactly `0`. Groups0 and2, every other group1 field,
and its exception records are identical. The builder retains target sigma while
zeroing QM-particle epsilon in this MM-only force. Therefore this sigma cannot
contribute to its LJ energy or force. Target–MM LJ is separately represented in
group3 and remains metal-dependent.

`PREPARE_V2.py` changes the comparison only. Before normalizing that one sigma
in a temporary XML comparison copy, it requires the target's QM membership,
exact zero charge/epsilon, no parameter offsets, no target exceptions, and the
actual group3 QM/MM interaction membership and target exclusion checks. Every
other common-force XML field and exception must still agree exactly. Neither
the OpenMM systems nor their serialized output is modified. A further gate
requires every new actual La/Dy MM XML hash to equal its corresponding failed
12950 artifact, including the active metal-dependent group3 parameters.

The molecular protocol, region selection, native parameters, boundary/charge
policy, geometric tolerances and four-fixture/eight-state denominator are
unchanged. `admission_revision` records the comparison repair separately.
This is construction validation; no MM or hybrid force qualification is claimed.

Source AST parsing passed. No imports, numerical tests or molecular evaluations
were performed on login. Root runs the same four actual preparations in a fresh
Slurm workspace using `CONFIG_V2.json` and the unchanged `RUN.sbatch` interface:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python PREPARE_V2.py prepare --config CONFIG_V2.json --output FRESH
```
