# Independent static review

**PASS for the declared saved-only, one-CPU analysis. No remaining static blocker found.** Actual arithmetic replay is performed by the Slurm job, not this review.

## Frozen pins

| Artifact | SHA256 |
| --- | --- |
| LOCALIZE.py | `415027352268cbabb8efdda6e3ea414c420f58e1ff17ba2140b58a9b966f6ee7` |
| CONFIG.json | `7297b790f09b36d64fa6f480856cae003aa23b048dcb25d73e5f25f11b328514` |
| SPEC.md | `588a423de2d30f1179270562810e353eefa2b03e6e4500b3960c348bcee8f2bf` |
| Parent RESULT | `f42cffcf5717a988ed569db925b840f02dfb8bd3b4c04750b77e33bc231c4790` |
| Inventory | `e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957` |
| RUN.sbatch | `febc525488f6ce15079c979b66ec7e32283c948f16c250e888471a6d738fbc86` |

I read the complete source/specification/configuration and launcher, inspected actual atom metadata, parsed source syntax without importing it, and verified all four direct CONFIG artifact pins and the launcher hash.

## Coverage and physical identity

The source requires the exact complete parent with all eight available endpoints and all four matched pairs. It verifies report, original inventory, physical/position and actual MACE receipt/task pins; joins exact force-report source metadata and physical IDs/order; requires finite 256×3 vectors and matched La/Dy coordinates and atom metadata before localization.

The initial draft rejected existing `terminal_completion` atoms. The frozen source explicitly admits that recorded kind alongside `protein_source`. Actual inventories contain terminal H2/H3 and OXT in every cell. Classification remains exactly as declared: selected metal first, every H second, heavy N/CA/C/O/OXT third, remaining protein heavy atoms last. Thus terminal hydrogens are retained with all hydrogens and OXT is backbone heavy. No covalent assignment, distance criterion or state alteration is introduced.

Each of the 256 atom records receives one exhaustive group and one residue prefix; group/residue summaries select from that complete record list. All residues and the separate metal are retained. Top ten atoms are a deterministic display sorted by norm and physical-ID tie break; they do not replace the full tables. Eight endpoint and four pair slots exist before admission, with partial results and failed/unattempted slots retained if a later check fails.

## Arithmetic and interpretation

The source uses saved vacuum errors, defines their common part as the two-metal arithmetic mean and their difference as Dy minus La, with the correct signs. It replays all saved difference components, endpoint and matched maximum/RMS metrics, frozen-direction projections and the declared squared-norm identity.

Group RMS divides by that group's atom count. Squared-norm fractions divide by the corresponding full common or difference squared norm; zero denominators produce null fractions. These are distributions of squared force error, not energy contributions, stiffnesses or unique interaction causes.

The hardcoded/configured absolute `1e-10` and relative `1e-12` tolerances govern this saved floating-point bookkeeping only. Original endpoint/matched metrics and native caveats remain copied, with physical-state, force-accuracy and preference claims false. No new physical agreement gate or fitting path is introduced.

Only standard-library modules are imported. The launcher uses `/usr/bin/python3`, one CPU, `standard_expand`, no GPU and no memory/time override. The source requires Slurm and contains no scientific engine, model load, subprocess, geometry trial, optimization or molecular call path.

No numerical force calculation, scientific import, replay, job, source/config edit or production change was performed by this review.
