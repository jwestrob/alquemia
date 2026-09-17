# AMOEBA framework capability — 2026-09-17

**Three real framework preparations pass. No energy or force was evaluated.**
This establishes protein/water parameter coverage, not a working metal model
or a qualified MACE/AMOEBA hybrid. The production baseline is unchanged.

| Frozen physical system | Framework atoms | Retained waters | Framework charge | Result |
|---|---:|---:|---:|---|
| GGR 1GLG | 4697 | 0 | −6 | Pass |
| Alpha-lactalbumin 1F6S | 1931 | 2 | −7 | Pass |
| Alpha-lactalbumin 6IP9 | 1897 | 3 | −7 | Pass |

Each full physical system also retains its one metal in the source ledger;
that metal is explicitly **unparameterized**. It was not silently dropped
from a scored system. All protein and water coordinates and source identities
remain exact. Both alpha structures retain their four source disulfides.
Template matching identifies HID/HIE, terminal residues and CYX from the actual
graph and existing hydrogens. No atom renaming, hydrogen addition, geometry
generation, protonation change or structural optimization was necessary.

Protocol: `frozen_physical_AMOEBA2018_framework_capability_v1`.
OpenMM 8.5.1, installed `amoeba2018.xml` + `amoeba2018_gk.xml`.
The framework inventory records all force components, local multipoles,
polarizabilities, eight covalent/polarization maps, local-axis atom indices,
GK parameters and exact physical/system mappings. Serialized systems are
preserved. The preparation defaults include dielectric 1/78.3 and the GK
cavity term; these were inventoried, not evaluated or selected as a new score.

## What prevents an immediate scoring pilot

The [published La model](https://pmc.ncbi.nlm.nih.gov/articles/PMC9957963/)
uses pair-specific damping and van der Waals parameters. Table 1 gives default
La damping 0.250 but La–acetamide oxygen values 0.349 (AMOEBA09) or 0.299
(AMOEBA22), and distinct acetate values. AMOEBA22 also changes ligand
polarizabilities. These are not automatically AMOEBA2018 parameters.

The installed OpenMM API supplies one Thole value per atom, and its tagged
reference implementation uses the minimum of the two values. Keeping the
published 0.250 La value therefore cannot yield those larger amide pair values.
There is no pair-override API in the inspected installed header. Setting
different protein atom values would alter other interactions and would not
implement the published model. A frozen distributed QM source also differs
from the paper's formal +3 polarizable La ion; borrowing its parameters alone
would not validate that hybrid.

The next engineering investigation is Tinker's maintained POLPAIR route and
its handling of frozen source sites. No new energy pilot is declared until
the boundary/damping and subtraction definitions in
[AMOEBA_ACCOUNTING.md](AMOEBA_ACCOUNTING.md) are resolved.

## Verification and cost

Six tests passed in 3.659 seconds: actual saved coordinates/charges and
parameter artifacts, corrupted duplicate atom, changed heavy coordinate,
unsupported cofactor, missing real bond and modified cache settings. Tests
performed no additional force-field parameterizations or scientific energies.
Manifest dry-run passed. No scientific executable was mocked.

The three actual preparations took **69.9076944924891 wall seconds** and
**69.65033699600001 process CPU seconds** in total. Peak process RSS was
**480652 KiB**; **zero GPU time** and no scheduler reservation. These timings
exclude imports, source reading and tests. Solver/inference/production cost
is unmeasured. Scientific usefulness and numerical energy credibility remain
untested, separately from the successful framework preparation.

Workspace: `workspaces/mace_omol_20260917/amoeba_capability_v1/`.
Manifest SHA256:
`4bf2e37dad96c9eba1fb9be45bffae2f45aad318a0e8efa8678d5c6e838c9339`.
The authoritative exact hash and all receipts are in
[AMOEBA_CAPABILITY_RESULT.json](AMOEBA_CAPABILITY_RESULT.json).
Source files and verification are under `amoeba_sources_v1/`; the direct PMC
download returned a CAPTCHA page and is explicitly not treated as paper text.
The paper's full text/Table 1 were inspected through the web tool instead.

Recommendation: retain the baseline; pursue the supported polarization backend
without claiming this preparation check improved classification.
