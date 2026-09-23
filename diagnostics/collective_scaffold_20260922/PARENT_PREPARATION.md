# Exact scaffold parents: preparation complete

Jacob approved this next round. The [round agreement](../scaffold_restart_round_20260922/PLAN.md)
and [role target resolution](../scaffold_restart_round_20260922/TARGET_MAPPING.md)
are pinned in the prepared manifest. This record covers **preparation only**;
[EXECUTION_PLAN.md](EXECUTION_PLAN.md) separately governs mechanics.

## Result

All eight declared consumed sources have exact ff19SB protein parents and all
24 archived original/Ca-adaptive/La-adaptive starting geometries. No source,
protonation, bond inventory, hydrogen coordinate or water inventory was changed.
The two crystals reuse the existing serialized OpenMM Systems after verifying
exact parameter equivalence and source/cap maps. Their previously documented,
distant OXT completions remain explicit; no new repair was made. The other six
parents matched ff19SB directly. No metal, PQQ or water force-field parameters
were assigned. All eight sources are dry.

| Source | Parent atoms | Mobile atoms / residues | Role target atoms | Caps |
|---|---:|---:|---:|---:|
| 1H4I | 9060 | 1965 / 138 | 3 | 12 |
| 4MAE | 8826 | 2099 / 146 | 5 | 11 |
| Q9Z4J7 canonical | 9431 | 1959 / 136 | 3 | 8 |
| Q88JH5 canonical | 9603 | 1941 / 135 | 3 | 8 |
| A0A3F2YLY8 Ca sample1 | 9299 | 1972 / 138 | 5 | 12 |
| A0A3F2YLY8 Ca sample3 | 9299 | 1972 / 138 | 5 | 15 |
| A0ACD6B9F2 Ca sample4 | 9504 | 2090 / 142 | 5 | 10 |
| A0ACD6B9F2 La sample4 | 9504 | 2080 / 141 | 5 | 10 |

The mobile set is fixed from q0: complete residues having a heavy atom within
8 Å of any canonical-core heavy atom, all core-contributing residues, and one
layer of directly bonded peptide neighbors. Actual topology bonds determine
neighbors; no residue-number inference or recursive shell expansion is used.
The set is identical across the three targets and scored metals for each source.

Role-defined coordination atoms are Glu OE1/OE2, Asn OD1, and actual extra-acid
Asp OD1/OD2 where present. These are pose constraints, not experimentally asserted
bond assignments. Membership has no distance cutoff; actual source distances and
legacy coordination evidence are retained. Catalytic Asp and the cation partner
are not independently fixed as donors.

## Interchange and verification

Authoritative manifest:
`workspaces/collective_scaffold_20260922/parents_v1/manifest.json`.
Each case has a separately readable `CASE.json` and pins:

- `parent`: serialized System, full topology/atom/bond metadata, exact positions
  in Å, stable source IDs, unchanged proton/template inventory and force-field XML.
- `context_maps`, `context_parent_mapping`, `origins`: physical/core/parent mapping,
  complete source and boundary links, endpoint charges and multiplicities.
- `mobile_set`: canonical heavy points, core parent indices, source-distance
  selection, actual peptide bonds and complete mobile/frozen atom indices.
- `donors`, `nonprotein_inventory`: exact role IDs/positions and fixed source
  PQQ/metal plus local-only PQQ H. The metal symbol changes only in paired contexts.
- `targets`: original/Ca-adaptive/La-adaptive full parent coordinates, full q,
  exact archived context pins and reconstructed paired contexts.

The legacy cap field `q0_xyz_A` denotes the original **ideal** link-H position.
Reconstruction retains the serialized cap offset:
`new_cap = original_context_cap + new_ideal_cap - original_ideal_cap`.
All 24 targets reproduce actual archived contexts with maximum discrepancy
1.4210854715202004e−14 Å (declared tolerance 1e−12 Å). All original parent bond
lengths and frozen exterior coordinates are preserved by those target mappings.
No cap becomes an independent physical degree of freedom.

Preparation took **58.859577963128686 wall seconds / 55.462222 process CPU-seconds**
on the existing CPU environment, with **zero molecular energy calls, force calls,
searches or cluster allocations**. OpenMM 8.5.1; ff19SB XML and implementation
hashes are recorded. The force field uses its original NoCutoff terms/exceptions
and no constraints; mechanics separately enforces source bond lengths. Parent
FF energy is for proposal generation only and is not added to electronic scores.

Seven real-artifact tests cover full population, source H/atom preservation,
archived parameter reuse, every target/bond/cap replay, independent mobile-shell
construction, role targets and explicitly corrupted real maps. Two initial test-only assertions assumed JSON key order and identical file paths
for byte-identical Ca/La mapping contents. Both assertions now check the actual
contract; no prepared artifact or molecule changed. Earlier logs remain.

## Commands

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/scaffold_parent_prepare.py \
  --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json \
  --agreement diagnostics/scaffold_restart_round_20260922/PLAN.md \
  --archive diagnostics/scaffold_environment_20260922/RESULT.json \
  --output workspaces/collective_scaffold_20260922/parents_replay_v1

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p test_scaffold_parent_prepare.py -v
```

The existing `parents_v1` is immutable; the replay command uses a new output.
Unsupported source chemistry/templates are explicit unavailable case rows.
There is no preparation blocker on this common8 panel. Mechanical feasibility,
force-field usefulness and electronic discriminatory benefit remain separate,
yet-unmeasured questions for this experiment. Production is unchanged.
