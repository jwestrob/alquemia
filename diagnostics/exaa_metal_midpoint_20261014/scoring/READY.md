# Finite midpoint scoring runner

Root bound the actual preparation and submitted job4137. No numerical execution, model loading, scientific imports or tests occurred on the login host; actual integration runs on Slurm. See ../CURRENT.md before acting.

The runner consumes the actual two-source midpoint preparation, revalidates the 20 saved original/transplanted cells from job4133, and uses the existing warm CPU MACE transport followed by the qualified standalone xTB/tblite component-v2 and native fixed-charge ALPB path. Maximum new work is eight MACE forwards, 32 SCC calls and 16 component evaluations, with no optimization or DFT. Failed cells remain explicit; complete historical3 and previous5 baselines remain separately available if new7 is incomplete.

Selection uses the existing deterministic order and strict work < −0.1 kcal/mol operational rule, with mathematical minima reported separately. Both metals use the same seven geometries. No class bands are applied. Lower energy alone does not demonstrate improved accuracy.

## Binding and execution

Metadata-only binding, after actual preparation:

```bash
/usr/bin/python3 diagnostics/exaa_metal_midpoint_20261014/scoring/RUN.py prepare \
  --preparation workspaces/exaa_metal_midpoint_20261014/preparation_v1/execution/PREPARATION.json \
  --output workspaces/exaa_metal_midpoint_20261014/scoring_v1
```

The command prints the absolute CONFIG path and SHA256. Root submits its frozen wrapper with exactly three arguments:

Already submitted as4137; this is the exact executable invocation, not a request to repeat it:

```bash
sbatch --parsable --nodelist=node-64-768g-5 \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/exaa_metal_midpoint_20261014/scoring_v1/RUN.sbatch \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/exaa_metal_midpoint_20261014/scoring_v1/CONFIG.json \
  f0edcf081a019a5b5bffb5278050fb1ce6052361837064527a339d213b2c7147 execute
```

`preflight` can replace `execute` for an allocation-bound check without model calls; `execute` includes the same preflight. A new immutable execution directory and workspace lock prevent accidental repeated dispatch. The wrapper requests one node and seven CPU tasks, using the existing lanm_qmmm parent Python and pinned model child environment. It adds no memory or time overrides. The model closes before the seven one-thread solvent workers start.

Prior eight-cell MACE job4130 took118.38 s; the corresponding solvent job4132 took166.64 s. These are observed references, not a runtime guarantee. The prior qualified warm model peaked at8931864 KiB, below the default shared seven-CPU allocation of approximately28 GB. Actual placement is checked inside Slurm. Phase and execution receipts retain actual call counts and pin measured child CPU/RSS records.

## Frozen implementation

- RUN.py SHA256: `575f143e1193fccf52274e8d6f7c1d5147e78ae0e6430213826e33574df53d32`
- RUN.sbatch SHA256: `d4c8e9cb6084f586380a6c6b974fe21208ed845a464c59ad6cd91464e676262c`

Independent static review found no scientific launch blocker. The final narrow reporting change ensures that a collection exception still records failed execution and its traceback. Actual Slurm integration remains the execution test.
