# Crystal-only terminal completion recovery ready

V2 job4088 is preserved: both sources completed protonation but failed the old
zero-terminal-completion assertion. Saved manifests report internal0/terminal1;
actual PDB text identifies A:TRP616/OXT. Numerical distance and invariance checks
are deliberately reserved for Slurm.

Frozen manifest: `workspaces/exaa_source_packing_20261014/crystal_v3/manifest.json`
SHA256 `360455889a1063f1ea577f08ae0f71fbf74215220be95ccc1d9e5d652e96894a`.

```bash
sbatch --nodelist=node-64-768g-5 \
  --export=ALL,MANIFEST="$PWD/workspaces/exaa_source_packing_20261014/crystal_v3/manifest.json" \
  --output="$PWD/workspaces/exaa_source_packing_20261014/crystal_v3/slurm-%j.out" \
  --error="$PWD/workspaces/exaa_source_packing_20261014/crystal_v3/slurm-%j.err" \
  diagnostics/exaa_source_packing_20261014/crystal/v3/RUN.sbatch
```

One CPU, no memory/time overrides. Root owns submission and completion wake.
No scientific import/test/preparation has run for v3 yet. The metadata-only
freeze reads and pins both real4088 outputs and truthful protonation manifests.
The job reuses all generated H and OXT coordinates without rerunning protonation.

Actual result will be `crystal_v3/execution/FINAL_PREPARATION.json`. Each admitted
row retains existing maps/state/source/context/spectator fields and adds a pinned
`terminal_completion` receipt proving exact sole OXT addition, deposited-coordinate
preservation and exclusion/distance from the final context. No PQQ repair,
threshold change, new heavy atom or redox assignment is permitted.
