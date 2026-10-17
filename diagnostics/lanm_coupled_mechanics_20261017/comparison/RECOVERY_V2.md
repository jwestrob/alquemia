# Saved comparison stage-dispatch correction

The frozen whole-LBT paired-pool CLI records stage `run`; the Mex paired-pool
CLI records stage `pool`. The original comparator 062f1167 incorrectly required
`pool` for both and would leave the completed LBT menus unavailable. Root confirmed
that no numerical comparison had been run with that source.

`COMPARE_V2.py` changes exactly one condition: accept `run` for LBT and `pool`
for Mex. Physics, thresholds, source/state admission, arithmetic, costs and all
other behavior are unchanged. Original COMPARE.py and immutable comparison_v1
staging remain intact. Bind a new config with COMPARE_V2.py from a new copy of the
staged inputs, then run the same one-CPU Slurm command with the V2 filename.
