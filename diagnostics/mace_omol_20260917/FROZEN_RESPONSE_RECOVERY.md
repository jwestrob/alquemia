# Native wrapper initialization repair

The first four operator replays, job 1201309, completed but failed the energy
and residual checks. All vacuum mutual-field entries were exactly zero.
Inspection found a concrete wrapper defect: the isolated frontend omitted
`switch('MPOLE')` before invoking `ufield0d` and `dfield0d`. Native `induce0c`
initializes this state before either call. The uninitialized distance cutoff
excluded non-self interactions. The later `egk` calls initialize their own
cutoff and their static permanent energy checks passed.

Retain the first build, inputs, logs, failed checks and cost: 29 allocated
seconds, 1,856 allocated core-seconds, 103.659 reported CPU-seconds. This is an
implementation failure, not evidence that the archived response is unstable.

Build a new isolated wrapper that restores the exact native initialization
and prints/asserts the squared cutoff before evaluating either field. Replay
the same four actual endpoints. No native kernel, physical source, tolerance,
energy expression, or selection rule changes. Use new software and execution
directories; no overwrite or silent repair of the first experiment. Existing
pilot authorization covers this implementation retry.
