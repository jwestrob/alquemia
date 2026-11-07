# CPU continuation of the declared periodic checks

Declared 7 November 2026, before these controls are evaluated. Job 52877 proved
that the unchanged full periodic L calculation fits the default 24-CPU memory
allocation. The H200 job remains pending. Complete the already specified
translation and remote-water checks for start 0 on CPU while waiting.

Reuse the two completed endpoint predictions from job 52877. Make exactly eight
new evaluations: translated A/B and three individually removed waters, A/B each.
Use precisely the same model, precision, PME, coordinates and geometric water
selection as the planned GPU controls. No new model or threshold is introduced.

Allocate 24 CPUs with default scheduler memory and time, `standard_expand`.
An internal 330-second timeout plus five-second termination grace bounds the
additional CPU cost to approximately 8,040 allocated seconds. Together with
4,953 already used, this remains within the four-CPU-hour preparation/static
CPU envelope. Preserve partial results if it cannot finish. Do not retry it.

The pending GPU calculation is still needed to measure actual H200 throughput
and verify CPU/GPU agreement. Adding the explicit CPU-controls command-line
branch to `periodic.py` does not change that GPU execution path. Each job saves
its actual source at startup.
