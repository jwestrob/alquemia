# ExaA target: independent DFT supports the large metal-dependent motion response

**The target's La-directed accommodation response is reproduced by independent
vacuum DFT.** On the same saved physical path, DFT gives +17.135466 kcal/mol
versus MACE +16.403481. This weakens a MACE-only local-response artifact as the
explanation for this ambiguous ExaA source. It does not establish a biological
metal preference or improve a classifier by itself.

## Actual matched result

| Same-path contribution (kcal/mol) | MACE vacuum | r2SCAN-3c vacuum | DFT minus MACE |
|---|---:|---:|---:|
| Ca work, moved minus origin | −10.379323 | −9.567475 | +0.811848 |
| La work, moved minus origin | −26.782804 | −26.702941 | +0.079862 |
| Differential, Ca work minus La work | +16.403481 | +17.135466 | +0.731985 |

Both metals lower their energy; La lowers it more, increasing E_Ca−E_La.
These are within-method differences, not cross-element absolute-energy fits.
The discrepancy lies below the predeclared1kcal/mol qualitative resolution floor;
that floor is neither an established error bar nor a biological acceptance test.
Unrounded values and exact joined output pins are in TARGET_RESPONSE_RESULT.json.
Slurm4125 used the unchanged comparator and four completed original target cells.
It keeps the eight-cell denominator and all four unjoined control cells explicit;
two control origins have separately finished, giving six completed endpoints overall.

## What this narrows—and what it does not

The preceding fixed-pool packing expansion changed the differential by only about
0.375kcal/mol, while the conditional crystal preparation had a much smaller
response. Together, these findings shift attention toward source-dependent local
strain and the shared chemical/environment representation. They do not prove AF3
is wrong, that every MACE response is reliable, or that dry oxidized-PQQ3− is the
correct experimental state. Both treatments share the same finite vacuum region.

This common origin→adaptive_La path also differs from the operational three-state
pool, which can select a different Ca proposal. No reference band, class call,
solvent correction, biological label or production default changes here.

## Numerical and state limits

Target Ca origin4117 ends in TRAH/NR with residual9.438462e−7 below1e−5;
both target La endpoints also pass their final TRAH gates. Ca proposal4118 ends
on native SOSCF energy convergence with final orbital gradient1.0425e−5 above
1e−5. Its16 applicable checks implement the frozen rule; it has weaker gradient
qualification. Late gap/Hessian/level-shift warnings remain on target origins/La
cells, and no independent state-stability qualification or error bound exists.
The observed agreement is therefore conditional local electronic evidence.
See method_review/TARGET_RESPONSE_REVIEW.md and the raw collected residuals.

## Next decision and actual work

Finish the two already-running control proposals4120/4122. Full collector4123
will compare the target's response with that of the declared sample0 control.
No further chemistry is added before that comparison. A supported source contrast
would favor examining the physical preparation/environment shared by both methods;
a discrepant control would instead require resolving transfer/state dependence.
Neither outcome licenses selecting a convenient fold or moving calibration bands.

Latest origin4117 took5716s on12CPUs (68592allocatedCPU-s; actual TotalCPU18:46:13,
batchMaxRSS12757568K). Target join4125 reports0integer wallseconds/0.087CPU-seconds.
Six completed molecular cells total379092allocatedCPU-s, zeroGPU; running work and
preflight/reporting are separate. These are development-reference costs, not
routine Nikasha latency. All successful outputs/GBWs and original failures remain.
