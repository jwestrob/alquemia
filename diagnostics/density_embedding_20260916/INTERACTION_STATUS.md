# Native interaction diagnostic: running, reference check failed

The two-task Asp303 EDA diagnostic remains **development only**. No complete
decomposition or new affinity score is available at this checkpoint. The
baseline/default remains unchanged.

Initial job 1199985 failed before any SCF because ORCA requires the geometry
before `%Frag`. Exactly the same scientific inputs, reordered, run as 1199986
on 32 allocated CPUs with two 16-rank workers. Adduct and core-fragment SCFs
converged. Asp303 fragment SCFs are unstable in the full adduct ghost basis;
ORCA automatically switched to TRAH. No manual functional, basis, state,
geometry, convergence threshold or time/CPU budget change was made.

## A definite reference limitation

Native EDA generates ghost basis functions from the partner fragment. The
completed core-fragment energies differ from the archived own-basis energies:

| Metal | Native minus archived core energy, kcal/mol | Frozen 0.01 check |
|---|---:|---|
| La | −0.250407924799 | failed |
| Ca | −0.316792779640 | failed |

This is a measured basis-reference discrepancy, not a changed molecular atom
inventory. The 47/7 physical atom mapping and 7/47 ghost centers are retained
separately. Native coordinate serialization changes are below 10⁻⁷ Å after
removing a common sub-microangstrom origin shift. It is not appropriate to
declare this decomposition an exact attribution of the original partition
error even if its remaining SCFs converge. The original gate stays failed.

The SCF instability has no uniquely established cause. Huge intermediate
TRAH model steps are not physical interaction energies and are never scored.
These results neither validate a cheap polarization model nor disprove all
environmental models.

## Automatic completion and recovery

The task-owned completion process uses the existing `affordable_watch.py` to
wait for **this job only**. It launches no new scientific calculation. Its
frozen collector and dependencies passed an immediate preflight against the
actual failed first attempt. At terminal status it writes:

`workspaces/density_embedding_20260916/eda_completion_v1/`

- `collection.json`: native outputs, electronic states, components if available,
  source mappings, all internal SCF/atomic-reference outputs and failed checks.
- `REPORT.md`: explicit completed, failed or unavailable result.
- `cost.json`: actual scheduler allocation including 1199985; no inferred zeros.
- `completion.json`: report/result/accounting hashes.

The same process appends the outcome to the existing vault note. It changes
no shared watcher, production setting, score or prior output. Launcher/PID and
log are `eda_completion_watch_v1.json` and `.log` in the parent workspace.
Consult `completion.json` when it exists; this checkpoint is not a permanent
claim that the job is still running.

Thirteen real-artifact software tests passed, zero skips: seven density/field
tests and six interaction input, reference-accounting and failure-path tests.
**The complete native EDA component parser has not yet been exercised on a
successful decomposition.** Scientific execution and parser validation are
separate. See [the machine record](INTERACTION_STATUS.json) and
[runnable collection/preparation commands](RUNBOOK.md).

Primary interface evidence: [native EDA](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/nocv.html),
[ghost basis notation and counterpoise](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/counterpoise.html),
[native SCF/TRAH behavior](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html).
