# One memory-only recovery of the unfinished compatibility pilot

## Observed failure

21618 ended OUT_OF_MEMORY after126scheduler seconds (126.611890s worker wall),
24CPUs,96,000MiB granted; batch peak98,304,464KiB. Slurm records two OOM kills.
The launcher exited9. The requested isolated charged-fragment SCF and ALMO
R/F/P/U chain were not completed. The auxiliary neutral-atom guess energies
printed during initialization are not the requested M3+ fragment result.
All original input, raw output, failed parser result and accounting are preserved.

This failed resource allocation does not establish that CP2K is electronically
incompatible. It produced no full-system comparison or CT result to retune.

## Narrow recovery and authorization

Correct the resource request once to **256GiB (262144MiB)** on the same shared
standard partition with standard_expand,24CPUs,12MPI×2OMP, no exclusivity.
Actual QoS permits24CPUs/job and100T memory/job; the prior96,000MiB was the
default request, not the QoS memory maximum. Current AGENTS.md explicitly
heads its older OMol/default-memory notes as historical, while the current
23October direction authorizes use of standard_expand resources. This recovery
uses the active bounded diagnostic's existing compute authorization.

The original RUN_SPEC's no-retry/default-memory execution choice is amended
only for this initialization OOM. There is no automatic retry loop, numerical
reoptimization, alternative electronic protocol or scope expansion. The one
scientific La-aquo compatibility pilot is still unfinished. Its physical and
numerical admission gates are unchanged.

No input line changes: exact four geometry decks, basis/ECP assets, native
binary/container, fragment definitions, PBE0/GAPW/cell/grid, convergence
controls, full-range exchange, D4 bookkeeping and parser are all the original
pinned artifacts. New config/runner/task/output paths preserve the failed run.
The new runner additionally checks the actual262144MiB memory grant.

The original two-hour native budget is **not increased**. Conservatively charge
ceil(126.611890)=127seconds for the first worker and limit the new native call
to7073seconds. For simplicity this shorter7073second cap also applies to each
remaining system if the pilot is admitted. No walltime override is requested.
Memory256GiB is a finite allowance, not a proved peak-memory estimate; parent
and current fragment grid environments coexist. No further automatic memory
escalation or restart is authorized by this recovery specification.

## Release and outcome

Root releases only aquo8_central_La after source/config review. On actual pilot
admission, the other three exact systems run concurrently under the same new
resource request and unchanged electronic protocol. No completed geometry is
repeated. All state/energy/decomposition checks remain required for every cell.

A repeated resource failure or a failed scientific admission yields an explicit
unavailable/incompatible/inconclusive report as appropriate; it does not yield
a CT estimate. No fitting or additional molecular campaign follows.
