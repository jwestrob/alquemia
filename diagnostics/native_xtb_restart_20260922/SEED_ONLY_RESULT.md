# Seed-only activation attempt: native restart unconfirmed

Job1210179 completed eight native GFN2 calls, but **zero confirmed restarts**.
Every output explicitly reports `INITIAL GUESS: SAD`; none reports XTBRESTART.
The immutable source seeds were present at the correct runtime basename, and
NoAutostart was removed consistently. No GBW was copied, as the initial plan required.
Both seed choices reproduce each destination's prior unseeded energy exactly at
the parser's full recorded precision.

| Destination | Medium | Self-seeded observed Eh | Cross-seeded observed Eh |
|---|---|---:|---:|
| old adaptive | vacuum | −252.116197615034 | −252.116197615034 |
| new template | vacuum | −252.108540025437 | −252.108540025437 |
| old adaptive | ALPB | −252.570227444230 | −252.570227444230 |
| new template | ALPB | −252.570227099694 | −252.570227099694 |

These are observed SAD calculations, not evidence that cross-initialized electronic
branches persist. All qualified restart energies and transfer fields are null,
with status `restart_not_confirmed`. Existing convergence diagnostics, charges,
seed-before/after files and executor receipts remain available in
`workspaces/native_xtb_restart_20260922/run_v2/collection_1210179.json`.

Cost: **26 scheduler seconds ×64 CPU =1664 allocated core-seconds**, zero GPU.
All eight attempts count. `run_v2/COSTS.json` includes inner worker allocation;
peak memory was not available from scheduler accounting. Local setup/reporting
is additional and unmetered. Five real-artifact tests pass after execution, with
zero skips; this does not convert unsuccessful restart activation into a success.

## Documented technical recovery distinction

The initial-guess documentation describes AutoStart checking an existing
same-basename GBW before setting MORead. Native documentation separately says
xtbw takes precedence for a restart. Thus omitting the GBW removes the documented
AutoStart trigger; the SAD result is consistent with this activation explanation.
It does not prove the native restart is unsupported. Explicit MORead with a GBW,
or a matched same-basename GBW for AutoStart, provides a documented activation
route; a native XTBRESTART marker must still be verified.
[Initial-guess documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/initialguess.html),
[native xTB documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb).

The exact matching GBW files are already pinned with the four source seeds. Using
them would revise the initial no-GBW preparation rule without changing any chemical
state or Hamiltonian. Root requested this documented investigation as technical
recovery of the same eight logical tests. No recovery call has been made when
this note is written; preserve this failed activation attempt regardless of any
subsequent recovery. No ordinary-SCF or standalone-backend experiment is added.
