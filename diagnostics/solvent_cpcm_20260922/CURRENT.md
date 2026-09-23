# Final checkpoint — matched CPCM route unqualified

Jobs 1209970/1209980 are terminal; the finite observer completed successfully and
exited. Do not relaunch it or recollect the immutable result. All 16 attempts ran:
four Ca CPCM and one Ca vacuum converged; 11 SCFs failed at the printed 499 cycles
under MaxIter500. **Zero complete Ca/La pairs; zero new classifications.**

The actual backend contains the expected CPCM field/cavity but forces ordinary
ORCA SCF. Matched seeded vacuum controls also fail 7/8 despite confirmed MOREAD.
The one complete 4MAE-Ca vacuum agrees with native to +0.0007018304 kcal/mol; this
cannot qualify missing La results. Preserve the released ALPB default.

Final [report](REPORT.md) and [artifact pins](FINAL_ARTIFACTS.json) distinguish
backend implementation from numerical feasibility and untested discriminatory
utility. Scheduler allocation totals 290048 core-seconds; inner executors
289550.752817 core-seconds; 0 GPU. No retries or wider submissions are planned.

Authoritative workspace result:
`workspaces/solvent_cpcm_20260922/matched_vacuum_v1/final_collection.json`.
Its SHA is `60a80f536d13a0745fee1e0cd535e3f2b1e1422631140e465b5e4dbaedc952b2`.
Original inputs, outputs, failed attempts and frozen code remain unchanged.
Vault note:
`/home/jwestrob/jwestrob/obsidian-vault/agent-captures/2026-09-22_Nikasha-matched-CPCM-solvent-pilot.md`.
