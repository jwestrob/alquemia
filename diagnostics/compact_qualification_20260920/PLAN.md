# Same-state numerical qualification of the compact solvent challenger

Jacob's 2026-09-20 direction: "pursue!" follows the recommendation to resolve
numerical reproducibility and measure the usable scanner path. Existing scientific
autonomy and parallel-work authorization apply. Preserve every previous result
and the promoted baseline; no score coefficient, label, threshold or chemistry
changes in this continuation.

## Question and actual failure diagnosis

Are the successful native GFN2 transfer corrections stable to tighter electronic
convergence and a different solver initialized from the same saved solution?
The original ordinary-SCF checks oscillate by up to several hartrees, with large
density changes, rather than approaching the stopping tolerance. Simply extending
their iteration cap would not answer this efficiently. Native default convergence
uses printed TolE=1e-6Eh and the xTB-specific charge/multipole mixer.

## Fixed numerical work:32 new low-level calls

1. **native_tight_fresh:**24 calls. Same complete-context structures for1H4I,
   Q9Z4J7,4MAE,A0ACD6B9F2,1F6S,1GLG; each Ca/La and vacuum/ALPBWater.
   These cover the original charged/neutral pilot, La-family PQQ contexts and
   the two non-PQQ development proteins. Keep UseXTBMixer=true and request
   `%scf Convergence Tight`, retaining300K smearing and all installed parameters.
   Verify actual printed TolE=1e-8Eh; an ignored keyword is not a passed check.
2. **ordinary_tight_seeded:**8 calls on the original1H4I/Q9Z4J7 pilot contexts.
   UseXTBMixer=false, Convergence Tight, same default125-iteration safeguard.
   Initialize from the actual saved native `.xtbw` atomic charges/multipoles,
   copied into fresh directories. Check the output for actual restart use.
   This tests the same electronic basin with another update algorithm; it is
   not an independent search for a global electronic minimum.

No new MACE/DFT, nuclear motion, water/protonation change, parameter refit, solvent
change or basis substitution. Native and seeded ordinary outcomes are separate;
never select whichever gives the desired classification. If a solver still
fails, retain the failure. No automatic full-panel tightening/rescore is defined.
The separate integrated-scanner agent times four fixed sites with8 fresh MACE
and16 unchanged-primary GFN2 calls; its independent manifest counts that cost.

## Acceptance and interpretation

Keep the prior tolerances: endpoint solvation-transfer difference<=0.10kcal/mol,
Ca-minus-La correction difference<=0.20kcal/mol; charge closure<=5e-4e and
max|atomic charge|<=4e. These thresholds are unchanged after the predictive results.
Compare both new variants directly with the original default primary outputs;
record raw endpoint, transfer and pair differences, charge differences, all SCF
cycles and actual requested/observed tolerances. Verify exact source coordinates,
charge/spin/electron counts and identical parameter exports.

Passing tighter convergence supports numerical adequacy on these six cases.
Passing seeded ordinary SCF supports same-basin solver agreement on two cases.
Neither establishes broad biological validity, unique electronic ground states,
or exact continuum physics. Failed original unseeded checks remain in the record.
The successful composite's calibration and default production workflow stay fixed.

## Resources and execution

Existing ORCA6.1.1 binary and manifest runner;64CPUs/128GiB, eight concurrent
eight-rank tasks, no GPU. Previous primary calls took roughly10–17s each; original
ordinary failures roughly60s per eight-call batch. This predicts minutes of
allocation for this finite32-call experiment; actual costs, including failures
and collection, will be recorded. No project time/compute-budget stop condition.

References verified against actual outputs and official ORCA6.1 manual:
- [Native xTB controls and .xtbw restarts](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb)
- [SCF convergence controls](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html)

Root owns compact_solvation_qualification.py and this directory. Scanner agent
owns its separate new interface/receipts. Existing primary code, sources, bands,
calculation records, active CC job and other agents' edits remain unchanged.
