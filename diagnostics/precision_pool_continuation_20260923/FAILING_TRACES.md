# Read-only diagnosis of the five settling failures

The traces favor transient early energy plateaus over a monotonic march toward
a common solution. They do not establish a stable competing electronic state,
and three endpoints per cell do not establish a persistent cycle.

| Source | Cold→pass1 energy,kcal/mol | Cold→pass2 | Interpretation from actual trace |
|---|---:|---:|---|
| C5B120 |+0.631954|+0.000087|Pass1 departs; pass2 returns close to cold.|
| BBL57595.1 |+0.636661|−0.000652|Same excursion/return pattern.|
| P15279 |+0.352082|−0.001468|Same excursion/return pattern.|
| P38539 |+0.000959|+4.111933|Pass2 stops much earlier with a different, still changing density.|
| Q4W6G0 |−0.010536|+0.111100|Small unresolved drift; no demonstrated stable new state.|

For the first three cells, the pass2 atomic-charge vectors also return close
to cold: maximum atomic-charge differences0.00155/0.00215/0.00239e. Pass1
deviations are0.04574/0.03264/0.02454e. Thus the apparent stage1 “improvement”
was not a consistent physical trend.

P38539 is the clearest stopping concern. Pass1 lasts22cycles and agrees with
cold energy to0.000959kcal/mol. Pass2 stops at11cycles after an energy step
of about2.7e−7Eh, while its last five iterations span45.03kcal/mol.
The final MaxDP is0.039706 and RMSDP0.00047481. La Mulliken charge changes
0.81221→0.74028e and the largest atomic-charge change is0.11409e. This is
evidence of a different transient density, not proof of a distinct stationary
electronic solution.

All three runs per cell print `ConvCheckMode Total+1el-Energy` and
`Energy Check signals convergence`. Parent's ORCA6.1§3.5.3.3 verification
states that the native mixer overrides other SCF settings. Therefore the
printed density limits alone do **not** prove that the native solver violated
its own implemented stop condition. The directly demonstrated failure is of
our independent0.1kcal repeated-energy test. Actual density and iteration
histories explain why a single small energy step is weak evidence here.

Recommendation: test a tighter, explicitly verified native **energy** stopping
condition from two already available self-seeds, preserving the native mixer.
This targets the observed plateau rather than adding another blind continuation.
The previous native1e−8Eh six-context test is acknowledged; it did not include
these exact cells. The failed ordinary-SCF restart is not proposed again.
Root separately authorized a new20-call1e−10Eh diagnostic with exact Ca
partners and cold/pass2 seeds. That future execution belongs to its own plan,
not to the completed32-pool results above. No additional call was made for this
trace inspection.

Exact15 traces, charges, energy shifts and output pins:
`workspaces/precision_pool_continuation_20260923/run_v1/FAILING_TRACE_AUDIT_v1.json`.
