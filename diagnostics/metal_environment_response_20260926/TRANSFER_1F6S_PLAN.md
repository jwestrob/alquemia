# Consumed non-PQQ response: explicit normalized52-atom1F6S state

Declared before new1F6Sresponse energies,27September2026. Third handoff development
core is the previouslyidentified alpha-lactalbumin1F6S chainA strongsite, chosen
for realpreparation quality before inspecting response, not its score. Use the
existing normalized52atom fullymappedamide representation, NOT the original40atom
repair or a splice of their fields. Pin this distinctpreparation identity.

Reuse exact normalizedsource coordinates,1880ff19SBMMcharges, Ca-1/La0 singlets,
11mappedcaps, whole peptideamide chemistry and two complete coordinatingwaters
211/212. Verify actual Lys79-Phe80 andAsp84-Leu85 peptideC-Nbonds; physicalatoms,
sourceboundaries and electron/ECP accounting must match. Noforcefield rebuild.

Same fixedstructuralrule selectsThr86HG1; +10degree orientation aroundCB->OG1,
PCindex1292. Allcorecoords fixed acrossA/B andpairedCa/La; onlyselectedfieldHmoves.
OH0.96A inheritednormalizedstate; O-to-core distance4.64009A makes this a more remote
perturbation thanPQQexamples. Do not enlarge angle or changegroup if response is
small. Preserve state/watercounts; no labelselection or inferredsubstrate.

Sixnew nativeendpoints: Ca/La x A/B/isolated; originalr2SCAN-3c NoAutostart DefGrid3
TightSCF EnGrad, DoEQfalseembedding, nativeLa46ECP. Sameenergyexpression, permetal
B-A anddelta=(LaB-LaA)-(CaB-CaA). No absoluteenergy comparison acrosscores or
crossprotocolbands. Sourceassociation, directaffinity andfunction remain separate;
these calls test electronicresponse, not whether1F6S is physiologicallyLa/Ca.

Known1H4Inumericalrigidfailures remainexplicit. No1F6SFD/grid/rigidqualification
is inferred; smallresponse is unresolved, not a classification. ExactMLweights
are unavailable; noMLcomparison. Mechanicalsolventmodel/optimization remain gated.

6workers x57MPI=342working/344allocatedslots; exclusivefullRealMemory/--mem=0,
standardruntimeheadroom,normalpriority,noPQQdependency. Existingrunner,collector,
actualCodexwake andMondaycutoff. Complete finite3core reference matrix; anynext
molecular experiment must address actualfindings, not expandlibraryblindly.
