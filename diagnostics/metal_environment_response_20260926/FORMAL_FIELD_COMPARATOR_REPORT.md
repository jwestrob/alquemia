# Formal-metal field comparator: incomplete PQQ response

Retrospective comparison on already-inspected actualA/B fields; zero newmolecular
calls. Modelpredicted delta=(qLa-qCa)*[phi_B(rmetal)-phi_A(rmetal)], with formal
chargedifference+1e andsame metalcoordinates. No ligandcharges are invented; their
metal-dependent response is deliberately omitted by this comparator.

| Core | Pointmetal,kcal/mol | Native,kcal/mol | Native minus pointmetal |
|---|---:|---:|---:|
|1H4I|0.186901735|0.352795817|0.165894082|
|4MAE|0.129209904|0.254914852|0.125704948|
|1F6S|0.308733177|0.269010220|-0.039722957|

For these PQQ motions the baremetal approximation accounts for about half the
computed differentialresponse. For1F6S the absolute difference is0.0397kcal/mol,
comparable to the original intended referenceuncertainty0.05; this is not a proven
accuracy bound or equivalence result. All are consumeddevelopment geometries.

This supports retaining the full electronicregion when testing environmentresponse
on PQQ instead of claiming thatmetalpotentialalone is complete. Differences may
include ligandcharge redistribution, polarization and representation effects;
no unique attribution is established. There is no classifiertraining, fittedscale,
biologicalsuccess claim, or modelpromotion. Native rigidqualification remainsfailed.

Exactfield/input/output hashes and values are in FORMAL_FIELD_COMPARATOR.json;
atomicunits implementation is scripts/metal_environment_fields.py:potential_field.
Eachrecord pins theactualcollection, which pins themanifest andXYZ/pointcharges.
ReadCa_A/Ca_B (same metalcoordinate), convertAngstromtobohr, evaluatephi atfirst
XYZatom, subtractB-A, multiplyoncebyHA_TO_KCAL andformalchargedifference1.
