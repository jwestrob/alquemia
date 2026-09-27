# Exact additive classical components: finite real-fixture checks

Declared before energy/force calls27September2026. Use pinnedledger_v1 and actual
1H4I physicalA/B source. This validates classicalimplementation only, not nativeQM
numerics, metalLJ biologicalaccuracy, solvation or wholeproteinrelaxation.

Use maintained OpenMM Reference platform (doubleprecision,singlethread) and
native bonded/CMAP objects; separate forcegroups for retainedbonded,MM-MM LJ,
MM-MM Coulomb(qtilde) andrealQM-MM LJ. Nativeexceptions replacepairs exactly.
No capFFparticle,QMQMclassical,C4,QM-MMclassicalCoulomb,periodicreactionfield or
secondsolventterm. Ledger andmetalpure12-6TIP3Pcandidate remain unchanged.

Finite32configurations: permetal A,B,exactrepeatA,jointrigidA; plus3physical
directions with +/-0.001 and+/-0.0005A displacements. Directions: originalThr159
hydroxyl arc andGlu177 CA boundary direction fromactualsource; metal translation
along actualCa-to-Glu177OE1 unitvector fixedatA andsharedbothmetals. Otherphysical
atoms unchanged. RebuildsourcecapJacobians if electronicmapping is reported;
noQM is rerun or combinedqualification inferred. Keep all failed/missingcells.

Acceptance declared before results (energieskcal/mol, forceskcal/mol/A): repeat
energy<=1e-8; jointrigidenergy<=1e-6 andmaxrotatedforce<=1e-5; eachgroupdirectional
FDandstepconvergence residual<=.001+.0001*abs(analytic). Reference doubleprecision
and simple analytical classicalforms justify tighter checks thanDFTquadrature;
report rawresiduals andlarge-energycancellation limits. Do not relaxgates onoutput.
Count explicitallMM/cross/omittedterms andeveryexception; check metal-independent
MMcomponents equal forCa/La atpairedgeometry. ReportcrossmetalLJ separately.

Execution: firstsingleCa_A task on1sharedCPU, --mem=0 (noartificialRAMcap), measure
actualtime/memory/error. Only if executable/admissible, run remaining31 independent
Reference tasks in32CPUsharedallocation, no repeatedCa_A. No GPU or highlevelcalls.
Tasks are intrinsicallysinglethread; do not reserve344CPUexclusive foronecall.
Normalpriorities, noPQQdependencies, no testqueue. Explicitmanifest/execution
receipts andactualcompletionwake forbothstages. Finitework, nooptimization/MD.
No fixed-timeallocation limit beyond clusterpolicy andMondayprestartshutdownrule.
