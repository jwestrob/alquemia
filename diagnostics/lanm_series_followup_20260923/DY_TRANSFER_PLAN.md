# Frozen structural transfer: Dy-conditioned Hans LanM

## Authorization and question

Jacob authorized contained overnight experiments on September23. Root selected
this next test after the completed original Hans/Mex La/Dy vector. Does its
relativeLa signal survive an actual Dy-conditioned Hans structure, or depend on
using the La-conditioned crystal? This is structural robustness/mechanism within
one consumed biological comparison, not independent site-affinity validation.

Use only deposited **8FNR author chainA, Dy201/202/203**, in EF1/2/3 order.
SourceSHA256:`b7da288458075fe269820a8dfdd6bc371bfadd9c9f07a59e03a096cc46677c88`.
ChainA contains all three sites; no replacement chain or score-selected site.
EF4 is not part of the established three-primary-site vector.

Read-only inventory before scoring shows9O contacts within3.3Å per selected
8FNR site, versus10 for8DQ2. E9Glu42/66/91 has onlyOE1 in that shell; itsOE2 is
outside. No deposited water is inside3.3Å at the three primary sites. Whole
sidechain membership remains determined by the existing carver, not contact
count tuning. Outer waters are not newly added to this frozen first-shell model.

## Fixed preparation and energy policy

Use the same existing occupancy-aware source-conformer selection, standard
pH5protonation, no missing residues,3.3Ågeneric selection and peptide-amide-v3
repair as the first experiment. Preserve every deposited heavy coordinate.
Unsupported/missing source chemistry remains explicit. Extend only the private
source selector to recognize depositedDy; production selection is unchanged.

At each site, La andDy use identical coordinates, caps, ligands, protonation and
water inventory. Native float64 OMOL uses physical multiplicities1/6, respectively.
Native GFN2 uses its explicitly recorded f-in-core effective singlet, unchanged
300K/native mixer/MaxIter500, initialNoAutostart plus exactly one same-cell,
same-medium GBW+xtbw continuation with positiveXTBRESTART confirmation. No DFT,
optimization, new electronic model, extra initial states or favorable root choice.

**Maximum6freshMACE +24GFN calls:** threeHans sites×two metals;12GFNcells×initial
and one continuation. Reuse all six immutable Mex endpoints from the completed
first experiment. NativeLa/Dy capability was already demonstrated there; no extra
capability calculation is needed. Retain partial failures and every declared site.
Use the same oneH200/32CPU/200000MiB warmMACE and64CPU/128GiB native8-rank×8worker
allocation policy. Prior complete12MACE/48GFN cost6480core-s/21GPU-s; this is a
smaller fixed-site continuation, not a new cost regime or project budget cap.

## Comparisons declared before energies

1. Report ordered `D_i=[HansDy−HansLa]−[MexDy−MexLa]` for each Hans source, with
   the exact same Mex endpoints. Preserve all native/solvent components and
   initial→continued sensitivity. No pooling/minimum selection between crystals.
2. Check actual atom identities, elements, charge, proton inventory, source
   fragments and cap recipes across8DQ2/8FNR. Only if they match, report each
   metal's work `E_Hans,8FNR(M)−E_Hans,8DQ2(M)` and itsDy-minusLa difference.
   Coordinate differences and generatedH differences remain explicit; this is
   not an isolated torsion intervention. If mappings/compositions differ, retain
   distinct conditional comparisons and leave that work unavailable.
3. Compare the full ordered vectors with the same qualitative protein-level
   apparent-affinity evidence. No site labels, fitted weights, thresholds,
   populations, averaged Kd or favorable-site selection. A source-sensitive sign
   is a limitation to report, not a reason to optimize it away.

New code/workspace versions preserve the completed original pilot, its source
snapshot and every numerical receipt. Root receives the prepared-source/mapping
finding before execution; no further user permission round is required.
