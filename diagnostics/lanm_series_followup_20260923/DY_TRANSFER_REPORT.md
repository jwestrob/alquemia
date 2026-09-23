# LanM structural transfer fails: the initial signal is not robust discrimination

**The apparent Hans/Mex La/Dy signal does not survive an actual Dy-conditioned
Hans crystal.** The original positive EF1/EF3 values reverse, and all three
Dy-source values are negative. Same composition, formal charges, mapped atoms
and cap recipes are retained. These results withdraw any interpretation of the
initial pilot as robust within-series discrimination.

All6freshMACE and24nativeGFN2 calls succeeded;6immutable Mex endpoints were
reused. No geometry optimization, DFT, changed water inventory, favorable-source
selection or default change occurred. The completed original result remains
preserved in REPORT.md and its immutable collection; this is the follow-on result.

## Ordered relative-selectivity result

`D_i = [HansDy−HansLa] − [MexDy−MexLa]`, positive for stronger relativeLa
preference inHans. Every E uses native OMOL plus native GFN2(ALPB−vacuum) with
the declared self-continuation. The same Mex source/endpoints occur in both rows.

| Site | Hans8DQ2,La-conditioned | Hans8FNR,Dy-conditioned | Source change |
|---|---:|---:|---:|
|EF1|+4.866368|**−9.568461**|−14.434829|
|EF2|−0.264440|**−28.187049**|−27.922610|
|EF3|+12.499116|**−20.673113**|−33.172229|

Units are modelkcal/mol, not measured affinity/free energy. The same protein-level
qualitative apparent-affinity comparison applies; no site labels or fitted
weights were introduced. Both source vectors are retained, not averaged or
selected for agreement with the biological direction.

| Site | Dy-source native MACE component ofD | Dy-source solvent component ofD |
|---|---:|---:|
|EF1|−9.119257|−0.449204|
|EF2|−29.013876|+0.826827|
|EF3|−19.814857|−0.858256|

The reversal already exists in native MACE. The solvent transfer does not remove
it. This identifies a serious source-geometry dependence of the current
descriptor; it does not prove a unique microscopic cause or that only MACE is
responsible for the missing affinity physics.

## Full four-cell Hans matrix at every site

These are absolute **composite model energies in kcal/mol**. Different elements
carry different atomic-energy offsets; their raw totals are not comparable
affinities. The balanced exchange above cancels those offsets. The full native
MACE eV and native vacuum/ALPB Hartree terms are also retained in
DY_TRANSFER_RESULT.json, without rounding.

| Site | Crystal source | La endpoint | Dy endpoint |
|---|---|---:|---:|
|EF1|8DQ2,La-conditioned|−880572.435086290|−1416113.845125990|
|EF1|8FNR,Dy-conditioned|−880644.042880789|−1416199.887749609|
|EF2|8DQ2,La-conditioned|−880611.334927121|−1416151.070208625|
|EF2|8FNR,Dy-conditioned|−880746.946523477|−1416314.604414634|
|EF3|8DQ2,La-conditioned|−880648.238856529|−1416182.184414610|
|EF3|8FNR,Dy-conditioned|−880663.070875450|−1416230.188662690|

The two structures have matching chemical compositions and mapped atoms, so
within-metal crystal-source energy differences are reportable:

| Site | La work,8FNR−8DQ2 | Dy work,8FNR−8DQ2 |
|---|---:|---:|
|EF1|−71.607794|−86.042624|
|EF2|−135.611596|−163.534206|
|EF3|−14.832019|−48.004248|

**Both metals prefer the Dy-source core at every site; Dy benefits more.** Thus
even selecting the lower endpoint over these same two source geometries for
each metal would select8FNR for both and preserve the negative pattern. This is
an algebraic consequence of the completed matrix, not a newly executed ensemble
or a rescued affinity estimate. Mex still has only its original Nd-conditioned
source, so such a search would also have asymmetric structural coverage.

Large source works include all represented coordinate changes, including
generatedH differences and cap geometry. They are not pure E9torsion energies or
protein deformation free energies. No mechanical scaffold cost, conformational
entropy or complete metal-binding thermodynamic cycle was computed.

## Physical source checks

The fixed source is8FNR author chainA, Dy201/202/203 mapped toEF1/2/3. All three
sites are supported. Each retains50atoms, charge−1 for both endpoints and no
deposited water inside the unchanged3.3Å rule. EF4 and its water remain outside
the declared primary vector; no alternate chain was substituted.

Every source-heavy coordinate is preserved exactly. The existing pH5 protonation
and amide repair produce the same atom identities/elements, H inventory and
sigma-cap identities/lengths as8DQ2. Both oxygens of E9Glu42/66/91 remain in each
whole sidechain fragment, even though onlyOE1 lies inside3.3Å in8FNR. Each source
has its own crystal geometry; matching means chemical/source-atom identity, not
identical coordinates or an artificially aligned/rebuilt donor arrangement.

The actual OMOL charge/spin batches and native effective electron counts pass:
La physical multiplicity1, Dy6; native GFN f-in-core effective singlet with150
valence electrons at all six new endpoints. All12 new continuations positively
reportXTBRESTART with unchanged parameters and charge closure. The largest new
individual initial→continued change is0.004388kcal/mol. The new Hans continuations
change the three reportedD values by−0.000021,+0.000536,−0.008451kcal/mol,
respectively. This is much smaller than the14–33kcal/mol source dependence.
It does not establish universal native-SCF convergence or a verified Dy
electronic ground state; those limits remain as in the first pilot.

## Execution, tests and disposition

| Job | Actual calls | Allocated time | CPUs | Core-seconds |
|---|---|---:|---:|---:|
|1210367|6MACE|10s|32|320|
|1210368|12GFN initial|26s|64|1664|
|1210372|12GFN self-continuation|26s|64|1664|
|Total|6MACE/24GFN| | |**3648**|

OneH200 was requested for10s: **10requestedGPU-seconds**. Native jobs requested
noGPU. Scheduler peak-memory values are unavailable; do not report them as zero.
Local preparation/reporting is additional. Combined with the original pilot:
18MACE/72GFN,10128allocatedcore-seconds and31requestedGPU-seconds. Zero molecular
failures and zero DFT calls in either pilot.

Five new real-fixture tests pass in0.702s, zero skips: exact finite source scope,
source/cap/composition mapping, retention of the secondE9oxygen despite its
changed distance, exact Mex reuse and exchange-cycle algebra, and actual physical
spin/effective-valence/restart receipts. All30 new molecular calls completed.

No more lanthanide calculations are authorized by this completed branch. Keep
the baseline unchanged. **Do not use this current fixed-core score as a robust
within-series affinity discriminator.** A future method would need to address
state/geometry thermodynamics and structural coverage; this result supplies no
validated correction for those omissions.

Final collection:
`workspaces/lanm_series_followup_20260923/dy_transfer_v1/final_collection.json`.
Costs: same directory,`COSTS.json`. Full four-cell/component matrix and sensitivity:
`DY_TRANSFER_RESULT.json`. Exact commands and provenance:
`DY_TRANSFER_COMMANDS.md`,`DY_TRANSFER_ARTIFACTS.json`.
