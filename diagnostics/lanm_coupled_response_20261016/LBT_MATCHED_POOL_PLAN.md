# Finite matched-FF rescore of the complete consumed LBT pool

Declared before new energies,16 October2026. Root stages/submits after review.
Question: do the existing actual accommodations improve the matched-model
comparison, including source sensitivity, when reselected with fixed native FF
ALPB? These candidates were searched under MACE plus self-consistent GFN solvent.
This is transfer rescoring of that fixed consumed menu, not optimization or
stationarity under the new Hamiltonian.

## Frozen matrix

Use complete Slurm12560 `lbt3_basis_pool_v1/execution/RESULT.json`, SHA
`b3d1b4dba88040549d69499af78a6179b911d3649b353933f02eb5e3e58066cf`.
It contains18 unique geometries and36 accepted La/Dy states:9 geometries/source.
All raw MACE records and XYZ/physical pins were found in a read-only inventory:
4 origin records in12500,16 own-metal POINT records,16 opposite-metal CELL records.
All carry whole256,q−2,La1/Dy6,omol metadata. Inventory is pinned separately.

|Menu|7CCO geometry indices|7CCN geometry indices|Global size|
|---|---|---|---:|
|Origin|00|01|2|
|Prior range pool|00,02,03,06,07|01,04,05,08,09|10|
|Matched7-coordinate control|prior+10,12|prior+14,16|14|
|10-coordinate candidate|prior+11,13|prior+15,17|14|
|Full union|00,02,03,06,07,10,11,12,13|01,04,05,08,09,14,15,16,17|18|

Preserve original12560 menus verbatim and add explicitly named origin/full-union
summaries. Both metals get each complete menu. All36 states remain in the result;
failures retain their slots and withhold every dependent menu. No fallback,
replacement, retry, further search, changed bounds, additional trial candidates,
water changes or excluded source. The original source0.40Å/target0.20Å/torsion
0.20rad search bounds and accepted source/contact receipts remain historical
admission evidence; this task changes no coordinates or guard.

## Calls and reuse

Reuse four exact origin FF-ALPB cells from Slurm12715 and all36 existing MACE
energies. Make at most32 new native256-atom FF-ALPB calls, sequentially in
geometry00→17 order, La then Dy, skipping the four origins. No new MACE, GFN SCC,
MM, DFT or optimization. Use the exact12715 binary and its pinned build/library
lineage. One CPU, zeroGPU. Native gradients are archived by the same wrapper but
are not used for additional response work. Count actual attempts and preserve
all failures. Root records scheduler costs; do not extrapolate prior wall time
into a promised runtime.

Reuse `LBT_MATCHED_MODEL_V2.inputs` and `solvent` for unchanged native FF charge,
terminal/state, MACE-model, Mex12685 and native-functional admission. The solvent
vector remains the original255 ff19SB charges plus formal+3; no charge response,
normalization or whole-peptide FF vacuum term. Whole256 is Q; field/MM terms
are partition-defined zero. Reuse old MACE only; exclude all old GFN solvent
energies from the new total.

Before the first call, verify the exact12560 snapshot,18 aliases/coordinate hashes,
36 accepted raw MACE receipts, each state/model, and all XYZ/physical pins.
Each candidate must retain the source's ordered256 IDs/elements/residue identity;
its coordinates may differ. Require saved physical/XYZ agreement at1e−12Å and
the original canonical float64 coordinate hash. For origins require exact12715
input/charge/binary pins and raw component replay. Do not apply the deposited
source-coordinate equality gate to moved candidates. Existing search/admission
records are preserved, not rerun with a different geometry rule.

## Reporting

For every source-local and global menu report selected La/Dy geometry IDs, exact
ties, matched native component totals, selected Dy−La contrast, each metal's
component work from the matching origin menu, and balanced exchange against
unchanged Mex12685:

R = (E_Mex,Dy−E_Mex,La) − (E_LBT,Dy−E_LBT,La).

Positive R favors Mex–La+LBT–Dy. Compute differences in native units before
conversion using the exact V2 constants. Report component R and change from the
origin menu. Keep both source contexts visible even if the global union chooses
one. Preserve old12560 summary and selections alongside new menus as separately
named evidence; never insert incompatible old solvent offsets into new scores.
The source-local/control/candidate/global menus are all retained; no favorable
menu becomes the sole result. A later saved-only comparison may combine the
new matrix with Mex12717 accommodation, but this scorer pins Mex12685 baseline.

The finite global menus are not populations. Existing source magnitude spread,
terminal-charge hypothesis, dry LBT versus hydrated Mex, whole-Q versus capped
local-Q partition, frozen-f and spin limitations remain. Improvement within this
menu is not a matched-FF optimizer result, equilibrium affinity, protein Kd,
independent validation or successful classifier.
