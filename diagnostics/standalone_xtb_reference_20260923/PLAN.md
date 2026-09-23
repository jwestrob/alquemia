# Standalone xTB reference fidelity: static origin and minimal adaptive pool

Root explicitly authorized this finite follow-on under Jacob's Sept23 overnight
instruction after the224-cell backend diagnostic passed4/4 energy and14/14 tight
force checks with useful known-class ordering. No optimization, new MACE, DFT,
225-fold execution or default change is included.

## Frozen selection and calls

Use the original completed30-case adaptive collection, retain exactly its25
canonical calibration sources plus three consumed crystals1H4I/4MAE/1KB0.
Exclude only the two unlabeled PLM cases. Canonical membership/labels come from
the original released reference's context/calibration rows, not current scores.
Both metals receive exactly origin, adaptive_Ca and adaptive_La, unchanged.
28 sources ×3 geometries ×2 metals ×2 media =336 logical solver cells. The actual
completed accuracy0.02 pilot supplies36 exact coordinate/state cells for1H4I,
4MAE andQ88JH5. **300 new standalone calls**, zero MACE/DFT/optimization.
Source geometries/charges/waters/caps/protonation and native MACE values are
unchanged and pinned. Missing cells remain unavailable; no native-solvent or
origin-only adaptive fallback. The earlier source pool remains preserved.

## Backend and numerical policy

Exactly the tested xTB6.7.1 backend and installed pinned GFN2 parameter file,
accuracy0.02 only, fresh/no-restart, electronic300K, singlet/unpaired0,500maximum
SCC iterations. Both media are standalone; ALPB water gsolv, normal230-point
surface, default solvent parameter298.15K/dielectric80.2/P16/GBOBC, H-bond term,
no ionic screening. Eight workers×eight threads on64CPUs/128GiB, noGPU. Exact
control/defaults and parser state checks are reused from the completed pilot.
Eight reused pilot q0 cells used gradients instead of SP; their exact same energy
functional/printed total is reused without duplicating calculations.

This is a new descriptor E = nativeOMOL(vac) + standaloneGFN2(ALPB−vacuum),
R=E_Ca−E_La. Do not borrow the old native calibration, use a mixed backend
solvation pair or remove ALPB surface/H-bond/reference components. The prior
small numerical diagnostic is supporting evidence, not proof that every source
has a unique SCC root. Report actual failures and any state/charge anomaly.

## Separate calibration and reporting

Compute two explicit variants at fixed geometry/state:
1. Static origin: both endpoint energies atq0.
2. Minimal adaptive pool: both metals see the same3candidates; record mathematical
   minima separately, report the existing0.1kcal/mol origin-retention rule.

Calibrate EACH variant independently on only the designated25 canonical cases,
using the unchanged class-extrema/minimum-gap rule from the existing comparator.
Require all25; insufficient separation/missing cases yield unavailable bands,
not threshold adjustment. Freeze separate reference records BEFORE any225-fold
transfer. Apply each variant's own bands to the three consumed crystals. Keep
old native scores/decisions visible and old-band transfers explicit. These are
consumed calibration/transfer structures, not new blind biological evidence.

Report all28 denominator statuses,25calibration calls/gap,3crystal calls, raw
static/adaptive differences and candidate choices. Group biological labels
appropriately; no affinity/free-energy/probability or universalzero claim.
A useful result supports a later frozen225-fold test; no such launch belongs to
this manifest. No third accuracy or automatic failed-cell retries.

Record exact336 logical cells,36 actual reused receipts,300 new attempts,
finite task manifest, executable/source pins, actual wall/CPU/memory cost,
regression/mapping/real-output tests, commands and vault note. Technical parser
fixes use preserved raw output without additional chemistry.
