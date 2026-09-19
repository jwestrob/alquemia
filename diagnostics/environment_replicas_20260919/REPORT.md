# Static quantum second-shell gains survive both GGR replicas

**All six native DFT alpha/GGR margins improve under the unchanged context rule.**
The weakest increases from 0.654840 to 4.564286 kcal/mol. Native OMOL applied to
the same complete contexts improves 2/6→6/6 directions. These are three GGR
structures and two alpha structures representing one consumed biological
comparison, not six independent validations. Production/default remain unchanged.

| Alpha / GGR | Native DFT core | Native DFT context | Native OMOL core | Native OMOL context |
|---|---:|---:|---:|---:|
| 1F6S / 1GLG | 10.923654 | 14.814843 | 8.303640 | 35.437320 |
| 1F6S / 2FW0 | 0.654840 | 6.473524 | −2.763855 | 25.319605 |
| 1F6S / 2FVY | 1.835024 | 4.564286 | −4.405836 | 20.391748 |
| 6IP9 / 1GLG | 14.572699 | 18.421517 | 9.613261 | 28.416934 |
| 6IP9 / 2FW0 | 4.303885 | 10.080199 | −1.454233 | 18.299219 |
| 6IP9 / 2FVY | 5.484069 | 8.170960 | −3.096215 | 13.371362 |

Entries are `R_alpha−R_GGR`, R=E_Ca−E_La. DFT uses kcal/mol and native OMOL the
converted model-kcal scale. Quantum CPCM and learned vacuum energies are distinct
Hamiltonians; do not add them or interpret equality of magnitudes. OMOL's effects
are much larger. Core→context GGR shifts are −6.635734/−8.563229/−5.473806 kcal/mol
in DFT, versus +10.835687/+9.885907/+13.171782 in OMOL. Relative directions improve
in both even though their individual environmental shifts have opposite signs.
This is a useful predictive development result, not a validated numerical DFT
correction or proof of a unique physical cause.

The reused static PQQ gap remains positive but narrows: DFT30.571889→28.1792,
OMOL102.222180→89.825205. Only two PQQ development structures have been evaluated
with this context rule. All25 reference fidelity and protocol-specific absolute
bands remain unestablished for it. A full reference test is the relevant next
step before advocating broad scanner use.

## Frozen chemistry and what changed

Four expanded native r2SCAN-3c/CPCM(Water)/DefGrid3 SPs were executed on 2FW0/2FVY;
all four corresponding actual original-core DFT energies were reused. All prior
alpha/1GLG/PQQ static inputs/results were reused. There was no geometry search,
new water, changed label, protonation change, microstate selection or threshold
fit. This is separate from the completed, unfavorable coordination-proposal
experiment; no geometry from that experiment entered this result.

Both new replicas expand52→115atoms and add formal charge−1. 1GLG expanded52→125
with the same charge change. The identical3.5Å complete polar-neighbor rule and
3.3Å original donor anchors select these different source contexts. Original
source atoms/surviving caps remain exact. Overlap completion replaces cut-bond
caps with their actual source atoms. Unsupported nearby species remain errors;
outer waters retain the explicitly excluded/fixed-inventory policy. The full
graphs, added residues, charges and coordinates are recorded in the manifests.

The effect includes composition, cavity, electrostatics and quantum interactions.
Do not call the gain uniquely a hydrogen-bond or polarization contribution.
GGR structural spread remains substantial; these results strengthen margins
without eliminating all representation/structural sensitivity. The smaller
six-gap denominator cannot establish broad biological affinity prediction.

## Actual calculations and cost

- DFT1202453: all4endpoints normal/converged,555s ×64CPUs=35520allocatedCPU-s.
- OMOL1202454: all8core/context endpoints computed,20s ×16CPUs=320CPU-s;
  20GPU-s. Eight-call worker time19.298328s. No failed attempts.
- Total **35840allocatedCPU-seconds,20GPU-seconds**. DFT batch MaxRSS26581792KiB.
  Local preparation/tests/reporting are additional, not assumed free.

Three real-source preflight tests pass, zero skips. They check exact regeneration,
all original coordinates/energy receipts, unchanged charges/paired geometry and
explicitly unavailable references. Actual scientific outputs above are separate
from parser/geometry checks. This is development cost, not a matched end-to-end
production speed benchmark.

Unrounded DFT result: `workspaces/environment_replicas_20260919/prepared_v1/collection_1202453.json`.
Actual OMOL endpoints: `.../mace_collection_1202454.json`; compact full matrix:
`MACE_RESULT.json`. `SUBMISSIONS.json` pins both finite manifests and runners.
No missing result was replaced with a baseline score. Old absolute bands were
not transferred. No default change, push or production rescore.

**Recommendation:** pursue the compact-context native-OMOL candidate on the
complete PQQ reference panel. Its accuracy/coverage must earn promotion; the
larger raw margin alone is not evidence of generalization.
