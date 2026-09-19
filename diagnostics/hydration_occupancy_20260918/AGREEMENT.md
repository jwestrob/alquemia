# Joint water occupancy development — approved 2026-09-18

Jacob: “I fully approve.” This follows the explicit proposal to enumerate all
four 1F6S and eight 6IP9 water arrangements for Ca and La, use MACE for geometry
and native DFT for energy, reuse four compatible full-water endpoints (24 total,
20 new), and report the missing bound-water free-energy contribution required
to change occupancy/discrimination. Production/default remains unchanged.

## Executed scope

Two consumed alpha-lactalbumin structures, one biological group. The fixed
source-defined 70/76-atom hydrogen-bond contexts and candidate water identities
come from hydration_network_20260918/CONFIG.json. All outer waters, protein,
metal and retained water oxygens remain fixed. Enumerate every subset of the
two/three variable waters; never select a deletion by classification. Neutral
H2O deletion preserves formal charges (La−1/Ca−2) and singlet parity.

For eight nonempty, nonfull patterns x two metals, run two fixed starting
orientations (source/radial-away) through the existing exact rigid-water MACE
optimizer. Sixteen tasks /32 searches; identical model/settings to the completed
full-state pilot (BFGS, gradient tolerance .001 eV/radian, maxiter200 is a
numerical safeguard, not a compute-time budget). Require both starts to converge;
choose lower MACE energy only within the same composition/metal. Do not use MACE
absolute energies to rank water counts. Empty patterns require no optimization.

DFT:20 native r2SCAN-3c/NoAutostart/CPCM(Water)/DefGrid3/TightSCF EnGrad endpoints,
four concurrent16-rank tasks on an existing64CPU allocation. EnGrad supplies
analytic force information for later water-motion checks at little extra cost.
Reuse the four actual, matched full-water endpoints from proposal_dft_v1.
No second composite correction, numerical gradients, Hessians or new functional.
GPU preparation uses the established16CPU/oneGPU/one-eighth host-RAM allocation.
No project compute/time cap; record actual receipts. Do not disturb native
orientation comparator jobs1201824/1201825.

## Outputs and interpretation

Complete pattern/metal DFT energies, water addition edges and nonadditivity,
conditional fixed-water-count minima, raw R=E_Ca−E_La, and changes from the full
pattern. Tie the exchange ledger to the computed gas/liquid-water reference.
Report the missing bound contribution for neutral exchange, not invented free
energies or probabilities. Bound-water motion/vibration and non-electrostatic
solvent terms remain unavailable. Fixed oxygen states are conditional electronic
comparisons, not equilibrated occupancies. No old threshold transfers.

The broader approved direction includes locally estimating bound-water motion
with MACE checked against DFT and proposing genuinely missing water sites from
experimental alignments/coordination space. These require concrete follow-on
manifests. This manifest implements the first arrangement experiment; it does
not silently include new sites, oxygen relaxation or an assumed entropy penalty.
