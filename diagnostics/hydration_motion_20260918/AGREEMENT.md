# Local water-motion check — approved 2026-09-18

Jacob: “proceed.” This follows completion of the joint arrangement experiment
and the explicit next step of bound-water motion/free-energy corrections.
The concrete scope is the local check in
hydration_occupancy_20260918/NEXT_MOTION.md, not a production occupancy score.

## Question and finite experiment

Does inexpensive native MACE-OMOL curvature approximate the native DFT/CPCM
restoring force for collective radial translation of bound water? Use four
completed real endpoints:1F6S pattern11 Ca/La and6IP9 pattern110 Ca/La.
The latter includes A322's poor MACE/DFT translation direction, selected by
physical gradient disagreement before these displacement energies are inspected.

Each retained variable water moves as a rigid H2O along its own original
metal–oxygen radial unit vector. A common scalar h moves every variable water
by h Angstrom. Metal, protein, outer waters, water orientation and internal shape
remain fixed. Directions are frozen at the center, identical between metals.
Positive h moves outward. Use h=±0.05A for8 native DFT energy/analytic-gradient
checks; native MACE evaluates those and h=±0.025A (16 new cheap calls).
Reuse all four center energies/gradients, including actual MACE force arrays.
No water insertion/deletion, protonation change or geometry optimization.

DFT: unchanged native r2SCAN-3c/NoAutostart/CPCM(Water)/DefGrid3/TightSCF EnGrad.
MACE: unchanged native unmasked OMOL0-100M checkpoint, float64, vacuum, fixed
charge/multiplicity per endpoint. No extra composite correction or solvent fit.

## Prediction and frozen acceptance

DeltaE_pred(h) = h*g_DFT(0) + E_MACE(h)−E_MACE(0)−h*g_MACE(0).
Projected g is dE/dh in kcal/mol/A, never a force with the wrong sign.
This anchors the target CPCM energy/force at the center and approximates force
changes with vacuum MACE. It is not self-consistent MACE/CPCM or an entropy model.

Retain all endpoint results and the Ca−La difference separately. For each metal:
- Anchored energy error ≤max(0.02kcal/mol,25% of |DFT even energy|), at each sign.
- Even-energy error ≤max(0.005kcal/mol,25% of |DFT even energy|), reusing the
  earlier local-curvature criterion. Even = [E(+h)+E(−h)−2E(0)]/2.
- MACE coarse/fine even energies, with the fine value scaled by4, agree within
  max(0.005kcal/mol,5% of |coarse even energy|).
- MACE central-difference gradient residual times h ≤max(0.02kcal/mol,
  5% of |h*g_MACE(0)|), reusing the previous analytic-gradient consistency check.
- The anchored projected gradient g_DFT(0)+g_MACE(h)−g_MACE(0) differs from
  actual g_DFT(h) by ≤max(0.4kcal/mol/A,25% of |g_DFT(h)−g_DFT(0)|).
  The absolute floor corresponds to0.02kcal/mol across the0.05A displacement.

The unit/sign/geometry checks must pass independently. Negative curvature is
reported without clamping. A passed single collective direction does not qualify
all6n water coordinates, basin entropy, oxygen relaxation or equilibrium occupancy.
Score/entropy corrections remain disabled in this experiment regardless of outcome.

## Execution

Existing finite task executors only. DFT:64CPU, four16-rank tasks,256GiB.
MACE:one standard GPU,16CPU,64474MiB (one eighth of host RAM), existing pinned
venv/checkpoint. No project compute/time cap. Estimate roughly two recent DFT
batches (~12min allocation wall) and a short GPU run, then record actual cost.
Preserve baseline/default, previous results, concurrent work and native orientation
jobs1201824/1201825. No automatic promotion, new threshold or production rescore.
