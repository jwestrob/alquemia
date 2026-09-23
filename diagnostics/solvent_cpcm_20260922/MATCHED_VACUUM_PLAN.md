# Actual solver-switch qualification

The initial eight CPCM calculations explicitly print “Setting regular ORCA SCF
and ShellPairData” despite `UseXTBMixer true`. This invokes the numerical
contingency declared in PLAN.md. Parent acknowledged the material switch and
directed retaining compatible converged CPCM cells rather than rerunning them.

This next finite submission is **eight fresh vacuum calculations only**, one for
each already declared physical endpoint. It changes no geometry, chemical state,
GFN2 parameters, electronic temperature, convergence tolerance or MACE energy.
Explicit ordinary SCF (`UseXTBMixer false`), genuine block `Guess MORead` from
each matching archived native-vacuum GBW, default TolE=1e-6 hartree, MaxIter500.
Both scalar parameters and actual initial-guess behavior must be checked.

The collector pairs these with actual converged ordinary-SCF CPCM results from
the initial eight. Initial failed/unfinished CPCM remains unavailable; this
submission contains **zero CPCM retries**, no DFT and no MACE calls. The starting
guess is a numerical aid, not a changed Hamiltonian. A separately recorded seeded
CPCM retry can be considered only after initial failures are known.

Report ordinary-minus-native vacuum differences separately. The new transfer is
the ordinary-CPCM minus ordinary-vacuum total-energy difference, never the printed
CPCM dielectric component alone. Then add native MACE vacuum. All new decisions
remain unavailable until a complete compatible canonical reference exists.

Use 64 CPU/128 GiB, eight concurrent eight-rank calls on the explicitly assigned
GPU-partition CPU host; no GPU, arbitrary project budget or new time limit.
Initial pilot and this qualification keep separate receipts and costs.
