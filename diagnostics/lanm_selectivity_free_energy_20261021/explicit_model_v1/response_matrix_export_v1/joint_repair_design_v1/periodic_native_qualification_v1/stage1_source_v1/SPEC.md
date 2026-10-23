# Isolated Stage1 source contract

This source implements only Stage1 of the frozen `../NUMERICAL_CONTRACT.md` and root intent `e4dc6021`. Build and native execution require root release after source review. `PROPOSAL.md` describes the separate prospective energy/force stage. No such entry is implemented here.

## Implementation and binding

`PARITY_SETUP.f` is derived mechanically from the admitted source's complete checked `periodicresponsesetup`. `SETUP.diff` records only its new symbol and removal of the final backend-guard call from this isolated entry. The original library, normal setup, backend guard, stability guard, geometry kernel, native field routines and all energy/force paths remain unchanged. `PARITY.f90` links this entry externally against the admitted `libtinker.a`; it never calls normal induction or energy/gradient entry points. Default applications remain disabled for active periodic response.

`PREPARE.py` binds only the old e41f La candidate, exact composition and job20053 static input. It retains the whole composed key and appends the exact previously qualified PME settings. No force-field value is recomposed. It writes literal atom/type/xyz/alpha/physical-coefficient expectations and the exact ten frozen probe vectors. The fixed expected normalized metal beta comes from the saved static manifest. Input preparation performs no operator or geometry-response calculation. The prepared key must not be used for energy or dynamics under this Stage1 release.

The harness validates the actual native map/alpha/coordinates/coefficients exactly, positive alpha, active sites, original domain/topology, simple-loop orthogonal TINFOIL/FFTW/MUTUAL mode, 49.324 Å cell, 9 Å cutoff, .4 inverse-Å Ewald split and 64³/order5 grid. It checks native alpha*beta against the admitted static beta before applications. It calculates h/H using the unchanged native compact kernel; H is not copied from a saved result. Runtime SHA verification of all prepared files, executable, immutable library/build and release is the launcher responsibility. Source readiness is not a runnable release.

## Operations, state and raw output

One `dfield0c` call initializes the native direct fields and PME state without an induced solve. Raw Ed/Ep are reported separately, along with the h that the unchanged common `dfield0` wrapper adds to both. Stage1 does not claim to test the common RHS addition during an actual induced solve; the separately declared Stage2 owns that check.

Thirty calls to `ufield0` supply complete normalized actions `z-C*field`; thirty separate primitive-chain calls supply positive normalized reciprocal actions. The latter chain matches the qualified harness and uses the native unsymmetrized qfac, grids, transforms and interpolation. The total includes sixty reciprocal applications. Modes are `(probe,0)`, `(0,probe)` and `(probe,next)` in that exact order; next is cyclic. The two channels are preserved independently. There are no induced solves, energy/gradient calls, sparse assemblies or factorizations.

`parity_state.dat` contains actual atom/type/coordinates/alpha/pdamp/Thole, both coefficient records, all h/H and raw direct fields, native metal–O POLPAIR gamma, and current grid/splines/qfac. `pme_norm_export_001.dat` must be byte-identical to the admitted central snapshot. The state exporter is tied to the same initialization as all actions. qfac, splines, alpha, h/H and grid map are checked after every application; u/v and FFT scratch are restored between calls. The receipt counts complete actions and reference applications separately. Raw channels have the same electronic/Å units as native AMOEBA; normalized action/probes are dimensionless.

## Saved comparison

`REPLAY.py` never calls native or a solver. It verifies exact prepared-input hashes and literal current snapshot/site/state equality, uses the saved admitted lower sparse CSC real+self+H matrix, and constructs a sparse grid scatter map from the pinned snapshot reader. Full complex FFT replay retains cross-channel leakage from unsymmetrized qfac; it never builds a dense global response matrix. For every mode/probe/channel it compares native reciprocal versus independent reciprocal, complete native versus sparse-real plus reciprocal, and complete-minus-native-reciprocal versus sparse-real. Each comparison uses the unchanged sampled bound `1e-11*max(1,max(abs(predicted)))`.

Additional implementation correspondence checks use that same scaled arithmetic bound for ordered independently reconstructed compact h/H, native damping-length readback, and normalized beta versus job20053. Thole and POLPAIR tokens, alpha/site/xyz/coefficient identities are exact. These added checks are prospective numerical correspondence checks, not electrical acceptance limits. The old initialized snapshot's existing norm admission is reused; its large norm sums are not recalculated here. Failed comparisons retain the completed check records and all raw files.

Passing provides finite ordinary-double action evidence at this one La state. `e_app` remains unknown; no uniform bound, force qualification, Dy qualification, trajectory, MD throughput or physical-accuracy claim follows. Root must review the actual result before any Stage2 guard implementation or evaluation release.
