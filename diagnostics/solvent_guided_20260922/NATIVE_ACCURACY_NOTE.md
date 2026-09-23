# Native GFN2 numerical controls — read-only finding, fallback proposed only

**No separately documented native-mixer accuracy control was found.** A working
standalone xTB installation does expose one. No molecular calculation, installation,
parameter change or backend switch was performed for this note. Existing numerical
replays and restart investigation belong to the separate source-start track.

## Native ORCA 6.1.1: what is and is not established

The version-specific manual states that the native mixer replaces the ordinary
SCF settings. Its native keyword table supplies `UseXTBMixer` and `XTBFOD`, with
`SmearTemp` described separately; it gives no accuracy multiplier or explicit
charge/multipole residual tolerance. `%xtb Accuracy` is documented for the
external `otool_xtb` interface. It is not evidence of a native control.
[ORCA 6.1.1 manual](https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/semiempirical.html#scf-with-native-xtb).

Actual installed-version outputs prevent a stronger claim that generic settings
have absolutely no effect. The archived native `Convergence Tight` 1H4I-Ca vacuum
run prints TolE=1e−8 Eh and stops at |ΔE|=1.7610e−9 Eh. However, its maximum and
RMS density changes are 8.1311e−5 and 2.5347e−6, versus printed generic tolerances
1e−7 and 5e−9. It explicitly reports the special xTB mixer and energy convergence.
The primary Q88JH5-Ca ALPB probe similarly stops at |ΔE|=2.1000e−7 Eh with
TolE=1e−6; maximum/RMS density changes exceed the displayed generic tolerances.

These records show that the displayed density thresholds are not jointly enforced
as ordinary density-convergence tests here. They do not reveal the native internal
charge/multipole stopping rule or prove that changing TolE controls it. Energy-only
agreement cannot certify a unique converged electronic solution. No undocumented
keyword is proposed. Exact input/output pins are in `NATIVE_ACCURACY_EVIDENCE.json`.

## Already installed standalone fallback

`/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb` runs successfully for `--version`
and `--help`: **6.7.1 (edcfbbe)**. Its help includes `--acc`, `--alpb`, `--gfn`,
explicit charge/unpaired-electron settings and electronic temperature. No model
was evaluated. The installed default GFN2 parameter file explicitly contains
Z=20 and Z=57 entries, so Ca and La are covered; this is parameter availability,
not evidence that their relative affinity is accurate.
[Versioned upstream parameter source](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/param_gfn2-xtb.txt).

Standalone `--acc` changes integral screening and SCC convergence; lower values
are tighter. The documented examples include accuracy1 and0.2, corresponding to
energy criteria1e−6 and2e−7 Eh. GFN2 uses a tighter wavefunction criterion than the
generic example table. This is a documented numerical control, unlike inferring
native behavior from ordinary ORCA printouts.
[xTB accuracy documentation](https://xtb-docs.readthedocs.io/en/latest/sp.html#accuracy).

GFN2 ALPB water is supported. The standalone solvent includes electrostatic,
surface/H-bond and reference-state terms; surface-grid and reference-state choices
must also be pinned. Thus shared names do not establish equality with native ORCA.
[xTB solvation documentation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html).

**Proposed next step only:** if the separately owned native replay remains
discontinuous, a fixed-coordinate comparison using the same problematic archived
geometries and both documented accuracy levels could test the standalone backend.
Both vacuum and ALPB must come from that same backend/settings, with identical
charge, spin, 300K temperature and coordinates. Do not combine standalone ALPB
with an old native vacuum energy, inherit a reference, or change production from
this documentation finding. Parent must freeze the exact finite qualification
manifest before any new calls. No ordinary-SCF rerun or larger benchmark follows
automatically.
