# Occupancy free-energy accounting

For a frozen protein/outer-water context and explicit variable-water state i:

Omega_i - Omega_empty = (E_DFT,i - E_DFT,empty - n_i E_gas_water)
                + DeltaE_relax,i + F_rigid_water,i
                + F_internal_water,i + DeltaG_nel,i
                - n_i (mu_liquid_water - E_gas_water).

E terms use the same native Hamiltonian, composition ledger and electronic units.
The liquid reference is the existing actual gas-water ORCA energy plus recorded
thermochemical terms at298.15K/1bar. No CPCM-water energy is silently substituted
for that chemical potential. The water-reference approximation remains explicit.

F_rigid_water describes the coupled translations/librations of bound waters
conditional on the fixed protein and water shapes. Its harmonic version requires
supported stationary, stable basins and physically valid thermal extent. For
mass-weighted frequencies omega_j, F=sum_j[hbar*omega_j/2 + kT*log(1-exp(-hbar*
omega_j/kT))]. No artificial spring, deleted unstable mode, frequency floor or
quasi-RRHO cutoff is introduced to obtain a favorable finite value. Actual masses
and inertia determine frequencies. Alternative physical orientation basins must
be counted once; H-label permutations are not separate physical minima. Site
subsets already enumerate distinct placements of indistinguishable waters.

F_internal_water is not supplied by a rigid-body calculation. It includes bound
water internal vibration relative to the electronic geometry used above; the gas
ZPE cannot simply be relabeled as a measured bound correction. Changes in core
vibrations/polarization beyond the frozen-context model also are not established.

DeltaG_nel is the change in non-electrostatic solvent free energy relative to the
empty state. ORCA's CPCM energy/derivatives omit this contribution. SMD supplies a
CDS component, but also changes electrostatic radii; replacing CPCM with SMD would
be a different target. The documented alternative GVDW_nel/GSES_nel parameter
sets coverH/C/N/O, not the present Ca/La solutes. Do not invent metal parameters
or call an undocumented zero missing-term a complete occupancy calculation.
See the primary [ORCA6.1 solvation accounting](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html#calculation-of-the-free-energy-of-solvation-within-the-c-pcm).

Omega is the grand free energy, G_i - n_i*mu; the expression above already
subtracts the reservoir term exactly once. If all required state terms and
competing basins are supported, use weights exp[-(Omega_i-Omega_empty)/RT]
with one common reference, normalize with log-sum-exp,
and retain the full denominator. Otherwise occupancy remains unavailable; report
known components and the actual missing terms. Baseline decision bands never
transfer automatically to this changed model.

The local vibrational formula follows harmonic statistical mechanics, as described
in the primary [ORCA thermochemistry manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html).
The available pinned MACE environment contains ASE's maintained HarmonicThermo;
its implementation can evaluate a supplied, supported spectrum without omitting
imaginary modes. That machinery does not itself validate the supplied curvature.
