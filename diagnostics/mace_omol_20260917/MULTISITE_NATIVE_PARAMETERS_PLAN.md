# Multisite native GK parameter preparation recovery

The five real multisite preparations reached AMOEBA matching but all failed in
OpenMM's unused GK-system builder: its default Bondi radius table has no Ca
(atomic number20). The first earlier preflight also exposed the real selected
metal source ID versus the single-metal adapter's `metal` alias. Both failures
are preserved in multisite_amoeba_capability_v1/v2; no energies ran.

Use the same physical five cases, protein/ACE/water/backgroundCa inventory and
coordinates from MULTISITE_AMOEBA_PREPARATION_PLAN. New preparation only:
parameterize permanent multipoles/polarization with pinned AMOEBA2018 XML,
explicitly leave OpenMM GK parameters unavailable, then read the actual GK
parameters using the already qualified native Tinker26.2 SOLUTE backend.
No dependency replacement, radius guess, monkey patch or changed solvent model.
Installed amoebabio18.prm type358 is Ca+2, charge2, polarizability0.5500A^3,
GK diameter3.6497A (cavity3.6670, neck0.1350). Retain its native parameters;
all other native settings are those of the established framework preflight.

Five parameterizations plus five native parameter-only reads. Zero energy,
force, DFT, MACE or response calls. Selected QM metal stays in physical ledger
and absent only from this framework; source/core boundary remains pending.
All background Ca are present/active with actual native polarizability/radii.
Compare native permanent parameters/coordinates against OpenMM/source records;
validate ACE bond, water inventory, charge closure and every physical atom.
No score or model-quality inference follows from this preparation test.

New products multisite_amoeba_capability_v3; failure/retry receipts retained.
This is a compatibility recovery within the authorized multisite preparation
scope, recorded before new native calls. Existing single-metal defaults and
all frozen GGR science remain unchanged.
