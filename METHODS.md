# Methods — theory and chemistry of the vertical-swap Ca²⁺/Ln³⁺ discriminator

This document explains what the discriminator does, why each design choice was
made, and which assumptions are baked into the numbers. For pipeline mechanics
see [PIPELINE.md](PIPELINE.md); for classification labels see
[TIERS.md](TIERS.md); for validation evidence see [VALIDATION.md](VALIDATION.md).

---

## 1. The discriminator equation

For a candidate metal-binding pocket modelled with La placed in the site:

```
ΔΔE = ( E[Ca-cluster] − E[La-cluster] ) − DIFF_AQUO

DIFF_AQUO = E[Ca·6H₂O]²⁺ − E[La·6H₂O]³⁺ = −646.064555551745 Ha
            (at r²SCAN-3c / CPCM(Water) / DefGrid3)
```

`ΔΔE` is reported in kcal/mol after multiplication by `HA2KCAL = 627.5095`.

**Sign convention.** `ΔΔE > 0` ⇒ the pocket prefers Ln³⁺ over Ca²⁺.
`ΔΔE < 0` ⇒ the pocket prefers Ca²⁺.

**Why a vertical swap.** The Ca-form energy is computed on the *same* nuclear
coordinates as the La-form: only the metal element is changed (and the basis
set / charge bookkeeping adjusted to match). No re-optimisation, no
re-folding. This is the *Born–Haber* style cycle that lets the aquo reference
cancel several systematic errors cleanly:

- **PDBFixer non-determinism.** PDBFixer picks a tautomer / hydrogen placement
  on each run, and absolute SP energies can drift by 15–66 kcal/mol between
  protonation rounds on the same input. Because both metals see identical
  geometry, this noise cancels in `E[Ca] − E[La]`. We measured the residual
  at <0.5 kcal/mol on identical structures.
- **Basis-set incompleteness.** r²SCAN-3c is composite — TZVP-quality with
  basis-set superposition correction baked in — but no finite basis is
  complete. Identical-coordinate, identical-basis Ca and La SPs share the
  same incompleteness error, which subtracts out.
- **ECP vs all-electron offset.** La uses a 46-core (def2-ECP) plus def2-TZVP
  small-core valence basis; Ca is all-electron. The two metals therefore
  carry different absolute zero references. Subtracting DIFF_AQUO — which is
  the same ECP/AE offset evaluated on idealised aquo ions — removes it to
  the same accuracy as the aquo reference itself.

The result is a difference-of-differences that is **insensitive to absolute
energy noise** but **sensitive to the chemistry of the carboxylate cluster
the metal sees**.

See [VALIDATION.md](VALIDATION.md) for cross-run reproducibility (Marco family
paralogs <0.05 kcal/mol; Colin's 00_Lav rerun 0.2–3 kcal/mol).

---

## 2. The DIFF_AQUO reference

`DIFF_AQUO = −646.064555551745 Ha` is the difference of two single-point
energies on idealised six-coordinate aquo complexes:

- **[Ca(H₂O)₆]²⁺** at r²SCAN-3c / CPCM(Water) / DefGrid3 →
  `E_aquo_Ca = −1288.921041427522 Ha`
- **[La(H₂O)₆]³⁺** at r²SCAN-3c / CPCM(Water) / DefGrid3 (La with
  def2-ECP + def2-TZVP large-core, mult = 1) →
  `E_aquo_La = −642.856485875777 Ha`

These constants are hard-coded in
[`scripts/update_results_jsonl.py`](scripts/update_results_jsonl.py):

```python
E_AQUO_CA = -1288.921041427522
E_AQUO_LA =  -642.856485875777
DIFF_AQUO = E_AQUO_CA - E_AQUO_LA  # -646.064555551745 Ha
HA2KCAL   = 627.5095
```

**Do not change these without re-running the aquo references at the same
level of theory and updating every downstream consumer.** Class-threshold
calibration ([TIERS.md](TIERS.md)) is anchored on this number.

A fourth SP on a single bulk water molecule (`bulk_water.xyz`,
`* xyzfile 0 1`) is run for every candidate as a **sanity check**: it should
recover the canonical r²SCAN-3c/CPCM(Water) energy for H₂O to <0.001 Ha. If
it doesn't, the ORCA installation, basis library, or solvation parameters have
drifted and the panel needs a re-baseline.

---

## 3. Cluster carving rules

The QM cluster is carved by [`scripts/carve_generic.py`](scripts/carve_generic.py)
(or `carve_with_pqq.py` when a PQQ ligand is present). Inputs are a protonated
PDB (PDBFixer at pH 7); outputs are three XYZ files (`{stem}_La_qm.xyz`,
`{stem}_Ca_qm.xyz`, `{stem}_apo_qm.xyz`) plus a fourth `bulk_water.xyz`.

### Geometric cutoff

```
FIRST_SHELL_CUT = 3.0 Å  (tightened from 3.2 Å on 2026-05-11)
```

Any non-metal heavy atom (O, N, or S) within 3.0 Å of the metal counts as a
first-shell donor and pulls its parent residue into the QM region.

### Sidechain dict (allowed donor residues)

```python
SIDECHAIN_QM_ATOMS = {
    "ASP": ["CB", "HB2", "HB3", "CG", "OD1", "OD2"],
    "GLU": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2"],
    "ASN": ["CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22"],
    "GLN": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22"],
    "SER": ["CB", "HB2", "HB3", "OG", "HG"],
    "TYR": ["CB", "HB2", "HB3", "CG", "CD1", "HD1", "CE1", "HE1",
            "CZ", "OH", "HH", "CE2", "HE2", "CD2", "HD2"],
}
```

Donors are restricted to **ASP, GLU, ASN, GLN, SER, TYR** — six residue types
covering all the hard-O donors the discriminator was calibrated against.
TYR was added on 2026-05-11 to capture Tyr-OH coordination.

For each first-shell residue, the listed sidechain atoms are pulled into the
QM region and a hydrogen is inserted at the Cβ position (1.09 Å along the
Cβ→Cα vector) to cap the broken Cβ-Cα bond. This is a standard "link-H"
treatment for QM cluster carves.

### Backbone-O carving

When a residue's backbone carbonyl O lies within `FIRST_SHELL_CUT` of the
metal, the C=O group is carved separately (independent of whether the
residue's sidechain qualifies):

- The backbone C and O atoms are added to the QM region.
- Two link-H atoms cap the broken peptide bonds — one along the C→Cα
  direction, one along the bisector of the (−O, −Cα) directions (the
  approximate "where the next N would be" direction).

This captures GLY, PRO, and other non-sidechain-coordinated sites that the
original dict-only carver missed.

### Inner-shell waters

Water oxygens within `FIRST_SHELL_CUT` of the metal are included as full
H₂O groups (O + 2 H atoms identified by O-H bond length ≤ 1.2 Å).

### Carve outputs

Three carves are written, identical except for the metal:

- **La form.** `[La, sidechain atoms…, water atoms…, link Hs…]`
  Charge = `+3 − n_carbox`.
- **Ca form** (vertical metal swap). Same atoms, La replaced by Ca.
  Charge = (La charge) − 1.
- **Apo form.** Metal and inner waters removed; sidechain atoms + link Hs only.
  Charge = (La charge) − 3.

Plus the standalone `bulk_water.xyz` (1 O + 2 H, charge 0).

### Policy: what is *not* in the dict

CYS, HIS, MET, LYS, ARG are **deliberately excluded** from
`SIDECHAIN_QM_ATOMS`. This is a methodology choice, not an oversight:

- The discriminator was calibrated on hard-O carboxylate / amide / hydroxyl
  donors (the same residue set that biology uses for high-affinity Ca and Ln
  binding). The XoxF/MxaF gold standard, LanM EF1, tannase, and SilE all bind
  through this set.
- Adding S-donors (CYS SG, MET SD) or N-aromatic donors (HIS NE2/ND1) would
  change the *chemistry being measured*. They favour transition metals over
  hard f-block ions, and including them in the QM region would mix
  Ln-vs-transition-metal selectivity into the Ca-vs-Ln signal we want.
- Long-chain charged sidechains (LYS NZ, ARG NH1/NH2) are similarly excluded:
  they form salt bridges rather than first-shell metal donors at biologically
  realistic geometries.

When a re-carve still produces a 1-atom cluster after the TYR + backbone-O
fixes, it is flagged as `SOLVENT_EXCLUDE` (closest donor > 4 Å — genuinely
empty pocket) or `CARVE_AMBIGUOUS` (closest donor is an unsupported residue)
rather than auto-expanded. See [TIERS.md](TIERS.md).

An attempt to add THR/CYS/HIS/MET/LYS/ARG to the dict was made on 2026-05-11
and reverted out per Jacob's explicit direction. This is recorded in the
script docstring of `carve_generic.py` and in [CHANGELOG.md](CHANGELOG.md).

---

## 4. PQQ-aware carving

PQQ (pyrroloquinoline quinone, 14 C + 8 O + 2 N = 24 heavy atoms) is the
primary metal coordinator in PQQ-MDH/ADH active sites: PQQ-O5 and PQQ-O4 form
a tight bidentate clamp on the metal. A protein-only carve of a PQQ-binder
systematically underestimates ΔΔE by ≈ 10 kcal/mol because the carve drops
the dominant donor pair.

### Detection

[`scripts/process_inbox.sh`](scripts/process_inbox.sh) inspects each CIF for a
residue named `LIG_*` with elemental composition:

```
≥ 12 C   AND   ≥ 6 O   AND   ≥ 1 N
```

(canonical PQQ is 14 C + 8 O + 2 N; the lenient cut-offs accommodate
Protenix/AF3 ligand-naming variation and missing-H representations).

When detected, the CIF is routed through the **tier-1 pipeline**:

```
CIF
  → normalize_af3_cif.py     (LIG_* → "LA"/"PQQ" residue renames)
  → protonate_cif.py         (PDBFixer at pH 7)
  → carve_with_pqq.py        (24 PQQ atoms enter the QM region)
```

### PQQ charge bookkeeping

The 24 PQQ heavy atoms enter the QM cluster with **formal charge −2**
(`PQQ_CHARGE = -2` in
[`scripts/carve_with_pqq.py`](scripts/carve_with_pqq.py)). This was set
empirically by an electron-parity calibration — at −2 PQQ, both the XoxF
carve (3 ASP/GLU + 1 ASN protein donors) and the MxaF carve (1 GLU + 1 ASN
protein donors) give even electron counts with La³⁺ or Ca²⁺. Chemically, −2
corresponds to PQQ with one of its three carboxylates protonated (or an
equivalent N-H tautomer). Without explicit ligand H atoms in the Protenix
output this is the cleanest closed-shell approximation available.

The Ca/Ln *selectivity* (ΔΔE) is **insensitive to this choice** because both
metal forms see the same PQQ. Absolute energies shift, but the difference
does not.

### Typical PQQ-aware cluster

55–58 atoms: metal + 2–3 carboxylates + 1 ASN + 24-atom PQQ + 3–4 link Hs.

---

## 5. Electronic structure

```orca
! r²SCAN-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
%basis
  NewECP La "def2-ECP" end
  NewGTO La "def2-TZVP" end
end
* xyzfile {charge} 1 {stem}_{kind}_qm.xyz
```

### Why r²SCAN-3c

[r²SCAN-3c](https://pubs.aip.org/aip/jcp/article/154/6/064103) (Grimme group,
2021) is a composite DFT method built around the meta-GGA r²SCAN functional,
a custom small-basis-set adapted for r²SCAN, a Becke–Johnson D4 dispersion
correction, and a geometric counterpoise correction. It is **basis-set-corrected
out of the box**, so we get TZVP-quality energetics on charged metal clusters
without paying for an explicit TZVP basis on every job. r²SCAN itself is
SCAN-fixed (numerically stable, no SCAN's grid-sensitivity pathology), and
the composite gives reasonable transition-metal and lanthanide geometries
without re-parameterisation. For our purposes — comparing the same nuclear
configuration with two different metals — it delivers reproducibility at
<0.05 kcal/mol while keeping wall-clock at 5–15 min per SP on 64 OpenMP
threads. We considered pure r²SCAN (more expensive, no basis correction),
PBE-D3 (cheaper but ECP-La geometries less stable), and B97-3c (next-cheapest
composite; cross-validated in [VALIDATION.md §7](VALIDATION.md)). r²SCAN-3c
won on the cost/accuracy curve.

### La basis set and ECP

```
NewECP La "def2-ECP" end           # 46-electron core (Karlsruhe def2-ECP)
NewGTO La "def2-TZVP" end          # TZVP valence basis on La
```

This is a **large-core ECP**: the 46-core treatment puts the 4f electrons in
the core (4f⁰ for La³⁺ is closed-shell singlet, multiplicity 1). Small-core
ECPs (28-core, 4f explicit) are available but would require re-validation
against the Tier-1 controls and are not needed at the discriminator-level
question we are asking. We did test small-core in a pre-discriminator pilot;
it added 4× cost without changing the ΔΔE sign on the canonical pair.

### Ca treatment

Ca is **all-electron** (def2-TZVP via the composite default). Closed-shell
singlet. Multiplicity 1.

### Solvation

```
CPCM(Water)
```

Implicit conductor-like polarizable continuum, water dielectric.
**Required** — gas-phase ΔΔE on charged clusters is wildly overstabilized for
the higher-charge state and gives unphysically large positive ΔΔE values
across the panel. CPCM cancels most of this in the difference but
the dielectric is still needed inside the SCF for stable convergence.

### Integration grid

```
DefGrid3
```

ORCA's "tight" default integration grid. Lower grids (1, 2) introduce ~1
kcal/mol noise that does not cancel cleanly in ΔΔE. We did not see further
ΔΔE convergence going to higher grids.

### NoAutostart

```
NoAutostart
```

Forces a fresh SCF on every run, preventing ORCA from picking up a stale
`.gbw` from a different system. (Workspace dirs accumulate Ca/La/apo/water
runs; we do not want La's wavefunction guess seeded by the prior Ca run on
the previous candidate.)

### Multiplicity

```
* xyzfile {charge} 1 ...     # mult = 1 throughout
```

La³⁺ (f⁰) and Ca²⁺ (d⁰) are both closed-shell singlets. We do not run open-shell
mid-Ln panels (Nd, Sm, Eu, Tb, Dy, …) in the production discriminator — those
were explored in earlier work (`qmmm/cluster_panel/`) and are not part of the
Ca-vs-Ln decision.

---

## 6. Class thresholds

The class assignment in
[`scripts/update_results_jsonl.py`](scripts/update_results_jsonl.py) is:

```python
def classify(ddE_kcal, n_atoms):
    if ddE_kcal is None:                       return "pending"
    if n_atoms < 15:                           return "EXCLUDE (under-carved)"
    if abs(ddE_kcal) > 30 and n_atoms < 50:    return "OUTLIER (...)"
    if ddE_kcal >= 20:   return "Ln-evolved"
    if ddE_kcal >= 5:    return "Ln-preferring"
    if ddE_kcal >= 0:    return "marginal"
    if ddE_kcal >= -5:   return "ambiguous"
    return "Ca-evolved"
```

The thresholds were locked at the 2026-05-07 panel review against the
validation set ([VALIDATION.md](VALIDATION.md)):

- **±5 kcal/mol noise band** spans `marginal` (0..+5) and `ambiguous` (−5..0).
  Anchored on the PDBFixer-noise-floor measurement and on cross-run agreement
  between independent reruns.
- **+20 / Ln-evolved threshold** chosen so the XoxF +24 anchor lands in the
  tier above the bulk of Ca-binding proteins that score Ln-preferring at the
  chemistry level (calmodulin +18, parvalbumin +11, calexcitin +11).
- **OUTLIER tier (|ΔΔE| > 30 on a < 50-atom cluster)** is the Mn²⁺/Zn²⁺
  small-acidic-pocket signature: SCF gets cornered into a state that can't
  accommodate La's larger ionic radius geometrically, and the energy gap
  blows up. Treat as biology, not numerical failure.
- **SEVERE asymmetric carve (|ΔΔE| > 100)** is a manual flag for atom-count
  mismatch between the Ca and La carves — typically a PQQ inclusion mismatch.
  Not auto-labelled; diagnose by diffing the two XYZ files.

See [TIERS.md](TIERS.md) for tier semantics and example stems.

---

## 7. What the discriminator does NOT measure

This is the most important paragraph in this file. ΔΔE is a **coordination-chemistry
score** at fixed nuclear geometry. It is *necessary but not sufficient* for
calling a pocket an in-vivo Ln-binder.

The discriminator does **not** measure:

- **In-vivo selectivity.** Cytoplasmic [Ca²⁺] ≈ 1 mM; environmental [Ln³⁺] is
  trace nanomolar at best. Most canonical Ca-binding proteins (calmodulin,
  calexcitin, calbindin, parvalbumin, calcineurin B, RTX β-roll) score
  Ln-preferring at the chemistry level (+5 to +20 band) but bind Ca²⁺ in
  vivo solely because of concentration ratios. This is the *same chemistry*
  that powers Tb³⁺ / Eu³⁺ luminescence probes — those probes outcompete Ca²⁺
  at canonical Ca-binding sites in vitro because the carboxylate cluster
  intrinsically prefers the harder f-block cation. The probe biochemistry is
  the chemistry-level analogue of what the discriminator reports.
- **Kinetics.** ΔΔE is a static energy difference. Binding on-rates, off-rates,
  conformational barriers, and induced-fit are all outside the scope. A
  positive ΔΔE pocket may not be kinetically accessible to Ln in vivo even
  if it is thermodynamically preferred at that snapshot.
- **Absolute binding affinity.** ΔΔE is a *difference* of binding energies,
  not an absolute. A +20 kcal/mol Ln-evolved pocket and a +5 kcal/mol
  Ln-preferring pocket might both bind Ln tightly in absolute terms; the
  discriminator says one is selective over Ca²⁺ by ~30× more energy than the
  other, not that one binds 10⁴× tighter.
- **Pocket function.** Whether a Ln-binding site is catalytic (XoxF-like
  PQQ-MDH activity), structural (LanM-like sensor), regulatory (DUF4394
  hypothesized chaperone), or accessory (uptake / storage / efflux). All of
  those distinctions live in the protein context — Pfam, KO, operon, genome
  neighborhood — not in the cluster energetics.

A positive ΔΔE filters out the negative space (Ca-evolved pockets) cleanly,
and lights up the chemistry-positive space. Biology then has to do the
selection.

---

## 8. Source material

The methodology arguments above are condensed from:

- **2026-05-07 vault note** `~/jwestrob/obsidian-vault/agent-captures/2026-05-07_ca-ln-dft-discriminator.md` — the original validation-and-thresholds writeup.
- **2026-05-11 vault note** `~/jwestrob/obsidian-vault/agent-captures/2026-05-11_dft-discriminator-state-of-pipeline.md` — compaction-resilient bootstrap; carve-dict policy and PQQ-aware routing rationale.
- [CONTEXT.md](CONTEXT.md), [HOWTO.md](HOWTO.md), and the docstrings of the
  scripts in [scripts/](scripts/).
