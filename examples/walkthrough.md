# Walkthrough — `a0a423vtm4_domain1` (PQQ-ADH, Ln-preferring)

A worked example through the full pipeline on a single protein, using real
files from the workspace. All numbers below are pulled from
`03_Euk_a0a423vtm4_domain1chaina-pqq-la_model_qm/` and
`results/colin_euk_pqq_adh_results.tsv`.

The candidate: **A0A423VTM4 domain1**, a eukaryotic PQQ-ADH from Colin's
`03_Euk` batch. Fold by Protenix with iPTM(La) = 0.98, ptm(protein) = 0.95,
ranking score 0.97 — high model confidence. ΔΔE = **+18.70 kcal/mol**, class
**Ln-preferring**. Sits just below the Ln-evolved threshold (+20). Picked for
this walkthrough because it is a clean PQQ-bearing case, the full output set is
intact, and the fold confidence is rock-solid.

For the methodology see [METHODS.md](../METHODS.md); for class semantics
[TIERS.md](../TIERS.md); for the full pipeline reference
[PIPELINE.md](../PIPELINE.md).

---

## Step 1 — Source CIF

```
/groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone/Euk/a0a423vtm4_domain1chaina-pqq-la_model.cif
```

This is a Protenix model (La placed in the candidate pocket plus the PQQ
cofactor encoded as a `LIG_*` residue). 55 atoms total in the QM region
after carving.

---

## Step 2 — Drop in the inbox

```bash
ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs

# Symlink-style staging — no disk duplication
ln -s /groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone/Euk/a0a423vtm4_domain1chaina-pqq-la_model.cif "$ALCH/inbox/"
```

(In the production run, the file was dropped through Colin's larger Euk
batch staging — the symlink is shown here as a single-CIF equivalent.)

---

## Step 3 — Run `process_inbox.sh`

```bash
bash $ALCH/scripts/process_inbox.sh
```

The driver:

1. **Iterates** `inbox/*.cif`, sanitises the stem to
   `03_Euk_a0a423vtm4_domain1chaina-pqq-la_model`.
2. **Checks for prior completion** — no `sp_<stem>_La.out` exists yet, so it
   proceeds.
3. **Filters for La**. The CIF contains `LA1` HETATM records, so the relaxed
   regex (`HETATM[[:space:]]+[0-9]+[[:space:]]+La`) matches.
4. **Detects PQQ**. The `LIG_C` residue in the CIF contains
   **14 C + 8 O + 2 N** = 24 heavy atoms — passes the
   ≥12 C + ≥6 O + ≥1 N test. **Routes through the tier-1 (PQQ-aware) pipeline.**

Inline output:
```
=== 03_Euk_a0a423vtm4_domain1chaina-pqq-la_model from a0a423vtm4_domain1chaina-pqq-la_model.cif ===
  PQQ detected → tier1 pipeline
```

---

## Step 4 — Normalise + protonate + carve (tier-1 pipeline)

### Normalise (`normalize_af3_cif.py`)

```
normalized: 03_Euk_a0a423vtm4_domain1chaina-pqq-la_model_normalized.pdb
            La=1  PQQ=1
```

Renames `LIG_A` (containing La) → residue `LA` with atom `LA`, and `LIG_C`
(PQQ-like composition) → residue `PQQ`. Output: `<stem>_normalized.pdb`.

### Protonate (`protonate_cif.py` in the `fep` env)

PDBFixer adds hydrogens at pH 7. Output: `<stem>_protonated.pdb`.

### Carve (`carve_with_pqq.py`)

```
metal: chain B LA  pos=(-3.80,-7.35,13.29)
First-shell donors (within 3.2 Å):
  PQQ /O5  d≈2.4
  PQQ /O4  d≈2.6
  ASP /OD1, OD2  (×3 carboxylates from the protein)
  ASN /OD1  (×1)
First-shell residues (protein): [...]
PQQ present: True
PQQ atoms included: 24
QM region: 55 atoms (3 carboxylates, 1 PQQ, 4 link Hs)
Charge: +3 (metal) -3 (carb) -2 (PQQ) = -2
```

Three XYZ files written:

- `<stem>_La_qm.xyz` — 55 atoms (1 La + 3 carboxylates + 1 Asn + 24-atom PQQ + 4 link Hs). Charge `-2`.
- `<stem>_Ca_qm.xyz` — same atoms, La replaced by Ca. Charge `-3`.
- `<stem>_apo_qm.xyz` — metal removed; 54 atoms. Charge `-5`.

Plus `bulk_water.xyz` (1 O + 2 H).

Element composition of the La carve (from `awk 'NR>2 {print $1}' ...xyz | sort | uniq -c`):

```
22 C    14 H    1 La    3 N    15 O
```

---

## Step 5 — Four ORCA single-points

`process_inbox.sh` adds the water SP (if not already there), patches every
`sp_*.inp` with `NoAutostart` (defensive), and submits the SLURM script.

Each input looks like:

```orca
! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
%basis
  NewECP La "def2-ECP" end
  NewGTO La "def2-TZVP" end
end
* xyzfile -2 1 03_Euk_a0a423vtm4_domain1chaina-pqq-la_model_La_qm.xyz
```

(For Ca, charge becomes `-3` and the `%basis` block is dropped — Ca is all-electron at the composite default.)

The SLURM job runs four SPs sequentially on a 64-core memory-partition node
with `--exclusive`. Total wallclock for this candidate: roughly 30–60 min.

### Energies (final, from `sp_<stem>_*.out`)

| SP    | charge | mult | E (Ha)                  |
|-------|-------:|-----:|------------------------:|
| La    | −2     | 1    | −2170.861424067753      |
| Ca    | −3     | 1    | −2816.896179035692      |
| apo   | −5     | 1    | −2139.318831017993      |
| water | 0      | 1    | −76.428549251187        |

All four converged with `ORCA TERMINATED NORMALLY`. The bulk-water energy
recovers the canonical r²SCAN-3c / CPCM(Water) reference for H₂O — sanity
check passed.

---

## Step 6 — ΔΔE computation

```
ΔΔE  =  ((E_Ca − E_La) − DIFF_AQUO) × 627.5095
      = ((−2816.896179035692 − (−2170.861424067753)) − (−646.064555551745)) × 627.5095
      = (−646.034754967939 − (−646.064555551745)) × 627.5095
      = +0.029800583806 Ha × 627.5095
      = +18.7001 kcal/mol
```

The aggregator records this as `ddE_kcal: 18.70` for stem
`03_Euk_a0a423vtm4_domain1chaina-pqq-la_model` in
`results/all_results.jsonl`.

---

## Step 7 — Classification

ΔΔE = **+18.70 kcal/mol**, n_atoms = 55. Following the rules in
[METHODS.md §6](../METHODS.md):

- `n_atoms = 55 ≥ 15` → not EXCLUDE.
- `|ΔΔE| = 18.70 ≤ 30` → not OUTLIER.
- `18.70 < +20` → not Ln-evolved.
- `18.70 ≥ +5` → **Ln-preferring**.

So the entry lands in the **Ln-preferring** tier, very near the Ln-evolved
threshold. Recorded in `results/colin_euk_pqq_adh_results.tsv` at rank 3 of
the Euk batch:

```
rank=3  uniprot=A0A423VTM4  domain=domain1  ddE_kcal=+18.70  class=Ln-preferring
iptm_La_fold=0.980  ptm_protein=0.950  ranking_score=0.970  n_atoms_QM=55
```

---

## Step 8 — Sanity check

A consistent Ln-preferring PQQ-ADH active site should show:

- **CN ≈ 7** for La: in this carve, 2 PQQ oxygens (O5 + O4 bidentate to the
  metal) + 3 protein carboxylate oxygens (one bidentate ASP, plus one OD1
  from each of two other ASPs) + 1 ASN amide O = **7 first-shell donors**.
- **Fold confidence:** iPTM(La) = 0.98 says Protenix is confident in the La
  placement. iPTM(protein) = 0.95 says the fold itself is solid.
- **n_atoms = 55** is in the typical PQQ-aware carve range (55–58 atoms) —
  metal + 3 carboxylates + 1 ASN + 24 PQQ + 4 link Hs.
- **Charge balance:** +3 (La) − 3 (3 carboxylates, each −1) − 2 (PQQ) = **−2**,
  which matches the `* xyzfile -2 1` line in the La input. The Ca form
  subtracts 1 from the metal charge → `-3`, matching its input.
- **Water SP sanity:** −76.428549 Ha. The expected r²SCAN-3c / CPCM(Water)
  energy for a single H₂O is within milli-Hartree of this; no SCF drift.

All consistent. The ΔΔE of +18.70 is a real Ln-preferring chemistry signal
from a clean carve on a high-confidence fold.

---

## Looking it up in the JSONL

```bash
grep -F '"03_Euk_a0a423vtm4_domain1chaina-pqq-la_model"' \
  results/all_results.jsonl \
  | jq .
```

Returns the full row:

```json
{
  "stem": "03_Euk_a0a423vtm4_domain1chaina-pqq-la_model",
  "source": "unknown",
  "n_atoms_la": 55,
  "n_carbox": 3,
  "charge_la": -2,
  "E_La_Ha":  -2170.861424067753,
  "E_Ca_Ha":  -2816.896179035692,
  "E_apo_Ha": -2139.318831017993,
  "E_water_Ha": -76.428549251187,
  "ddE_kcal": 18.7001,
  "class": "Ln-preferring",
  "term_la": true,
  "term_ca": true,
  "term_apo": true,
  "term_water": true,
  "cluster_panel_path": "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/03_Euk_a0a423vtm4_domain1chaina-pqq-la_model_qm"
}
```

(`source` shows `unknown` because the `03_Euk_*` naming predates the
source-detection regex in `update_results_jsonl.py`. The biology source is
"Colin Euk batch.")

---

## Takeaway

A0A423VTM4 domain1 is a textbook walk through the production pipeline:
PQQ-aware routing fires, all four SPs converge cleanly, the carve has the
expected size and composition, the ΔΔE math matches by hand, and the class
assignment is unambiguous. It sits just below the Ln-evolved threshold,
which is consistent with its biology — a eukaryotic PQQ-ADH that
*could* be a Ln-using oxidoreductase given environmental Ln availability,
worth following up biochemically (see [RESULTS.md §7](../RESULTS.md)).
