#!/bin/bash
# Batch overnight DFT discriminator panel.
# For each candidate: protonate CIF (PDBFixer) → carve → submit SLURM job.
set -uo pipefail

ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
CIF_BASE=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/results/overnight_v2/foldseek_cifs
FEP_PYTHON=/home/jwestrob/miniconda3/envs/fep/bin/python
SCAN_PYTHON=/home/jwestrob/miniconda3/envs/lanm_qmmm/bin/python
[ -x "$SCAN_PYTHON" ] || SCAN_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
[ -x "$SCAN_PYTHON" ] || SCAN_PYTHON=python3

# id  subdir  basename
declare -a CANDS=(
  "agmatinase_A0A2P8QC45  phase1_v2  A0A2P8QC45.cif"
  "arginase_A0A6N9AP46    phase1_v2  A0A6N9AP46.cif"
  "cytoc_A0A1Q7R3F2       phase1_v2  A0A1Q7R3F2.cif"
  "midas_A0A233S232       phase1_v2  A0A233S232.cif"
  "ferncalmod_A0A8T2TXN1  ferns      A0A8T2TXN1.cif"
  "m24pep_A0A2E7SM25      phase1_v2  A0A2E7SM25.cif"
  "adendeam_A0A950E6K5    phase1_v2  A0A950E6K5.cif"
  "aminopepP_A0A432JR12   phase1_v2  A0A432JR12.cif"
  "fernperox2_A0A9D4U0H5  ferns      A0A9D4U0H5.cif"
  "prolidase_A0A2V5ZP60   phase1_v2  A0A2V5ZP60.cif"
  "novelfold_A0A945JYC2   phase1_v2  A0A945JYC2.cif"
  "fernparv_A0A8T2SI74    ferns      A0A8T2SI74.cif"
  "rtxbroll_A0A3M2BIM9    phase2_ted AF-A0A3M2BIM9-F1-model_v4_TED02.cif"
)

for line in "${CANDS[@]}"; do
  set -- $line
  stem=$1; subdir=$2; cif_name=$3
  src_cif="$CIF_BASE/$subdir/$cif_name"
  out_dir="$ALCH/${stem}_qm"

  if [ ! -f "$src_cif" ]; then
    echo "[SKIP] $stem: source CIF missing: $src_cif"
    continue
  fi
  if [ -f "$out_dir/sp_${stem}_La.out" ] && grep -q "FINAL SINGLE POINT ENERGY" "$out_dir/sp_${stem}_La.out" 2>/dev/null; then
    echo "[DONE] $stem: already has La SP, skipping"
    continue
  fi

  mkdir -p "$out_dir"
  echo ""
  echo "=========================================="
  echo "[$stem] from $src_cif"
  echo "=========================================="

  # 1. Protonate (CIF → PDB → PDBFixer at pH 7)
  pdb_proton="$out_dir/${stem}_protonated.pdb"
  if [ ! -f "$pdb_proton" ]; then
    "$FEP_PYTHON" "$ALCH/scripts/protonate_cif.py" "$src_cif" "$pdb_proton" || {
      echo "[FAIL] protonation"; continue
    }
  fi

  # 2. Carve QM cluster
  if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
    "$SCAN_PYTHON" "$ALCH/scripts/carve_generic.py" "$pdb_proton" "$out_dir" --stem "$stem" || {
      echo "[FAIL] carve"; continue
    }
  fi

  # 3. Add water SP input (carve_generic doesn't include this)
  if [ ! -f "$out_dir/sp_${stem}_water.inp" ]; then
    cp "$ALCH/scripts/bulk_water.xyz" "$out_dir/bulk_water.xyz" 2>/dev/null || cat > "$out_dir/bulk_water.xyz" <<EOF
3
H2O reference
O   0.000000  0.000000  0.000000
H   0.756950  0.000000  0.585822
H  -0.756950  0.000000  0.585822
EOF
    cat > "$out_dir/sp_${stem}_water.inp" <<EOF
! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
* xyzfile 0 1 bulk_water.xyz
EOF
  fi

  # 4. Patch submit script to include water SP and use NoAutostart
  submit="$out_dir/submit_${stem}.sh"
  if [ -f "$submit" ]; then
    # Make sure water is in the kind list
    if ! grep -q "for kind in La Ca apo water" "$submit"; then
      sed -i 's/for kind in La Ca apo;/for kind in La Ca apo water;/' "$submit"
    fi
    # Add NoAutostart to inp files if not present
    for inp in "$out_dir"/sp_${stem}_*.inp; do
      [ -f "$inp" ] || continue
      if ! grep -q "NoAutostart" "$inp"; then
        sed -i 's/r2SCAN-3c CPCM/r2SCAN-3c NoAutostart CPCM/' "$inp"
      fi
    done
  fi

  # 5. Submit
  jid=$(sbatch --parsable "$submit" 2>&1)
  if [[ "$jid" =~ ^[0-9]+$ ]]; then
    echo "[OK] $stem submitted as job $jid"
  else
    echo "[FAIL] sbatch: $jid"
  fi
done

echo ""
echo "=========================================="
echo "Final queue state:"
squeue -u jwestrob -o "%.10i %.20j %.10T %.10M %R" 2>/dev/null | head -25
