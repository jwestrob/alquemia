#!/bin/bash
# Batch v3: Alphaproteo elite + Marco's β-roll family.
# 8 new candidates to add to the 39-protein panel.
set -uo pipefail

ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
FEP_PYTHON=/home/jwestrob/miniconda3/envs/fep/bin/python
SCAN_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python

# Format: stem  cif  [site-chain]  [site-resnum]
# For ETR2 with 2 La, we carve at each site separately (chain B = EF1, chain C = EF2)
declare -a CANDS=(
  "alphap_etr2_ef1     /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/elite_la_folds_v2/output_etr2_monomer_la2/A0A127ETR2_mono_la2/seed_101/predictions/A0A127ETR2_mono_la2_rank_0.cif         B 1"
  "alphap_etr2_ef2     /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/elite_la_folds_v2/output_etr2_monomer_la2/A0A127ETR2_mono_la2/seed_101/predictions/A0A127ETR2_mono_la2_rank_0.cif         C 1"
  "alphap_ivl5_excalibur  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/elite_la_folds/output/A0A2E2IVL5/seed_101/predictions/A0A2E2IVL5_rank_0.cif"
  "alphap_arginase    /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/elite_la_folds/output/A0A328AMX2/seed_101/predictions/A0A328AMX2_rank_0.cif"
  "alphap_aminopep_p  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/elite_la_folds/output/A0A246SZW5/seed_101/predictions/A0A246SZW5_rank_0.cif"
  "marco_akiy_44      /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/candidate_bundle/marco_beta_roll/NZ_AKIY01000241.1_44_BradyrhizSp000282615.cif"
  "marco_bsox_25      /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/candidate_bundle/marco_beta_roll/NZ_BSOX01000478.1_25_BliaoNingense.cif"
  "marco_cp030053_2053 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/candidate_bundle/marco_beta_roll/NZ_CP030053.1_2053_4thHomolog.cif"
)

submitted=0
for line in "${CANDS[@]}"; do
  read -ra parts <<< "$line"
  stem="${parts[0]}"
  src_cif="${parts[1]}"
  site_chain="${parts[2]:-}"
  site_resnum="${parts[3]:-}"

  if [ ! -f "$src_cif" ]; then
    echo "[SKIP] $stem: missing CIF"
    continue
  fi
  out_dir="$ALCH/${stem}_qm"
  if [ -f "$out_dir/sp_${stem}_La.out" ] && grep -q "FINAL SINGLE POINT ENERGY" "$out_dir/sp_${stem}_La.out" 2>/dev/null; then
    echo "[DONE] $stem: already has La SP, skipping"
    continue
  fi
  mkdir -p "$out_dir"
  echo ""
  echo "=========================================="
  echo "[$stem] from $(basename $src_cif)"
  if [ -n "$site_chain" ]; then echo "  site: chain $site_chain res $site_resnum"; fi
  echo "=========================================="

  # Need a fresh protonated PDB per stem (PDBFixer is shared, but each stem gets its own)
  pdb_proton="$out_dir/${stem}_protonated.pdb"
  if [ ! -f "$pdb_proton" ]; then
    "$FEP_PYTHON" "$ALCH/scripts/protonate_cif.py" "$src_cif" "$pdb_proton" 2>&1 | tail -3 || {
      echo "[FAIL] protonation"; continue
    }
  fi

  # Carve — pass site-chain/resnum if provided
  if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
    if [ -n "$site_chain" ]; then
      "$SCAN_PYTHON" "$ALCH/scripts/carve_generic.py" "$pdb_proton" "$out_dir" --stem "$stem" --site-chain "$site_chain" --site-resnum "$site_resnum" 2>&1 | tail -10 || {
        echo "[FAIL] carve"; continue
      }
    else
      "$SCAN_PYTHON" "$ALCH/scripts/carve_generic.py" "$pdb_proton" "$out_dir" --stem "$stem" 2>&1 | tail -10 || {
        echo "[FAIL] carve"; continue
      }
    fi
  fi

  # Add water SP
  if [ ! -f "$out_dir/sp_${stem}_water.inp" ]; then
    cat > "$out_dir/bulk_water.xyz" <<EOF
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

  submit="$out_dir/submit_${stem}.sh"
  if [ -f "$submit" ]; then
    if ! grep -q "for kind in La Ca apo water" "$submit"; then
      sed -i 's/for kind in La Ca apo;/for kind in La Ca apo water;/' "$submit"
    fi
    for inp in "$out_dir"/sp_${stem}_*.inp; do
      [ -f "$inp" ] || continue
      if ! grep -q "NoAutostart" "$inp"; then
        sed -i 's/r2SCAN-3c CPCM/r2SCAN-3c NoAutostart CPCM/' "$inp"
      fi
    done
  fi

  jid=$(sbatch --parsable "$submit" 2>&1)
  if [[ "$jid" =~ ^[0-9]+$ ]]; then
    echo "[OK] $stem submitted as job $jid"
    submitted=$((submitted+1))
  else
    echo "[FAIL] sbatch: $jid"
  fi
done

echo ""
echo "=========================================="
echo "Submitted: $submitted"
squeue -u jwestrob -o "%.10i %.30j %.10T %.10M" 2>/dev/null | head -15
