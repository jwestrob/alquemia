#!/bin/bash
# Batch v2 overnight DFT discriminator panel — round 2 (15 more candidates).
set -uo pipefail

ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
CB=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/candidate_bundle
PHASE1=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/results/overnight_v2/foldseek_cifs/phase1_v2
FEP_PYTHON=/home/jwestrob/miniconda3/envs/fep/bin/python
SCAN_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python

# stem  full_path
declare -a CANDS=(
  "hyphomyst_NZ_WMBQ      $CB/webserver_revalidated/NZ_WMBQ01000002.1_957_HyphomicrobiumMystery_1La_iPTM0p984.cif"
  "typeAorphan_NZ_JADPKR  $CB/webserver_revalidated/NZ_JADPKR010000108.1_381_NovelTypeAorphan_disorderToOrder_iPTM0p948.cif"
  "marcobetar_NZ_CP019948 $CB/marco_beta_roll/NZ_CP019948.1_943_BradyrhizUnnamed.cif"
  "duf882_NZ_FRXO         $CB/spici_lams_other/NZ_FRXO01000001.1_31_DUF882_TonBLocus.cif"
  "lanthhyd_NZ_JAKSYI     $CB/spici_lams_other/NZ_JAKSYI010000031.1_14_PG_NlpC_lanthanophoreHydrolase.cif"
  "cbb3oxid_NZ_RXIZ       $CB/spici_lams_other/NZ_RXIZ01000013.1_160_551aa_cbb3OxidaseLocus.cif"
  "lanmfusion_NZ_LFLZ     $CB/webserver_revalidated/NZ_LFLZ01000011.1_14_LanMfusion571aa_3La_chainIptm0p94.cif"
  "terbdiv_NZ_OBQD        $CB/spici_lams_other/NZ_OBQD01000021.1_7_TerB_divergentParalog.cif"
  "tannasepara_NZ_JYMT    $CB/tannase_paralog_controls/NZ_JYMT01000017.1_60_86pctID_LanMless_Ln3plus_likely.cif"
  "opaque_A0A4Q2Y9J7      $PHASE1/A0A4Q2Y9J7.cif"
  "phase1_A0A965Z5K8      $PHASE1/A0A965Z5K8.cif"
  "phase1_A0A933FB28      $PHASE1/A0A933FB28.cif"
  "phase1_A0A5M3W1H2      $PHASE1/A0A5M3W1H2.cif"
  "ferncellul_A0A8T2TYS0  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/results/overnight_v2/foldseek_cifs/ferns/A0A8T2TYS0.cif"
  "ferncellul2_A0A9D4V3A3 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/results/overnight_v2/foldseek_cifs/ferns/A0A9D4V3A3.cif"
)

submitted=0
for line in "${CANDS[@]}"; do
  set -- $line
  stem=$1; src_cif=$2

  if [ ! -f "$src_cif" ]; then
    echo "[SKIP] $stem: source CIF missing: $src_cif"
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
  echo "=========================================="

  # Verify CIF has La/Ce
  if ! grep -qE "^(ATOM|HETATM).*[ ]LA[ ]|[ ]CE[ ]|[ ]Y[ ]" "$src_cif"; then
    if ! grep -qE "LA |CE |Y " "$src_cif"; then
      echo "[SKIP] $stem: no La/Ce/Y in CIF"
      continue
    fi
  fi

  pdb_proton="$out_dir/${stem}_protonated.pdb"
  if [ ! -f "$pdb_proton" ]; then
    "$FEP_PYTHON" "$ALCH/scripts/protonate_cif.py" "$src_cif" "$pdb_proton" 2>&1 | tail -3 || {
      echo "[FAIL] protonation"; continue
    }
  fi

  if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
    "$SCAN_PYTHON" "$ALCH/scripts/carve_generic.py" "$pdb_proton" "$out_dir" --stem "$stem" 2>&1 | tail -10 || {
      echo "[FAIL] carve"; continue
    }
  fi

  # Add water SP input
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
echo "Final queue state:"
squeue -u jwestrob -o "%.10i %.20j %.10T %.10M %R" 2>/dev/null | head -25
