#!/bin/bash
# Tier 1: Colin's PQQ-La/Ca controls. Normalize AF3 CIF → PDB, protonate, carve_with_pqq, submit.
set -uo pipefail

ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
COLIN=/groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone
FEP_PYTHON=/home/jwestrob/miniconda3/envs/fep/bin/python
SCAN_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python

submitted=0; skipped=0
for tier in La-verified Ca-verified; do
  prefix="colinpqq_la"
  [ "$tier" = "Ca-verified" ] && prefix="colinpqq_ca"
  for cif in "$COLIN/$tier"/*.cif; do
    [ -f "$cif" ] || continue
    base=$(basename "$cif" .cif)
    # Sanitize: replace dots, hyphens with underscore
    safe=$(echo "$base" | tr -c '[:alnum:]_' '_' | sed 's/__*/_/g' | sed 's/_$//')
    stem="${prefix}_${safe}"
    out_dir="$ALCH/${stem}_qm"

    if [ -f "$out_dir/sp_${stem}_La.out" ] && grep -q "FINAL SINGLE POINT ENERGY" "$out_dir/sp_${stem}_La.out" 2>/dev/null; then
      echo "[DONE] $stem already processed"
      skipped=$((skipped+1))
      continue
    fi

    mkdir -p "$out_dir"
    echo ""
    echo "=== $stem ==="
    echo "  src: $cif"

    # 1. Normalize AF3 CIF → PDB
    pdb_norm="$out_dir/${stem}_normalized.pdb"
    if [ ! -f "$pdb_norm" ]; then
      "$SCAN_PYTHON" "$ALCH/scripts/normalize_af3_cif.py" "$cif" "$pdb_norm" 2>&1 | tail -2 || {
        echo "[FAIL] normalize"; continue
      }
    fi

    # 2. Protonate
    pdb_proton="$out_dir/${stem}_protonated.pdb"
    if [ ! -f "$pdb_proton" ]; then
      "$FEP_PYTHON" "$ALCH/scripts/protonate_cif.py" "$pdb_norm" "$pdb_proton" 2>&1 | tail -2 || {
        echo "[FAIL] protonation"; continue
      }
    fi

    # 3. Carve with PQQ.  Find what chain La actually ended up in.
    la_chain=$(grep "^HETATM" "$pdb_proton" | awk '$4=="LA"{print $5; exit}')
    if [ -z "$la_chain" ]; then
      la_chain="B"  # default
    fi

    if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
      "$SCAN_PYTHON" "$ALCH/scripts/carve_with_pqq.py" "$pdb_proton" "$out_dir" \
        --stem "$stem" --metal-chain "$la_chain" --metal-resname "LA" 2>&1 | tail -10 || {
        echo "[FAIL] carve"; continue
      }
    fi

    # 4. Patch inputs (NoAutostart already in carve_with_pqq output)
    submit="$out_dir/submit_${stem}.sh"
    if [ ! -f "$submit" ]; then
      echo "[FAIL] no submit script"; continue
    fi

    jid=$(sbatch --parsable "$submit" 2>&1)
    if [[ "$jid" =~ ^[0-9]+$ ]]; then
      echo "[OK] $stem -> job $jid"
      submitted=$((submitted+1))
    else
      echo "[FAIL] sbatch: $jid"
    fi
  done
done

echo ""
echo "Submitted: $submitted | Skipped: $skipped"
squeue -u jwestrob -h -t RUNNING -o "%i %j %M" 2>/dev/null | grep -vE "alphap|taxon_scan" | wc -l
echo "  ^our running jobs after submit"
