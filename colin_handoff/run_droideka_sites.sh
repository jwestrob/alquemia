#!/bin/bash
# Special handler for the Acidobacterium "droideka" structure with 6 La sites and NO PQQ.
# Carves at each La site separately (chains B/C/D/E/F/G in the AF3 mmCIF).
#
# Usage:
#   ./run_droideka_sites.sh <droideka_cif> <workspace>
#
# Outputs 6 subdirs (one per La site), each with full discriminator pipeline.

set -uo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <droideka_cif> <workspace_dir>"
  exit 1
fi

DROI_CIF="$1"
WORKSPACE="$2"
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)/scripts"

SHARED_ENV=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm
PYTHON="$SHARED_ENV/bin/python"

mkdir -p "$WORKSPACE"

# 1. Normalize + protonate ONCE (shared across all 6 sites)
NORM=$WORKSPACE/droideka_normalized.pdb
PROTON=$WORKSPACE/droideka_protonated.pdb
"$PYTHON" "$SCRIPTS_DIR/normalize_af3_cif.py" "$DROI_CIF" "$NORM" 2>&1 | tail -2
"$PYTHON" "$SCRIPTS_DIR/protonate_cif.py" "$NORM" "$PROTON" 2>&1 | tail -2

# 2. For each of 6 La chains, carve at that specific site
submitted=0
for site in B C D E F G; do
  stem="droideka_site${site}"
  out_dir="$WORKSPACE/${stem}_qm"
  if [ -f "$out_dir/sp_${stem}_La.out" ] && grep -q "FINAL SINGLE POINT ENERGY" "$out_dir/sp_${stem}_La.out" 2>/dev/null; then
    echo "[DONE] $stem"
    continue
  fi
  mkdir -p "$out_dir"
  cp "$PROTON" "$out_dir/${stem}_protonated.pdb"
  echo ""
  echo "=== $stem (chain $site) ==="

  "$PYTHON" "$SCRIPTS_DIR/carve_generic.py" "$out_dir/${stem}_protonated.pdb" "$out_dir" \
    --stem "$stem" --site-chain "$site" --site-resnum 1 2>&1 | tail -8

  # Add water SP input
  if [ ! -f "$out_dir/sp_${stem}_water.inp" ]; then
    cat > "$out_dir/bulk_water.xyz" <<'EOF'
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

  # Patch submit script: add water to the loop, add NoAutostart to inps
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

    jid=$(sbatch --parsable "$submit" 2>&1)
    if [[ "$jid" =~ ^[0-9]+$ ]]; then
      echo "[OK] $stem -> job $jid"
      submitted=$((submitted+1))
    fi
  fi
done

echo ""
echo "Submitted $submitted of 6 droideka sites."
echo "Aggregate when done: $PYTHON $SCRIPTS_DIR/analyze_results.py $WORKSPACE"
