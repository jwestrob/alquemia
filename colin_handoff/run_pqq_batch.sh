#!/bin/bash
# Run the DFT Ca²⁺/Ln³⁺ discriminator on a batch of PQQ-bound (MDH-like) AlphaFold3 CIFs.
#
# Usage:
#   ./run_pqq_batch.sh <input_dir> <workspace>
#
# <input_dir>:  directory containing AF3 .cif files (any number)
# <workspace>:  output directory; one subdir per candidate; results.tsv at the end
#
# Each CIF goes through: normalize → protonate → carve_with_pqq → submit 4 SPs to SLURM.
# Re-running is idempotent (skips already-completed candidates).
#
# The pipeline is the SAME r²SCAN-3c CPCM(Water) DefGrid3 used by Jacob's panel,
# so results are directly comparable. Sign convention: ΔΔE(Ca - La) > 0 → Ln-preferring.

set -uo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <input_dir_with_cifs> <workspace_dir>"
  echo ""
  echo "Example:"
  echo "  $0 /groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone/Euk \\"
  echo "     /groups/banfield/users/colinr/discriminator_run_euk"
  exit 1
fi

INPUT_DIR="$1"
WORKSPACE="$2"
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)/scripts"

# Universal env paths (shared, accessible by both Jacob and Colin)
SHARED_ENV=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm
PYTHON="$SHARED_ENV/bin/python"
ORCA_PATH=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg

# Sanity checks
if [ ! -f "$PYTHON" ]; then
  echo "ERROR: shared python env not found at $SHARED_ENV"
  echo "Tell Jacob; the env should be in /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/"
  exit 1
fi
if [ ! -d "$INPUT_DIR" ]; then
  echo "ERROR: input dir not found: $INPUT_DIR"
  exit 1
fi
if [ ! -x "$ORCA_PATH/orca" ]; then
  echo "ERROR: ORCA not found at $ORCA_PATH"
  exit 1
fi

mkdir -p "$WORKSPACE"

submitted=0
skipped=0
failed=0

for cif in "$INPUT_DIR"/*.cif; do
  [ -f "$cif" ] || continue
  base=$(basename "$cif" .cif)
  # Sanitize: replace non-alphanum with underscore
  stem=$(echo "$base" | tr -c '[:alnum:]_' '_' | sed 's/__*/_/g' | sed 's/_$//')
  out_dir="$WORKSPACE/${stem}_qm"

  # Idempotent skip
  if [ -f "$out_dir/sp_${stem}_La.out" ] && grep -q "FINAL SINGLE POINT ENERGY" "$out_dir/sp_${stem}_La.out" 2>/dev/null; then
    echo "[DONE] $stem"
    skipped=$((skipped+1))
    continue
  fi

  mkdir -p "$out_dir"
  echo ""
  echo "=== $stem ==="

  # 1. Normalize AF3 CIF (LIG_B → LA, LIG_C → PQQ etc.)
  pdb_norm="$out_dir/${stem}_normalized.pdb"
  if [ ! -f "$pdb_norm" ]; then
    "$PYTHON" "$SCRIPTS_DIR/normalize_af3_cif.py" "$cif" "$pdb_norm" 2>&1 | tail -2
  fi

  # 2. Protonate (PDBFixer @ pH 7)
  pdb_proton="$out_dir/${stem}_protonated.pdb"
  if [ ! -f "$pdb_proton" ]; then
    "$PYTHON" "$SCRIPTS_DIR/protonate_cif.py" "$pdb_norm" "$pdb_proton" 2>&1 | tail -2
    if [ ! -f "$pdb_proton" ]; then
      echo "[FAIL] protonation"; failed=$((failed+1)); continue
    fi
  fi

  # 3. Find La chain after protonation (PDBFixer can renumber)
  la_chain=$(grep "^HETATM" "$pdb_proton" | awk '$4=="LA"{print $5; exit}')
  if [ -z "$la_chain" ]; then
    la_chain="B"
  fi

  # 4. Carve QM cluster (PQQ-aware)
  if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
    "$PYTHON" "$SCRIPTS_DIR/carve_with_pqq.py" "$pdb_proton" "$out_dir" \
      --stem "$stem" --metal-chain "$la_chain" --metal-resname "LA" 2>&1 | tail -8
    if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
      echo "[FAIL] carve"; failed=$((failed+1)); continue
    fi
  fi

  # 5. Submit SLURM job
  submit="$out_dir/submit_${stem}.sh"
  if [ ! -f "$submit" ]; then
    echo "[FAIL] no submit script"; failed=$((failed+1)); continue
  fi

  jid=$(sbatch --parsable "$submit" 2>&1)
  if [[ "$jid" =~ ^[0-9]+$ ]]; then
    echo "[OK] $stem -> job $jid"
    submitted=$((submitted+1))
  else
    echo "[FAIL sbatch] $jid"
    failed=$((failed+1))
  fi
done

echo ""
echo "=========================================="
echo "Submitted: $submitted | Skipped (already done): $skipped | Failed: $failed"
echo ""
echo "Once jobs complete, aggregate results:"
echo "  $PYTHON $SCRIPTS_DIR/analyze_results.py $WORKSPACE"
echo ""
echo "Live monitoring:"
echo "  squeue -u \$USER -o \"%.10i %.40j %.10T %.10M\""
