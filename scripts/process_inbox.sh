#!/bin/bash
# Process all CIFs in the inbox: protonate, carve, submit.
# Drop CIFs in alchemical_bvs/inbox/<stem>.cif and run this.
# Guarded by flock so the watcher + a manual run don't double-submit.
set -uo pipefail

# Acquire exclusive lock; if another instance is mid-run, exit cleanly.
exec 200>/tmp/process_inbox.lock
if ! flock -n 200; then
    echo "[INFO] process_inbox.sh: another instance holds the lock — skipping"
    exit 0
fi

ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
INBOX=$ALCH/inbox
PROCESSED=$INBOX/processed
FEP_PYTHON=/home/jwestrob/miniconda3/envs/fep/bin/python
SCAN_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python

mkdir -p "$PROCESSED"
submitted=0
skipped=0

# Per-invocation submit cap (prevents QOSMaxSubmit hits when inbox is huge).
# Override with MAX_SUBMIT=N bash process_inbox.sh
MAX_SUBMIT=${MAX_SUBMIT:-100}

for cif in "$INBOX"/*.cif; do
  [ -f "$cif" ] || continue
  if [ "$submitted" -ge "$MAX_SUBMIT" ]; then
    echo "[CAP] hit MAX_SUBMIT=$MAX_SUBMIT, leaving rest of inbox for next invocation"
    break
  fi
  base=$(basename "$cif" .cif)
  # Sanitize: replace '.' and other special chars in stem
  stem=$(echo "$base" | tr -c '[:alnum:]_-' '_' | sed 's/__*/_/g' | sed 's/_$//')

  out_dir="$ALCH/workspaces/${stem}_qm"
  if [ -f "$out_dir/sp_${stem}_La.out" ] && grep -q "FINAL SINGLE POINT ENERGY" "$out_dir/sp_${stem}_La.out" 2>/dev/null; then
    echo "[DONE] $stem already processed"
    mv "$cif" "$PROCESSED/"
    skipped=$((skipped+1))
    continue
  fi

  # 2026-05-11: relaxed from ` LA ` (which only matched mmCIF chem_comp_dict lines
  # for the bare residue name "LA") to also accept HETATM records with element La,
  # regardless of residue name. Colin's PQQ-ADH drops use comp_id "LA1" instead of
  # "LA", which made the old filter reject all 192 of them as "no La in CIF".
  if ! grep -qE '(^| )LA |HETATM[[:space:]]+[0-9]+[[:space:]]+La[[:space:]]' "$cif" 2>/dev/null; then
    echo "[SKIP] $stem: no La in CIF"
    skipped=$((skipped+1))
    continue
  fi

  mkdir -p "$out_dir"
  echo ""
  echo "=== $stem from $base.cif ==="

  # 2026-05-11: detect PQQ-bearing CIFs (Colin's PQQ-MDH/ADH drops and the
  # Euk PQQ-ADH bolus all encode PQQ as a 24-atom LIG_* residue). If we see
  # a ligand residue with PQQ-like elemental composition (≥12 C + ≥6 O + ≥1 N),
  # route through normalize_af3_cif → protonate → carve_with_pqq so the PQQ
  # cofactor (the primary Ln coordinator at PQQ-O5) is in the QM cluster.
  # Otherwise: protonate → carve_generic (original behaviour).
  has_pqq=$("$SCAN_PYTHON" - "$cif" <<'PYEOF' 2>/dev/null
import sys, gemmi
try:
    st = gemmi.read_structure(sys.argv[1])
    for model in st:
        for chain in model:
            for res in chain:
                if not res.name.startswith("LIG"):
                    continue
                c = sum(1 for a in res if a.element.name == "C")
                o = sum(1 for a in res if a.element.name == "O")
                n = sum(1 for a in res if a.element.name == "N")
                if c >= 12 and o >= 6 and n >= 1:
                    print("yes"); sys.exit(0)
    print("no")
except Exception:
    print("no")
PYEOF
)

  if [ "$has_pqq" = "yes" ]; then
    # ── PQQ pipeline: normalize → protonate → carve_with_pqq ──
    echo "  PQQ detected → tier1 pipeline"
    pdb_norm="$out_dir/${stem}_normalized.pdb"
    if [ ! -f "$pdb_norm" ]; then
      "$SCAN_PYTHON" "$ALCH/scripts/normalize_af3_cif.py" "$cif" "$pdb_norm" 2>&1 | tail -2 || {
        echo "[FAIL] normalize"; continue
      }
    fi
    pdb_proton="$out_dir/${stem}_protonated.pdb"
    if [ ! -f "$pdb_proton" ]; then
      "$FEP_PYTHON" "$ALCH/scripts/protonate_cif.py" "$pdb_norm" "$pdb_proton" 2>&1 | tail -2 || {
        echo "[FAIL] protonation"; continue
      }
    fi
    la_chain=$(grep "^HETATM" "$pdb_proton" | awk '$4=="LA"{print $5; exit}')
    [ -z "$la_chain" ] && la_chain="B"
    if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
      "$SCAN_PYTHON" "$ALCH/scripts/carve_with_pqq.py" "$pdb_proton" "$out_dir" \
        --stem "$stem" --metal-chain "$la_chain" --metal-resname "LA" 2>&1 | tail -8 || {
        echo "[FAIL] carve_with_pqq"; continue
      }
    fi
  else
    # ── Generic pipeline: protonate → carve_generic ──
    pdb_proton="$out_dir/${stem}_protonated.pdb"
    if [ ! -f "$pdb_proton" ]; then
      "$FEP_PYTHON" "$ALCH/scripts/protonate_cif.py" "$cif" "$pdb_proton" 2>&1 | tail -2 || {
        echo "[FAIL] protonation"; continue
      }
    fi
    if [ ! -f "$out_dir/${stem}_La_qm.xyz" ]; then
      "$SCAN_PYTHON" "$ALCH/scripts/carve_generic.py" "$pdb_proton" "$out_dir" --stem "$stem" 2>&1 | tail -8 || {
        echo "[FAIL] carve"; continue
      }
    fi
  fi

  # Add water SP
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
    echo "[OK] $stem -> job $jid"
    mv "$cif" "$PROCESSED/"
    submitted=$((submitted+1))
  else
    echo "[FAIL] sbatch: $jid"
  fi
done

echo ""
echo "Submitted: $submitted | Skipped: $skipped"
