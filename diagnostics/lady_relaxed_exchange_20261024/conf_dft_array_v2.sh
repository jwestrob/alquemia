#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -c 1
#SBATCH --array=1-5
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
# v2: like conf_dft_array.sh, but task i takes the i-th lowest CREST conformer whose
# metal coordination (O/N within 2.95 A of atom 1) matches seed.xyz, because
# --noreftopo sampling can move a water out of the first shell.
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=1
C="$1"; i=$SLURM_ARRAY_TASK_ID
[ -s "$C/crest_conformers.xyz" ] || { echo "no conformers"; exit 0; }
R="$C/dft_c$i"; mkdir -p "$R"
python3 - "$C" $i "$R/start.xyz" <<'PY' || exit 0
import sys, numpy as np
C, i, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
def frames(p):
    L = open(p).read().splitlines(); k = 0
    while k < len(L) and L[k].strip():
        n = int(L[k]); yield L[k:k + n + 2]; k += n + 2
def cn(fr):
    a = [l.split() for l in fr[2:]]; X = np.array([[float(v) for v in r[1:4]] for r in a])
    return sum(1 for j in range(1, len(a)) if a[j][0] in ("O", "N") and np.linalg.norm(X[j] - X[0]) < 2.95)
ref = cn(next(frames(f"{C}/seed.xyz")))
ok = [f for f in frames(f"{C}/crest_conformers.xyz") if cn(f) == ref]
print(f"seed CN {ref}; {len(ok)} CREST conformers match")
if i > len(ok): sys.exit(1)
open(out, "w").write("\n".join(ok[i - 1]) + "\n")
PY
cd "$R"
sed -e 's/ TightOpt / Opt /' -e "s/nprocs [0-9@A-Z]* end/nprocs $SLURM_NTASKS end/" "$2" > job.inp
grep -q "start.xyz" job.inp || { echo "bad base input"; exit 1; }
$ORCA/orca job.inp > job.out 2> job.err
echo "orca rc=$?"
