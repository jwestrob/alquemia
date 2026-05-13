#!/bin/bash
# Set up B97-3c functional robustness panel.
# Re-runs LanM, XoxF, MxaF, tannase at B97-3c CPCM(Water).
# Plus aquo-Ca / aquo-La / bulk-water at B97-3c for the discriminator cycle.
set -uo pipefail
ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ORCA=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg

mkdir -p $ALCH/b97_3c_panel

# Aquo references (use existing CN8 xyz files)
mkdir -p $ALCH/b97_3c_panel/aquo
cp $ALCH/qmmm/aquo_panel/Ca/aquo_Ca_cn8.xyz $ALCH/b97_3c_panel/aquo/ 2>/dev/null
cp $ALCH/qmmm/aquo_test/aquo_La_cn8.xyz $ALCH/b97_3c_panel/aquo/ 2>/dev/null

cat > $ALCH/b97_3c_panel/aquo/sp_aquo_Ca_b97.inp <<EOF
! B97-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
* xyzfile 2 1 aquo_Ca_cn8.xyz
EOF
cat > $ALCH/b97_3c_panel/aquo/sp_aquo_La_b97.inp <<EOF
! B97-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
%basis
  NewECP La "def2-ECP" end
  NewGTO La "def2-TZVP" end
end
* xyzfile 3 1 aquo_La_cn8.xyz
EOF
cat > $ALCH/b97_3c_panel/aquo/bulk_water.xyz <<EOF
3
H2O reference
O   0.000000  0.000000  0.000000
H   0.756950  0.000000  0.585822
H  -0.756950  0.000000  0.585822
EOF
cat > $ALCH/b97_3c_panel/aquo/sp_water_b97.inp <<EOF
! B97-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
* xyzfile 0 1 bulk_water.xyz
EOF

cat > $ALCH/b97_3c_panel/aquo/submit_aquo.sh <<EOF
#!/bin/bash
#SBATCH -p memory
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J b97_aquo
#SBATCH -o $ALCH/b97_3c_panel/aquo/slurm_%j.out
#SBATCH -e $ALCH/b97_3c_panel/aquo/slurm_%j.err
set -uo pipefail
export PATH=$ORCA:\$PATH
export LD_LIBRARY_PATH=$ORCA:\${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=\$SLURM_CPUS_ON_NODE
cd $ALCH/b97_3c_panel/aquo
export ORCA_TMPDIR=\$PWD; export TMPDIR=\$PWD
echo "started=\$(date -Iseconds) host=\$(hostname)"
for k in Ca La water; do
  inp=sp_\${k}_b97.inp
  [ "\$k" = "water" ] && inp=sp_water_b97.inp
  echo "=== \$k ==="
  $ORCA/orca \$inp > sp_\${k}_b97.out 2>&1
  grep "FINAL SINGLE POINT" sp_\${k}_b97.out | tail -1
done
echo "finished=\$(date -Iseconds)"
EOF
chmod +x $ALCH/b97_3c_panel/aquo/submit_aquo.sh

# Per-candidate B97-3c: just La + Ca (skip apo since we're comparing ΔΔE not absolute)
for stem in xoxf mxaf tannase calexcitin; do
  out=$ALCH/b97_3c_panel/${stem}
  mkdir -p $out
  src=$ALCH/${stem}_qm
  cp $src/${stem}_La_qm.xyz $out/ 2>/dev/null
  cp $src/${stem}_Ca_qm.xyz $out/ 2>/dev/null

  charge_la=$(grep "xyzfile" $src/sp_${stem}_La.inp | awk '{print $3}')
  charge_ca=$(grep "xyzfile" $src/sp_${stem}_Ca.inp | awk '{print $3}')

  # Do La need ECP block?
  needs_la_ecp=""
  if grep -q "La " $out/${stem}_La_qm.xyz 2>/dev/null; then
    needs_la_ecp='%basis
  NewECP La "def2-ECP" end
  NewGTO La "def2-TZVP" end
end'
  fi

  cat > $out/sp_${stem}_La_b97.inp <<EOF
! B97-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
$needs_la_ecp
* xyzfile $charge_la 1 ${stem}_La_qm.xyz
EOF
  cat > $out/sp_${stem}_Ca_b97.inp <<EOF
! B97-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
* xyzfile $charge_ca 1 ${stem}_Ca_qm.xyz
EOF
  cat > $out/submit_b97.sh <<EOF
#!/bin/bash
#SBATCH -p memory
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J b97_${stem}
#SBATCH -o $out/slurm_%j.out
#SBATCH -e $out/slurm_%j.err
set -uo pipefail
export PATH=$ORCA:\$PATH
export LD_LIBRARY_PATH=$ORCA:\${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=\$SLURM_CPUS_ON_NODE
cd $out
export ORCA_TMPDIR=\$PWD; export TMPDIR=\$PWD
echo "started=\$(date -Iseconds) host=\$(hostname)"
for k in La Ca; do
  echo "=== \$k ==="
  $ORCA/orca sp_${stem}_\${k}_b97.inp > sp_${stem}_\${k}_b97.out 2>&1
  grep "FINAL SINGLE POINT" sp_${stem}_\${k}_b97.out | tail -1
done
echo "finished=\$(date -Iseconds)"
EOF
  chmod +x $out/submit_b97.sh
done

# Submit all
echo "=== Submitting B97-3c jobs ==="
for f in $ALCH/b97_3c_panel/aquo/submit_aquo.sh $ALCH/b97_3c_panel/xoxf/submit_b97.sh $ALCH/b97_3c_panel/mxaf/submit_b97.sh $ALCH/b97_3c_panel/tannase/submit_b97.sh $ALCH/b97_3c_panel/calexcitin/submit_b97.sh; do
  jid=$(sbatch --parsable $f 2>&1)
  echo "  $f → $jid"
done
echo ""
squeue -u jwestrob -o "%.10i %.20j %.10T %.10M %R" 2>/dev/null | head -15
