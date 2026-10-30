#!/bin/bash
# La3+ / Dy3+ in a cubic OPC box (no counter-ions: uniform Ewald background, same charge both ends).
# Writes <ws>/aq/{La,Dy}.prmtop (12-6-4 via parmed add12_6_4) and a shared rst7.
set -euo pipefail
P=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm; export AMBERHOME=$P PATH=$P/bin:$PATH
WS=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lady_explicit_solvent_20261026/aq
mkdir -p $WS && cd $WS
for M in La Dy; do U=$([ $M = La ] && echo LA || echo Dy); R=$([ $M = La ] && echo LA || echo Dy)
cat > tleap_$M.in <<L
source leaprc.water.opc
loadAmberParams $P/dat/leap/parm/frcmod.ionslm_1264_opc
s = copy $U
solvateBox s OPCBOX 15.5 iso
saveAmberParm s ${M}_126.prmtop ${M}.rst7
quit
L
tleap -f tleap_$M.in > tleap_$M.log
cat > p_$M.in <<L
add12_6_4 :$R watermodel OPC polfile $P/dat/leap/parm/lj_1264_pol.dat tunfactor 1.0
outparm $M.prmtop
quit
L
parmed -p ${M}_126.prmtop -i p_$M.in > parmed_$M.log
done
cmp La.rst7 Dy.rst7 && echo "same coordinates"
grep -c WAT La.rst7 >/dev/null; grep "residues\|Added" tleap_La.log | tail -3
