# Fixed-coordinate donor second shell — approved 2026-09-19

Jacob requested parallel second-shell modeling and full compute utilization before
biotite downtime. Root reviewed this concrete five-case scope and approved it:
1H4I/MxaF, 4MAE/XoxF, 1GLG/GGR, 1F6S and 6IP9/alpha-lactalbumin.
The temporary goal discussion pause ended with explicit instruction to resume.

Question: does explicit local donor-neighbor chemistry alter Ca/La discrimination
usefully, and does cheap native MACE reproduce its direction? Fixed rule: select
nonwater core O/N within3.3A of metal, complete their carbonyl/carboxylate polar
functional partners, and include complete standard peptide/sidechain fragments
with N/O/S within3.5A of those atoms. Covalent connectivity closes overlaps.
No new waters, optimization, protonation changes, thresholds or fitted parameters.
Keep every parent source coordinate and surviving cap exact; overlapping fragments
replace existing cut-bond caps with actual source connectivity, recorded explicitly.
PQQ cofactor remains complete and unchanged. Alpha starts from archived improved
water-H preparation; identical H coordinates are used before/after expansion.

Exactly10new native r2SCAN-3c/CPCM(Water)/DefGrid3 single points (Ca/La at5expanded
contexts), using the parent nativeSP recipe. No gradient requirement or baseline
rerun. Up to20native unmasked OMOL-0 100M float64 endpoint evaluations cover core
and expanded contexts. Reuse compatible archived results when available.
CPU execution:64allocatedCPUs, four concurrent16rank endpoints,256GiB host memory.
GPU execution: oneGPU,16CPUs, scheduler-approved RAM perGPU. Expected context sizes
roughly120–220atoms; anticipatedCPU scale tens ofminutes, GPU scale minutes.
These are estimates; measured receipts determine actual cost. No arbitrary budget
or time cutoff. Finite manifests retain failures and existing allocation policies.

Primary outputs: same-site deltaR and changes in XoxF-minus-MxaF and alpha-minus-GGR
contrasts, R=E_Ca-E_La. LargerR is moreLa-like on each scale. Expanded charges,
composition, CPCM boundary and fragment changes are explicit. No old threshold or
universalzero is inherited. Added-residue electronic/cavity effects are not called
a uniquely decomposed hydrogen-bond energy. MACE andDFT retain distinct scales.
Allcases are consumed development examples. Alpha replicas are onebiologicalgroup;
PQQ functional association and condition-qualified affinity evidence stay separate.
No claim of broad accuracy or production promotion follows from this smallpanel.

Pipeline default, archivedresults and canonical25casePQQpanel remain untouched.
