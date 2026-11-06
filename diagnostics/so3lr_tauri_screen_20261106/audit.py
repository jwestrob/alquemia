import os,re,csv,time
import numpy as np
from common import *
assert os.environ.get('SLURM_JOB_ID')
start=time.time(); panel=read(W/'PANEL.json')['rows'];by={r['id']:r for r in panel}
old=ROOT/'workspaces/lady_ingredient_tests_20261026';results={}
for kind in ('gas','sc_gas'):
 for state in ['La_aq8','La_aq9','Dy_aq8','Dy_aq9','Dy_ion','La_ion','La_acmono_w7','Dy_acmono_w7']+(['water'] if kind=='gas' else []):
  d=old/kind/state;p=d/'job.out'
  if not p.exists():continue
  s=p.read_text(errors='replace');inp=(d/'job.inp').read_text();z,x=xyzread(d/'start.xyz')
  en=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)',s);s2=re.findall(r'Expectation value of <S\*\*2>\s*:\s*([-\d.]+)',s)
  vals={}
  for k,pattern in [('charge',r'Total Charge\s+Charge\s+\.\.\.\.\s+([-\d]+)'),('multiplicity',r'Multiplicity\s+Mult\s+\.\.\.\.\s+(\d+)'),('electrons',r'Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)'),('nbasis',r'Number of basis functions\s+.*?(\d+)')]:
   m=re.findall(pattern,s,re.I);vals[k]=m[-1] if m else None
  match=by.get(('hydrated_acetate_'+state[:2]) if 'acmono' in state else ('aquo'+state[-1]+'_central_'+state[:2]))
  geom_equal=bool(match and z==match['Z'] and np.array_equal(np.array(x),np.array(match['positions_A'])))
  results[kind+'/'+state]=dict(energy_Eh=float(en[-1]) if en else None,S2=float(s2[-1]) if s2 else None,
    normal_termination='ORCA TERMINATED NORMALLY' in s,scf_cycles=[int(v) for v in re.findall(r'SCF CONVERGED AFTER\s+(\d+) CYCLES',s)],
    ecp_lines=[line.strip() for line in s.splitlines() if 'replacing ' in line and 'core electrons' in line],
    stability_lines=[line.strip() for line in s.splitlines() if 'stability' in line.lower()][-12:],
    metadata=vals,geometry_matches_fixed_panel=geom_equal,input=pin(d/'job.inp'),xyz=pin(d/'start.xyz'),output=pin(p),
    orbital_occupation_excerpts='\n'.join(s[s.rfind('ORBITAL ENERGIES'):].splitlines()[:150]),
    qualified_for_new_panel=False,reason='PBE0-D4 vs requested OMol-style method; no matched analytic gradient; geometry check explicit')
water=results['gas/water']['energy_Eh'];cyc={}
for kind in ('gas','sc_gas'):
 a=results[kind+'/Dy_aq8']['energy_Eh'];b=results[kind+'/Dy_aq9']['energy_Eh']
 cyc[kind]=(a+water-b)*EH_KCAL
rows=list(csv.DictReader(open(ROOT/'workspaces/ln_structure_increment_20261030/evaluation_records.tsv'),delimiter='\t'))
def select(lid,I):return next(r for r in rows if r['lig_id']==lid and float(r['I_M'])==I)
pairs=[]
for name,l1,l2,I in [('cis_trans','599','600',1.),('fallback_EDTA_cyclopent','576','597',.1)]:
 rr=[select(l1,I),select(l2,I)]
 raw=[]
 v=(ROOT/'workspaces/ln_selectivity_calibration_20261027/nist46_dat/F8558_VALUES.DAT').read_bytes()
 for r in rr:
  for met in ('La','Dy'):
   off=int(r[met+'_offset']);b=v[off:off+15]
   assert int.from_bytes(b[:2],'little')==int(r['lig_id'])
   assert float(b[6:13].decode().strip())==float(r['logK_'+met])
   raw.append(dict(ligand=r['lig_id'],metal=met,offset=off,bytes_hex=b.hex()))
 d=-1.3642470130342386*((float(rr[1]['logK_Dy'])-float(rr[1]['logK_La']))-(float(rr[0]['logK_Dy'])-float(rr[0]['logK_La'])))
 pairs.append(dict(name=name,records=rr,raw_records=raw,D_kcal_conditional=d,status='unresolved_not_admitted',
   reason='Specific background anion is not identified by native annotations; compilation/source rows do not establish a common electrolyte. Fallback also has adjusted vs unadjusted annotations. No salt is guessed.',
   row_conditions_matched=all(r['T_K']=='298.15' and r['protonation']=='H4L' and r['equilibrium']=='ML=M.L' for r in rr)))
save(W/'OLD_REFERENCE_AUDIT.json',dict(results=results,Dy_W9_to_W8_plus_water_kcal=cyc,large_minus_small_kcal=cyc['gas']-cyc['sc_gas'],compatible_reference_available=False,new_calls=0,cost_wall_seconds=time.time()-start,job=os.environ['SLURM_JOB_ID']))
save(W/'EXPERIMENTAL_PAIR_AUDIT.json',dict(pairs=pairs,admitted_pair=None,no_additional_pair_search=True))
print('old Dy contrasts',cyc,'neither experimental pair admitted; missing explicit electrolyte/source medium')
