"""Audit saved native outputs. Never launches ORCA or changes root criteria."""
import os,re,argparse,datetime,subprocess
import numpy as np
from common import *
assert os.environ.get('SLURM_JOB_ID')
ap=argparse.ArgumentParser();ap.add_argument('--partial',action='store_true');args=ap.parse_args()
manifest=read(W/'REFERENCE_MANIFEST.json');panel={r['id']:r for r in read(W/'PANEL.json')['rows']}
jobids=[read(Path(r['directory'])/'RUNTIME.json')['job'] for r in manifest['rows']]
accounting=subprocess.check_output(['sacct','-n','-P','-j',','.join(jobids),'--format=JobIDRaw,State,ExitCode'],text=True)
states={row[0]:dict(state=row[1],exit_code=row[2]) for line in accounting.splitlines() if len(row:=line.split('|'))>=3 and row[0] in jobids}
# ORCA's documented native coordinate convention (manual 6.0, General 1.2),
# also verified against the 12-decimal Bohr coordinates in 6.1.1 property files.
ORCA_BOHR_A=.5291772083

def last_number(text,pattern,cast=float):
 a=re.findall(pattern,text,re.M);return cast(a[-1]) if a else None

def parse_native(text):
 vals={}
 for m in re.finditer(r'^\s*Last (.*?)\s+\.\.\.\s+([-+\d.eE]+)\s+Tolerance\s*:\s*([-+\d.eE]+)',text,re.M):
  v,t=float(m[2]),float(m[3]);vals[m[1].strip()]=dict(value=v,tolerance=t,passes=abs(v)<=t)
 meta={k:last_number(text,p,int) for k,p in dict(charge=r'Total Charge\s+Charge\s+\.\.\.\.\s+([-\d]+)',
      multiplicity=r'Multiplicity\s+Mult\s+\.\.\.\.\s+(\d+)',electrons=r'Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)',
      basis_dimension=r'Basis Dimension\s+Dim\s+\.\.\.\.\s+(\d+)').items()}
 stab=text[text.rfind('WAVEFUNCTION STABILITY ANALYSIS'):] if 'WAVEFUNCTION STABILITY ANALYSIS' in text else ''
 stab=stab.split('TOTAL SCF ENERGY')[0]
 eig=[float(x) for x in re.findall(r'E\(\s*\d+\)\s*=\s*([-+\d.]+)\s+Eh',stab)]
 occupations={}
 orb=text[text.rfind('\nORBITAL ENERGIES\n'):] if '\nORBITAL ENERGIES\n' in text else ''
 for name,end in [('UP','SPIN DOWN ORBITALS'),('DOWN','MULLIKEN')]:
  start='SPIN '+name+' ORBITALS'
  section=orb.split(start,1)[1].split(end,1)[0] if start in orb else ''
  rows=re.findall(r'^\s*(\d+)\s+(\d+\.\d+)\s+([-+\d.]+)\s+([-+\d.]+)\s*$',section,re.M)
  n=meta['basis_dimension']
  rows=rows[:n] if n else rows
  occ=[(int(i),float(o),float(e)) for i,o,e,_ in rows]
  occupied=[e for _,o,e in occ if o>.5];virtual=[e for _,o,e in occ if o<.5]
  occupations[name]=dict(n_orbitals=len(occ),sum_occupations=sum(o for _,o,_ in occ),
    HOMO_Eh=max(occupied) if occupied else None,LUMO_Eh=min(virtual) if virtual else None,
    frontier=[dict(index=i,occupation=o,energy_Eh=e) for i,o,e in occ if i in range(len(occupied)-7,len(occupied)+8)])
 return dict(metadata=meta,normal_termination='ORCA TERMINATED NORMALLY' in text,
     energy_Eh=last_number(text,r'FINAL SINGLE POINT ENERGY\s+([-+\d.]+)'),
     S2=last_number(text,r'Expectation value of <S\*\*2>\s*:\s*([-+\d.]+)'),
     scf_converged_cycles=[int(x) for x in re.findall(r'SCF CONVERGED AFTER\s+(\d+) CYCLES',text)],
     printed_scf_criteria=vals,all_printed_scf_criteria_pass=len(vals)>=4 and all(v['passes'] for v in vals.values()),
     native_ecp_lines=[l.strip() for l in text.splitlines() if 'replacing ' in l and 'core electrons' in l],
     native_basis_groups=[l.strip() for l in text.splitlines() if re.match(r'\s*Group\s+\d+ Type\s+\w+\s+:',l)],
     stability=dict(eigenvalues_Eh=eig,stable_statement='wavefunction is stable' in stab,
       converged='CONVERGENCE OF RESIDUAL NORM REACHED' in stab,excerpt=stab[-3000:]),
     occupations=occupations,native_method=dict(functional='Functional name                       .... wB97M-V' in text,
       UHF='Hartree-Fock type      HFTyp           .... UHF' in text,
       automatic_SCNL='Switching to SCNL' in text,automatic_AutoTRAH_off='Switching AutoTRAH off' in text,
       NL_energy_Eh=last_number(text,r'NL Energy, E\(C,NL\)\s*:\s*([-+\d.]+)'),
       engrad_module='ORCA SCF GRADIENT CALCULATION' in text),
     relevant_warnings=[l.strip() for l in text.splitlines() if 'WARNING' in l or 'Warning' in l][:35])

def engrad(path):
 tokens=[]
 for line in Path(path).read_text().splitlines():
  if line.strip() and not line.lstrip().startswith('#'):tokens.extend(line.split())
 n=int(tokens[0]);assert len(tokens)==2+7*n,(path,n,len(tokens))
 energy=float(tokens[1]);g=np.array([float(t) for t in tokens[2:2+3*n]]).reshape(n,3)
 rows=np.array([float(t) for t in tokens[2+3*n:]]).reshape(n,4)
 assert np.isfinite(g).all() and np.isfinite(rows).all() and np.isfinite(energy)
 # ORCA prints engrad coordinates to seven decimals in Bohr. Compare within
 # the file's own last-digit rounding plus the Bohr conversion convention.
 ct=np.array(tokens[2+3*n:]).reshape(n,4)[:,1:]
 rounding=np.array([[.5*10**(-len(t.split('.')[-1]))*ORCA_BOHR_A for t in row] for row in ct])
 return energy,g,rows[:,0].astype(int).tolist(),rows[:,1:]*ORCA_BOHR_A,rounding+1e-12

results=[]
for r in manifest['rows']:
 d=Path(r['directory']);p=panel[r['id']];text=(d/'job.out').read_text(errors='replace')
 runtime=read(d/'RUNTIME.json');scheduler=states.get(runtime['job'],{})
 terminal=scheduler.get('state','').split(' ')[0] in ('COMPLETED','FAILED','CANCELLED','OUT_OF_MEMORY','TIMEOUT','NODE_FAIL')
 v=parse_native(text);v.update(id=r['id'],terminal=terminal,scheduler=scheduler,output=pin(d/'job.out'),input=pin(d/'job.inp'),
                             xyz=pin(d/'start.xyz'),basis=pin(d/'def2-tzvpd.bas'))
 if not v['terminal']:
  assert args.partial,'Wait for all actual terminal results before final collection.'
  v.update(qualified=False,reason='native_job_still_running');results.append(v);continue
 v['exit']=read(d/'EXIT.json') if (d/'EXIT.json').exists() else dict(returncode=None,native_exit_record_missing=True)
 v['runtime']=runtime
 expected=(d/'prepared.inp').read_text().replace('@NPROCS@',str(v['runtime']['nprocs'])).replace('@MAXCORE@',str(v['runtime']['maxcore_MB']))
 assert expected==(d/'job.inp').read_text(),'Native input differs from frozen preparation.'
 assert pin(d/'def2-tzvpd.bas')==r['basis']
 z,x=xyzread(d/'start.xyz');assert z==p['Z'] and x==p['positions_A']
 v['engrad_complete']=(d/'job.engrad').exists()
 if v['engrad_complete']:
  energy,g,z,positions,tolerance=engrad(d/'job.engrad')
  prop=(d/'job.property.txt').read_text();section=prop.split('&CartesianCoordinates',1)[1].split('$End',1)[0]
  prows=re.findall(r'^\s*(\w+)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s*$',section,re.M)
  assert len(prows)==len(p['Z'])
  ppos=np.array([[float(t) for t in row[1:]] for row in prows])*ORCA_BOHR_A
  pz=[next(z for z,label in SYMBOLS.items() if label==row[0]) for row in prows]
  property_match=pz==p['Z'] and bool(np.max(abs(ppos-np.array(p['positions_A'])))<1e-11)
  v.update(engrad=pin(d/'job.engrad'),native_property=pin(d/'job.property.txt'),gradient_energy_Eh=energy,
           forces_kcal_mol_A=(-g*EH_KCAL/ORCA_BOHR_A).tolist(),native_property_coordinates_match=property_match,
           native_property_coordinate_max_error_A=float(np.max(abs(ppos-np.array(p['positions_A'])))),
           engrad_coordinates_match=property_match and z==p['Z'] and bool(np.all(abs(positions-ppos)<=tolerance)),
           engrad_coordinate_max_rounding_error_A=float(np.max(abs(positions-np.array(p['positions_A'])))),
           engrad_coordinate_tolerance_A=float(tolerance.max()),
           energy_gradient_match=v['energy_Eh'] is not None and abs(energy-v['energy_Eh'])<5e-8)
 expected_spin=(r['multiplicity']-1)/2;S2=expected_spin*(expected_spin+1)
 v['S2_target']=S2;v['S2_gate']=v['S2'] is not None and abs(v['S2']-S2)<=.5
 v['metadata_gate']=all(v['metadata'][k]==r['expected_'+k] if 'expected_'+k in r else v['metadata'][k]==r[k] for k in ('charge','multiplicity','electrons'))
 met='Dy' if r['multiplicity']==6 else 'La'
 v['ecp_gate']=any('Type '+met+' ECP' in l and f'replacing {r["expected_ecp_core"]} core electrons' in l for l in v['native_ecp_lines'])
 st=v['stability'];v['stability_gate']=bool(st['eigenvalues_Eh']) and st['stable_statement'] and st['converged'] and min(st['eigenvalues_Eh'])>=-1e-4
 oc=v['occupations'];v['occupation_gate']=all(oc[spin]['n_orbitals']==v['metadata']['basis_dimension'] and abs(oc[spin]['sum_occupations']-(r['expected_electrons']+sign*(r['multiplicity']-1))/2)<1e-6 for spin,sign in [('UP',1),('DOWN',-1)])
 gates=dict(normal_termination=v['normal_termination'],zero_exit=v['exit']['returncode']==0,
            SCF=bool(v['scf_converged_cycles']),all_printed_Tight_residuals=v['all_printed_scf_criteria_pass'],
            physical_spin=v['S2_gate'],stability=v['stability_gate'],state_metadata=v['metadata_gate'],core=v['ecp_gate'],occupations=v['occupation_gate'],
            method=all(v['native_method'][k] for k in ['functional','UHF','automatic_SCNL','engrad_module']),
            gradient=v.get('engrad_complete',False) and v.get('engrad_coordinates_match',False) and v.get('energy_gradient_match',False))
 v['gates']=gates;v['qualified']=all(gates.values());v['reason']=[k for k,b in gates.items() if not b]
 results.append(v)

# Audit the two historical Dy core cycles with exactly the same residual parser.
old=read(W/'OLD_REFERENCE_AUDIT.json');old_details={}
for kind in ('gas','sc_gas'):
 for state in ('Dy_aq8','Dy_aq9'):
  path=Path(old['results'][kind+'/'+state]['output']['path']);v=parse_native(path.read_text(errors='replace'))
  old_details[kind+'/'+state]=dict(output=pin(path),printed_scf_criteria=v['printed_scf_criteria'],
     all_printed_scf_criteria_pass=v['all_printed_scf_criteria_pass'],stability=v['stability'],occupations=v['occupations'])

finished=all(r['terminal'] for r in results)
out=dict(rows=results,all_terminal=finished,all_qualified=all(r['qualified'] for r in results),
         plan=pin(D/'PLAN.md'),manifest=pin(W/'REFERENCE_MANIFEST.json'),old_reference_residual_audit=old_details,
         collector=pin(__file__),job=os.environ['SLURM_JOB_ID'],new_molecular_calls=0)
if finished:dest=RESULTS_DIR/'REFERENCE_RESULTS.json'
else:dest=W/('REFERENCE_PARTIAL_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'.json')
save(dest,out)
for r in results:print(r['id'],'terminal',r['terminal'],'qualified',r['qualified'],'reasons',r['reason'])
print('saved',dest)
