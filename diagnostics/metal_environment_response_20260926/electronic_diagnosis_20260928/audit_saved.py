"""Read only the five declared failed Dy outputs; never launch molecular work."""
import argparse,hashlib,json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
from slurm_native_health_watch import inspect_text,F,number
CASES={'embedded_A':'lanm_ef3_hans_scout_v1/Dy_A','embedded_B':'lanm_ef3_hans_scout_v1/Dy_B','isolated_native_PModel':'dy_small_guess_v3/PModel','isolated_native_HCore':'dy_small_guess_v3/HCore','isolated_PBE0_PModel':'dy_pbe0_capability_v1/PModel'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();result={}
 for name,rel in CASES.items():
  f=ROOT/'workspaces/metal_environment_response_20260926'/rel/'endpoint.out';t=f.read_text();r=inspect_text(t)
  rows=[{'iteration':int(i),'energy_hartree_diagnostic_only':number(e),'residual':number(g),'solver':s} for i,e,g,s in re.findall(r'^\s*(\d+)\s+('+F+r')\s+('+F+r').*\((TRAH|NR)\s+MAcro\)',t,re.M)]
  matches=re.findall(r'Smallest eigenvalue\s+\.{2,}\s+('+F+')',t)
  result[name]={'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'health':r,'macro_history':rows,'smallest_overlap_eigenvalue':number(matches[-1]) if matches else None,'accepted_energy':None,'occupation_assignment':'unavailable: no qualified converged molecular state','gbw_is_qualified_reference':False}
 out=Path(a.output)
 with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
 for name,r in result.items(): print(name,'macros',len(r['macro_history']),'overlap',r['smallest_overlap_eigenvalue'],'last',r['macro_history'][-1] if r['macro_history'] else None)
if __name__=='__main__':main()
