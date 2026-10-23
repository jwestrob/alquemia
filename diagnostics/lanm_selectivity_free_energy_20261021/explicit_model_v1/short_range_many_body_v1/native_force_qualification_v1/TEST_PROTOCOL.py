#!/usr/bin/env python3
"""Tiny manufactured grammar and fixed-gate checks; no native or torch import."""
import json,tempfile,time,platform
from pathlib import Path
from PREPARE import pin,save
import PROTOCOL
from COLLECT import slope

def main(out):
 t=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False);checks=[]
 slot={'geometry':{'xyz_A':[[0.,0.,0.]],'types':[500]}};fixture={'native_operator':{'K_kcal_A_per_mol_e2':332.0637133,'site_atom_indices':[0],'alpha_A3':[1.]}}
 text='''SRMB_BASELINE_V1
COUNTS 1 1 3
CONSTANTS 332.0637133 1 4.803206802 0.529177210544
OPERATOR_CHECKS 0 0
PERMANENT_ORIGIN_A 0 0 0
PERMANENT_MU_EA 0 0 0
CHARGE 3
MODE VACUUM MUTUAL THOLE DIELECTRIC_1 ZERO_FIELD ALL_ACTIVE NO_VIRTUAL_SITES
FLAGS_DPEQUAL_PCGGUESS_TCGGUESS F T F
CALLS 1 3 0
SITE 1 1 500 1 0 0 0
VECTOR 1 0.1 0.2
A 1 1 0 0
VECTOR 2 0.3 0.4
A 2 0 1 0
VECTOR 3 0.5 0.6
A 3 0 0 1
END_SRMB_BASELINE_V1
'''
 with tempfile.TemporaryDirectory() as tmp:
  p=Path(tmp)/'baseline';p.write_text(text);b=PROTOCOL.baseline(p,slot,fixture);assert b['Ed']!=b['Ep'];checks.append('unequal dual RHS preserved')
  for label,bad in [('wrong atom',text.replace('SITE 1 1 500','SITE 1 2 500')),('nonfinite field',text.replace('VECTOR 1 0.1','VECTOR 1 nan')),('trailing record',text+'EXTRA\n')]:
   p.write_text(bad)
   try:PROTOCOL.baseline(p,slot,fixture)
   except (AssertionError,ValueError,StopIteration):checks.append(label+' rejected')
   else:raise AssertionError(label)
  uv=Path(tmp)/'uv';u=[.1,.3,.5];v=[.2,.4,.6];PROTOCOL.write_uv(uv,b,u,v);assert 'ED 0.10000000000000001' in uv.read_text() and 'EP 0.20000000000000001' in uv.read_text();checks.append('17-digit dual-state writer')
  raw='''SRMB_FIXED_UV_GRADIENT_V1
COUNTS 1 1 3
CALLS 1 3 1 0
OTHER_ENERGY 1
RAW_REDUCED_POLARIZATION -2
RAW_NATIVE_ENERGY -1
SUPPLIED_UV_UNCHANGED T
UV 1 1 0.1 0.3 0.5 0.2 0.4 0.6
BASELINE_GRADIENT 1 2 3 4
POLAR_CONTRACTION_GRADIENT 1 1 1 1
RAW_NATIVE_VIRIAL 1 1 0 0
RAW_NATIVE_VIRIAL 2 0 1 0
RAW_NATIVE_VIRIAL 3 0 0 1
END_SRMB_FIXED_UV_GRADIENT_V1
''';p.write_text(raw);r=PROTOCOL.result(p,b,u,v);assert r['calls']['induced_solves']==0 and r['baseline_gradient']==[[2.,3.,4.]];checks.append('full gradient and zero SCF counts parsed')
  p.write_text(raw.replace('0.1 0.3 0.5 0.2','0.11 0.3 0.5 0.2'))
  try:PROTOCOL.result(p,b,u,v)
  except AssertionError:checks.append('changed supplied u rejected')
  else:raise AssertionError('changed supplied u accepted')
  assert slope(.0001,-.0001,.0001,[[1.,0.,0.]],[[1.,0.,0.]])['passed'];assert not slope(.0002,-.0002,.0001,[[1.,0.,0.]],[[1.,0.,0.]])['passed'];checks.append('declared FD gate passes and fails without padding')
 r={'status':'passed_manufactured_protocol_controls','implementation':pin(__file__),'protocol':pin(Path(__file__).with_name('PROTOCOL.py')),'collector':pin(Path(__file__).with_name('COLLECT.py')),'checks':checks,'count':len(checks),'native_calls':0,'operator_solves':0,'actual_fixture_evaluations':0,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1,'host':platform.node()}};save(out/'RESULT.json',r);print(json.dumps(r['cost']))
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.out)
