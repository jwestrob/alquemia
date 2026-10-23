#!/usr/bin/env python3
"""Repair only constructed cap hydrogens; four parameter-only invocations."""
import argparse,copy,hashlib,importlib.util,json,math,os,shutil,time
from pathlib import Path
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):assert pin(p['path'])==p;return Path(p['path'])
def sub(a,b):return [x-y for x,y in zip(a,b)]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def norm(v):return math.sqrt(dot(v,v))
def unit(v):return [x/norm(v) for x in v]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def angle(a,b):return math.degrees(math.acos(max(-1,min(1,dot(unit(a),unit(b))))))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();start=time.monotonic();old=json.loads(Path(a.parent).read_text());assert pin(a.parent)['sha256']=='42db5a66e95edf2af5a5fd1eb77c779dc18cac4bcfeaadde39d6beba9fe71b67'
 assert old['energy_calls']==old['force_calls']==0
 helper=verify(old['helper']);sp=importlib.util.spec_from_file_location('protein_saved_helper',helper);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);result={k:copy.deepcopy(old[k]) for k in ('source_cif','parameters','declared_fragment_ids','declared_geometry_ids','executables')}
 result.update(status='unavailable',implementation=pin(__file__),specification=pin(Path(__file__).with_name('CAP_REPAIR_SPEC_V3.md')),helper=old['helper'],parent=pin(a.parent),receipts=[],fragments=[],geometries=[],energy_calls=0,force_calls=0,optimization_steps=0,model_selected=False)
 withdrawal={'parent':pin(a.parent),'status':'withdrawn_from_reference_and_scoring_use','reason':'NME constructed methyl H overlaps retained H1 despite distinct coordinates','native_calls':0,'fragments':[]}
 try:
  assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
  for p in old['executables'].values():verify(p)
  for f in old['fragments']:
   typed=h.xyz_atoms(verify(f['xyz']));original=copy.deepcopy(typed);joins=f['source_joins']
   carbon=next(j['index'] for j in joins if j['fragment_residue']=='NME' and j['fragment_name']=='C')
   nitrogen=next(j['index'] for j in joins if j['fragment_residue']=='NME' and j['fragment_name']=='N')
   h1=next(j['index'] for j in joins if j['fragment_residue']=='NME' and j['fragment_name']=='H1')
   changed=[r['index'] for r in f['constructed_atoms'] if typed[r['index']]['bonds']==[carbon]];assert len(changed)==2
   hs=[h1]+changed;withdrawal['fragments'].append({'id':f['id'],'methyl_H_indices':hs,'minimum_methyl_HH_A':min(norm(sub(typed[i]['xyz'],typed[j]['xyz'])) for k,i in enumerate(hs) for j in hs[k+1:])})
   axis=unit(sub(typed[nitrogen]['xyz'],typed[carbon]['xyz']));u=unit(sub(typed[h1]['xyz'],typed[carbon]['xyz']));ux=cross(axis,u);projection=dot(axis,u)
   for i,phi in zip(changed,(2*math.pi/3,-2*math.pi/3)):
    v=[u[k]*math.cos(phi)+ux[k]*math.sin(phi)+axis[k]*projection*(1-math.cos(phi)) for k in range(3)]
    typed[i]['xyz']=[typed[carbon]['xyz'][k]+1.10*v[k] for k in range(3)]
   d=out/f['id'];d.mkdir();h.write_xyz(d/'fragment.xyz',typed,f['id']+' cap-onlyV3; no energy');typed=h.xyz_atoms(d/'fragment.xyz')
   for i,(x,y) in enumerate(zip(typed,original)):
    assert x['type']==y['type'] and x['bonds']==y['bonds']
    if i not in changed:assert x['xyz']==y['xyz']
   lengths=[norm(sub(typed[i]['xyz'],typed[carbon]['xyz'])) for i in changed];assert all(abs(x-1.10)<=1e-6 for x in lengths)
   hh=[norm(sub(typed[i]['xyz'],typed[j]['xyz'])) for k,i in enumerate(hs) for j in hs[k+1:]];assert min(hh)>1.5
   vectors=[sub(typed[i]['xyz'],typed[carbon]['xyz']) for i in [nitrogen]+hs];angles=[angle(x,y) for k,x in enumerate(vectors) for y in vectors[k+1:]];assert all(90<x<130 for x in angles)
   # The complete cap-hydrogen screen includes unchanged ACE hydrogens too.
   parent_g=next(g for g in old['geometries'] if g['geometry_id']==f['id']+'_La');symbols=parent_g['symbols'][1:];contacts=[]
   for ca in f['constructed_atoms']:
    i=ca['index']
    for j in range(len(typed)):
     if i==j or j in typed[i]['bonds']:continue
     dist=norm(sub(typed[i]['xyz'],typed[j]['xyz']));limit=.8 if symbols[j]=='H' else 1.1
     contacts.append((dist,i,j,limit));assert dist>=limit,(f['id'],dist,i,j,limit)
   shutil.copyfile(verify(old['parameters']),d/'amoebabio09.prm');shutil.copyfile(verify(f['key']),d/'fragment.key');shutil.copyfile(verify(f['sequence']),d/'fragment.seq')
   h.run([old['executables']['analyze']['path'],'fragment.xyz','GP'],d,d/'parameters.out',result['receipts'])
   text=(d/'parameters.out').read_text();previous=verify(f['parameter_output']).read_text();assert 'Total Potential Energy' not in text and 'undefined' not in text.lower()
   def multipoles(s):return s.split('Atomic Multipole Parameters :',1)[1].split('Dipole Polarizability Parameters :',1)[0]
   assert multipoles(text)==multipoles(previous)
   updated={**f,'xyz':pin(d/'fragment.xyz'),'key':pin(d/'fragment.key'),'parameter_output':pin(d/'parameters.out'),'sequence':pin(d/'fragment.seq'),'v2_xyz':f['xyz'],'changed_constructed_indices':changed,'methyl_HH_A':hh,'methyl_angles_deg':angles,'minimum_constructed_H_nonbonded_contact':list(min(contacts)),'exact_native_local_multipoles_equal_v2':True}
   updated['constructed_cap_bonds']=[{'atom_index':r['index'],'parent_index':typed[r['index']]['bonds'][0],'length_A':norm(sub(typed[r['index']]['xyz'],typed[typed[r['index']]['bonds'][0]]['xyz']))} for r in f['constructed_atoms']]
   result['fragments'].append(updated)
   for metal in ('La','Dy'):
    pg=next(g for g in old['geometries'] if g['geometry_id']==f['id']+'_'+metal);rows=[{'label':metal,'xyz':pg['positions_A'][0],'type':pg['source_tinker_types'][0],'bonds':[]}]+[{**t,'bonds':[j+1 for j in t['bonds']]} for t in typed]
    h.write_xyz(d/(metal+'.txyz'),rows,pg['geometry_id']+' cap-onlyV3');xp=d/(metal+'.xyz');xp.write_text(f'{len(rows)}\n{pg["geometry_id"]}; only constructed NME H repaired\n'+''.join(s+' '+' '.join(f'{v:.10f}' for v in row['xyz'])+'\n' for s,row in zip(pg['symbols'],rows)))
    positions=[[float(x) for x in line.split()[1:]] for line in xp.read_text().splitlines()[2:]]
    assert all(positions[i]==pg['positions_A'][i] for i in range(len(positions)) if i not in [j+1 for j in changed])
    result['geometries'].append({**pg,'positions_A':positions,'xyz':pin(xp),'tinker_xyz':pin(d/(metal+'.txyz')),'parent_fragment':updated['xyz'],'v2_geometry':pg['xyz'],'geometry_convention':'exact deposited source; only two constructed NME capH repaired geometrically'})
  assert len(result['fragments'])==4 and len(result['geometries'])==8
  result['status']='prepared_no_scoring_release'
 except Exception as e:result['error']=repr(e)
 result['elapsed_seconds']=time.monotonic()-start;h.save(out/'V2_WITHDRAWAL.json',withdrawal);result['withdrawal']=pin(out/'V2_WITHDRAWAL.json');h.save(out/'PANEL.json',result);print(json.dumps({'status':result['status'],'panel':pin(out/'PANEL.json'),'error':result.get('error')}))
 return 0 if result['status']=='prepared_no_scoring_release' else 1
if __name__=='__main__':raise SystemExit(main())
