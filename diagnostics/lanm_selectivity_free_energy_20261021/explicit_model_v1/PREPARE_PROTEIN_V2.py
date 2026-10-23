#!/usr/bin/env python3
"""Exact 9CEQ model-1 typing and unequilibrated explicit-water geometry; no energy."""
import argparse, datetime, hashlib, json, math, os, re, resource, shlex, shutil, subprocess, time
from pathlib import Path

def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def save(p,x): Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def cif_atoms(path):
    lines=Path(path).read_text().splitlines(); i=next(i for i,l in enumerate(lines) if l.strip()=='_atom_site.group_PDB')
    keys=[]
    while lines[i].startswith('_atom_site.'):
        keys.append(lines[i].strip().split('.')[1]); i+=1
    rows=[]
    while i<len(lines) and not lines[i].startswith('#'):
        t=shlex.split(lines[i]); assert len(t)==len(keys); i+=1
        r=dict(zip(keys,t))
        if r['pdbx_PDB_model_num']=='1':
            assert r['label_alt_id']=='.' and float(r['occupancy'])==1
            r['xyz']=[float(r[k]) for k in ('Cartn_x','Cartn_y','Cartn_z')]; rows.append(r)
    return rows

def xyz_atoms(path):
    lines=Path(path).read_text().splitlines(); n=int(lines[0].split()[0]); rows=[]
    assert len(lines)==n+1
    for i,line in enumerate(lines[1:]):
        t=line.split(); assert int(t[0])==i+1
        rows.append({'label':t[1],'xyz':list(map(float,t[2:5])),'type':int(t[5]),'bonds':[int(v)-1 for v in t[6:]]})
    return rows

def write_xyz(path,rows,title,box=None):
    text=f'{len(rows)}  {title}\n'
    if box is not None: text+=f'{box:16.8f} {box:16.8f} {box:16.8f} 90.0 90.0 90.0\n'
    text+=''.join(f'{i+1:6d} {r["label"]:3s} '+''.join(f'{v:18.10f}' for v in r['xyz'])+f' {r["type"]:5d}'+''.join(f' {v+1:6d}' for v in r['bonds'])+'\n' for i,r in enumerate(rows))
    Path(path).write_text(text)

def run(argv,cwd,output,receipts):
    start=time.monotonic(); before=resource.getrusage(resource.RUSAGE_CHILDREN)
    r={'argv':list(map(str,argv)),'cwd':str(cwd),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'threads':1,'energy_calls':0}
    receipts.append(r)
    try:
        with open(output,'w') as f:
            p=subprocess.run(argv,cwd=cwd,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,timeout=60)
        r['returncode']=p.returncode; assert p.returncode==0
    finally:
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        r.update(elapsed_seconds=time.monotonic()-start,user_cpu_seconds=after.ru_utime-before.ru_utime,
                 system_cpu_seconds=after.ru_stime-before.ru_stime,output=pin(output))

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--cif',required=True);p.add_argument('--tinker-source',required=True);p.add_argument('--tinker-build',required=True);p.add_argument('--previous',required=True);a=p.parse_args()
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False); src=Path(a.tinker_source).resolve();build=Path(a.tinker_build).resolve()
    start=time.monotonic(); result={'status':'unavailable','implementation':pin(__file__),'source_cif':pin(a.cif),'receipts':[],
        'energy_calls':0,'force_calls':0,'minimization_steps':0,'dynamics_steps':0,'production_ready':False}
    try:
        assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
        assert result['source_cif']['sha256']=='956ed5bf63bc4e4c05c223d9e34912a0f57c0d1b757a1648c7c385cc4fab967b'
        atoms=cif_atoms(a.cif); metals=[r for r in atoms if r['type_symbol']=='LA']; assert len(metals)==1
        prot=[r for r in atoms if r['label_asym_id']=='A']; assert len(prot)==257 and len(atoms)==258
        expected=['PHE','ILE','ASP','THR','ASN','ASN','ASP','GLY','TRP','ILE','GLU','GLY','ASP','GLU','LEU','LEU','ALA','NH2']
        assert all({r['label_comp_id'] for r in prot if int(r['label_seq_id'])==i+1}=={name} for i,name in enumerate(expected))
        assert {r['label_atom_id'] for r in prot if r['label_comp_id']=='NH2'}=={'N','HN1','HN2'}
        prm=src/'params/amoebabio09.prm';shutil.copyfile(prm,out/'amoebabio09.prm')
        result['protein_parameters']=pin(prm);result['executables']={n:pin(build/n) for n in ('pdbxyz','analyze')}
        lines=['REMARK 9CEQ model 1; complete amidated LBT3, metal separately retained']
        rename=[]
        for i,r in enumerate(prot):
            name=r['label_atom_id'];res=r['label_comp_id']
            if res=='NH2' and name in ('HN1','HN2'):
                original=name;name={'HN1':'H1','HN2':'H2'}[name];rename.append({'source_id':r['id'],'from':original,'to':name})
            atom_name=(' '+name).ljust(4) if len(name)<4 else name
            x,y,z=r['xyz'];lines.append(f'ATOM  {i+1:5d} {atom_name} {res:3s} A{int(r["label_seq_id"]):4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {r["type_symbol"]:>2s}')
        lines+=['TER','END'];(out/'protein.pdb').write_text('\n'.join(lines)+'\n')
        (out/'protein.key').write_text('parameters amoebabio09.prm\nopenmp-threads 1\ndigits 10\n')
        previous=Path(a.previous).resolve(); old=json.loads((previous/'RESULT.json').read_text())
        assert pin(previous/'RESULT.json')['sha256']=='ff567c1eede1a1310e7d878f9725d48f0c28159fc3d529cfe82471a31d82c033'
        assert old['status']=='unavailable' and old['error']=="ValueError('No closing quotation')"
        assert old['source_cif']==result['source_cif'] and old['protein_parameters']==result['protein_parameters']
        assert old['executables']==result['executables']
        for name in ('protein.pdb','protein.key','amoebabio09.prm'):
            assert (out/name).read_bytes()==(previous/name).read_bytes()
        for name in ('protein.xyz','protein.seq','pdbxyz.out'):shutil.copyfile(previous/name,out/name)
        result['reused_native_topology']={'parent_result':pin(previous/'RESULT.json'),'xyz':pin(previous/'protein.xyz'),'sequence':pin(previous/'protein.seq'),'original_receipts':old['receipts']}
        xyz=xyz_atoms(out/'protein.xyz')
        defs={}
        for line in prm.read_text().splitlines():
            if not re.match(r'^\s*atom\s+',line):continue
            t=shlex.split(line)
            if t and t[0]=='atom':defs[int(t[1])]={'class':int(t[2]),'Z':int(t[5]),'description':t[4]}
        assert all(r['type'] in defs and r['type']>0 for r in xyz)
        for i,r in enumerate(xyz):
            assert all(0<=j<len(xyz) and i in xyz[j]['bonds'] for j in r['bonds'])
        joins=[];used=set()
        for r in prot:
            matches=[i for i,x in enumerate(xyz) if all(abs(v-w)<=1e-8 for v,w in zip(r['xyz'],x['xyz']))]
            if r['type_symbol']!='H': assert len(matches)==1, f'Heavy atom not preserved: {r}'
            if matches:
                assert len(matches)==1 and matches[0] not in used;used.add(matches[0])
                joins.append({'source_id':r['id'],'residue':r['label_comp_id'],'residue_number':int(r['label_seq_id']),'atom':r['label_atom_id'],'element':r['type_symbol'],'tinker_index':matches[0]+1,'type':xyz[matches[0]]['type'],'definition':defs[xyz[matches[0]]['type']]})
            else:
                joins.append({'source_id':r['id'],'residue':r['label_comp_id'],'residue_number':int(r['label_seq_id']),'atom':r['label_atom_id'],'element':r['type_symbol'],'tinker_index':None,'status':'deposited_hydrogen_replaced_by_native_geometry'})
        assert len(xyz)==257
        run([build/'analyze','protein.xyz','GP'],out,out/'protein_parameters.out',result['receipts'])
        parameter_text=(out/'protein_parameters.out').read_text()
        assert 'Total Potential Energy' not in parameter_text and 'Undefined' not in parameter_text
        assert 'Atomic Multipole Parameters' in parameter_text
        result.update(protein_atoms=len(xyz),deposited_atoms_retained_by_coordinates=len(used),generated_hydrogens=len(xyz)-len(used),
                      source_joins=joins,cap_name_changes=rename,protein_xyz=pin(out/'protein.xyz'),sequence=pin(out/'protein.seq'))
        # Only coordinates are used from the native TIP3P example; no TIP3P parameters are imported.
        water_source=src/'example/waterbox.xyz';waters=xyz_atoms(water_source);assert len(waters)==1500
        L0=24.662;L=2*L0
        center=[(max(r['xyz'][k] for r in xyz)+min(r['xyz'][k] for r in xyz))/2 for k in range(3)]
        solute=[]
        for r in xyz:solute.append({**r,'xyz':[r['xyz'][k]-center[k] for k in range(3)]})
        metal_xyz=[metals[0]['xyz'][k]-center[k] for k in range(3)]
        heavy=[r['xyz'] for r in solute if defs[r['type']]['Z']!=1]+[metal_xyz]
        hydrogen=[r['xyz'] for r in solute if defs[r['type']]['Z']==1]
        def d2(a,b):return sum((x-y-L*round((x-y)/L))**2 for x,y in zip(a,b))
        candidates=[]
        for ix in (0,1):
          for iy in (0,1):
           for iz in (0,1):
            shift=[(v-.5)*L0 for v in (ix,iy,iz)]
            for j in range(0,len(waters),3):
                group=waters[j:j+3];assert [r['type'] for r in group]==[1,2,2]
                o=group[0]['xyz'];pts=[[o[k]+shift[k] for k in range(3)]]
                for h in group[1:]:pts.append([pts[0][k]+h['xyz'][k]-o[k]-L0*round((h['xyz'][k]-o[k])/L0) for k in range(3)])
                if any(d2(pts[0],v)<2.4**2 for v in heavy):continue
                if any(d2(h,v)<1.5**2 for h in pts[1:] for v in heavy):continue
                if any(d2(h,v)<1.0 for h in pts[1:] for v in hydrogen):continue
                nearest=min(d2(pts[0],v) for v in heavy)
                candidates.append({'xyz':pts,'min_heavy_distance_squared':nearest,'source_water_index':j//3,'tile':[ix,iy,iz]})
        selected=[]
        for idx in sorted(range(len(candidates)),key=lambda i:(-candidates[i]['min_heavy_distance_squared'],i)):
            if all(d2(candidates[idx]['xyz'][0],candidates[k]['xyz'][0])>=25 for k in selected):selected.append(idx)
            if len(selected)==15:break
        assert len(selected)==15
        common=[]
        for i,w in enumerate(candidates):
            if i in selected:continue
            n=len(solute)+1+len(common)
            common.extend([{'label':'O','xyz':w['xyz'][0],'type':402,'bonds':[n+1,n+2]},
                           {'label':'H','xyz':w['xyz'][1],'type':403,'bonds':[n]},
                           {'label':'H','xyz':w['xyz'][2],'type':403,'bonds':[n]}])
        for i,k in enumerate(selected):common.append({'label':'Na' if i<8 else 'Cl','xyz':candidates[k]['xyz'][0],'type':405 if i<8 else 414,'bonds':[]})
        boxes={}
        for metal,typ in [('La',900),('Dy',901)]:
            rows=solute+[{'label':metal,'xyz':metal_xyz,'type':typ,'bonds':[]}]+common
            for i,r in enumerate(rows):assert all(0<=j<len(rows) and i in rows[j]['bonds'] for j in r['bonds'])
            p=out/(metal+'_explicit_water.xyz');write_xyz(p,rows,'PREPARATION ONLY; donor/metal transfer and sampling unqualified',L);boxes[metal]=pin(p)
        result.update(status='geometry_prepared_model_unqualified',boxes=boxes,box_length_A=L,water_source=pin(water_source),
            water_count=len(candidates)-15,Na_count=8,Cl_count=7,total_atoms=len(solute)+1+len(common),
            declared_peptide_charge=-4,declared_complex_charge=-1,declared_box_charge=0,
            salt_pairs=7,nominal_NaCl_M=7/(6.02214076e23*(L**3)*1e-27),
            padding_A=[L/2-max(abs(r['xyz'][k]) for r in solute) for k in range(3)],
            solute_translation_A=[-v for v in center],ion_replacement_sites=[candidates[k] for k in selected],
            paired_geometry_identical_except_metal_type=True,metal_parameters_selected=False,
            caveats=['No protein Hamiltonian selected: donor correction transfer must be qualified.',
                     'TIP3P packing coordinates are unequilibrated AMOEBA initial water positions.',
                     'No MES buffer represented; seven salt pairs and one counterion are an initial finite-box choice.',
                     'No energy, force, equilibration, free energy or accuracy result.'])
    except Exception as e:result['error']=repr(e)
    result['elapsed_seconds']=time.monotonic()-start;save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'error':result.get('error')}))
    return 0 if result['status']=='geometry_prepared_model_unqualified' else 1

if __name__=='__main__':raise SystemExit(main())
