#!/usr/bin/env python3
"""Four exact-source capped peptide contexts, no energy or force evaluations."""
import argparse, hashlib, importlib.util, json, os, re, shlex, shutil, time
from decimal import Decimal
from pathlib import Path

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--cif',required=True);p.add_argument('--parameters',required=True);p.add_argument('--build',required=True);a=p.parse_args()
    leaf=Path(__file__).resolve().parent;helper=leaf/'PREPARE_PROTEIN_V2.py'
    sp=importlib.util.spec_from_file_location('protein_geometry_helper',helper);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);build=Path(a.build).resolve();prm=Path(a.parameters).resolve();start=time.monotonic()
    result={'status':'unavailable','implementation':pin(__file__),'helper':pin(helper),'specification':pin(leaf/'PEPTIDE_FRAGMENT_SPEC.md'),
            'source_cif':pin(a.cif),'parameters':pin(prm),'receipts':[],'fragments':[],'geometries':[],
            'energy_calls':0,'force_calls':0,'optimization_steps':0,'model_selected':False}
    try:
        assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
        assert result['source_cif']['sha256']=='956ed5bf63bc4e4c05c223d9e34912a0f57c0d1b757a1648c7c385cc4fab967b'
        atoms=h.cif_atoms(a.cif);metal=next(r for r in atoms if r['type_symbol']=='LA')
        definitions={}
        for line in prm.read_text().splitlines():
            if re.match(r'^\s*atom\s+',line):
                t=shlex.split(line);definitions[int(t[1])]={'Z':int(t[5]),'description':t[4],'class':int(t[2])}
        symbol={1:'H',6:'C',7:'N',8:'O'}
        for residue,name,charge,donor_names,donor_type in [(3,'ASP',-1,['OD1','OD2'],143),(11,'GLU',-1,['OE1','OE2'],161),(5,'ASN',0,['OD1'],153),(9,'TRP',0,['O'],11)]:
            cid=f'ACE_{name}{residue}_NME';d=out/cid;d.mkdir()
            chain=[r for r in atoms if r['label_asym_id']=='A'];central=[r for r in chain if int(r['label_seq_id'])==residue]
            before={r['label_atom_id']:r for r in chain if int(r['label_seq_id'])==residue-1};after={r['label_atom_id']:r for r in chain if int(r['label_seq_id'])==residue+1}
            assert central and all(r['label_comp_id']==name for r in central)
            entries=[(before['CA'],'ACE',1,'CH3'),(before['C'],'ACE',1,'C'),(before['O'],'ACE',1,'O')]
            entries += [(r,name,2,r['label_atom_id']) for r in central]
            entries += [(after['N'],'NME',3,'N'),(after['CA'],'NME',3,'C')]
            if 'H' in after:entries.append((after['H'],'NME',3,'H'))
            lines=[f'REMARK {cid}; geometry-only qualification candidate from exact 9CEQ model1']
            for i,(r,res,number,atom_name) in enumerate(entries):
                an=(' '+atom_name).ljust(4) if len(atom_name)<4 else atom_name;x,y,z=r['xyz']
                lines.append(f'ATOM  {i+1:5d} {an} {res:3s} A{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {r["type_symbol"]:>2s}')
            (d/'fragment.pdb').write_text('\n'.join(lines+['TER','END'])+'\n');shutil.copyfile(prm,d/'amoebabio09.prm')
            (d/'fragment.key').write_text('parameters amoebabio09.prm\nopenmp-threads 1\ndigits 10\n')
            h.run([build/'pdbxyz','fragment.pdb'],d,d/'pdbxyz.out',result['receipts'])
            typed=h.xyz_atoms(d/'fragment.xyz');assert all(r['type'] in definitions for r in typed)
            joins=[];used=set();donors=[]
            for r,res,n,an in entries:
                matches=[i for i,t in enumerate(typed) if t['xyz']==r['xyz']];assert len(matches)==1
                idx=matches[0];assert idx not in used;used.add(idx)
                joins.append({'source_id':r['id'],'source_residue':r['label_comp_id'],'source_residue_number':int(r['label_seq_id']),
                              'source_name':r['label_atom_id'],'fragment_residue':res,'fragment_name':an,'index':idx,'type':typed[idx]['type']})
                if n==2 and an in donor_names:donors.append(idx);assert typed[idx]['type']==donor_type
            assert len(donors)==len(donor_names)
            h.run([build/'analyze','fragment.xyz','GP'],d,d/'parameters.out',result['receipts'])
            text=(d/'parameters.out').read_text();assert 'Total Potential Energy' not in text and 'Undefined' not in text
            section=text.split('Atomic Multipole Parameters :',1)[1].split('Dipole Polarizability Parameters :',1)[0]
            charges=[]
            for line in section.splitlines():
                t=line.split()
                if len(t)>=5 and t[0].isdigit() and t[1].isdigit() and any(x in line for x in ['Z-then-X','Bisector','Z-Bisect','3-Fold','Z-Only','None']):charges.append((int(t[1]),Decimal(t[-1])))
            assert [i for i,q in charges]==list(range(1,len(typed)+1))
            qsum=sum(q for i,q in charges);envelope=Decimal(len(charges))*Decimal('0.000005');assert abs(qsum-Decimal(charge))<=envelope
            f={'id':cid,'charge':charge,'atoms':len(typed),'source_joins':joins,'constructed_atoms':[{'index':i,'type':t['type'],'element':symbol[definitions[t['type']]['Z']]} for i,t in enumerate(typed) if i not in used],
               'donor_indices':donors,'donor_type':donor_type,'native_printed_monopole_sum':str(qsum),'printed_rounding_envelope':str(envelope),
               'xyz':pin(d/'fragment.xyz'),'key':pin(d/'fragment.key'),'parameter_output':pin(d/'parameters.out'),'sequence':pin(d/'fragment.seq')}
            assert all(x['element']=='H' for x in f['constructed_atoms']);result['fragments'].append(f)
            for m,typ,mult in [('La',900,1),('Dy',901,6)]:
                rows=[{'label':m,'xyz':metal['xyz'],'type':typ,'bonds':[]}]+[{**r,'bonds':[j+1 for j in r['bonds']]} for r in typed]
                path=d/(m+'.txyz');h.write_xyz(path,rows,cid+'_'+m+'; preparation only')
                symbols=[m]+[symbol[definitions[r['type']]['Z']] for r in typed]
                xyz=d/(m+'.xyz');xyz.write_text(f'{len(rows)}\n{cid}_{m}; no optimization\n'+''.join(s+' '+ ' '.join(f'{v:.10f}' for v in r['xyz'])+'\n' for s,r in zip(symbols,rows)))
                result['geometries'].append({'geometry_id':cid+'_'+m,'case_id':cid,'metal':m,'metal_index':0,'charge':charge+3,
                    'physical_multiplicity':mult,'classical_spin_degree_of_freedom':False,'atoms':len(rows),'symbols':symbols,'positions_A':[r['xyz'] for r in rows],
                    'source_tinker_types':[r['type'] for r in rows],'source_bonds_zero_based':[r['bonds'] for r in rows],'donor_indices':[i+1 for i in donors],
                    'xyz':pin(xyz),'tinker_xyz':pin(path),'parent_fragment':f['xyz'],'geometry_convention':'exact deposited center residue and metal; neighbor-defined neutral caps, native constructed methyl H'})
        assert len(result['fragments'])==4 and len(result['geometries'])==8
        result['status']='prepared_no_scoring_release';result['executables']={n:pin(build/n) for n in ('pdbxyz','analyze')}
    except Exception as e:result['error']=repr(e)
    result['elapsed_seconds']=time.monotonic()-start;h.save(out/'PANEL.json',result);print(json.dumps({'status':result['status'],'panel':pin(out/'PANEL.json'),'error':result.get('error')}))
    return 0 if result['status']=='prepared_no_scoring_release' else 1

if __name__=='__main__':raise SystemExit(main())
