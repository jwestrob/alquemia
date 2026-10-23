#!/usr/bin/env python3
"""Extract exact published Tinker assets and shared fixed physical fixtures; no energy calls."""
import argparse, datetime, hashlib, json, re, shutil, socket, time
from pathlib import Path

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--panel',required=True); p.add_argument('--si',required=True)
    p.add_argument('--dy-parameters',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    start=time.monotonic(); out=Path(a.output).resolve(); out.mkdir(parents=True,exist_ok=False)
    panel=json.loads(Path(a.panel).read_text()); si=Path(a.si).read_text(); dy=Path(a.dy_parameters).read_text().splitlines()
    assert len(panel['geometries'])==28 and len(panel['free_fragments'])==4
    assert pin(a.si)==panel['published_coordinate_source']
    assert '# AMOEBA22 *.key file' in si and '# AMOEBA09.prm' in si
    prm=si.split('# AMOEBA09.prm',1)[1]
    key=si.split('# AMOEBA22 *.key file',1)[1].split('###############################################################################',1)[0]
    key=key.replace('parameters    ./amoeba09.prm','parameters amoeba09.prm')
    assert re.search(r'^atom\s+36\s+34\s+O\s+"Water O"',prm,re.M)
    assert 'polpair   105  500     0.2987' in key and 'polpair   403  500     0.1518' in key
    ion=[]
    for i,line in enumerate(dy):
        tok=line.split()
        if len(tok)>1 and tok[1]=='10' and tok[0] in ('atom','vdw','polarize','multipole'):
            line=re.sub(r'^(\w+)\s+10\b',r'\1 501',line)
            if tok[0]=='atom': line=re.sub(r'^(atom\s+501)\s+10\b',r'\1 501',line)
            ion.append(line)
            if tok[0]=='multipole': ion+=dy[i+1:i+5]
    assert len(ion)==8 and '0.7280' in '\n'.join(ion) and '15.0000' in '\n'.join(ion)
    key+='\n# Dy from author-contributed amoebail22; unqualified common-water/donor transfer\n'+'\n'.join(ion)+'\n'
    key+='\npolar-eps 1.0e-8\ndigits 10\nopenmp-threads 1\n'
    (out/'amoeba09.prm').write_text(prm); (out/'candidate.key').write_text(key)
    rows=[]
    for g in panel['geometries']:
        assert g['metal_index']==0 and g['metal'] in ('La','Dy') and g['atoms']==len(g['symbols'])==len(g['positions_A'])
        assert pin(g['xyz']['path'])==g['xyz']
        types=list(g.get('source_tinker_types',[])); bonds=g.get('source_bonds_zero_based')
        if not types:
            assert g['kind']=='constructed_aquo' and g['atoms']==1+3*g['water_count']
            types=[500]+[v for _ in range(g['water_count']) for v in (36,37,37)]
            bonds=[[]]
            for k in range(g['water_count']):
                o=1+3*k; bonds += [[o+1,o+2],[o],[o]]
        types[0]=500 if g['metal']=='La' else 501
        d=out/g['geometry_id']; d.mkdir()
        (d/'system.xyz').write_text(str(g['atoms'])+'  '+g['geometry_id']+'\n'+''.join(
            f'{i+1:5d} {s:3s} '+''.join(f'{v:20.12f}' for v in x)+f' {t:5d}'+''.join(f' {b+1:5d}' for b in bs)+'\n'
            for i,(s,x,t,bs) in enumerate(zip(g['symbols'],g['positions_A'],types,bonds))))
        shutil.copyfile(out/'amoeba09.prm',d/'amoeba09.prm'); shutil.copyfile(out/'candidate.key',d/'system.key')
        rows.append(dict(id=g['geometry_id'],case_id=g['case_id'],metal=g['metal'],kind=g['kind'],charge=g['charge'],
            atoms=g['atoms'],source_geometry=g,xyz=pin(d/'system.xyz'),key=pin(d/'system.key'),parameters=pin(d/'amoeba09.prm'),
            types=types,bonds_zero_based=bonds,fragment_id=g['case_id']+'_free' if g['kind']=='published_dimer' else None,
            derivative_test=g['distortion']=='central' and g['case_id'] in ('water','acetate_bi')))
    for g in panel['free_fragments']:
        src=next(r for r in rows if r['case_id']==g['case_id'] and r['metal']=='La' and r['source_geometry']['distortion']=='central')
        sg=src['source_geometry']; d=out/g['geometry_id']; d.mkdir(); assert pin(g['xyz']['path'])==g['xyz']
        types=src['types'][1:]; bonds=[[j-1 for j in b] for b in src['bonds_zero_based'][1:]]
        (d/'system.xyz').write_text(str(g['atoms'])+'  '+g['geometry_id']+'\n'+''.join(
            f'{i+1:5d} {s:3s} '+''.join(f'{v:20.12f}' for v in x)+f' {t:5d}'+''.join(f' {b+1:5d}' for b in bs)+'\n'
            for i,(s,x,t,bs) in enumerate(zip(sg['symbols'][1:],sg['positions_A'][1:],types,bonds))))
        shutil.copyfile(out/'amoeba09.prm',d/'amoeba09.prm'); shutil.copyfile(out/'candidate.key',d/'system.key')
        rows.append(dict(id=g['geometry_id'],case_id=g['case_id'],metal=None,kind='free_fragment',charge=g['charge'],atoms=g['atoms'],
            source_geometry=g,xyz=pin(d/'system.xyz'),key=pin(d/'system.key'),parameters=pin(d/'amoeba09.prm'),types=types,
            bonds_zero_based=bonds,fragment_id=None,derivative_test=False))
    manifest=dict(protocol='explicit_model_v1',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        implementation=pin(__file__),panel=pin(a.panel),published_La_SI=pin(a.si),Dy_author_parameters=pin(a.dy_parameters),
        model='La2023_AMOEBA22_and_Dy_IL_ion_common_SI_water_transfer_candidate',transfer_qualified=False,
        implicit_solvent=False,ALPB=False,polarization='mutual induced dipoles',polar_eps=1e-8,
        rows=rows,declared_geometry_cells=28,free_fragments=4,local_preparation=dict(host=socket.gethostname(),
        elapsed_seconds=time.monotonic()-start,threads=1,molecular_calls=0))
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps(pin(out/'MANIFEST.json')))

if __name__=='__main__': main()
