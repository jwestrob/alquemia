#!/usr/bin/env python3
"""Bind exactly12 existing noncentral geometries; zero engine calls."""
import argparse,copy,hashlib,importlib.util,json,os,shutil,time
from pathlib import Path
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();t=time.monotonic();b=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
 source=b/'diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/REFERENCE.py';assert pin(source)['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129'
 sp=importlib.util.spec_from_file_location('frozen_reference',source);ref=importlib.util.module_from_spec(sp);sp.loader.exec_module(ref)
 parent=b/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/bound_v1/NATIVE_PANEL.json';old=ref.load_panel(parent)
 p=b/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/prepared_v2/PANEL.json';assert pin(p)['sha256']=='1191e3a00e3c02aa861faf83aae22eb250c94064e7ac780bf2fbebb00eb6c21a';physical=json.loads(p.read_text());assert physical['status']=='prepared' and physical['contact_failures']==[]
 gs=[copy.deepcopy(g) for g in physical['geometries'] if g['distortion'] in ('minus','plus')];assert len(gs)==12
 assert len({g['geometry_id'] for g in gs})==12 and all(g['geometry_id'] not in {c['geometry_id'] for c in old['cells']} for g in gs)
 (out/'xyz').mkdir();cells=[]
 for g in gs:
  src=ref.verify(g['xyz']);dest=out/'xyz'/src.name;shutil.copyfile(src,dest);assert pin(dest)['sha256']==g['xyz']['sha256'];g['source_xyz']=g['xyz'];g['xyz']=pin(dest)
  core=46 if g['metal']=='La' else 55;ne=sum(ref.Z[s] for s in g['symbols'])-g['charge']-core;assert ne==(100 if g['donor_family']=='acetate_bi' else 110)
  cells.append({'cell_id':'ecp_pbe0__'+g['geometry_id'],'recipe':'ecp_pbe0','geometry_id':g['geometry_id'],'charge':g['charge'],'multiplicity':1,'physical_multiplicity':g['physical_multiplicity'],'metal':g['metal'],'ecp_core':core,'explicit_electrons':ne,'geometry':g,'pilot':True})
 native={k:copy.deepcopy(old[k]) for k in ('protocol','adapter','published_coordinate_source','archived_recipe_source','orca','assets')}
 native.update(created_utc=ref.now(),specification=pin(Path(__file__).with_name('TRAINING_RADIAL_SPEC.md')),parent_panel=pin(parent),source_inventory=pin(p),geometries=gs,free_fragments=[],cells=cells,pilot_cell_ids=[c['cell_id'] for c in cells],interpretation='Consumed development training only; exact12existing noncentral radial geometries; no native calls in preparation')
 np=ref.write(out/'NATIVE_PANEL.json',native);vp=ref.write(out/'VALIDATION.json',ref.validate(np['path']))
 receipt={'status':'prepared_no_molecular_calls','implementation':pin(__file__),'native_panel':np,'saved_validation':vp,'native_cells':12,'new_molecular_calls':0,'local_threads':1,'elapsed_seconds':time.monotonic()-t,'native_command_template':['/usr/bin/python3',str(source),'run','--panel',np['path'],'--cell','CELL_ID','--out','EXECUTION_DIRECTORY'],'cell_ids':[c['cell_id'] for c in cells]}
 ref.write(out/'READY.json',receipt);print(json.dumps(receipt))
if __name__=='__main__':main()
