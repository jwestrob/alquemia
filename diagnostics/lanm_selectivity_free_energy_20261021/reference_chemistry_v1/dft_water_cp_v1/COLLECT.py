#!/usr/bin/env python3
"""Pinned saved-only final comparison; no native or model invocation."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('dft_cp',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
b=m.b

def canonical(value):return json.loads(json.dumps(value))

def same(actual,expected,label):b.check(canonical(actual)==canonical(expected),label)

def convergence(text):
    native=text.split('****END OF INPUT****')[-1]
    return {'mode':re.findall(r'Convergence Check Mode\s+ConvCheckMode\s*\.{2,}\s*([^\n]+)',native),
        'stop_signals':re.findall(r'\*+\s*(.*?)\s+signals convergence\s*\*+',native),
        'numeric_printed_criteria':[{'name':n.strip(),'value':m.cp.number(v),'tolerance':m.cp.number(t),'passes_printed_tolerance':abs(m.cp.number(v))<=m.cp.number(t)} for n,v,t in re.findall(r'Last\s+(.+?)\s+\.{3}\s+('+m.F+r')\s+Tolerance\s*:\s*('+m.F+r')',native)]}

def admit(path,panel_path,panel,adapter):
    r=json.loads(Path(path).read_text());b.check(r['status']=='complete' and all(r['checks'].values()),'Failed native admission '+str(path))
    same(r['panel'],b.pin(panel_path),'Result panel changed');same(r['adapter'],panel['adapter'],'Native source changed')
    c=next(x for x in panel['cells'] if x['cell_id']==r['cell_id']);same(r['cell'],c,'Result physical cell changed')
    for key in ('output','stderr','input'):b.verify(r[key])
    startpath=Path(path).parent/'START.json';start=json.loads(startpath.read_text());same(start['cell'],c,'Launch cell changed');same(start['panel'],r['panel'],'Launch panel changed')
    b.check(start['ranks']==4 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4),'Actual resource record changed')
    same(b.verify(r['input']).read_text(),adapter.deck(panel,c,4,start['maxcore_MiB']),'Executed input deck differs')
    text=b.verify(r['output']).read_text();replay=adapter.parse(text,c,r['returncode'])
    for key,val in replay.items():same(r[key],val,'Saved/native replay differs: '+key)
    return r,{'result':b.pin(path),'start':b.pin(startpath),'cell_id':r['cell_id'],'checks':r['checks'],'raw_printed_scf_criteria':r['printed_scf_criteria'],'convergence':convergence(text),'allocated_cpu_seconds':r['allocated_cpu_seconds']}

def prepare(path):
    cfg={'protocol':'matched_water_cp_final_v1','created_utc':b.now(),'source':b.pin(__file__),'specification':b.pin(HERE/'COLLECTION_SPEC.md'),
        'dft_panel':b.pin(m.WS/'dft_water_cp_v1/prepared_v3/PANEL.json'),
        'canonical_collection':b.pin(m.WS/'canonical_water_cp_v2/collection_v1/RESULT.json'),
        'scf_interpretation':b.pin(HERE.parent/'scf_mode_interpretation_v1/REPORT.md'),
        'scf_interpretation_data':b.pin(HERE.parent/'scf_mode_interpretation_v1/RESULT.json'),
        'decision_resolution_kcal_mol':0.25,'primary_method':'PBE0_D4_CP',
        'task_root':str(m.WS/'dft_water_cp_v1/run_v1/tasks'),'expected_new_cells':['La__ion_ghost_water','La__water_ghost_ion','Dy__ion_ghost_water','Dy__water_ghost_ion','real_water_D4']}
    m.load(cfg['dft_panel']['path']);b.write(path,cfg);print(json.dumps({'config':b.pin(path),'molecular_calls':0}))

def run(config,out):
    tick=time.monotonic();cfg=json.loads(Path(config).read_text());b.check(cfg['protocol']=='matched_water_cp_final_v1','Protocol changed')
    for key in ('source','specification','dft_panel','canonical_collection','scf_interpretation','scf_interpretation_data'):b.verify(cfg[key])
    same(cfg['source'],b.pin(__file__),'Collector source changed');b.check(cfg['decision_resolution_kcal_mol']==.25 and cfg['primary_method']=='PBE0_D4_CP','Decision changed')
    panelpath=cfg['dft_panel']['path'];p=m.load(panelpath);cp_path=p['canonical_panel']['path'];cp=m.cp.load(cp_path)
    b.check(set(cfg['expected_new_cells'])=={x['cell_id'] for x in p['cells']},'Inventory changed')
    newpaths=[str(Path(cfg['task_root'])/cell/'execution/RESULT.json') for cell in cfg['expected_new_cells']]
    native=[];terms={};seen=set()
    for path in newpaths:
        r,e=admit(path,panelpath,p,m);b.check(r['cell_id'] not in seen,'Duplicate new cell');seen.add(r['cell_id']);native.append(e)
        terms[r['cell_id']]={'electronic_Eh':r['electronic_energy_hartree'],'D4_Eh':r['D4_energy_hartree'],'final_Eh':r['energy_hartree']}
    cc=json.loads(b.verify(cfg['canonical_collection']).read_text());ccpaths=[str(b.verify(x)) for x in cc['result_pins']];seen=set()
    for path in ccpaths:
        r,e=admit(path,cp_path,cp,m.cp);b.check(r['cell_id'] not in seen,'Duplicate canonical cell');seen.add(r['cell_id']);native.append(e)
    b.check(seen=={x['cell_id'] for x in cp['cells']},'Canonical denominator changed')
    same(m.cp.collect(cp_path,ccpaths),cc,'Canonical collection replay differs')
    old=b.load_panel(p['parent_panel']['path']);original=[]
    for metal,saved in p['saved_complexes'].items():
        path=b.verify(saved['result']);r=json.loads(path.read_text());same(r['panel'],p['parent_panel'],'Original panel changed');same(r['adapter'],old['adapter'],'Original source changed')
        c=next(x for x in old['cells'] if x['cell_id']==r['cell_id']);replay=b.collect_one(old,c,path.parent,r['returncode'])
        for key,value in replay.items():same(r[key],value,'Original complex native replay differs '+key)
        g=next(x for x in cp['cells'] if x['cell_id']==metal+'__complex');g=dict(g,d4=True);extra=m.parse(b.verify(r['output']).read_text(),g,r['returncode']);b.check(extra['status']=='complete','Extra original complex admission failed')
        same(extra['electronic_energy_hartree'],saved['electronic_energy_hartree'],'Original electronic term changed');same(extra['D4_energy_hartree'],saved['D4_energy_hartree'],'Original D4 term changed')
        original.append({'result':b.pin(path),'cell_id':r['cell_id'],'checks':r['checks'],'additional_checks':extra['checks'],'raw_printed_scf_criteria':r['printed_scf_criteria'],'convergence':convergence(b.verify(r['output']).read_text()),'allocated_cpu_seconds':r['allocated_cpu_seconds']})
        terms[metal+'__saved_complex']={key:saved[key] for key in ('electronic_energy_hartree','D4_energy_hartree','final_energy_hartree')}
    result=m.collect(panelpath,newpaths,cfg['canonical_collection']['path']);delta=result['Dy_minus_La'];ref=delta['CCSD_T_CP_kcal_mol']
    differences={key:delta[key]-ref for key in ('PBE0_CP_kcal_mol','PBE0_D4_CP_kcal_mol')};primary=differences['PBE0_D4_CP_kcal_mol'];within=abs(primary)<=cfg['decision_resolution_kcal_mol']
    result.update(config=b.pin(config),source=b.pin(__file__),specification=cfg['specification'],electronic_dispersion_terms=terms,
        paired_method_minus_CC_kcal_mol=differences,
        individual_method_minus_CC_kcal_mol={metal:{key:values[key]-values['CCSD_T_CP_kcal_mol'] for key in ('PBE0_CP_kcal_mol','PBE0_D4_CP_kcal_mol')} for metal,values in result['per_metal'].items()},
        decision={'primary_method':'PBE0_D4_CP','resolution_kcal_mol':.25,'within_diagnostic_resolution':within,'conclusion':('agrees_with_canonical_CC_at_this_resolution' if within else 'not_interchangeable_with_canonical_CC_at_this_resolution'),'physical_accuracy_qualified':False},
        native_replay_reviews=native,original_complex_reviews=original,scf_interpretation=cfg['scf_interpretation'],scf_interpretation_data=cfg['scf_interpretation_data'],
        local_threads=1,elapsed_seconds=time.monotonic()-tick,command=sys.argv,utc=b.now(),
        cost_accounting={'new_DFT_allocated_CPU_seconds':sum(x['allocated_cpu_seconds'] for x in native[:5]),'reused_CC_allocated_CPU_seconds':sum(x['allocated_cpu_seconds'] for x in native[5:]),'reused_original_DFT_allocated_CPU_seconds':sum(x['allocated_cpu_seconds'] for x in original),'collection_new_molecular_calls':0})
    out.mkdir(parents=True,exist_ok=False);b.write(out/'RESULT.json',result)
    lines=['# Matched water counterpoise comparison','',f"**PBE0-D4 {'agrees with' if within else 'is not interchangeable with'} canonical CCSD(T) at the declared0.25kcal/mol resolution on this paired fixed-water fixture.**",'',
        '| Method | La interaction | Dy interaction | Dy minus La | Paired difference from CC |','|---|---:|---:|---:|---:|']
    for key,label in [('CCSD_T_CP_kcal_mol','CCSD(T), CP'),('PBE0_CP_kcal_mol','PBE0 electronic, CP'),('PBE0_D4_CP_kcal_mol','PBE0-D4, CP')]:
        lines.append(f"| {label} | {result['per_metal']['La'][key]:.6f} | {result['per_metal']['Dy'][key]:.6f} | {delta[key]:.6f} | {delta[key]-ref:+.6f} |")
    lines+=['','All energies are kcal/mol. Every value uses the same fixed coordinates, orbital bases and frozen-f effective-core representation. Individual-metal errors remain in RESULT.json; cancellation in the paired difference is not individual accuracy.',
        '',f"Real-water D4 is {terms['real_water_D4']['D4_Eh']:.12f}Eh. It is subtracted in each interaction; isolated-ion D4 is the explicit zero empty-pair/triple sum. No D4 was evaluated on ghosts. Electronic-only and dispersion-added results remain separate.",
        '', 'All five new fragments, two reused complexes and six canonical native results were replayed against their frozen sources and inputs. Raw SCF residuals and actual stopping signals are retained in RESULT.json. The printed thresholds are not retroactively treated as an all-six stopping gate; the pinned SCF interpretation explains the configured mode and stale DIIS fields in prior outputs.',
        '', 'This finite decision compares two approximate electronic methods at one geometry/basis. It does not establish basis convergence, frozen-f or spin-orbit accuracy, force accuracy, donor selectivity, or protein/water exchange free energy. The current numerical stopping rules are retained; observable precision can be tested separately without relabeling these completed calculations as failed.',
        '', 'The practical next chemistry question remains the paired water-to-amide contrast and its change with hydration, using the already declared correlated qualification design. A water-only agreement would not validate donor transfer; a discrepancy here would prevent adopting PBE0-D4 as an interchangeable correlated comparator without a bounded follow-up.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'status':'complete','result':b.pin(out/'RESULT.json'),'decision':result['decision'],'paired_differences':differences}))

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True);a=sub.add_parser('prepare');a.add_argument('--out',type=Path,required=True);a=sub.add_parser('run');a.add_argument('--config',required=True);a.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.action=='prepare':prepare(a.out)
    else:run(a.config,a.out)

if __name__=='__main__':main()
