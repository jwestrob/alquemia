#!/usr/bin/env python3
"""Copy the sixteen admitted donor radial/rotation endpoints; no new scoring."""
import argparse
import time
from pathlib import Path
import BUILD_DATA as b

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--baseline',required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args();t=time.monotonic();cpu=time.process_time();out=Path(a.out).resolve()
    b.check(not out.exists(),'fresh displaced materialization');out.mkdir(parents=True)
    sources=b.read(b.pin(b.HERE/'INPUTS.json'));ds=b.read(sources['design_sources']);auth=ds['donor_authorities']
    model=b.read(auth['model_derivative_result']);reference=b.read(auth['reference_derivative_result'])
    manifest=b.read(auth['model_manifest']);base=b.read(b.pin(a.baseline));resolver=b.load(sources['resolver'])
    key=b.verify(base['key']).read_text();prm=b.verify(base['parameter_file']).read_text();shared=resolver.parse_parameters(prm,key)
    b.check(model['complete_cells']==16 and reference['admitted_cells']==16,'exact admitted sixteen')
    b.check(reference['scoring']['all_numerical_checks_passed'],'original two-step qualifications retained')
    rows={}
    for cid,r in model['rows'].items():
        slot=reference['slots'][cid];q=b.read(slot['result']);m=manifest['rows'][cid]
        b.check(r['status']==slot['status']==q['status']=='complete','complete endpoint '+cid)
        g,p,k=b.typed_geometry(r['input'],'acetate',resolver)
        c=b.compatibility(p,k,g,shared,key,r['input']['parameters']['sha256']==base['parameter_file']['sha256'])
        b.check(c['compatible'],'same shared baseline '+cid)
        rows[cid]=dict(geometry=g,inputs=r['input'],baseline_compatibility=c,
            coordinate=dict(direction=m['direction'],q_A=m['q_A'],shift_label=m['shift_label'],source_geometry=m['source_geometry']),
            model=dict(energy_kcal_mol=r['parsed']['native']['total_energy_kcal_mol'],
                       components=r['parsed']['native']['components'],literal_tokens=r['parsed']['tokens'],
                       print_bounds=r['parsed']['half_print_units_kcal_mol'],source=auth['model_derivative_result'],
                       forces=None,force_unavailable_reason='These are saved analyze energy calls; no endpoint gradient was requested.'),
            reference=dict(energy_hartree=q['energy_hartree'],electronic_energy_hartree=q['electronic_energy_hartree'],
                           D4_energy_hartree=q['D4_energy_hartree'],method='Original ECP PBE0-D4 protocol',
                           printed_scf_criteria=q['printed_scf_criteria'],native_checks=q['checks'],
                           source=slot['result'],forces=None,force_unavailable_reason='Saved single-point references; center derivatives qualified separately.'),
            baseline_operator=None,consumed_development=True)
    b.save(out/'RESULT.json',dict(status='materialized_sixteen_saved_donor_displacements',implementation=b.pin(__file__),
        baseline=b.pin(a.baseline),rows=rows,original_model_derivative_comparison=model['comparison'],
        original_reference_derivative_scoring=reference['scoring'],sources=auth,
        new_native_calls=0,new_quantum_calls=0,operator_solves=0,fits=0,
        cost=dict(wall_seconds=time.monotonic()-t,CPU_seconds=time.process_time()-cpu,threads=1)))
    print(b.pin(out/'RESULT.json'))

if __name__=='__main__':main()
