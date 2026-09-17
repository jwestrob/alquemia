"""Compare actual intact-chain native and batched OMOL energies, without scoring."""
import argparse
import json
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_hybrid import EV_TO_KCAL
from mace_omol_intact import qualified as native_qualified
from mace_omol_edge_run import validate,collect
from mace_omol import TOL
from mace_file_checks import cached_file_checks


@cached_file_checks
def compare(native_collection,edge_collection):
    ns=read_json(native_collection);np=verify(ns['manifest']);nm=read_json(np)
    if nm['stage']=='cpu_intact':
        validate(np);native=collect(np)
        if ns!=native or not native['numerical_gate_pass']:raise InvalidArtifact('passing actual native CPU intact reference required')
        native_inputs=nm['intact_manifest']
    else:
        native,nm=native_qualified(native_collection,'intact_qualification');native_inputs=native['manifest']
    saved=read_json(edge_collection);mp=verify(saved['manifest']);m=read_json(mp)
    validate(mp);edge=collect(mp)
    if saved!=edge or m['stage']!='edge_intact' or not edge['numerical_gate_pass']:
        raise InvalidArtifact('actual passing batched intact qualification required')
    if m['intact_manifest']!=native_inputs or set(edge['rows'])!=set(native['rows']):
        raise InvalidArtifact('native/adapter inputs do not match')
    errors={};checks=[]
    for key,row in edge['rows'].items():
        error=(row['energy_eV']-native['rows'][key]['energy_eV'])*EV_TO_KCAL;errors[key]=error
        checks.append({'name':key,'error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
    for variant in ('primary','repeat','rotate','farther'):
        pair={}
        for position in (('detached',) if variant=='farther' else ('bound','detached')):
            e=errors[f'ALPHA_1F6S_Ca_{position}_{variant}']-errors[f'ALPHA_1F6S_La_{position}_{variant}'];pair[position]=e
            checks.append({'name':'Ca_minus_La_'+position+'_'+variant,'error_kcal_mol':e,'pass':abs(e)<=TOL['energy_kcal_mol']})
        if variant!='farther':
            e=pair['bound']-pair['detached'];checks.append({'name':'R_coord_'+variant,'error_kcal_mol':e,'pass':abs(e)<=TOL['energy_kcal_mol']})
    result={'status':'complete','native_collection':record(native_collection),'edge_collection':record(edge_collection),
            'tolerances':TOL,'checks':checks,'native_equivalence_pass':all(c['pass'] for c in checks),
            'native_peak_cuda_allocated_bytes':max((r['peak_cuda_allocated_bytes'] for r in native['rows'].values() if r['peak_cuda_allocated_bytes'] is not None),default=None),
            'native_peak_host_RSS_KiB':max(r['peak_host_RSS_KiB'] for r in native['rows'].values()),
            'edge_peak_cuda_allocated_bytes':max(r['peak_cuda_allocated_bytes'] for r in edge['rows'].values()),
            'native_evaluation_seconds':sum(r['evaluation_seconds'] for r in native['rows'].values()),
            'edge_evaluation_seconds':sum(r['evaluation_seconds'] for r in edge['rows'].values()),
            'native_device':next(iter(native['rows'].values()))['device'],
            'edge_device':next(iter(edge['rows'].values()))['device'],
            'runtime_comparison_caveat':'Different devices; not a matched-hardware speedup claim.',
            'force_validation_status':'not_requested_energy_only','predictive_validation':None,'baseline_changed':False}
    return result


@cached_file_checks
def verified(path):
    saved=read_json(path);verify(saved['implementation'])
    result=compare(verify(saved['native_collection']),verify(saved['edge_collection']))
    if {k:v for k,v in saved.items() if k!='implementation'}!=result or not result['native_equivalence_pass']:
        raise InvalidArtifact('actual native/intact edge equivalence required')
    return result


def report(native_collection,edge_collection,output):
    result=compare(native_collection,edge_collection);checks=result['checks']
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/Path(__file__).name)
    result['implementation']=record(out/Path(__file__).name);write_new(out/'result.json',result)
    (out/'REPORT.md').write_text('# Exact OMOL interaction batching\n\n'
        f'Native intact equivalence: {result["native_equivalence_pass"]};{len(checks)}actual energy/paired checks.\n\n'
        f'Maximum error: {max(abs(c["error_kcal_mol"]) for c in checks):.12g}kcal/mol.\n\n'
        f'Native reference device: {result["native_device"]["name"]}; peak host RSS: {result["native_peak_host_RSS_KiB"]} KiB.\n\n'
        f'Batched GPU peak allocation: {result["edge_peak_cuda_allocated_bytes"]} bytes. Native GPU allocation: {result["native_peak_cuda_allocated_bytes"] if result["native_peak_cuda_allocated_bytes"] is not None else "not applicable (CPU reference)"}.\n\n'
        'The physical graph, weights, charge/spin, precision and energy expression are unchanged. No force, prediction or production promotion is claimed. Devices differ; timings do not establish a matched-hardware speedup.\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('native-collection','edge-collection','output'):p.add_argument('--'+key,required=True)
    r=report(**vars(p.parse_args()));print(json.dumps({k:r[k] for k in ('status','native_equivalence_pass','native_peak_cuda_allocated_bytes','edge_peak_cuda_allocated_bytes')},indent=2))
