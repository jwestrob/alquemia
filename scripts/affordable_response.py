"""Analytic-gradient adapter with physical cap chain rules; no enabled response score."""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import numpy as np

from affordable_common import (InvalidArtifact, BOHR_TO_A, HA_TO_KCAL, verify, read_json,
                               record, write_new, energy, xyz)
from site_mechanics import read_engrad


def cap_jacobians(retained_xyz, omitted_xyz, length_A):
    """d[x + l(y-x)/|y-x|]/dx, /dy; coordinates and l are in Angstrom."""
    d = np.asarray(omitted_xyz, dtype=float) - np.asarray(retained_xyz, dtype=float)
    r = np.linalg.norm(d)
    if not np.isfinite(r) or r <= length_A:
        raise InvalidArtifact('invalid physical link-atom geometry')
    e = d/r
    b = length_A/r * (np.eye(3)-np.outer(e,e))
    return np.eye(3)-b, b


def source_key(source):
    return '/'.join(str(source[k]) for k in ('chain_index','residue_index','atom'))


def map_to_source(gradient_kcal_A, repair):
    """Map link gradients to BOTH true endpoint atoms, not to free cap motions."""
    g = np.asarray(gradient_kcal_A, dtype=float)
    mapping = repair['atom_graph']['source_to_qm']
    if g.shape != (len(mapping)+1, 3) or not np.isfinite(g).all():
        raise InvalidArtifact('gradient atom count/nonfinite values')
    from affordable_peptide import SourceGraph
    graph = SourceGraph(verify(repair['source_structure']),verify(repair['topology_definition']))
    result={'metal':g[0].copy()}
    for row in mapping:
        local=g[row['qm_index']]
        if row['kind']=='source':
            k=source_key(row['source']); result[k]=result.get(k,np.zeros(3))+local
        elif row['kind']=='sigma_link_H':
            coords=[]
            for label in ('retained','omitted'):
                s=row[label]; a=graph.atoms[(s['chain_index'],s['residue_index'],s['atom'])]
                coords.append(np.array(tuple(a.pos)))
            ja,jb=cap_jacobians(*coords,row['length_A'])
            for label,j in [('retained',ja),('omitted',jb)]:
                k=source_key(row[label]); result[k]=result.get(k,np.zeros(3))+j.T@local
        else:
            raise InvalidArtifact('unknown coordinate mapping')
    return {k:v.tolist() for k,v in sorted(result.items())}


def extract(engrad_path, output, input_path, xyz_path, repair=None):
    text=Path(input_path).read_text()
    if re.search(r'\b(?:NumGrad|NumFreq|Opt|Freq)\b',text,re.I) or not re.search(r'\bEnGrad\b',text,re.I):
        raise InvalidArtifact('requires explicit analytic EnGrad single endpoint; no numerical derivatives')
    output_text=Path(output).read_text()
    if re.search(r'numerical gradient|numerical differentiation',output_text,re.I):
        raise InvalidArtifact('output indicates numerical differentiation')
    e=energy(output)
    grad=read_engrad(Path(engrad_path))
    atoms=xyz(xyz_path)
    if len(atoms)!=grad['atom_count'] or abs(e-grad['energy_Ha'])>1e-8:
        raise InvalidArtifact('gradient energy/atom mismatch')
    import gemmi
    if not np.array_equal(grad['atomic_numbers'],[gemmi.Element(a[0]).atomic_number for a in atoms]):
        raise InvalidArtifact('gradient element ordering mismatch')
    if not np.allclose(grad['coordinates_bohr']*BOHR_TO_A,[a[1:] for a in atoms],atol=1e-6,rtol=0):
        raise InvalidArtifact('gradient coordinates differ from frozen XYZ')
    g=grad['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A
    result={'schema_version':'alquemia.analytic_gradient.v1','energy_hartree':e,
            'gradient_kcal_mol_per_A':g.tolist(),'quantity':'gradient_not_force',
            'artifacts':{k:record(p) for k,p in [('engrad',engrad_path),('output',output),('input',input_path),('xyz',xyz_path)]},
            'energy_scope':'isolated_CPCM_endpoint', 'environment_gradient_included':False,
            'combined_environment_response_status':'unsupported_missing_environment_derivative',
            'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,
            'entropy_correction_kcal_mol':None}
    if repair is not None:
        result['mapped_source_gradient_kcal_mol_per_A']=map_to_source(g,read_json(repair))
        result['repair_mapping']=record(repair)
    return result


def bounded_response(gradient, curvature, *, validation=None):
    # There is deliberately no numeric bypass. No curvature source is validated
    # in this development version, regardless of a supplied matrix's eigenvalues.
    return {'status':'response_model_not_validated','relaxation_energy_kcal_mol':None,
            'reason':'No independently validated inexpensive curvature backend is registered.'}


def paired_sensitivity(la_path,ca_path):
    """Export grad(E_Ca-E_La) in the same physical coordinate measure."""
    la,ca=read_json(la_path),read_json(ca_path)
    for r in (la,ca):
        if r['energy_scope']!='isolated_CPCM_endpoint' or r['quantity']!='gradient_not_force':
            raise InvalidArtifact('unsupported gradient energy scope or quantity')
        for artifact in r['artifacts'].values(): verify(artifact)
    a,b=[xyz(verify(r['artifacts']['xyz'])) for r in (la,ca)]
    if len(a)!=len(b) or a[0][0]!='La' or b[0][0]!='Ca' or a[0][1:]!=b[0][1:] or a[1:]!=b[1:]:
        raise InvalidArtifact('gradient endpoint coordinates/order differ')
    inputs=[verify(r['artifacts']['input']).read_text() for r in (la,ca)]
    stripped=[re.sub(r'^\s*\*\s+xyzfile.*$', '',t,flags=re.M|re.I) for t in inputs]
    if stripped[0]!=stripped[1]: raise InvalidArtifact('gradient endpoint input methods differ')
    mapped='mapped_source_gradient_kcal_mol_per_A'
    if mapped not in la or mapped not in ca or la.get('repair_mapping')!=ca.get('repair_mapping'):
        raise InvalidArtifact('matching physical source mappings are required')
    if set(la[mapped])!=set(ca[mapped]): raise InvalidArtifact('physical coordinate measures differ')
    delta={k:(np.asarray(ca[mapped][k])-la[mapped][k]).tolist() for k in la[mapped]}
    if not all(np.isfinite(v).all() for v in map(np.asarray,delta.values())):
        raise InvalidArtifact('nonfinite paired gradient')
    return {'quantity':'grad_E_Ca_minus_E_La','units':'kcal_mol_per_A','physical_source_gradient':delta,
            'sources':{'La':record(la_path),'Ca':record(ca_path)},'energy_scope':'isolated_CPCM_endpoint',
            'environment_gradient_included':False,'uncertainty_covariance':None,
            'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('engrad','orca-output','input','xyz','la-gradient','ca-gradient'):
        p.add_argument('--'+name,type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repair-manifest',type=Path)
    a=p.parse_args()
    if a.la_gradient or a.ca_gradient:
        if not (a.la_gradient and a.ca_gradient): p.error('both gradient records are required')
        result=paired_sensitivity(a.la_gradient,a.ca_gradient)
    else:
        if not all((a.engrad,a.orca_output,a.input,a.xyz)): p.error('extraction requires --engrad, --orca-output, --input, --xyz')
        result=extract(a.engrad,a.orca_output,a.input,a.xyz,a.repair_manifest)
    write_new(a.output,result)
    print('analytic gradient extracted; response_model_not_validated')


if __name__=='__main__': main()
