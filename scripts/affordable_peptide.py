"""Opt-in peptide-amide repair from a pinned, already selected generic core.

This never changes the default carver. Standard-residue topology supplies
intraresidue connectivity, coordinates validate bonds, and geometrically bonded
C/N neighbours in the deposited polymer order supply peptide connectivity.
Residue numbers are identifiers, never arithmetic connectivity rules.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

import gemmi
import numpy as np

import carve_generic as old
from coordination_policy import selected_residue_atoms, atom_altloc
from affordable_common import InvalidArtifact, record, verify, read_json, write_new, paired

PROTOCOL = "generic_peptide_amide_vertical_native_r2scan3c_v3"
POLICY = "source_graph_amide_preserving_sigma_caps_v1"


class SourceGraph:
    def __init__(self, source, topology):
        self.source = Path(source)
        self.topology = Path(topology)
        self.structure = gemmi.read_structure(str(source))
        indexed, _ = old.index_residues(self.structure[0])
        self.residues = {r.key: r for r in indexed}
        self.atoms, self.meta, self.by_residue = {}, {}, {}
        self.edges, self.required, self.hparents = set(), defaultdict(set), {}
        self.amide_next, self.amide_prev = {}, {}
        self.supported_topology_residues = set()
        templates = {r.attrib['name']: r for r in ET.parse(topology).findall('.//Residues/Residue')}
        for r in indexed:
            selected = selected_residue_atoms(r.residue)
            names = {}
            for a in selected:
                name = a.name.strip()
                if name in names:
                    raise InvalidArtifact(f"ambiguous selected atom: {r.source_dict()} {name}")
                key = (*r.key, name)
                self.atoms[key] = a
                self.meta[key] = {**r.source_dict(), "chain_index": r.key[0],
                                  "residue_index": r.key[1], "atom": name,
                                  "altloc": atom_altloc(a), "element": a.element.name}
                names[name] = key
            self.by_residue[r.key] = names
            template = templates.get(r.resname) or templates.get(r.canonical_resname)
            # HIS aliases affect hydrogen names, not the heavy-atom graph.
            if template is None and r.canonical_resname == 'HIS':
                template = templates.get('HID')
            if template is not None:
                self.supported_topology_residues.add(r.key)
                for b in template.findall('Bond'):
                    n1, n2 = b.attrib['atomName1'], b.attrib['atomName2']
                    if n1.startswith('H') or n2.startswith('H'):
                        continue
                    self.require_edge(r.key, n1, n2)
            if 'OXT' in names:
                self.require_edge(r.key, 'C', 'OXT')
            heavy = [a for a in selected if not old.is_hydrogen(a)]
            for h in selected:
                if not old.is_hydrogen(h) or not heavy:
                    continue
                parent = min(heavy, key=lambda a: a.pos.dist(h.pos))
                if parent.pos.dist(h.pos) <= parent.element.covalent_r + h.element.covalent_r + .30:
                    self.hparents[names[h.name.strip()]] = names[parent.name.strip()]
        # Polymer order + covalent distance + label_seq metadata, never auth_seq+1.
        for chain_i, chain in enumerate(self.structure[0]):
            for residue_i in range(len(chain) - 1):
                akey, bkey = (chain_i, residue_i), (chain_i, residue_i + 1)
                an, bn = self.by_residue.get(akey, {}), self.by_residue.get(bkey, {})
                if 'C' not in an or 'N' not in bn:
                    continue
                a, b = self.atoms[an['C']], self.atoms[bn['N']]
                if not 1.15 <= a.pos.dist(b.pos) <= 1.65:
                    continue
                ar, br = chain[residue_i], chain[residue_i + 1]
                if ar.label_seq is not None and br.label_seq is not None and br.label_seq != ar.label_seq + 1:
                    continue
                if atom_altloc(a) and atom_altloc(b) and atom_altloc(a) != atom_altloc(b):
                    continue
                if 'OXT' in an:
                    raise InvalidArtifact("terminal OXT conflicts with covalent peptide neighbour")
                self.add_edge(an['C'], bn['N'])
                self.amide_next[akey], self.amide_prev[bkey] = bn['N'], an['C']

    def add_edge(self, a, b):
        self.edges.add(tuple(sorted((a, b))))
        self.required[a].add(b)
        self.required[b].add(a)

    def require_edge(self, residue, n1, n2):
        a, b = (*residue, n1), (*residue, n2)
        self.required[a].add(b)
        self.required[b].add(a)
        if a in self.atoms and b in self.atoms:
            self.edges.add(tuple(sorted((a, b))))

    def locate(self, selector):
        matches = [r for r in self.residues.values()
                   if r.chain == selector['chain'] and r.resnum == selector['resnum']
                   and r.insertion_code == selector.get('insertion_code', '')
                   and r.resname == selector['resname']]
        if len(matches) != 1:
            raise InvalidArtifact(f"ambiguous/missing source residue: {selector}")
        return matches[0]

    def key(self, residue, name):
        k = (*residue, name)
        if k not in self.atoms:
            raise InvalidArtifact(f"missing required source atom: {k}")
        return k

    def peptide(self, residue):
        selected = {self.key(residue, n) for n in ('C', 'O')}
        if 'OXT' in self.by_residue[residue]:
            selected.add(self.key(residue, 'OXT'))
        else:
            n = self.amide_next.get(residue)
            if n is None:
                raise InvalidArtifact(f"absent/ambiguous covalently bonded amide N for {residue}")
            selected.add(n)
            if self.residues[n[:2]].canonical_resname == 'PRO':
                selected.update(self.key(n[:2], x) for x in ('CA', 'CB', 'CG', 'CD'))
        return selected

    def close_overlaps(self, selected):
        """Join selected branches through CA; never cut a peptide C-N bond."""
        selected = set(selected)
        while True:
            previous = set(selected)
            for rkey, names in self.by_residue.items():
                branches = [names[n] for n in ('N', 'C', 'CB') if n in names and names[n] in selected]
                if len(branches) > 1:
                    selected.add(self.key(rkey, 'CA'))
                if names.get('N') in selected:
                    carbon = self.amide_prev.get(rkey)
                    if carbon is None:
                        raise InvalidArtifact(f"selected amide N has no peptide carbonyl: {rkey}")
                    selected.update(self.peptide(carbon[:2]))
            if selected == previous:
                return selected

    def materialize(self, selected):
        selected = self.close_overlaps(selected)
        for k in selected:
            if k[:2] not in self.supported_topology_residues:
                raise InvalidArtifact(f"unsupported residue topology for selected atom: {k}")
            for other in self.required[k]:
                if other not in self.atoms:
                    raise InvalidArtifact(f"missing covalent neighbour {other} of selected {k}")
                a, b = self.atoms[k], self.atoms[other]
                distance = a.pos.dist(b.pos)
                if distance < .7 or distance > a.element.covalent_r + b.element.covalent_r + .35:
                    raise InvalidArtifact(f"invalid source bond {k}--{other}: {distance:.3f} A")
        hs = {h for h, parent in self.hparents.items() if parent in selected}
        # Carbonyl and amide valence: missing protonation is an explicit failure.
        for k in selected:
            if k[2] == 'N':
                expected = 0 if self.residues[k[:2]].canonical_resname == 'PRO' else 1
                if sum(p == k for p in self.hparents.values()) != expected:
                    raise InvalidArtifact(f"unsupported amide-N protonation at {k}")
            if k[2] == 'CA':
                expected = 2 if self.residues[k[:2]].canonical_resname == 'GLY' else 1
                if sum(p == k for p in self.hparents.values()) != expected:
                    raise InvalidArtifact(f"missing/extra alpha hydrogen at {k}")
            if self.residues[k[:2]].canonical_resname == 'PRO' and k[2] in ('CB','CG','CD'):
                if sum(p == k for p in self.hparents.values()) != 2:
                    raise InvalidArtifact(f"missing/extra proline-ring hydrogen at {k}")
        order = sorted(selected | hs)
        atoms = [old.atom_to_qm(self.atoms[k]) for k in order]
        mapping = [{"qm_index": i + 1, "kind": "source", "source": self.meta[k],
                    "xyz_A": list(atoms[i][1:])} for i, k in enumerate(order)]
        cuts, retained = [], []
        allowed = {frozenset(p) for p in [('C','CA'), ('N','CA'), ('CA','CB')]}
        for a, b in sorted(self.edges):
            if a in selected and b in selected:
                retained.append({"a": self.meta[a], "b": self.meta[b]})
            elif (a in selected) != (b in selected):
                keep, omit = (a, b) if a in selected else (b, a)
                if keep[:2] != omit[:2] or frozenset((keep[2], omit[2])) not in allowed:
                    raise InvalidArtifact(f"unsupported boundary bond: {keep}--{omit}")
                x = np.array(old.atom_to_qm(self.atoms[keep])[1:])
                y = np.array(old.atom_to_qm(self.atoms[omit])[1:])
                length = 1.01 if self.atoms[keep].element.name == 'N' else 1.09
                pos = x + length * (y - x) / np.linalg.norm(y - x)
                atoms.append(('H', *pos.tolist()))
                cap = {"qm_index": len(atoms), "kind": "sigma_link_H", "retained": self.meta[keep],
                       "omitted": self.meta[omit], "length_A": length, "xyz_A": pos.tolist(),
                       "source_bond_length_A": float(np.linalg.norm(y - x))}
                mapping.append(cap)
                cuts.append(cap)
        # Duplicate source/cap positions are chemical failures, never deduplicated by distance.
        positions = [tuple(a[1:]) for a in atoms]
        if len(set(positions)) != len(positions):
            raise InvalidArtifact("overlapping atoms/caps")
        return atoms, {"source_to_qm": mapping, "retained_bonds": retained,
                       "cut_bonds_and_caps": cuts, "selected_source_heavy_atoms": len(selected)}


def repair(manifest_path, output, topology):
    manifest_path, output = Path(manifest_path), Path(output)
    baseline = read_json(manifest_path)
    if baseline['protocol_id'] != old.PROTOCOL_ID or baseline['status'] != 'ready_for_orca':
        raise InvalidArtifact('only prepared generic-v2 cores are accepted; PQQ stays unchanged')
    source = verify(baseline['source_structure'])
    for metal in ('La','Ca'):
        for kind in ('xyz','orca_input'):
            verify(baseline['outputs'][metal][kind])
    graph = SourceGraph(source, topology)
    selected, water_fragments, charge_entries = set(), [], []
    for fragment in baseline['qm_fragments']:
        r = graph.locate(fragment['source'])
        kind = fragment['fragment_type']
        if kind == 'sidechain':
            # Reuse the existing microstate/explicit-H checker, not inferred charges.
            _, checked = old.sidechain_fragment(r, [], 0)
            if checked['formal_charge'] != fragment['formal_charge']:
                raise InvalidArtifact('side-chain microstate changed')
            selected.update(graph.key(r.key, n) for n in old.SIDECHAIN_QM_HEAVY_ATOMS[r.canonical_resname])
        elif kind in ('backbone_carbonyl', 'c_terminal_carboxylate'):
            selected.update(graph.peptide(r.key))
        elif kind in ('water', 'explicit_water'):
            atoms, _ = old.water_fragment(r, [], 0)
            water_fragments.append((r, atoms))
        else:
            raise InvalidArtifact(f'unsupported fragment: {kind}')
        charge_entries.append({"source": fragment['source'], "kind": kind,
                               "formal_charge": fragment['formal_charge']})
    atoms, graph_record = graph.materialize(selected)
    for r, water in water_fragments:
        names = [a.name.strip() for a in selected_residue_atoms(r.residue)]
        # Map by exact position, since the legacy water helper uses its own order.
        for atom in water:
            matches = [k for k in graph.by_residue[r.key].values() if old.atom_to_qm(graph.atoms[k]) == atom]
            if len(matches) != 1:
                raise InvalidArtifact('unmapped explicit-water atom')
            atoms.append(atom)
            graph_record['source_to_qm'].append({"qm_index": len(atoms), "kind": "source",
                                                 "source": graph.meta[matches[0]], "xyz_A": list(atom[1:])})
    ligand_charge = sum(x['formal_charge'] for x in charge_entries)
    if ligand_charge != baseline['charge_ledger']['fragment_formal_charge_sum']:
        raise InvalidArtifact('charge ledger changed')
    if output.exists():
        raise InvalidArtifact(f'refusing to overwrite workspace {output}')
    output.mkdir(parents=True)
    stem = baseline['stem'] + '_amide_v3'
    result = {"schema_version": "alquemia.peptide_repair.v1", "protocol_id": PROTOCOL,
              "status": "prepared", "policy_id": POLICY, "baseline": record(manifest_path),
              "source_structure": record(source), "topology_definition": record(topology),
              "implementation": {"repair": record(__file__), "generic": record(old.__file__)},
              "selected_site": baseline['selected_site'], "coordination": baseline['coordination'],
              "protonation_manifest": baseline.get('protonation_manifest'),
              "explicit_water_inventory": [r.source_dict() for r, _ in water_fragments],
              "charge_ledger": charge_entries, "atom_graph": graph_record,
              "electronic_structure": baseline['electronic_structure'],
              "decision": "uncalibrated_protocol", "absolute_reference_status": "not_registered",
              "outputs": {}, "tasks": []}
    metal_pos = baseline['selected_site']['xyz_A']
    for metal, q in [('La', ligand_charge + 3), ('Ca', ligand_charge + 2)]:
        rows = [(metal, *metal_pos)] + atoms
        ne = old.require_closed_shell(metal, rows, q)
        xp = output / f'{stem}_{metal}.xyz'
        xp.write_text('\n'.join([str(len(rows)), PROTOCOL] +
                               [f'{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}' for a in rows]) + '\n')
        inp = output / f'{stem}_{metal}.inp'
        inp.write_text(f'! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n* xyzfile {q} 1 {xp.name}\n')
        result['outputs'][metal] = {"xyz": record(xp), "input": record(inp), "charge": q,
                                    "multiplicity": 1, "all_electron_count_for_parity": ne}
        result['tasks'].append({"task_id": metal, "input": record(inp), "xyz": record(xp),
                                "output_path": str(output.resolve() / f'{stem}_{metal}.out')})
    result['paired_invariants'] = paired(result['outputs']['La']['xyz']['path'], result['outputs']['Ca']['xyz']['path'],
                                         ligand_charge + 3, ligand_charge + 2)
    from affordable_common import xyz
    serialized = xyz(result['outputs']['La']['xyz']['path'])
    displacements=[]
    for item in graph_record['source_to_qm']:
        if item['kind']!='source' or item['source']['element'] in ('H','D'):
            continue
        s=item['source']; k=(s['chain_index'],s['residue_index'],s['atom'])
        displacements.append(float(np.linalg.norm(np.array(serialized[item['qm_index']][1:])-
                                                 np.array(tuple(graph.atoms[k].pos)))))
    result['heavy_coordinate_max_displacement_A'] = max(displacements,default=0.)
    if result['heavy_coordinate_max_displacement_A']>1e-8:
        raise InvalidArtifact('serialized source heavy coordinates moved')
    write_new(output / 'repair_manifest.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline-manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--topology', type=Path, required=True)
    a = p.parse_args()
    result = repair(a.baseline_manifest, a.output, a.topology)
    print(f"{result['protocol_id']}: {result['paired_invariants']['atom_count']} atoms, prepared")


if __name__ == '__main__':
    main()
