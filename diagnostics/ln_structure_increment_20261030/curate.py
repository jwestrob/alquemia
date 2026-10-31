"""Bounded real-record curation; metadata only, no selectivity fitting.

Reuse the historical parser's binary layout and donor SMARTS, preserving offsets
and annotations discarded by its original tabulation. Outputs stay private/local.
"""
import ast
import csv
import hashlib
import json
import os
import re
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize as std

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = Path(os.environ.get('NIKASHA_WORKSPACE', str(ROOT / 'workspaces/ln_structure_increment_20261030')))
RECORDS = Path(os.environ.get('NIKASHA_RECORD_DIR', str(HERE)))
OLD = ROOT / 'workspaces/ln_selectivity_calibration_20261027'
KC = 0.00198720425864083 * 298.15 * __import__('math').log(10)


def read(p):
    return list(csv.DictReader(open(p), delimiter='\t'))


def write(p, rows):
    assert rows, p
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(p, 'w') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def donor_patterns():
    tree = ast.parse((ROOT / 'diagnostics/ln_selectivity_calibration_20261027/fit_donors.py').read_text())
    return ast.literal_eval(next(n.value for n in tree.body if isinstance(n, ast.Assign)
                                and any(isinstance(t, ast.Name) and t.id == 'SMARTS' for t in n.targets)))


def identity(smi, stereo=True):
    m = Chem.MolFromSmiles(smi)
    m = std.Cleanup(m)
    m = std.Uncharger().uncharge(m)
    m = std.TautomerEnumerator().Canonicalize(m)
    if not stereo:
        Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m, isomericSmiles=stereo)


def family(m):
    p = {k: Chem.MolFromSmarts(v) for k, v in donor_patterns().items()}
    c = {k: len(m.GetSubstructMatches(q)) for k, q in p.items()}
    if any(len(r) >= 9 for r in m.GetRingInfo().AtomRings()):
        return 'aza_or_oxa_macrocycles'
    if c['amine_N'] and c['carboxyl']:
        return 'aminocarboxylates'  # includes every EDTA backbone, IDA/NTA/DTPA/TTHA
    if c['phosphonate']:
        return 'phosphates_nucleotides'
    if c['pyridine_N'] or m.HasSubstructMatch(Chem.MolFromSmarts('[n+]([O-])')):
        return 'aromatic_nitrogen_carboxylates'
    if any(a.GetSymbol() == 'N' for a in m.GetAtoms()) and c['phenol']:
        return 'nitroso_phenols'
    if c['carboxyl']:
        if any(a.GetIsAromatic() for a in m.GetAtoms()):
            return 'aromatic_carboxylates'
        if c['alcohol']:
            return 'hydroxy_carboxylates'
        if c['ether_O']:
            return 'ether_carboxylates'
        return 'aliphatic_carboxylates'
    if any(a.GetSymbol() == 'C' and a.IsInRing() for a in m.GetAtoms()):
        return 'cyclic_oxygen_chelators'
    if c['carbonyl_O']:
        return 'acyclic_carbonyl_chelators'
    return 'small_inorganic_anions'


def add_acetates(smi):
    m = Chem.MolFromSmiles(smi)
    ns = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'N']
    assert len(ns) == 2
    for n in ns:
        for _ in range(2):
            offset = m.GetNumAtoms()
            rw = Chem.RWMol(Chem.CombineMols(m, Chem.MolFromSmiles('CC(=O)O')))
            rw.AddBond(n, offset, Chem.BondType.SINGLE)
            m = rw.GetMol()
            Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m, isomericSmiles=True)


def raw_records():
    d = OLD / 'nist46_dat'
    def lines(fn):
        return (d / fn).read_bytes().decode('latin1').splitlines()
    metals = {}
    for line in lines('F8548_CODE2.DAT')[1:]:
        m = re.match(r'\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if m:
            metals[int(m[1])] = m[4]
    equilibria = {}
    for line in lines('F8549_CODE3.DAT')[1:]:
        m = re.match(r'\s*(\d+)\s+(\S+)\s+(\S+)', line)
        if m:
            equilibria[int(m[1])] = (m[2], m[3])
    conditions = {int(i): v.split() for i, v in re.findall(
        r'(\d+) ([KHS\*] [\dtv\*]+ [\d\.\*]+)', (d / 'F8550_CODE4.DAT').read_bytes().decode('latin1'))}
    notes = {}
    for line in lines('F8552_CODE6.DAT'):
        m = re.match(r'\s*(\d+)\s+(\S+)\s+(.*)', line)
        if m:
            notes[int(m[1])] = m[2] + ': ' + m[3].strip()
    v = (d / 'F8558_VALUES.DAT').read_bytes()
    out = []
    for offset in range(15, len(v) - 14, 15):
        r = v[offset:offset+15]
        metal = metals.get(r[2], '')
        eq = int.from_bytes(r[3:5], 'little')
        cond = conditions.get(r[5], [])
        if metal not in ('La/3+', 'Dy/3+') or equilibria.get(eq) != ('ML', '[ML]/[M][L]'):
            continue
        if len(cond) != 3 or cond[:2] != ['K', '25']:
            continue
        out.append(dict(lig_id=str(int.from_bytes(r[:2], 'little')), metal=metal[:2],
                        I=float(cond[2]), T_C=25, equilibrium='ML=M.L', stoichiometry='1:1',
                        value_raw=r[6:13].decode('latin1').strip(),
                        byte_offset=offset, condition_code=r[5], equilibrium_code=eq,
                        precision_code=r[13], annotation_code=r[14],
                        annotation=notes.get(r[14], '') if r[14] > 1 else ''))
    return out


def main():
    W.mkdir(exist_ok=True)
    RECORDS.mkdir(exist_ok=True,parents=True)
    pairs = read(OLD / 'pairs.tsv')
    smiles = {r['lig_id']: r for r in read(OLD / 'smiles.tsv')}
    recovered = {}
    cyclic = {'599': 'cis-cyclohexane-1,2-diamine', '600': 'trans-cyclohexane-1,2-diamine',
              '601': 'cis-cyclohexane-1,3-diamine', '602': 'trans-cyclohexane-1,3-diamine',
              '603': 'cis-cyclohexane-1,4-diamine', '604': 'trans-cyclohexane-1,4-diamine',
              '597': 'trans-cyclopentane-1,2-diamine', '598': 'trans-cyclopentane-1,3-diamine'}
    for lid, nm in cyclic.items():
        j = json.loads((W / 'sources' / (nm + '.json')).read_text())
        assert j['status'] == 'SUCCESS'
        recovered[lid] = (add_acetates(j['smiles']), 'OPSIN systematic relative-stereo diamine + two acetate arms per N; ' + nm)
    for lid, smi in {'614': 'N(CCCCCCN(CC(=O)O)CC(=O)O)(CC(=O)O)CC(=O)O',
                     '621': 'N(C(C(=O)O)CN(CC(=O)O)CC(=O)O)(CC(=O)O)CC(=O)O',
                     '632': 'N(CC(O)CN(CC(=O)O)CC(=O)O)(CC(=O)O)CC(=O)O'}.items():
        recovered[lid] = (smi, 'NIST systematic backbone name; explicit graph construction; DL/achiral as named')
    historical_ids = set(r['lig_id'] for r in read(OLD/'features.tsv'))
    ligands, exclusions = [], []
    for r in pairs:
        lid = r['lig_id']
        sr = smiles[lid]
        smi, source = recovered.get(lid, (sr['smiles'], 'historical PubChem formula-checked; query=' + sr['query']))
        reason = ''
        if lid not in recovered and sr['formula_ok'] != 'yes':
            reason = 'structure lookup/formula unresolved; outside bounded tetraacetate recovery'
        if not reason:
            m = Chem.MolFromSmiles(smi)
            formula = re.sub(r'([A-Z][a-z]?)1(?![0-9])', r'\1', r['formula'])
            if m is None or rdMolDescriptors.CalcMolFormula(m) != formula:
                reason = 'graph formula differs from NIST chemical identity'
            elif re.match(r'(?:D-|L-|cis-|trans-|\(1R)', r['name']) and lid not in recovered:
                # Formula matching alone did not verify these stereospecific identities.
                if '?' in [ch for _, ch in Chem.FindMolChiralCenters(m, includeUnassigned=True)]:
                    reason = 'historical SMILES lacks named stereochemistry'
                if r['name'].startswith('trans-') and not any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in m.GetBonds()) and not Chem.FindMolChiralCenters(m):
                    reason = 'historical SMILES lacks named ring stereochemistry'
        if reason:
            exclusions.append(dict(lig_id=lid, name=r['name'], stage='identity', reason=reason))
            continue
        m = Chem.MolFromSmiles(smi)
        # Charge belongs to the equilibrium L, not the protonated drawing.
        match = re.fullmatch(r'H(\d*)L', r['prot'])
        if not match:
            exclusions.append(dict(lig_id=lid, name=r['name'], stage='charge', reason='unresolved NIST L charge'))
            continue
        ligands.append(dict(lig_id=lid, name=r['name'], formula=r['formula'], protonation=r['prot'],
                            smiles=Chem.MolToSmiles(m, isomericSmiles=True),
                            chemical_identity=lid+'|'+r['name']+'|'+Chem.MolToSmiles(m, isomericSmiles=True),
                            normalized_identity=identity(smi, False),
                            stereo_identity=identity(smi, True), ligand_charge=-int(match[1] or 1),
                            family=family(m), tetraacetate_challenge=int('dinitrilotetraacetic acid' in r['name']),
                            historical_101=int(lid in historical_ids), structure_provenance=source))
    raw = raw_records()
    write(W / 'raw_direct_records.tsv', raw)
    by_key = defaultdict(lambda: defaultdict(list))
    ligmap = {r['lig_id']: r for r in ligands}
    for r in raw:
        if r['lig_id'] in ligmap:
            by_key[r['lig_id'], r['condition_code']][r['metal']].append(r)
    records = []
    for (lid, cond), metals in sorted(by_key.items()):
        if set(metals) != {'La', 'Dy'}:
            continue
        reason = ''
        if any(len(v) != 1 for v in metals.values()):
            reason = 'multiple values within condition; no unique paired measurement'
        if reason:
            exclusions.append(dict(lig_id=lid, stage='measurement', condition_code=cond, reason=reason))
            continue
        la, dy = metals['La'][0], metals['Dy'][0]
        if not all(re.fullmatch(r'[+-]?\d*\.?\d+', r['value_raw']) for r in (la, dy)):
            reason = 'flagged/estimated/censored native constant'
        if la['annotation_code'] != dy['annotation_code']:
            reason = 'different native condition/species annotations; unverified match'
        if reason:
            exclusions.append(dict(lig_id=lid, stage='measurement', condition_code=cond,
                                   La_annotation=la['annotation'], Dy_annotation=dy['annotation'], reason=reason))
            continue
        lk, dk = float(la['value_raw']), float(dy['value_raw'])
        records.append(dict(**ligmap[lid], record_id=f'nist46_{lid}_C{cond}', T_K=298.15,
                            I_M=la['I'], logK_La=lk, logK_Dy=dk, y_kcal=-KC*(dk-lk),
                            equilibrium='ML=M.L', metal_valence=3, stoichiometry='1:1',
                            La_offset=la['byte_offset'], Dy_offset=dy['byte_offset'],
                            native_annotation=la['annotation'], annotation_code=la['annotation_code'],
                            source='NIST SRD46 v8 F8558_VALUES.DAT; exact native byte offsets',
                            provenance_limit='critically selected compilation; row-level original-study link not decoded'))
    # Select one primary condition without looking at labels; retain every other condition.
    for lid in ligmap:
        rr = [r for r in records if r['lig_id'] == lid]
        if rr:
            selected = min(rr, key=lambda r: (abs(r['I_M']-.1), r['I_M'], r['record_id']))
            for r in rr:
                r['primary'] = int(r is selected)
        else:
            exclusions.append(dict(lig_id=lid, stage='measurement', reason='no admitted direct La/Dy condition'))
    primary = [r for r in records if r['primary']]
    write(W/'evaluation_records.tsv', records)
    write(W/'primary.tsv', primary)
    write(W/'exclusions.tsv', exclusions)
    write(RECORDS/'groups.tsv', [{k:r[k] for k in ('lig_id','name','family','tetraacetate_challenge','structure_provenance')}
                              for r in primary])
    split = {g: {'test_ligands':[r['lig_id'] for r in primary if r['family']==g],
                 'train_ligands':[r['lig_id'] for r in primary if r['family']!=g]}
             for g in sorted(set(r['family'] for r in primary))}
    (RECORDS/'splits.json').write_text(json.dumps(split,indent=2)+'\n')
    pins = [HERE/'PLAN.md', RECORDS/'groups.tsv', RECORDS/'splits.json', W/'evaluation_records.tsv', W/'primary.tsv', W/'exclusions.tsv',
            *[(OLD/'nist46_dat'/f) for f in ('F8558_VALUES.DAT','F8547_CODE1.DAT','F8548_CODE2.DAT','F8549_CODE3.DAT','F8550_CODE4.DAT','F8552_CODE6.DAT')]]
    manifest = {'frozen_before_new_fits':True,'grouping_version':2,
                'post_initial_results_chemistry_audit_repair':True,
                'kcal_per_log_unit':KC,'n_primary':len(primary),
                'n_conditions':len(records),'n_families':len(split),
                'tetraacetates':sum(r['tetraacetate_challenge'] for r in primary),
                'pins':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}}
    (RECORDS/'FREEZE.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k!='pins'},indent=2))


if __name__=='__main__':
    main()
