#!/usr/bin/env python3
"""
Rebuild results/ln_class_hits.tsv from results/all_results.jsonl,
enriched with spicy_lams Sharur DB classifications and FP annotations.

Standalone script — call from update_results_jsonl.py or directly.
"""
from __future__ import annotations
import json, csv, sys
from pathlib import Path

ALCH = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/"
            "lanthanide_binding/on_density_scanner/alchemical_bvs")
JSONL = ALCH / "results" / "all_results.jsonl"
FP_TSV = ALCH / "results" / "known_false_positives.tsv"
OUT_TSV = ALCH / "results" / "ln_class_hits.tsv"

SHARUR = Path("/groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/dft_la_candidates")
CLASS_TSVS = [SHARUR / "holo_low_CN" / "_classifications.tsv",
              SHARUR / "holo" / "_classifications.tsv"]

FIELDS = ['rank','stem','source','iptm','n_atoms','ddE_kcal','class',
          'spicy_class','pfam','ko','novel','length_aa','bin_id','gene_id',
          'fp_note','cif_path','workspace_path']


def main():
    # Load spicy_lams classifications, keyed by sanitized stem
    classmap = {}
    for path in CLASS_TSVS:
        if not path.exists():
            print(f"warn: classifications TSV missing: {path}", file=sys.stderr)
            continue
        with open(path) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                stem = row['cif_name'].replace('.cif', '').replace('.', '_').replace('__', '_')
                classmap[stem] = row

    # Load Ln-class entries from JSONL (deduped by stem; first occurrence wins)
    hits = {}
    if not JSONL.exists():
        print(f"error: {JSONL} not found", file=sys.stderr)
        sys.exit(1)
    with open(JSONL) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get('class') in ('Ln-preferring', 'Ln-evolved'):
                s = r['stem']
                if s not in hits:
                    hits[s] = r

    # Load FP catalogue
    fp = {}
    if FP_TSV.exists():
        with open(FP_TSV) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    fp[parts[0]] = parts[1].strip()

    # Assemble rows
    rows = []
    for stem, r in hits.items():
        cands = [stem, stem.removeprefix('NZ_')]
        sl = {}
        for c in cands:
            if c in classmap:
                sl = classmap[c]
                break
        iptm = r.get('iptm')
        iptm_s = f'{iptm:.3f}' if isinstance(iptm, (int, float)) else '-'
        cls = r['class']
        fp_note = fp.get(stem, '')
        if fp_note and 'FALSE_POSITIVE' not in cls:
            cls = f'FALSE_POSITIVE ({cls})'
        rows.append({
            'stem': stem,
            'source': r.get('source', 'spicylams'),
            'iptm': iptm_s,
            'n_atoms': r.get('n_atoms', ''),
            'ddE_kcal': f'{r["ddE_kcal"]:+.2f}',
            'class': cls,
            'spicy_class': sl.get('class', '-'),
            'pfam': sl.get('pfam', '-'),
            'ko': sl.get('ko', '-'),
            'novel': sl.get('novel', '-'),
            'length_aa': sl.get('length_aa', '-'),
            'bin_id': sl.get('bin_id', '-'),
            'gene_id': sl.get('gene_id_prod', '-'),
            'fp_note': fp_note,
            'cif_path': r.get('cif_path', ''),
            'workspace_path': r.get('workspace_path', ''),
        })

    # Sort by ΔΔE descending
    rows.sort(key=lambda x: -float(x['ddE_kcal']))

    # Write
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, 'w') as f:
        f.write('# DFT Ca²⁺/Ln³⁺ discriminator — Ln-class hits (ΔΔE ≥ +5 kcal/mol), '
                'enriched with spicy_lams Sharur DB annotations\n')
        f.write('\t'.join(FIELDS) + '\n')
        for i, r in enumerate(rows, 1):
            r['rank'] = i
            f.write('\t'.join(str(r[k]) for k in FIELDS) + '\n')

    print(f"Tracker rebuilt: {len(rows)} Ln-class rows to {OUT_TSV}", file=sys.stderr)


if __name__ == "__main__":
    main()
