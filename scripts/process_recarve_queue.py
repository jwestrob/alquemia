#!/usr/bin/env python3
"""
process_recarve_queue.py — re-carve and re-submit all entries in
alchemical_bvs/results/recarve_queue.tsv after the 2026-05-11 carve_generic.py
fixes (S in donor scan; CYS/HIS/THR/MET/LYS/ARG in SIDECHAIN_QM_ATOMS;
--first-shell-cut CLI flag).

Per entry:
  1. Verify the protonated PDB exists.
  2. Skip if a SLURM job for this stem is already queued/running.
  3. Move stale carve/ORCA outputs into <qm_dir>/_pre_recarve_<date>/ (audit trail).
  4. Choose carver:
       - PQQ atom within 4 Å of the metal → carve_with_pqq.py
       - Else                            → carve_generic.py at 3.0 Å,
                                            retry at 3.5 Å if N_atoms == 1
  5. If new La_qm.xyz has > 1 atom → emit water SP, sbatch the panel.
  6. Otherwise classify and log:
       - closest O/N/S > 4 Å  → SOLVENT_EXCLUDE  (don't submit; truly empty pocket)
       - else                 → CARVE_AMBIGUOUS  (don't submit; needs manual look)
  7. Write per-stem outcome to recarve_log_<UTC>.tsv next to the queue.

Idempotent: re-running re-uses the protonated PDB and won't double-submit if
the stem already has a job in squeue or already has a finished La SP.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import gemmi  # type: ignore

ALCH = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
QUEUE_TSV = ALCH / 'results' / 'recarve_queue.tsv'
SCAN_PY = '/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'
CARVE_GENERIC = ALCH / 'scripts' / 'carve_generic.py'
CARVE_PQQ = ALCH / 'scripts' / 'carve_with_pqq.py'

SOLVENT_EXCLUDE_THRESHOLD = 4.0   # Å — if closest donor > this, call it solvent
PQQ_DETECT_RADIUS = 4.0           # Å — if any PQQ atom within this, use PQQ carver
RETRY_WIDE_CUTOFF = 3.5           # Å — second-try cutoff for borderline cases

ARTIFACT_GLOBS = [
    '*_La_qm.xyz', '*_Ca_qm.xyz', '*_apo_qm.xyz', 'bulk_water.xyz',
    'sp_*.inp', 'sp_*.out', 'sp_*.gbw', 'sp_*.densities*',
    'sp_*.cpcm*', 'sp_*.bibtex', 'sp_*.property.txt',
    'submit_*.sh', 'slurm_*.out', 'slurm_*.err',
]


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')


def find_metal(model) -> tuple[str, int, str, gemmi.Position] | None:
    """Return (chain, resnum, resname, position) for the first La/Ce in the model."""
    for chain in model:
        for res in chain:
            for at in res:
                if at.element.name in ('La', 'Ce', 'Y'):
                    return (chain.name, res.seqid.num, res.name, at.pos)
    return None


def closest_donor_dist(model, la_pos: gemmi.Position) -> tuple[float, str]:
    """Return (distance_Å, descriptor) of the closest non-metal O/N/S atom."""
    best = (1e9, '-')
    for chain in model:
        for res in chain:
            if res.name in ('LA', 'CA', 'CE', 'HOH', 'WAT'):
                continue
            for at in res:
                if at.element.name not in ('O', 'N', 'S'):
                    continue
                d = la_pos.dist(at.pos)
                if d < best[0]:
                    best = (d, f'{res.name}{res.seqid.num}/{at.name}')
    return best


def has_pqq_contact(model, la_pos: gemmi.Position, radius: float) -> bool:
    for chain in model:
        for res in chain:
            if res.name != 'PQQ':
                continue
            for at in res:
                if la_pos.dist(at.pos) <= radius:
                    return True
    return False


def stem_in_squeue(stem: str) -> bool:
    try:
        out = subprocess.check_output(
            ['squeue', '-u', os.environ.get('USER', 'jwestrob'),
             '-h', '-o', '%j'],
            text=True, timeout=30,
        )
    except Exception:
        return False
    name = f'orca_{stem}'
    for line in out.splitlines():
        if line.strip() == name:
            return True
    return False


def la_panel_already_done(qm_dir: Path, stem: str) -> bool:
    """Check whether a working La/Ca/apo/water panel already exists."""
    needed = [f'sp_{stem}_{k}.out' for k in ('La', 'Ca', 'apo', 'water')]
    for n in needed:
        f = qm_dir / n
        if not f.exists():
            return False
        try:
            tail = subprocess.check_output(['tail', '-50', str(f)], text=True, timeout=10)
        except Exception:
            return False
        if 'FINAL SINGLE POINT ENERGY' not in tail:
            return False
    return True


def stash_artifacts(qm_dir: Path, backup_root_name: str) -> Path:
    backup = qm_dir / backup_root_name
    backup.mkdir(exist_ok=True)
    for pattern in ARTIFACT_GLOBS:
        for f in qm_dir.glob(pattern):
            # Don't move the protonated.pdb (we need it for re-carve)
            if f.name.endswith('_protonated.pdb'):
                continue
            shutil.move(str(f), str(backup / f.name))
    return backup


def run_carve_generic(pdb: Path, qm_dir: Path, stem: str, first_shell_cut: float | None) -> int:
    """Run carve_generic.py and return the new La_qm.xyz atom count (0 if missing)."""
    cmd = [SCAN_PY, str(CARVE_GENERIC), str(pdb), str(qm_dir), '--stem', stem]
    if first_shell_cut is not None:
        cmd += ['--first-shell-cut', str(first_shell_cut)]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if res.returncode != 0:
        print(f'    [carve_generic FAILED] rc={res.returncode}: {res.stderr.strip()[-400:]}')
        return -1
    xyz = qm_dir / f'{stem}_La_qm.xyz'
    if not xyz.exists():
        return -1
    try:
        with xyz.open() as fh:
            return int(fh.readline().strip())
    except Exception:
        return -1


def run_carve_pqq(pdb: Path, qm_dir: Path, stem: str, metal_chain: str, metal_resname: str) -> int:
    cmd = [SCAN_PY, str(CARVE_PQQ), str(pdb), str(qm_dir), '--stem', stem,
           '--metal-chain', metal_chain, '--metal-resname', metal_resname]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if res.returncode != 0:
        print(f'    [carve_with_pqq FAILED] rc={res.returncode}: {res.stderr.strip()[-400:]}')
        return -1
    xyz = qm_dir / f'{stem}_La_qm.xyz'
    if not xyz.exists():
        return -1
    try:
        with xyz.open() as fh:
            return int(fh.readline().strip())
    except Exception:
        return -1


def emit_water_sp(qm_dir: Path, stem: str) -> None:
    """Write bulk_water.xyz and sp_<stem>_water.inp matching process_inbox.sh."""
    (qm_dir / 'bulk_water.xyz').write_text(
        '3\nH2O reference\n'
        'O   0.000000  0.000000  0.000000\n'
        'H   0.756950  0.000000  0.585822\n'
        'H  -0.756950  0.000000  0.585822\n'
    )
    (qm_dir / f'sp_{stem}_water.inp').write_text(
        '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n'
        '%maxcore 8000\n'
        '* xyzfile 0 1 bulk_water.xyz\n'
    )
    # Ensure NoAutostart is on every sp_*.inp (matches process_inbox.sh patch)
    for inp in qm_dir.glob(f'sp_{stem}_*.inp'):
        txt = inp.read_text()
        if 'NoAutostart' not in txt:
            inp.write_text(txt.replace('r2SCAN-3c CPCM', 'r2SCAN-3c NoAutostart CPCM', 1))
    # Patch the submit script's for-loop to include water if it doesn't already
    sub = qm_dir / f'submit_{stem}.sh'
    if sub.exists():
        txt = sub.read_text()
        if 'for kind in La Ca apo water' not in txt:
            sub.write_text(txt.replace('for kind in La Ca apo;', 'for kind in La Ca apo water;'))


def sbatch_submit(qm_dir: Path, stem: str, dry_run: bool) -> str:
    submit_sh = qm_dir / f'submit_{stem}.sh'
    if not submit_sh.exists():
        return 'no_submit_script'
    if dry_run:
        return 'DRY_RUN'
    res = subprocess.run(['sbatch', '--parsable', str(submit_sh)],
                         capture_output=True, text=True, timeout=60)
    out = (res.stdout or '').strip()
    err = (res.stderr or '').strip()
    if res.returncode == 0 and out.split(';')[0].isdigit():
        return out
    return f'sbatch_fail:{err or out}'


def process_entry(row: dict, *, max_submit: int, submitted: int, dry_run: bool,
                  backup_name: str) -> dict:
    stem = row['stem']
    qm_dir = Path(row['qm_dir'])
    pdb = Path(row['protonated_pdb'])

    rec = {
        'stem': stem,
        'outcome': '',
        'n_atoms_la_new': None,
        'closest_donor': None,
        'cutoff_used': None,
        'carver': None,
        'job_id': None,
        'note': '',
    }

    if not qm_dir.is_dir():
        rec['outcome'] = 'SKIP_NO_QM_DIR'; rec['note'] = str(qm_dir); return rec
    if not pdb.is_file():
        rec['outcome'] = 'SKIP_NO_PROTONATED_PDB'; rec['note'] = str(pdb); return rec
    if stem_in_squeue(stem):
        rec['outcome'] = 'SKIP_ALREADY_QUEUED'; return rec
    if la_panel_already_done(qm_dir, stem):
        rec['outcome'] = 'SKIP_PANEL_DONE'; return rec

    # Diagnostic: read PDB once, find metal + closest donor + PQQ
    try:
        st = gemmi.read_structure(str(pdb))
        st.setup_entities()
        model = st[0]
    except Exception as e:
        rec['outcome'] = 'SKIP_PDB_PARSE_ERROR'; rec['note'] = str(e); return rec

    metal = find_metal(model)
    if metal is None:
        rec['outcome'] = 'SKIP_NO_METAL_IN_PDB'; return rec
    m_chain, m_resnum, m_resname, m_pos = metal
    closest_d, closest_atom = closest_donor_dist(model, m_pos)
    rec['closest_donor'] = f'{closest_atom}@{closest_d:.2f}A'
    pqq_near = has_pqq_contact(model, m_pos, PQQ_DETECT_RADIUS)

    # Stash existing artifacts (only after we've decided to act on this entry)
    stash_artifacts(qm_dir, backup_name)

    # Choose carver
    if pqq_near:
        rec['carver'] = 'carve_with_pqq'
        n_atoms = run_carve_pqq(pdb, qm_dir, stem, m_chain, m_resname)
        rec['cutoff_used'] = 'pqq_default'
    else:
        rec['carver'] = 'carve_generic'
        n_atoms = run_carve_generic(pdb, qm_dir, stem, first_shell_cut=None)  # use default 3.0
        rec['cutoff_used'] = '3.0'
        if n_atoms <= 1:
            # Retry at 3.5 Å if the donor sits in the borderline band
            if closest_d <= RETRY_WIDE_CUTOFF:
                n_atoms = run_carve_generic(pdb, qm_dir, stem,
                                            first_shell_cut=RETRY_WIDE_CUTOFF)
                rec['cutoff_used'] = '3.5_retry'
    rec['n_atoms_la_new'] = n_atoms

    if n_atoms <= 1:
        # No first-shell — triage by absolute closest distance
        if closest_d > SOLVENT_EXCLUDE_THRESHOLD:
            rec['outcome'] = 'SOLVENT_EXCLUDE'
        else:
            rec['outcome'] = 'CARVE_AMBIGUOUS'
        return rec

    # We have a real carve. Add water SP and submit.
    emit_water_sp(qm_dir, stem)
    if submitted >= max_submit:
        rec['outcome'] = 'CARVED_DEFERRED_CAP'
        rec['note'] = f'hit MAX_SUBMIT={max_submit}; submit later'
        return rec
    job_id = sbatch_submit(qm_dir, stem, dry_run)
    rec['job_id'] = job_id
    if job_id.isdigit() or job_id == 'DRY_RUN':
        rec['outcome'] = 'RECARVED_AND_SUBMITTED' if not dry_run else 'RECARVED_DRY'
    else:
        rec['outcome'] = 'SUBMIT_FAILED'
        rec['note'] = job_id
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--queue', type=Path, default=QUEUE_TSV)
    ap.add_argument('--max-submit', type=int, default=int(os.environ.get('MAX_SUBMIT', '100')))
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--limit', type=int, default=None, help='Process only first N entries.')
    ap.add_argument('--only-stem', default=None, help='Process only this stem (debug).')
    args = ap.parse_args()

    if not args.queue.is_file():
        print(f'[FATAL] queue not found: {args.queue}', file=sys.stderr)
        return 2

    rows: list[dict] = []
    with args.queue.open() as fh:
        rdr = csv.DictReader((L for L in fh if not L.startswith('#')), delimiter='\t')
        for r in rdr:
            if not r.get('stem'):
                continue
            rows.append(r)

    if args.only_stem:
        rows = [r for r in rows if r['stem'] == args.only_stem]
    if args.limit:
        rows = rows[:args.limit]

    print(f'[INFO] queue size: {len(rows)} entries (max_submit={args.max_submit}, dry_run={args.dry_run})')

    backup_name = f'_pre_recarve_{now_utc()}'
    log_path = args.queue.parent / f'recarve_log_{now_utc()}.tsv'
    submitted = 0
    counts: dict[str, int] = {}

    with log_path.open('w') as out:
        out.write('stem\toutcome\tn_atoms_la_new\tclosest_donor\tcutoff_used\tcarver\tjob_id\tnote\n')
        for i, row in enumerate(rows, 1):
            print(f'[{i}/{len(rows)}] {row["stem"]}')
            rec = process_entry(row,
                                max_submit=args.max_submit,
                                submitted=submitted,
                                dry_run=args.dry_run,
                                backup_name=backup_name)
            counts[rec['outcome']] = counts.get(rec['outcome'], 0) + 1
            if rec['outcome'] in ('RECARVED_AND_SUBMITTED',):
                submitted += 1
            print(f'    → {rec["outcome"]} '
                  f'n_atoms={rec["n_atoms_la_new"]} '
                  f'closest={rec["closest_donor"]} '
                  f'cutoff={rec["cutoff_used"]} '
                  f'carver={rec["carver"]} '
                  f'job={rec["job_id"]}')
            out.write('\t'.join([
                rec['stem'], rec['outcome'],
                str(rec['n_atoms_la_new']) if rec['n_atoms_la_new'] is not None else '',
                rec['closest_donor'] or '',
                rec['cutoff_used'] or '',
                rec['carver'] or '',
                rec['job_id'] or '',
                rec['note'],
            ]) + '\n')
            out.flush()

    print()
    print(f'[INFO] outcome histogram: {json.dumps(counts, indent=2)}')
    print(f'[INFO] {submitted} jobs sbatched')
    print(f'[INFO] log: {log_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
