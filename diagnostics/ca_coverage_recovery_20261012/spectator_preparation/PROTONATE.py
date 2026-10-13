"""Research-only one/two-Ca wrapper around the exact original protonator.

No scientific import occurs before the actual Slurm guard.  The only chemistry
policy extension is admission of one explicitly preserved spectator Ca residue.
"""
import os
from pathlib import Path
import random

POLICY = 'standard_residues_only_explicit_Ca_spectator_no_CCD_v1'
SUBPROTOCOL = 'pdbfixer_standard_only_explicit_Ca_spectator_rng20260914_CPU1_v1'


def protonate(source, output, *, pins_path, ph, legacy, spectator_count):
    assert os.environ.get('SLURM_JOB_ID'), 'scientific preparation requires Slurm'
    assert spectator_count in (0, 1) and ph == 7.0
    import gemmi
    from openmm import Platform, app
    import pdbfixer as pdbfixer_module
    from pdbfixer import PDBFixer
    pins = legacy.read_object(pins_path)
    assert pins['protonation_subprotocol'] == legacy.EXPECTED_SUBPROTOCOL
    canonical_pin = pins['canonical_helpers']['protonate_cif']
    assert legacy.verify_file_record(canonical_pin, 'canonical protonator') == Path(legacy.canonical.__file__).resolve()
    legacy.verify_file_record(pins['experiment_helpers']['protonate_standard_only'], 'historical wrapper')
    expected = {'LA': 1, 'PQQ': 1, **({'CA': 1} if spectator_count else {})}

    def inspect(path):
        st = gemmi.read_structure(str(path)); assert len(st) == 1
        counts = {k: 0 for k in expected}; hydrogen = {k: 0 for k in expected}
        calcium = 0
        for chain in st[0]:
            for residue in chain:
                name = residue.name.strip().upper()
                assert name in legacy.STANDARD_PROTEIN_RESIDUES or name in expected, ('unsupported residue', name)
                if name in expected:
                    counts[name] += 1
                    hydrogen[name] += sum(a.element.name in ('H', 'D') for a in residue)
                    if name in ('LA', 'CA'):
                        assert len(residue) == 1 and residue[0].element.name == 'Ca', 'source ions must both remain native Ca'
                        calcium += 1
        assert counts == expected and not any(hydrogen.values())
        assert calcium == spectator_count + 1
        return {'residue_counts': counts, 'nonstandard_hydrogens': hydrogen, 'native_Ca_count': calcium}

    before = inspect(source)
    cpu = Platform.getPlatformByName('CPU')
    calls = {'constructor': 0, 'add_missing_hydrogens': 0, 'ccd_suppression': 0}

    class StandardOnly(PDBFixer):
        def __init__(self, *args, **kwargs):
            requested = kwargs.get('platform')
            assert requested is None or requested.getName() == 'CPU'
            kwargs['platform'] = cpu; calls['constructor'] += 1
            super().__init__(*args, **kwargs)

        def _downloadNonstandardDefinitions(self):
            calls['ccd_suppression'] += 1
            app.Modeller._loadStandardHydrogenDefinitions()
            return {}

        def addMissingHydrogens(self, pH=7.0, forcefield=None):
            assert forcefield is None and pH == 7.0
            assert cpu.getPropertyDefaultValue('Threads') == '1'
            calls['add_missing_hydrogens'] += 1
            return super().addMissingHydrogens(pH=pH, forcefield=None)

    original = pdbfixer_module.PDBFixer
    rng = random.getstate(); threads = cpu.getPropertyDefaultValue('Threads')
    random.seed(legacy.RNG_SEED); cpu.setPropertyDefaultValue('Threads', '1')
    pdbfixer_module.PDBFixer = StandardOnly
    try:
        legacy.canonical.protonate(source, output, ph=ph, add_missing_residues=False)
    finally:
        pdbfixer_module.PDBFixer = original
        cpu.setPropertyDefaultValue('Threads', threads); random.setstate(rng)
    assert calls == {'constructor': 1, 'add_missing_hydrogens': 1, 'ccd_suppression': 1}
    after = inspect(output)
    path = legacy.canonical.protonation_manifest_path(output)
    manifest = legacy.read_object(path)
    assert manifest['repaired_missing_atom_count'] == manifest['repaired_missing_terminal_atom_count'] == 0
    manifest.update(canonical_protonator=canonical_pin,
        experiment_wrapper={'path': str(Path(__file__).resolve()), 'sha256': legacy.sha256_file(Path(__file__).resolve())},
        experiment_protonation_protocol_id=SUBPROTOCOL,
        historical_protonation_wrapper=pins['experiment_helpers']['protonate_standard_only'],
        nonstandard_definition_policy={'policy_id': POLICY, 'before': before, 'after': after,
            'suppressed_method': 'PDBFixer._downloadNonstandardDefinitions', 'forcefield': None,
            'python_random_seed': legacy.RNG_SEED, 'openmm_platform': 'CPU', 'openmm_CPU_Threads': 1,
            'wrapper_call_ledger': calls, 'PQQ_hydrogen_policy': 'none_added_by_PDBFixer_all_supplied_by_pinned_carver',
            'spectator_count': spectator_count, 'metal_hydrogens': 'none', 'source_occupancy_changed_by_protonator': False})
    legacy.write_json_atomic(path, manifest)
    return path
