"""Review checks on real prepared chains and archived scalar values; no model calls."""
import copy
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,xyz,InvalidArtifact,HA_TO_KCAL
import lanm_global_occupancy as scoring

PREP=ROOT/'workspaces/lanm_global_occupancy_20260923/prepared_v1/manifest.json'

class ScoringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states=scoring.states(PREP)
        cls.pin,cls.s,cls.mapping=cls.states[0]
        cls.q0=np.array([a[1:] for a in xyz(verify(cls.s['endpoints']['La']['xyz']))])

    def test_real_inline_state_schema_and_all_origins_admitted(self):
        self.assertEqual(len(self.states),9)
        for _,s,m in self.states:
            q=np.array([a[1:] for a in xyz(verify(s['endpoints']['La']['xyz']))])
            g=scoring.geometry_check(q,q,m)
            self.assertTrue(g['admitted'],s['state_id'])
            self.assertFalse(g['origin_severe_heavy_clashes'])
        one=scoring.states(PREP,[self.s['state_id']])
        self.assertEqual(len(one),1)
        with self.assertRaises(InvalidArtifact):scoring.states(PREP,['absent_real_state'])

    def test_malformed_actual_state_rejected_without_new_geometry(self):
        m=read_json(PREP)
        m['states'][0]['charge']+=1
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad_inline.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'inline'):scoring.states(p)
            s=copy.deepcopy(self.s);s['endpoints']['Dy']['native_effective_multiplicity']=11
            sp=Path(d)/'bad_state.json';write_new(sp,s)
            p2=Path(d)/'bad_native.json';write_new(p2,{'states':[record(sp)]})
            with self.assertRaisesRegex(InvalidArtifact,'native'):scoring.states(p2)

    def test_cartesian_box_and_nonfinite_are_rejected(self):
        translated=self.q0+.2
        self.assertTrue(scoring.geometry_check(self.q0,translated,self.mapping)['admitted'])
        g=scoring.geometry_check(self.q0,self.q0+.46,self.mapping)
        self.assertFalse(g['admitted']);self.assertTrue(g['exceeds_box'])
        malformed=self.q0.copy();malformed[0,0]=np.nan
        self.assertFalse(scoring.geometry_check(self.q0,malformed,self.mapping)['admitted'])

    def test_actual_water_bond_and_source_chirality_guards(self):
        b=next(x for x in self.mapping['bonds'] if x['kind']=='water_covalent')
        i,j=b['indices'];q=self.q0.copy();v=q[j]-q[i];q[j]+=v/np.linalg.norm(v)*.25
        g=scoring.geometry_check(self.q0,q,self.mapping)
        self.assertFalse(g['admitted']);self.assertTrue(g['bad_bonds'])
        residues={}
        for a in self.mapping['atoms']:
            if a['kind']=='protein_source':
                residues.setdefault(a['source']['resid'],{})[a['source']['name']]=a['index']
        r=next(v for v in residues.values() if {'CA','N','C','CB'}<=v.keys())
        ca,n,c,cb=[r[k] for k in ('CA','N','C','CB')]
        q=self.q0.copy();normal=np.cross(q[n]-q[ca],q[c]-q[ca]);normal/=np.linalg.norm(normal)
        q[cb]-=2*np.dot(q[cb]-q[ca],normal)*normal
        g=scoring.geometry_check(self.q0,q,self.mapping)
        self.assertFalse(g['admitted']);self.assertTrue(g['bad_CA_stereocentres'])

    def test_new_unbonded_heavy_overlap_is_explicit(self):
        water=next(a['index'] for a in self.mapping['atoms'] if a['kind']=='retained_crystal_water' and a['element']=='O')
        target=next(a['index'] for a in self.mapping['atoms'] if a['kind']=='protein_source' and a['element']=='C')
        q=self.q0.copy();q[water]=q[target]
        g=scoring.geometry_check(self.q0,q,self.mapping)
        self.assertFalse(g['admitted'])
        self.assertIn(sorted([water,target]),g['new_severe_heavy_clashes'])

    def test_selection_and_units_on_archived_real_core_values_only(self):
        # This is numeric regression using completed prior energies. It does not
        # claim those local cores are whole-chain candidates or create outputs.
        old=read_json(ROOT/'diagnostics/lanm_series_followup_20260923/DY_TRANSFER_RESULT.json')
        actual=[r for r in old['four_Hans_cells_per_site'] if r['site']=='EF1' and r['metal']=='La']
        bysource={r['source']:r for r in actual}
        origin=bysource['8DQ2_La_conditioned']['composite_kcal_mol']
        proposal=bysource['8FNR_Dy_conditioned']['composite_kcal_mol']
        chosen=scoring.choose({'origin':origin,'archived_other_source':proposal})
        self.assertEqual(chosen['selected'],'archived_other_source')
        self.assertAlmostEqual(chosen['accommodation_work_kcal_mol'],proposal-origin,places=10)
        self.assertEqual(scoring.choose({'origin':origin,'exact_duplicate':origin})['selected'],'origin')
        for row in actual:
            composite=row['MACE_eV']*scoring.EV_TO_KCAL+(row['GFN_ALPB_Ha']-row['GFN_vacuum_Ha'])*HA_TO_KCAL
            self.assertAlmostEqual(composite,row['composite_kcal_mol'],places=8)

    def test_required_real_origin_pair_cannot_lose_a_cell_or_state(self):
        # Staging-only metadata from the actual paired XYZ files: no energies,
        # forces or successful molecular statuses are manufactured.
        pins=[record(p) for p,_,_ in self.states]
        cells=[{'candidate':'origin','metal':z,'xyz':self.s['endpoints'][z]['xyz'],
                'charge':self.s['charge'],
                'physical_multiplicity':self.s['endpoints'][z]['physical_multiplicity'],
                'status':'not_executed','energy_eV':None} for z in scoring.METALS]
        row={'state_id':self.s['state_id'],'state':record(self.pin),'cells':cells,
             'proposals':{},'candidate_aliases':{},'status':'not_executed'}
        staged={'selected_states':[self.s['state_id']],'rows':[row],'capability_only':True}
        scoring.validate_mace_rows({'states':pins},staged)
        absent=copy.deepcopy(staged);absent['rows'][0]['cells'].pop()
        with self.assertRaisesRegex(InvalidArtifact,'pool'):
            scoring.validate_mace_rows({'states':pins},absent)
        absent=copy.deepcopy(staged);absent['selected_states'].append(self.states[1][1]['state_id'])
        with self.assertRaisesRegex(InvalidArtifact,'denominator'):
            scoring.validate_mace_rows({'states':pins},absent)
        absent=copy.deepcopy(staged);absent['capability_only']=False
        with self.assertRaisesRegex(InvalidArtifact,'search'):
            scoring.validate_mace_rows({'states':pins},absent)

if __name__=='__main__':unittest.main()
