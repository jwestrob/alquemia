"""Algebra/interface tests on archived 1H4I charges, NOT molecular qualification.

The field-squared probe below is deliberately a mathematical derivative probe,
not a surrogate energy or invented scientific observation. No executable model
or physical boundary force mapping is qualified by these tests.
"""
import copy
import hashlib
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from metal_environment_fields import (BOHR_ANGSTROM, potential_field,
    chain_rule_forces, charged_gauge_shift, evaluation_fingerprint)

FIXTURE = ROOT / 'workspaces/affordable_challenger_20260915/solver_recovery/1h4i_qm33_La/primary/calculation'
PINS = {'core': '9f82122e64d01d239d2cbb6cc032671fc546cfb17ed3967248a839084c67151e',
        'environment': '6a2f369619b143984c7bef2bd483f8718b4d52acf186c6b2e3ad234a54eeda89'}


def read_real(name):
    path = FIXTURE / (name + '.pqr')
    if not path.exists():
        raise unittest.SkipTest('Archived real 1H4I fixture unavailable: ' + str(path))
    if hashlib.sha256(path.read_bytes()).hexdigest() != PINS[name]:
        raise AssertionError('Real fixture content hash changed: ' + str(path))
    rows = [list(map(float, line.split()[-5:])) for line in path.read_text().splitlines()
            if line.startswith('ATOM ')]
    rows = np.asarray(rows)
    return rows[:, :3] / BOHR_ANGSTROM, rows[:, 3]


class FiniteFieldTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x, cls.qx = read_real('core')
        cls.y, cls.qy = read_real('environment')
        # Archived PB PQR retains zero-charge QM cavity centers. They are not
        # MM point charges: exclude precisely zero-charge rows for this test,
        # not charged overlapping atoms. This is not a new boundary recipe.
        nonzero = cls.qy != 0
        cls.y, cls.qy = cls.y[nonzero], cls.qy[nonzero]

    def probe(self, x, y):
        phi, field = potential_field(x, y, self.qy)
        # Pure chain-rule probe, not a physical force-field energy.
        return float(self.qx @ phi + 0.5*np.sum(field*field)), phi, field

    def test_directional_derivatives_and_action_reaction(self):
        _, _, field = self.probe(self.x, self.y)
        fx, fy = chain_rule_forces(self.x, self.y, self.qy, self.qx, field)
        np.testing.assert_allclose(fx.sum(0) + fy.sum(0), 0, atol=2e-13)
        h = 1e-4  # bohr; declared before residual inspection
        # Actual nearest-environment direction; both QM and MM centers tested.
        j = np.argmin(np.sum((self.y-self.x[0])**2, axis=1))
        direction = self.y[j]-self.x[0]
        direction /= np.linalg.norm(direction)
        for which, idx, force in [('qm', 0, fx[0]), ('mm', j, fy[j])]:
            xp, xm, yp, ym = (v.copy() for v in (self.x,self.x,self.y,self.y))
            (xp if which == 'qm' else yp)[idx] += h*direction
            (xm if which == 'qm' else ym)[idx] -= h*direction
            numerical = -(self.probe(xp,yp)[0]-self.probe(xm,ym)[0])/(2*h)
            self.assertLess(abs(numerical - force@direction), 2e-8)

    def test_joint_rigid_transform_and_repeatability(self):
        energy, phi, field = self.probe(self.x, self.y)
        fx, fy = chain_rule_forces(self.x,self.y,self.qy,self.qx,field)
        angle = 0.71
        c,s = np.cos(angle), np.sin(angle)
        rotation = np.array([[c,-s,0],[s,c,0],[0,0,1]])
        x, y = self.x@rotation.T+3.25, self.y@rotation.T+3.25
        e2,p2,f2 = self.probe(x,y)
        gx,gy = chain_rule_forces(x,y,self.qy,self.qx,f2)
        self.assertAlmostEqual(energy,e2,places=11)
        np.testing.assert_allclose(p2,phi,atol=1e-12)
        np.testing.assert_allclose(f2,field@rotation.T,atol=1e-12)
        np.testing.assert_allclose(gx,fx@rotation.T,atol=1e-12)
        np.testing.assert_allclose(gy,fy@rotation.T,atol=1e-12)
        self.assertEqual(energy,self.probe(self.x,self.y)[0])

    def test_gauge_shift_charge_accounting(self):
        p,f = potential_field(self.x,self.y,self.qy)
        shifted,fs = potential_field(self.x,self.y,self.qy,gauge=.007)
        np.testing.assert_array_equal(f,fs)
        measured = self.qx @ (shifted-p)
        self.assertAlmostEqual(measured,charged_gauge_shift(self.qx.sum(),.007),places=13)

    def test_chunk_equivalence_and_corrupted_overlap(self):
        p,f = potential_field(self.x,self.y,self.qy)
        p2,f2 = potential_field(self.x,self.y,self.qy,chunk_size=71)
        np.testing.assert_allclose(p,p2,atol=1e-13)
        np.testing.assert_allclose(f,f2,atol=1e-13)
        broken = self.y.copy(); broken[0] = self.x[0]
        with self.assertRaisesRegex(ValueError,'Overlapping'):
            potential_field(self.x,broken,self.qy)

    def test_fingerprint_invalidates_every_declared_input(self):
        scopes = dict(state={'charge':-1,'multiplicity':1,'element':'La'},
            geometry={'core_pqr_sha256':PINS['core']},
            environment={'pqr_sha256':PINS['environment'],'gauge':0.,'field_override':None},
            model={'implementation':'algebra_test_only'},protocol={'units':'atomic'})
        first = evaluation_fingerprint(**scopes)
        self.assertEqual(first,evaluation_fingerprint(**copy.deepcopy(scopes)))
        for scope in scopes:
            changed = copy.deepcopy(scopes)
            changed[scope]['changed_setting'] = True
            self.assertNotEqual(first,evaluation_fingerprint(**changed))
        for key,value in [('charge',0),('multiplicity',3)]:
            changed = copy.deepcopy(scopes); changed['state'][key]=value
            self.assertNotEqual(first,evaluation_fingerprint(**changed))
        with self.assertRaises(ValueError):
            evaluation_fingerprint(**dict(scopes,model={}))


if __name__ == '__main__':
    unittest.main()
