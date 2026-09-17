"""Sign/accounting regression on actual saved MACE outputs."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from mace_short_range import paired_components,accepted_source
RESULT=ROOT/'diagnostics/mace_short_range_20260916/result.json'


@unittest.skipUnless(RESULT.exists(),'actual short-component assessment required')
class ShortRangeTests(unittest.TestCase):
    def test_exact_saved_values_and_failed_alpha_order_retained(self):
        r=read_json(RESULT)
        for label in ('medium','large'):
            c,m=accepted_source(verify(r['saved_sources'][label]))
            for case,row in r['rows'].items():
                actual=paired_components(*(c['rows'][f'{case}_{metal}_primary'] for metal in ('La','Ca')))
                self.assertEqual(actual,row['full'][label]);self.assertIsNone(actual['calibrated_class'])
            self.assertEqual([x['expected_order'] for x in r['contrasts'][label]],[True,False,False])
        self.assertFalse(r['primary_development_ordering_screen_pass'])


if __name__=='__main__':unittest.main()
