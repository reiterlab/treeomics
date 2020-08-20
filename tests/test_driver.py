"""Unittests around driver class and files"""

import unittest
import os

from treeomics.utils.driver import Driver, read_driver_file

__author__ = 'Johannes REITER'
__date__ = 'May 15, 2018'


class TestDriver(unittest.TestCase):

    def setUp(self):
        self.driver_fp = os.path.join('tests', 'drivers_stub.csv')
        print(os.getcwd())

    def tearDown(self):
        pass

    def test_read_driver_file(self):

        drs = read_driver_file(self.driver_fp, cancer_type='colon')
        self.assertEqual(len(drs), 1)

        drs = read_driver_file(self.driver_fp, cancer_type='pancreatic')
        self.assertEqual(len(drs), 0)

        drs = read_driver_file(self.driver_fp)
        self.assertEqual(len(drs), 6)

        self.assertIn('KRAS', drs)
        self.assertIn('MBD6', drs)

        self.assertEqual(len(drs['KRAS'].sources), 3)
        self.assertEqual(len(drs['MBD6'].sources), 0)

        self.assertEqual(Driver.MaxSourceSupport, 3)

    def test_init(self):

        dr = Driver('KRAS')
        self.assertEqual(dr.gene_name, 'KRAS')
        self.assertEqual(len(dr.sources), 0)

        dr = Driver('KRAS', sources=['MethodA', 'MethodB', 'MethodC'])
        self.assertEqual(dr.gene_name, 'KRAS')
        self.assertEqual(len(dr.sources), 3)
        self.assertEqual(Driver.MaxSourceSupport, 3)


def main():
    unittest.main()

if __name__ == "__main__":
    main()