import unittest

import os

from treeomics.__main__ import main


class IntegrationTest(unittest.TestCase):

    def setUp(self):
        self.input_dir = os.path.join('input', 'Makohon2017')
        self.fp_pam01_mr = os.path.join(self.input_dir, 'Pam01_1-6_mutant_reads.txt')
        self.fp_pam01_cov = os.path.join(self.input_dir, 'Pam01_1-6_phredcoverage.txt')
        print(os.getcwd())

    def tearDown(self):
        pass

    def test_example1(self):

        example1_args = ['-r', self.fp_pam01_mr,
                         '-s', self.fp_pam01_cov,
                         '-n', 'Pam01N3',
                         '-e', '0.005']
        main(example1_args)

        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
