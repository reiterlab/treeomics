
import pytest
import os

from treeomics.__main__ import main

@pytest.fixture
def input_dir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'input', 'Makohon2017')

@pytest.fixture
def fp_pam01_mr(input_dir):
    return os.path.join(input_dir, 'Pam01_1-6_mutant_reads.txt')

@pytest.fixture
def fp_pam01_cov(input_dir):
    return os.path.join(input_dir, 'Pam01_1-6_phredcoverage.txt')


def test_example1(fp_pam01_mr, fp_pam01_cov):

    cplex = pytest.importorskip('cplex')

    example1_args = ['-r', fp_pam01_mr,
                     '-s', fp_pam01_cov,
                     '-n', 'Pam01N3',
                     '-e', '0.005']
    return_value = main(example1_args)

    assert return_value is True
