__author__ = 'jreiter'


OUTPUT_FOLDER = 'output'        # path to output folder for all files

# output file prefix
incomp_mps_plot_prefix = 'incomp_mps_table_'    # mutation table plot of incompatible mutation patterns
artifacts_plot_prefix = 'artifacts_table_'      # table plot of putative artifacts

# Sample sequencing data filtering settings
SAMPLE_COVERAGE_THRESHOLD = 0  # samples with a lower median coverage are discarded
MAF_THRESHOLD = 0.0             # samples with a lower median mutant allele frequency are discarded

# Sequencing data analysis settings
MIN_ABSENT_COVERAGE = 100        # minimum coverage to call a variant absent (if null hypothesis was accepted)
FPR = 0.0023       # calculated false-positive rate for targeted sequencing
#FPR = 0.01          # assumed false-positive rate in Illumina WGS data
#FPR = 0.005          # assumed false-positive rate in Gundem et al. 2015
FDR = 0.01          # targeted false-discovery rate
BI_E = 0.005         # sequencing error rate for bayesian inference
BI_C0 = 0.5         # prior mixture parameter of delta function and uniform distribution for bayesian inference

# minimum reliability score of a incompatible mutation pattern that a subclone of different origin is considered
MIN_SC_SCORE = 3.0

# Minimal reliability score of the MPs in the mutation pattern overview graph if there are many MPs
MIN_MP_SCORE = 1.0
MAX_NO_MPS = 20     # apply MIN_MP_WEIGHT if there are more than this number of MPs in the data
