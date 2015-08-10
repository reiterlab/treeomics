__author__ = 'jreiter'

# Sample sequencing data filtering settings
SAMPLE_COVERAGE_THRESHOLD = 0   # samples with a lower median coverage are discarded
MAF_THRESHOLD = 0.0             # samples with a lower median mutant allele frequency are discarded

# Bayesian sequencing data analysis settings                ######################
BI_E = 0.005            # sequencing error rate for bayesian inference
BI_C0 = 0.5             # prior mixture parameter of delta function and uniform distribution for bayesian inference

# minimum likelihood that at least one variant has an incompatible mp such that this mp is considered as a subclone
MIN_MP_LH = 1.0


# Conventional binary sequencing data analysis settings     ######################
FPR = 0.005                 # assumed false-positive rate for targeted sequencing (Gundem et al. 2015)
# FPR = 0.01                # assumed false-positive rate in Illumina WGS data
FDR = 0.05                  # targeted false-discovery rate
MIN_ABSENT_COVERAGE = 100   # minimum coverage to call a variant absent (if null hypothesis was accepted)


# ########################## OUTPUT CONFIGURATIONS #################################
OUTPUT_FOLDER = 'output'    # path to output folder for all files

# output file prefix
incomp_mps_plot_prefix = 'incomp_mps_table_'    # mutation table plot of incompatible mutation patterns
artifacts_plot_prefix = 'artifacts_table_'      # table plot of putative artifacts

# Minimal reliability score for the depicted of the MPs in the mutation pattern overview graph
MIN_MP_SCORE = 1.0  # show only mutation patterns in the circos plots with a greater reliability score
MAX_NO_MPS = 25     # show at most this number of mutation patterns in the circos plots
