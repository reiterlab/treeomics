#!/usr/bin/python
"""Basic configurations to run Treeomics"""
__author__ = 'Johannes REITER'

# ######################### DEFAULT PARAMETER VALUES #############################

# Explore the solution space to assess the support of the inferred branches
POOL_SIZE = 1000        # number of best solutions explored by ILP solver

# number of bootstrapping samples
NO_BOOTSTRAP_SAMPLES = 0

# subclone detection enabled
SUBCLONE_DETECTION = False

# Time limit for MILP solver in seconds
# For data sets with less than 10-15 samples, a time limit is typically not required (TIME_LIMIT = None)
TIME_LIMIT = None

# For efficiency, the number of explored mutation pattern per variant can be limited
# Note that if MAX_NO_MPS is not None, the optimal solution is no longer guaranteed
MAX_NO_MPS = None

# #################### BAYESIAN sequencing data analysis settings ######################
BI_E = 0.01            # sequencing error rate for bayesian inference
BI_C0 = 0.5            # prior mixture parameter of delta function and uniform distribution for bayesian inference
MAX_ABSENT_VAF = 0.05  # maximal absent VAF before considering estimated purity

LOH_FREQUENCY = 0.0   # probability that a somatic SNV along a lineage is lost due to loss of heterozygosity


# Conventional binary sequencing data analysis settings     ######################
FPR = 0.005                 # assumed false-positive rate for targeted sequencing (Gundem et al. 2015)
# FPR = 0.01                # assumed false-positive rate in Illumina WGS data
FDR = 0.05                  # targeted false-discovery rate
MIN_ABSENT_COVERAGE = 100   # minimum coverage to call a variant powered absent (if null hypothesis was accepted)


# ########################### OPTIONAL FILTERING ##################################

# Sample sequencing data filtering settings
SAMPLE_COVERAGE_THRESHOLD = 0   # samples with a lower median coverage are discarded
MAF_THRESHOLD = 0.0             # samples with a lower median mutant allele frequency are discarded

# Minimum VAF of a variant in at least one of the provided samples
MIN_VAF = 0.0
# Minimum number of reads reporting a variant in at least one of the provided samples
MIN_VAR_READS = 0

# ########################## INPUT CONFIGURATIONS #################################
# path to a list of all cancer gene census from COSMIC in a comma separated value file
# if not available, provide None
# hg18 is not available, hg19 corresponds to GRCh37, hg38 corresponds to GRCh38
# CGC_PATH = '../input/cancer_gene_census_grch37_v75.csv'
CGC_PATH = None


# ########################## OUTPUT CONFIGURATIONS #################################
OUTPUT_FOLDER = 'output'    # path to output folder for all files

# additional gene names list to be highlighted in mutation table independent of the COSMIC list
# Drivers specific to pancreatic cancer
DRIVERS = ('ACVR1B', 'ARID1A', 'ARID1B', 'ATM', 'BRCA1', 'BRCA2', 'BRAF', 'CDK6', 'CDKN2A', 'CHD4', 'ERBB2', 'FGFR2', 'GATA6',
           'KDM6A', 'KMT2D', 'KRAS', 'MAP2K4', 'MET', 'MLH1',  'MSH2', 'MLL2', 'MYC', 'NOV', 'PALB2', 'PBRM1', 'PIK3C2B',
           'PIK3CA', 'PIK3R3', 'PREX2', 'RNF43', 'RPA1', 'ROBO1', 'ROBO2', 'SF3B1', 'SLIT2', 'SMAD4', 'SMARCA2',
           'SMARCA4', 'SOX9', 'STK11', 'TGFBR2', 'TP53', 'U2AF1')
# DRIVERS = set()

# output file prefix
incomp_mps_plot_suffix = '_incomp_mps_table'    # mutation table plot of incompatible mutation patterns
artifacts_plot_suffix = '_artifacts_table'      # table plot of putative artifacts

# Minimal reliability score for the depicted of the MPs in the mutation pattern overview graph
MIN_MP_SCORE = 0.01  # show only mutation patterns in the circos plots with a greater reliability score
CIRCOS_MAX_NO_MPS = 30     # show at most this number of mutation patterns in the circos plots
