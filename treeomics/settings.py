#!/usr/bin/python
"""Basic configurations to run Treeomics"""
__author__ = 'Johannes REITER'

# ######################### DEFAULT PARAMETER VALUES #############################

# Explore the solution space to assess the support of the inferred branches
POOL_SIZE = 1000            # number of best solutions explored by ILP solver
NO_PLOTTED_SOLUTIONS = 5   # number of best ranked solution trees that will be plotted (cannot be larger than the pool)

# number of bootstrapping samples
NO_BOOTSTRAP_SAMPLES = 0

# subclone detection enabled
SUBCLONE_DETECTION = False

# Time limit for MILP solver in seconds
# For data sets with less than 10-15 samples, a time limit is typically not required (TIME_LIMIT = None)
TIME_LIMIT = None

# For efficiency, the number of explored mutation pattern (MAX_NO_MPS) per variant can be limited
# Note that if MAX_NO_MPS is not None, the optimal solution is no longer guaranteed
MAX_NO_MPS = None

# assign variants classified as absent by the Bayesian inference model to their next most likely mutation pattern
# i.e. also show variants in tiny subclones in the inferred phylogeny
SHOW_BI_ABSENT_MUTS = False

# #################### BAYESIAN sequencing data analysis settings ######################
BI_E = 0.01            # sequencing error rate for bayesian inference
BI_C0 = 0.5            # prior mixture parameter of delta function and uniform distribution for bayesian inference
MAX_ABSENT_VAF = 0.05  # maximal absent VAF before considering estimated purity

LOH_FREQUENCY = 0.0   # probability that a somatic SNV along a lineage is lost due to loss of heterozygosity

# DEPRECATED FROM VERSION 1.7.0 ONWARD
# Conventional binary sequencing data analysis settings     ######################
FPR = 0.005                 # assumed false-positive rate for targeted sequencing (Gundem et al. 2015)
# FPR = 0.01                # assumed false-positive rate in Illumina WGS data
FDR = 0.05                  # targeted false-discovery rate
MIN_ABSENT_COVERAGE = 100   # minimum coverage to call a variant powered absent (if null hypothesis was accepted)


# ########################### OPTIONAL FILTERING ##################################

# Sample sequencing data filtering settings
SAMPLE_COVERAGE_THRESHOLD = 0   # samples with a lower median coverage are discarded
MAF_THRESHOLD = 0.0             # samples with a lower median mutant allele frequency are discarded

# Minimum VAF of a variant in at least one of the provided samples with a minimum number of variant reads (DEFAULT: 0.0)
MIN_VAF = 0.0
MIN_VAR_READS = 0

# minimum coverage of a variant across all samples
MIN_VAR_COV = 0

# Variants are excluded if they reach the below number of mutant reads and the variant allele frequency in normal sample
MUT_READS_NORMAL_TH = 3     # DEFAULT: 3
VAF_NORMAL_TH = 0.02        # DEFAULT: 0.02

# ########################## INPUT CONFIGURATIONS #################################
# path to a list of all cancer gene census from COSMIC in a comma separated value file
# if not available, provide None
# hg18 (GRCh36) is not available, hg19 corresponds to GRCh37, hg38 corresponds to GRCh38
CGC_PATH = 'input/cancer_gene_census_grch37_v84.csv'
# CGC_PATH = 'input/cancer_gene_census_grch38_v84.csv'
# CGC_PATH = None

# additional gene names list to be highlighted in mutation table and inferred phylogeny independent of the COSMIC list
# drivers specific to pancreatic cancer
# DRIVER_PATH = 'input/PDAC_drivers.csv'

# TCGA consensus driver gene list from Bailey et al, Cell 2018
DRIVER_PATH = 'input/BaileyDing2018_driverconsensus.csv'
# union of driver lists inferred by 20/20plus, TUSON and MutsigCV in Tokheim et al, PNAS, 2016
# DRIVER_PATH = 'input/Tokheim_drivers_union.csv'
# DRIVER_PATH = 'input/mCRC/mCRC_drivers.csv'

# only necessary to predict the mutation type (silent, intronic, etc)
# supporting grch36 (for hg18), grch37 (for hg19) and grch38 (for hg38)
REF_GENOME = 'grch37'
# REF_GENOME = None

# path to CSV file with variants common in ExAC or other panel of normals
# These variants are likely sequencing artifacts and therefore excluded by Treeomics
# ExAC_FILE = '../input/ExAC/ExAC_sites_0_0001.csv'
COMMON_VARS_FILE = 'input/ExAC/ExAC_sites_0_001.csv'
# COMMON_VARS_FILE = None

# ########################## OUTPUT CONFIGURATIONS #################################
OUTPUT_FOLDER = 'output'    # path to output folder for all files

# depending on various display settings, one has to play with the zoom to get a PDF report with a decent font size
ZOOM = 1.0

# output file prefix
incomp_mps_plot_suffix = '_incomp_mps_table'    # mutation table plot of incompatible mutation patterns
artifacts_plot_suffix = '_artifacts_table'      # table plot of putative artifacts

# Minimal reliability score for the depicted of the MPs in the mutation pattern overview graph
MIN_MP_SCORE = 0.01  # show only mutation patterns in the circos plots with a greater reliability score
CIRCOS_MAX_NO_MPS = 30     # show at most this number of mutation patterns in the circos plots
