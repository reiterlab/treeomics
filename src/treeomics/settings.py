#!/usr/bin/python
"""Basic configurations to run Treeomics"""
__author__ = 'Johannes REITER'

# ######################### DEFAULT PARAMETER VALUES #############################
# Sample sequencing data filtering settings
SAMPLE_COVERAGE_THRESHOLD = 0   # samples with a lower median coverage are discarded
MAF_THRESHOLD = 0.0             # samples with a lower median mutant allele frequency are discarded

# Bayesian sequencing data analysis settings                ######################
BI_E = 0.005            # sequencing error rate for bayesian inference
BI_C0 = 0.5             # prior mixture parameter of delta function and uniform distribution for bayesian inference

# minimum likelihood that at least one variant has an incompatible mp
# such that this mp is considered as a subclone
SUBCLONE_DETECTION = False

# Conventional binary sequencing data analysis settings     ######################
FPR = 0.005                 # assumed false-positive rate for targeted sequencing (Gundem et al. 2015)
# FPR = 0.01                # assumed false-positive rate in Illumina WGS data
FDR = 0.05                  # targeted false-discovery rate
MIN_ABSENT_COVERAGE = 100   # minimum coverage to call a variant powered absent (if null hypothesis was accepted)

# #### OPTIONAL FILTERING ####
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
# DRIVERS = ('ACVR1B', 'ARID1A', 'ARID1B', 'ATM', 'BRCA1', 'BRCA2', 'BRAF', 'CDK6', 'CDKN2A', 'ERBB2', 'FGFR2', 'GATA6',
#            'KDM6A', 'KRAS', 'MAP2K4', 'MET', 'MLH1',  'MSH2', 'MLL2', 'MYC', 'NOV', 'PALB2', 'PBRM1', 'PIK3C2B',
#            'PIK3CA', 'PIK3R3', 'PREX2', 'RNF43', 'RPA1', 'ROBO1', 'ROBO2', 'SF3B1', 'SLIT2', 'SMAD4', 'SMARCA2',
#            'SMARCA4', 'SOX9', 'STK11', 'TGFBR2', 'TP53')
DRIVERS = set()

# output file prefix
incomp_mps_plot_prefix = 'incomp_mps_table_'    # mutation table plot of incompatible mutation patterns
artifacts_plot_prefix = 'artifacts_table_'      # table plot of putative artifacts

# Minimal reliability score for the depicted of the MPs in the mutation pattern overview graph
MIN_MP_SCORE = 0.005  # show only mutation patterns in the circos plots with a greater reliability score
CIRCOS_MAX_NO_MPS = 30     # show at most this number of mutation patterns in the circos plots
