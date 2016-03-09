__author__ = 'jreiter'


# Treeomics version number
VERSION = (1, 2, 1)

# DEFAULT PARAMETER VALUES
# maximal number variants to show in a mutation table plot
MAX_MUTS_TABLE_PLOT = 1000

# bayesian inference model
# assumed posterior p0 probability when no sequencing data was reported in some sample
NO_DATA_P0 = 0.9

# presence probability of a variant for calculating reliability score
# is upper bounded because the same variant could have been independently acquired twice
# number of detected coding mutations divided by the size of the exome (1.5% of the genome)
MAX_PRE_PROB = 1.0 - (100.0 / 45e06)

# absence probability of a variant for calculating reliability score
# should be upper bounded because the variant could have been lost by LOH
# for most sequencing depth this lower bound is irrelevant
MAX_ABS_PROB = 1.0 - 1e-04

# default prior when variant is believed to be present
PSEUDO_ALPHA = 1.0  # alpha parameter for the beta prior
PSEUDO_BETA = 1.5     # beta parameter for the beta prior

# Necessary constants: do NOT change
POS_UNKNOWN = -2
NEG_UNKNOWN = -1
