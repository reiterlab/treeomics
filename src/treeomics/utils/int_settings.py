__author__ = 'jreiter'


# Treeomics version number
VERSION = (1, 0, 0)

# DEFAULT PARAMETER VALUES
# maximal number variants to show in a mutation table plot
MAX_MUTS_TABLE_PLOT = 1000

# bayesian inference model
# assumed posterior p0 probability when no sequencing data was reported in some sample
NO_DATA_P0 = 0.99
PSEUDO_ALPHA = 1  # alpha parameter for the beta prior
PSEUDO_BETA = 3   # beta parameter for the beta prior

# Necessary constants: do NOT change
POS_UNKNOWN = -2
NEG_UNKNOWN = -1
