"""Filter variants based on their location or mutation effect """
import logging

try:    # check if varcode and pyensembl is available (necessary for Windows)
    from utils.mutation_effects import get_top_effect_name
    VARCODE = True

except ImportError:
    # mutation effect prediction will not be performed since VarCode is not avilable
    VARCODE = False
    get_top_effect_name = None

__author__ = 'Johannes REITER'
__date__ = 'March 30, 2017'


# get logger for application
logger = logging.getLogger('treeomics')


def is_intronic(variant):
    if VARCODE:
        mut_effect = get_top_effect_name(variant)
        if mut_effect == 'Intronic':
            return True
        else:
            return False
    else:
        return None


def is_intergenic(variant):
    if VARCODE:
        mut_effect = get_top_effect_name(variant)
        if mut_effect == 'Intergenic':
            return True
        else:
            return False
    else:
        return None


def is_incompletetranscript(variant):
    if VARCODE:
        mut_effect = get_top_effect_name(variant)
        if mut_effect == 'IncompleteTranscript':
            return True
        else:
            return False
    else:
        return None
