"""Filter variants based on their location or mutation effect """
import logging
from enum import Enum, auto

try:    # check if varcode and pyensembl is available (necessary for Windows)
    from treeomics.utils.mutation_effects import get_top_effect_name
    VARCODE = True

except ImportError:
    # mutation effect prediction will not be performed since VarCode is not avilable
    VARCODE = False
    get_top_effect_name = None

__author__ = 'Johannes REITER'
__date__ = 'March 30, 2017'


# get logger
logger = logging.getLogger(__name__)


class Filter(Enum):
    PASSED = auto()
    MIN_VAF = auto()
    MIN_VAR_READS = auto()
    MIN_VAR_COV = auto()
    NORMAL = auto()
    INTRONIC = auto()
    INTERGENIC = auto()
    INCOMPL_TRANSCRIPT = auto()
    SNP = auto()

    def __str__(self):
        if self.value == Filter.MIN_VAF.value:
            return 'No significant VAF was reached in at least one sample.'
        elif self.value == Filter.MIN_VAR_READS.value:
            return 'Minimum number of supporting variant reads was not reached in any sample.'
        elif self.value == Filter.MIN_VAR_COV.value:
            return 'Minimum coverage of variant was not reached across all samples.'
        elif self.value == Filter.NORMAL.value:
            return 'Variant was also observed in the normal sample.'
        elif self.value == Filter.INTRONIC.value:
            return 'Intronic variant in WES data.'
        elif self.value == Filter.INTERGENIC.value:
            return 'Intergenic variant in WES data.'
        elif self.value == Filter.INCOMPL_TRANSCRIPT.value:
            return 'Incomplete transcript annotation (likely intron) in WES data.'
        elif self.value == Filter.SNP.value:
            return 'Variant is a common SNP or sequencing artifact.'
        elif self.value == Filter.PASSED.value:
            return 'Variant passed all filters.'
        else:
            raise ValueError(f'String function not yet implemented for {self.name}')


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
