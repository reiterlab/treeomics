""" Prediction mutation effects by using VarCode """
import logging
from itertools import chain
import sys
import traceback

from treeomics.varcode.effects import predict_variant_effect_on_transcript
from treeomics.varcode.effects.effect_classes import *        # https://github.com/hammerlab/varcode

__author__ = 'Johannes REITER'
__date__ = 'March 2, 2017'


# get logger
logger = logging.getLogger(__name__)


class EffectTypes:

    types = None

    @staticmethod
    def types():
        if EffectTypes.types is None:
            mutation_effects = set()
            for i in chain(_get_all_subclasses(MutationEffect), _get_all_subclasses(SpliceSite)):
                mutation_effects.add(i.__name__)

            EffectTypes.types = mutation_effects

        return EffectTypes.types


def get_top_effect_name(variant):
    """
    Return the name of the most likely effect of this mutation
    :param variant: varcode variant, see https://github.com/hammerlab/varcode
    :return: effect name
    """

    try:
        return type(variant.effects().top_priority_effect()).__name__
        # return get_variant_effect_longest_transcript(variant)

    except BaseException as e:
        logger.error('Error: {}'.format(str(e)))
        if logger.getEffectiveLevel() == logging.DEBUG:
            traceback.print_tb(sys.exc_info()[2])
        logger.warning('Mutation effect for variant {} could not be inferred.'.format(variant.short_description))
        return 'unknown'


def is_top_substitution(variant):
    """
    Returns True if the top predicted effect name of the variant is a substitution, e.g. Missense
    :param variant: varcode variant, see https://github.com/hammerlab/varcode
    :return: True or False
    """

    if get_top_effect_name(variant) == 'Substitution':
        return True
    else:
        return False


def get_variant_effect_longest_transcript(variant):
    """
    Predict variant effect based on longest transcript and return string
    :param variant: varcode variant, see https://github.com/hammerlab/varcode
    :return: effect name
    """

    try:
        ltc = _get_longest_transcript(variant.transcripts)
        if ltc is None:
            ve = Intergenic(variant).__class__.__name__
        else:
            ve = predict_variant_effect_on_transcript(variant, ltc).__class__.__name__

        return ve

    except:
        logger.warning('Mutation effect for variant {} could not be inferred.'.format(variant.short_description))
        return 'unknown'


def _get_longest_transcript(transcripts):
    """
    Returns the longest transcript from a list of given transcripts
    """
    longest_tc = None
    max_length = 0
    for tc in transcripts:
        if len(tc) > max_length:
            longest_tc = tc
            max_length = len(tc)

    return longest_tc


def _get_all_subclasses(cls):
    """
    Generator of all a class's subclasses
    """
    try:
        for subclass1 in cls.__subclasses__():
            yield subclass1
            for subclass2 in _get_all_subclasses(subclass1):
                yield subclass2
    except TypeError:
        return
