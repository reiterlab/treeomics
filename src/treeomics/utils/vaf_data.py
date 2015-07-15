__author__ = 'jreiter'

import logging
from collections import defaultdict
from utils.statistics import calculate_present_pvalue

# get logger for application
logger = logging.getLogger('treeomics')


def calculate_p_values(mut_reads, phred_coverage, false_positive_rate):
    """
    Calculate median coverage and median MAF of all confirmed present mutations in each sample
    Calculate p-values for confidence in presence
    :param mut_reads: mut_key dictionary with the number of mutant reads in each sample
    :param phred_coverage: mut_key dictionary with the phred coverage in each sample
    :param false_positive_rate: false positive read of the used sequencing technology
    :return: dictionary of present p-values per sample per variant key
    """

    present_p_values = defaultdict(dict)

    # run through each distinct variant and determine its present and absent p-value
    for mut_key in mut_reads.keys():
        for sample_name in mut_reads[mut_key].keys():

            if mut_reads[mut_key][sample_name] >= 0 and phred_coverage[mut_key][sample_name] >= 0:
                present_p_values[sample_name][mut_key] = calculate_present_pvalue(
                    mut_reads[mut_key][sample_name], phred_coverage[mut_key][sample_name], false_positive_rate)

    return present_p_values


def filter_mutation_functions(data, mut_functions, filtered_functions):
    """
    Filter all mutations which have the desired function (given set of filtered_functions)
    """

    logger.debug('Filter mutations which have these functions: {}'.format(str(filtered_functions)))

    filtered_muts = set()
    for mut_idx in range(len(data)):

        # remove remaining mutations
        if mut_functions[mut_idx] not in filtered_functions:
            for sa_idx in range(len(data[mut_idx])):
                data[mut_idx][sa_idx] = 0
            filtered_muts.add(mut_idx)

    logger.debug('These mutations ({}) have been filtered out {} with functions: {}'.format(len(filtered_muts),
                 str(filtered_muts), ','.join(set(mut_functions[mut] for mut in filtered_muts))))

    return data


def get_present_mutations(data, include_unknowns=True):
    """
    Return list of mutations which are present in a least one sample
    :param data: 2d-array of the filtered mutant allele frequencies; MAF < 0 correspond to unknown positions
    :param include_unknowns: denotes if unknown positions count as present variants
    :return: list of indices of present variants
    """

    # indices list of the present mutations
    present_mutations = []
    for mut_idx in range(len(data)):

        for maf in data[mut_idx]:

            if maf > 0 or (include_unknowns and maf < 0):
                present_mutations.append(mut_idx)
                break

    return present_mutations


def get_shared_mutations(data, include_unknowns=True):
    """
    Return list of mutations which are shared among at least 2 but not all samples
    :param data: 2d-array of the filtered mutant allele frequencies; MAF < 0 correspond to unknown positions
    :param include_unknowns: denotes if unknown positions count as present variants
    :return: list of indices of shared variants
    """

    # indices list of the shared mutations
    shared_mutations = []
    for mut_idx in range(len(data)):

        p = 0
        a = False
        for maf in data[mut_idx]:

            if maf > 0:     # mutation is present
                p += 1
            elif maf == 0:  # mutation is absent
                a = True
            elif maf < 0 and include_unknowns:  # mutations is unknown
                p += 1

            if p > 1 and a:
                shared_mutations.append(mut_idx)
                break

    return shared_mutations
