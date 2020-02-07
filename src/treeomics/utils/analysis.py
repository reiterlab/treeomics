#!/usr/bin/python
"""Generates analysis file providing an overview of the results"""

import logging
import csv
import math
from collections import defaultdict
from utils.int_settings import NEG_UNKNOWN, POS_UNKNOWN
import utils.int_settings as def_sets
import numpy as np
from utils.similarity_analysis import calculate_genetic_similarity, calculate_bi_genetic_similarity
from utils.data_tables import write_posterior_table, write_vep_input, write_cravat_input
from phylogeny.simple_phylogeny import SimplePhylogeny
from phylogeny.max_lh_phylogeny import MaxLHPhylogeny


__author__ = 'Johannes REITER'

# get logger for application
logger = logging.getLogger('treeomics')


def analyze_data(patient, post_table_filepath=None, vep_filepath=None, cravat_filepath=None):
    """
    Process data and write posterior table
    Possibility of apply various filters
    :param patient: data structure around sequencing data of a subject
    :param post_table_filepath: file path for generating a table with posterior probabilities
    :param fp_jsc: output path to file where Jaccard similarity matrix will be stored
    :param fp_gendist: output path to file where genetic distance matrix will be stored
    :param vep_filepath: output file path to write all substitution variants in a format acceptable to Ensembl VEP
    :param cravat_filepath: output file path to write all substitution variants in a format acceptable to CHASM/CRAVAT
    """

    # Possibility to implement filters!
    # Look at the average/median frequency of founder mutations
    # Filter mutations out with less than half of the average founder mutation frequency
    # data_utils.remove_contradicting_mutations(patient.data)

    # determine in which samples each mutation is present (positives)
    for mut in range(len(patient.data)):
        for sa_idx, maf in enumerate(patient.data[mut]):
            if 0 < maf:
                patient.samples[sa_idx].add(mut)
                patient.mutations[mut].add(sa_idx)

    avg_mutations = 0
    for sa_idx, sa_name in enumerate(patient.sample_names):
        # logger.debug("Mutations present in sample {}: {}, {}".format(
        #     sa_name, len(patient.samples[sa_idx]),
        #     str(patient.samples[sa_idx]) if len(patient.samples[sa_idx]) < 200 else ''))
        avg_mutations += len(patient.samples[sa_idx])

    avg_mutations /= float(len(patient.sample_names))
    logger.info('The average number of mutations per sample in patient {} is {:.1f}.'.format(
        patient.name, avg_mutations))

    for mut_idx in range(len(patient.mut_keys)):
        patient.mutations[mut_idx] = frozenset(patient.mutations[mut_idx])
        # logger.debug("Mutation {} ({}) is present in: {}".format(patient.mut_names[mut_idx],
        #            patient.gene_names[mut_idx],
        #            (','.join(patient.sample_names[sa_idx] for sa_idx in patient.mutations[mut_idx]))))

    if patient.vc_variants is not None:
        if vep_filepath is not None:
            write_vep_input(vep_filepath, patient)
        if cravat_filepath is not None:
            write_cravat_input(cravat_filepath, patient)

    # calculate homogeneity index (Jaccard similarity coefficient)
    patient.gen_dis, patient.sim_coeff, patient.sim_coeff_ex = calculate_genetic_similarity(patient)
    patient.bi_gen_dis, patient.bi_sim_coeff = calculate_bi_genetic_similarity(patient)

    patient.sharing_status = determine_sharing_status(patient)
    assert (len(patient.mut_keys) == len(patient.sharing_status)), 'Error inferring sharing status with BI model'

    # keep list of mutations present in some of the used samples
    patient.present_mutations = get_present_mutations(patient.log_p01)

    if post_table_filepath is not None:
        # write file with posterior probabilities
        write_posterior_table(
            post_table_filepath, patient.sample_names, patient.purities, patient.sample_mafs,
            patient.mut_positions, patient.gene_names, patient.log_p01, patient.betas)

    # calculate mutation numbers per sample based on the Bayesian inference model
    patient.calculate_no_present_vars()


def determine_sharing_status(patient):
    """
    Determine which samples share the various mutations based on the Bayesian inference model
    :param patient: data structure around sequencing data of a subject
    """

    pres_lp = math.log(0.5)  # log probability of 50%
    sharing_status = list()
    for mut_idx, ps in patient.log_p01.items():

        sa_idxs = [sa_idx for sa_idx, (_, p1) in enumerate(ps) if p1 > pres_lp]

        if len(sa_idxs) == len(patient.sample_names):
            patient.founders.add(mut_idx)
            sharing_status.append('Trunk')
        else:
            patient.shared_muts[sum(1 for _, p1 in ps if p1 > pres_lp)].add(mut_idx)

            if len(sa_idxs) > 1:      # shared
                if (all(sa_idx in sa_idxs for sa_idx, sa_name in enumerate(patient.sample_names)
                        if 'M' in sa_name and 'TM' not in sa_name)
                        and any('M' in sa_name and 'TM' not in sa_name for sa_name in patient.sample_names)):
                    sharing_status.append('UntrMetTrunk')

                elif (all(sa_idx in sa_idxs for sa_idx, sa_name in enumerate(patient.sample_names) if 'M' in sa_name)
                      and any('M' in sa_name for sa_name in patient.sample_names)):
                    sharing_status.append('MetTrunk')

                elif any('M' in patient.sample_names[sa_idx] and 'TM' not in patient.sample_names[sa_idx]
                         for sa_idx in sa_idxs):
                    sharing_status.append('UntrMetShared')

                elif any('M' in patient.sample_names[sa_idx] for sa_idx in sa_idxs):
                    sharing_status.append('MetShared')

                elif (all(sa_idx in sa_idxs for sa_idx, sa_name in enumerate(patient.sample_names)
                          if 'PT' in sa_name or 'Primary' in sa_name)
                        and any('PT' in sa_name or 'Primary' in sa_name for sa_name in patient.sample_names)):
                    sharing_status.append('PTTrunk')

                else:
                    sharing_status.append('Shared')

            elif len(sa_idxs) == 1:                                               # private
                sa_name = patient.sample_names[sa_idxs[0]]
                if 'TM' in sa_name:
                    sharing_status.append('TrMetPrivate')
                elif 'M' in sa_name and 'TM' not in sa_name:
                    sharing_status.append('UntrMetPrivate')
                elif 'PT' in sa_name or 'Primary' in sa_name:
                    sharing_status.append('PTPrivate')
                else:
                    sharing_status.append('Private')

            else:
                sharing_status.append('Absent')
                # raise RuntimeError('Mutation {} is not assigned to any samples!'.format(patient.mut_keys[mut_idx]))

    logger.info("{:.1%} ({}/{}) of all distinct mutations are founders.".format(
        float(len(patient.founders)) / len(patient.mutations), len(patient.founders), len(patient.mutations)))
    logger.info('In average {:.1f} ({:.1%}) mutations are unique (private) per sample.'.format(
        float(len(patient.shared_muts[1])) / len(patient.sample_names),
        (float(len(patient.shared_muts[1])) / len(patient.sample_names)) /
        (sum(len(muts) for sa_idx, muts in patient.samples.items()) / len(patient.sample_names))))

    # for shared in range(len(patient.sample_names)-1, -1, -1):
    #     if len(patient.shared_muts[shared]) > 0:
    #         logger.debug('Mutations shared in {} samples ({} of {} = {:.3f}): {}'.format(shared,
    #                      len(patient.shared_muts[shared]), len(patient.mutations),
    #                      float(len(patient.shared_muts[shared])) / len(patient.mutations),
    #                      patient.shared_muts[shared] if len(patient.shared_muts[shared]) < 200 else ''))

    # compute the number of shared and additional mutations among the sample pairs
    # similar to a distance matrix
    patient.common_muts = [[set() for _ in range(patient.n)] for _ in range(patient.n)]
    patient.add_muts = [[0 for _ in range(patient.n)] for _ in range(patient.n)]

    # helper to save calculation time:
    cmuts = [[0 for _ in range(patient.n)] for _ in range(patient.n)]

    for s1 in range(patient.n):
        for s2 in range(patient.n):
            patient.common_muts[s1][s2] = patient.samples[s1].intersection(patient.samples[s2])
            cmuts[s1][s2] = len(patient.common_muts[s1][s2])
            patient.add_muts[s1][s2] = patient.samples[s1].difference(patient.samples[s2])

            # logger.debug('Sample {} has {} mutations in common with sample {} and {} in addition. '.format(
            #    patient.sample_names[s1], cmuts[s1][s2], patient.sample_names[s2], len(patient.add_muts[s1][s2])))

            # logger.debug(
            #     ('Sample {} has most ({}) mutations in common with sample(s) {}'
            #      ' and least ({}) with sample(s) {}. ').format(
            #      patient.sample_names[s1],
            #      max(cms for idx, cms in enumerate(cmuts[s1]) if idx != s1),
            #      str([patient.sample_names[idx] for idx, cms in enumerate(cmuts[s1])
            #          if cms == max(cms for idx, cms in enumerate(cmuts[s1]) if idx != s1)]),
            #      min(cmuts[s1]),
            #      str([patient.sample_names[idx] for idx, cms in enumerate(cmuts[s1]) if cms == min(cmuts[s1])])))

    # displays common mutations among the samples excluding founder mutations
    # similar to a distance matrix
    # logger.debug('Parsimony-informative mutations among the samples: ')
    # for s1 in range(patient.n):
    #     logger.debug(patient.sample_names[s1]+': '
    #                  + ' '.join((str((muts-len(patient.founders)) if muts != -1 else -1)) for muts in cmuts[s1]))

    # mutation patterns are sharing its mutations in exactly the same samples (basis for the weighting scheme)
    patient.mps = defaultdict(set)

    # Compute all clones sharing mutations in exactly the same samples
    # Clones with mutations in all or only one sample are uninformative for
    # the creation of the most parsimonious tree
    for mut_idx, samples in patient.mutations.items():
        # if 1 < len(samples) < len(patient.sample_names):
        if 0 < len(samples):
            patient.mps[samples].add(mut_idx)

    # show the 10 clones supported by the most mutations
    # for key, value in islice(sorted(patient.mps.items(), key=lambda x: len(x[1]), reverse=True), 10):
    #     logger.debug('Mutation pattern {} shares mutations: {} '.format(str(key), value))

    logger.info('Total number of distinct mutation patterns: {}'.format(len(patient.mps)))

    return sharing_status


def get_present_mutations(log_p01):
    """
    Return list of mutations which are present in a least one sample
    :param log_p01: list of tuple posterior (log probability that VAF = 0, log probability that VAF > 0)
                    for each variant
    :return: list of indices of present variants
    """

    # indices list of the present mutations
    present_mutations = []
    pres_lp = math.log(0.5)  # log probability of 50%
    for mut_idx, ps in log_p01.items():

        if any(p1 > pres_lp for _, p1 in ps):
            present_mutations.append(mut_idx)

    return present_mutations


def create_analysis_file(patient, min_sa_cov, analysis_filepath, phylogeny=None, comp_node_frequencies=None,
                         no_replications=0):
    """"
    Create a file with the main data analysis results
    :param patient: instance of the class around the sequencing data of a subject
    :param min_sa_cov: minimum median coverage per sample
    :param analysis_filepath: path to the output analysis file
    :param phylogeny: instance of the phylogeny class
    :param comp_node_frequencies: frequency of reproduced identified patterns (if down-sampling was performed)
    :param no_replications: number of replications per down-sampled fraction of variants (if down-sampling)
    """

    logger.debug('Create analysis output file: {}'.format(analysis_filepath))

    # write analysis to file
    with open(analysis_filepath, 'w') as analysis_file:

        analysis_file.write('# Analyzed data from patient {}.\n'.format(patient.name))
        analysis_file.write('# {} samples passed the filtering.\n'.format(len(patient.sample_names)))
        analysis_file.write('# Total number of detected variants: {} \n'.format(len(patient.mutations)))
        no_present_mutations = len(patient.present_mutations)
        analysis_file.write('# Variants present in one of the passed samples: {} \n'.format(no_present_mutations))

        if min_sa_cov > 0:
            analysis_file.write('# Sample median coverage threshold (otherwise discarded): {} \n\n'.format(
                min_sa_cov))

        # provide some analysis about the raw sequencing data
        for sample_name in patient.sample_names:
            analysis_file.write('# Median phred coverage in sample {}: {} (mean: {:.2f})\n'.format(
                sample_name, np.nanmedian(patient.sample_coverages[sample_name]),
                np.mean(patient.sample_coverages[sample_name])))

            analysis_file.write('# Median MAF in sample {}: {:.2%}\n'.format(
                                sample_name, np.nanmedian(patient.sample_mafs[sample_name])))

        # median and mean coverage
        coverages = []
        for mut_key in patient.mut_reads.keys():
            # for sample_name in patient.mut_reads[mut_key].keys():     # all samples
            for sample_name in patient.sample_names:
                if patient.coverage[mut_key][sample_name] >= 0:
                    coverages.append(patient.coverage[mut_key][sample_name])

        analysis_file.write('# Median coverage in the used samples of patient {}: {} (mean: {:.2f})\n'.format(
            patient.name, np.nanmedian(coverages), np.mean(coverages)))

        analysis_file.write('# The average number of mutations per sample in patient {} is {}.\n'.format(patient.name,
                            (float(sum(len(muts) for sa_idx, muts in patient.samples.items()))
                             / len(patient.sample_names))))

        analysis_file.write("# {:.2%} ({}/{}) of all distinct mutations are founders. \n".format(
            float(len(patient.founders))/no_present_mutations, len(patient.founders), no_present_mutations))
        analysis_file.write('# In average {:.2%} ({}) mutations are unique (private) per sample. \n'.format(
            (float(len(patient.shared_muts[1])) / len(patient.sample_names)) /
            (sum(len(muts) for sa_idx, muts in patient.samples.items()) / len(patient.sample_names)),
            float(len(patient.shared_muts[1])) / len(patient.sample_names)))

        for sample_name in patient.sample_names:
            analysis_file.write('# Sample {} classifications: '.format(sample_name)
                                + '{} positives; {} negatives; {} positive unknowns, {} negative unknowns;\n'.format(
                                patient.positives[sample_name], patient.negatives[sample_name],
                                patient.unknowns[0][sample_name], patient.unknowns[1][sample_name]))

        if patient.gene_names is not None and len(patient.gene_names) > 0:
            for sa_idx, sample_name in enumerate(patient.sample_names):
                analysis_file.write("# Gene names of mutations present in sample {} ({}): {}\n".format(
                    sample_name, len(patient.samples[sa_idx]),
                    ', '.join(str(patient.gene_names[mut]) for mut in patient.samples[sa_idx])))

        for shared in range(len(patient.sample_names)-1, -1, -1):
            if len(patient.shared_muts[shared]) > 0:
                analysis_file.write('# Mutations present in {} samples ({} of {} = {:.2%}): {} \n'.format(
                    shared, len(patient.shared_muts[shared]), len(patient.mutations),
                    float(len(patient.shared_muts[shared])) / len(patient.mutations),
                    ', '.join(patient.mut_keys[m] for m in patient.shared_muts[shared])))

        # if phylogeny is not None:
        if isinstance(phylogeny, SimplePhylogeny) and phylogeny.compatible_tree is not None:
            # how many mutations are are compatible on an evolutionary tree
            analysis_file.write(
                '# Phylogeny: {:.2%} ({} / {}) of all mutations are compatible on an evolutionary tree. \n'.format(
                    float(len(phylogeny.solutions[0].compatible_mutations)) / no_present_mutations,
                    len(phylogeny.solutions[0].compatible_mutations), no_present_mutations))

            # Percentage of conflicting mutations versus shared mutations (excluding unique and founder mutations)
            # evidence for contradictions in the current evolutionary theory of cancer???
            # single cell sequencing will be need to shade light into this puzzle
            no_shared_muts = sum(len(patient.shared_muts[shared]) for shared in range(2, len(patient.sample_names)))
            analysis_file.write(
                '# Phylogeny: {:.2%} ({} / {}) of all shared (excluding unique and founding) '
                + 'mutations are conflicting. \n'.format(
                    float(len(phylogeny.solutions[0].conflicting_mutations)) / no_shared_muts,
                    len(phylogeny.solutions[0].conflicting_mutations), no_shared_muts))

        if isinstance(phylogeny, MaxLHPhylogeny) and phylogeny.mlh_tree is not None:
            # how many positions are evolutionarily incompatible
            analysis_file.write(
                '# Maximum likelihood phylogeny: {} putative false-positives and {} putative false-negatives. \n'
                .format(len(phylogeny.solutions[0].false_positives), len(phylogeny.solutions[0].false_negatives)))

            # add information about false-positives
            for mut_idx, samples in phylogeny.solutions[0].false_positives.items():
                analysis_file.write('# Putative false-positive of variant {} in samples {}\n'.format(
                    patient.mut_keys[mut_idx], ', '.join(patient.sample_names[sa_idx] for sa_idx in samples)))

            # add information about false-negatives
            for mut_idx, samples in phylogeny.solutions[0].false_negatives.items():
                analysis_file.write('# Putative false-negative of variant {} in samples {}\n'.format(
                    patient.mut_keys[mut_idx], ', '.join(patient.sample_names[sa_idx] for sa_idx in samples)))

            # add information about false-negatives due to too low coverage (unknowns)
            for mut_idx, samples in phylogeny.solutions[0].false_negative_unknowns.items():
                analysis_file.write('# Putative present mutation of unknown variant {} in samples {}\n'.format(
                    patient.mut_keys[mut_idx], ', '.join(patient.sample_names[sa_idx] for sa_idx in samples)))

        analysis_file.write('id\tname\t'+('\t'.join(
            'pres'+str(shared) for shared in range(len(patient.sample_names), 0, -1))) + '\t\n')

        for sa_idx, sa_name in enumerate(patient.sample_names):
            analysis_file.write(str(sa_idx+1)+'\t'+str(sa_name)+'\t'
                                + '\t'.join((str(sum(1 for mut, samples in patient.mutations.items()
                                                 if mut in patient.samples[sa_idx] and len(samples) == shared)))
                                            for shared in range(len(patient.sample_names), 0, -1))+'\t\n')

        logger.info('Created analysis file for patient {}: {} \n'.format(patient.name, analysis_filepath))


def create_data_analysis_file(patient, analysis_filepath):
    """
    Create raw data analysis file in TSV format
    """

    logger.debug('Create data analysis output file: {}'.format(analysis_filepath))

    # write analysis to file
    with open(analysis_filepath, 'w') as analysis_file:

        csvwriter = csv.writer(analysis_file, delimiter='\t')

        analysis_file.write('# Analyzed data from patient {}.\n'.format(patient.name))

        # write header
        csvwriter.writerow(('Sample', 'MedianCoverage', 'EstPurity', 'MedianMAF',
                            'BayPresent', 'BayAbsent', 'Present', 'Absent', 'Unknown'))

        # build up output data sequentially
        pres_lp = math.log(0.5)
        for sa_idx, sample_name in enumerate(patient.sample_names):

            row = list()
            row.append(sample_name)

            row.append(np.nanmedian(patient.sample_coverages[sample_name]))
            if sample_name in patient.purities:
                row.append('{:.5f}'.format(patient.purities[sample_name]))
            else:
                row.append('-')

            row.append('{:.3f}'.format(np.nanmedian(patient.sample_mafs[sample_name])))

            # Bayesian inference model
            # present if probability to be present is greater than 50%
            row.append(sum(1 for ps in patient.log_p01.values() if ps[sa_idx][1] > pres_lp))
            row.append(sum(1 for ps in patient.log_p01.values() if ps[sa_idx][1] <= pres_lp))

            # conventional classification
            row.append(patient.positives[sample_name])
            row.append(patient.negatives[sample_name])
            row.append(patient.unknowns[0][sample_name]+patient.unknowns[1][sample_name])

            # write row to file
            csvwriter.writerow(row)


def print_genetic_distance_table(patient):
    """
    Print the genetic distance among the samples and the Jaccard similarity coefficient
    :param patient: instance of class around sequencing data of a subject
    """

    gds = [[0 for _ in range(patient.n)] for _ in range(patient.n)]
    homogeneity = [[0 for _ in range(patient.n)] for _ in range(patient.n)]

    for s1_idx in range(len(patient.data[0])):
        for s2_idx in range(len(patient.data[0])):

            disagree = 0
            present_agree = 0
            no_known_variants = 0

            for mut_idx in range(len(patient.data)):

                if patient.data[mut_idx][s1_idx] > 0:          # mutation present
                    if patient.data[mut_idx][s2_idx] == 0:     # disagree
                        no_known_variants += 1
                        disagree += 1

                    # mutation present in both
                    elif patient.data[mut_idx][s2_idx] > 0:    # agree
                        no_known_variants += 1
                        present_agree += 1

                    elif patient.data[mut_idx][s2_idx] == NEG_UNKNOWN:
                        # no indication of a mutation at this position but low coverage
                        pass

                elif patient.data[mut_idx][s1_idx] == 0:     # mutation absent
                    if patient.data[mut_idx][s2_idx] > 0:
                        no_known_variants += 1
                        disagree += 1

                    elif patient.data[mut_idx][s2_idx] == 0:   # agree
                        pass

                    elif patient.data[mut_idx][s2_idx] == POS_UNKNOWN:
                        # indication of a mutation at this position but insufficient mutant reads
                        # maybe penalize distance with two epsilon
                        pass

                # mutation believed to be present but insufficient mutant reads
                elif patient.data[mut_idx][s1_idx] == POS_UNKNOWN:
                    if patient.data[mut_idx][s2_idx] == 0:
                        pass

                # mutation believed to be absent but insufficient coverage
                elif patient.data[mut_idx][s1_idx] == NEG_UNKNOWN:
                    if patient.data[mut_idx][s2_idx] > 0:
                        pass

            gds[s1_idx][s2_idx] = disagree
            homogeneity[s1_idx][s2_idx] = float(present_agree) / no_known_variants

    # Produce latex table with the genetic distance between samples
    print('Genetic distance across the samples:')
    for s1_idx in range(patient.n):
        print('{} & '.format(patient.sample_names[s1_idx])
              + ' & '.join('${}$'.format(gds[s1_idx][s2_idx]) for s2_idx in range(patient.n))+' \\\\')

    print('Homogeneity index based on the fraction of shared mutations:')
    for s1_idx in range(patient.n):
        print('{} & '.format(patient.sample_names[s1_idx])
              + ' & '.join('${:.2f}$'.format(homogeneity[s1_idx][s2_idx]) for s2_idx in range(patient.n))+' \\\\')
