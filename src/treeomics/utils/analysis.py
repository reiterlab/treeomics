#!/usr/bin/python
"""Generates analysis file providing an overview of the results"""
__author__ = 'Johannes REITER'


import logging
import csv
import settings
from utils.int_settings import NEG_UNKNOWN, POS_UNKNOWN
import numpy as np
from phylogeny.simple_phylogeny import SimplePhylogeny
from phylogeny.max_lh_phylogeny import MaxLHPhylogeny

# get logger for application
logger = logging.getLogger('treeomics')


def create_analysis_file(patient, min_sa_cov, analysis_filepath, phylogeny=None, comp_node_frequencies=None,
                         no_replications=0):
    """"
    Create a file with the main data analysis results
    :param min_sa_cov: minimum median coverage per sample
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
                sample_name, np.median(patient.sample_phred_coverages[sample_name]),
                np.mean(patient.sample_phred_coverages[sample_name])))

            analysis_file.write('# Median MAF in sample {}: {:.2%}\n'.format(
                                sample_name, np.median(patient.sample_mafs[sample_name])))

        # median and mean coverage
        coverages = []
        for mut_key in patient.mut_reads.keys():
            # for sample_name in patient.mut_reads[mut_key].keys():     # all samples
            for sample_name in patient.sample_names:
                if patient.phred_coverage[mut_key][sample_name] >= 0:
                    coverages.append(patient.phred_coverage[mut_key][sample_name])

        analysis_file.write('# Median coverage in the used samples of patient {}: {} (mean: {:.2f})\n'.format(
            patient.name, np.median(coverages), np.mean(coverages)))

        # total number of exonic mutations if information is available
        if patient.mut_functions is not None and len(patient.mut_functions):
            no_exonic_muts = sum(1 for mut_func in patient.mut_functions if mut_func in settings.EXONIC_FILTER)
            analysis_file.write('# Number of exonic mutations {}/{} (={:.3f}). \n'.format(
                no_exonic_muts, no_present_mutations, (float(no_exonic_muts) / no_present_mutations)))

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
                    float(len(phylogeny.compatible_mutations)) / no_present_mutations,
                    len(phylogeny.compatible_mutations), no_present_mutations))

            # Percentage of conflicting mutations versus shared mutations (excluding unique and founder mutations)
            # evidence for contradictions in the current evolutionary theory of cancer???
            # single cell sequencing will be need to shade light into this puzzle
            no_shared_muts = sum(len(patient.shared_muts[shared]) for shared in range(2, len(patient.sample_names)))
            analysis_file.write(
                '# Phylogeny: {:.2%} ({} / {}) of all shared (excluding unique and founding) '
                + 'mutations are conflicting. \n'.format(
                    float(len(phylogeny.conflicting_mutations)) / no_shared_muts,
                    len(phylogeny.conflicting_mutations), no_shared_muts))

            # write robustness analysis to file
            if comp_node_frequencies is not None:
                analysis_file.write('# Robustness analysis through {} replications \n'.format(no_replications))
                analysis_file.write('# Mutation pattern  \t {} \n'.format(' \t '.join(
                                    '{}%'.format(fr) for fr in sorted(comp_node_frequencies.keys()))))
                for node in sorted(phylogeny.compatible_nodes.keys(),
                                   key=lambda k: -phylogeny.node_weights[k]):
                    # print only parsimony informative MPs
                    if 1 < len(node) < len(patient.sample_names):
                        analysis_file.write('# ({}) \t {} \n'.format(','.join(str(n+1) for n in sorted(node)), ' \t '.join(
                            '{:.3f}'.format(
                                comp_node_frequencies[fr][node]) for fr in sorted(comp_node_frequencies.keys()))))

        if isinstance(phylogeny, MaxLHPhylogeny) and phylogeny.mlh_tree is not None:
            # how many positions are evolutionarily incompatible
            analysis_file.write(
                '# Maximum likelihood phylogeny: {} putative false-positives and {} putative false-negatives. \n'
                .format(len(phylogeny.false_positives), len(phylogeny.false_negatives)))

            # add information about false-positives
            for mut_idx, samples in phylogeny.false_positives.items():
                analysis_file.write('# Putative false-positive of variant {} in samples {}\n'.format(
                    patient.mut_keys[mut_idx], ', '.join(patient.sample_names[sa_idx] for sa_idx in samples)))

            # add information about false-negatives
            for mut_idx, samples in phylogeny.false_negatives.items():
                analysis_file.write('# Putative false-negative of variant {} in samples {}\n'.format(
                    patient.mut_keys[mut_idx], ', '.join(patient.sample_names[sa_idx] for sa_idx in samples)))

            # add information about false-negatives due to too low coverage (unknowns)
            for mut_idx, samples in phylogeny.false_negative_unknowns.items():
                analysis_file.write('# Putative present mutation of unknown variant {} in samples {}\n'.format(
                    patient.mut_keys[mut_idx], ', '.join(patient.sample_names[sa_idx] for sa_idx in samples)))

        analysis_file.write('id\tname\t'+('\t'.join(
            'pres'+str(shared) for shared in range(len(patient.sample_names), 0, -1)))
            + '\t' + ('exonic\t' if patient.mut_functions is not None and len(patient.mut_functions) else '')+'\n')

        for sa_idx, sa_name in enumerate(patient.sample_names):
            analysis_file.write(str(sa_idx+1)+'\t'+str(sa_name)+'\t'
                                + '\t'.join((str(sum(1 for mut, samples in patient.mutations.items()
                                                 if mut in patient.samples[sa_idx] and len(samples) == shared)))
                                            for shared in range(len(patient.sample_names), 0, -1))+'\t\n')

            # number of exonic mutations per sample if information is available
            if patient.mut_functions is not None and len(patient.mut_functions):
                analysis_file.write(str(sum(1 for mut in patient.samples[sa_idx]
                                            if patient.mut_functions[mut] in settings.EXONIC_FILTER)) + '\t\n')

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
        csvwriter.writerow(('Sample', 'Med. coverage', 'Med. dis. coverage', 'Med. MAF', 'Present', 'Absent',
                            'Unknown (present)', 'Unknown (absent)'))

        # build up output data sequentially
        for sample_name in patient.sample_names:

            row = list()
            row.append(sample_name)

            row.append(np.median(patient.sample_phred_coverages[sample_name]))
            if patient.sample_dis_phred_coverages is not None:
                row.append(np.median(patient.sample_dis_phred_coverages[sample_name]))
            else:
                row.append('n/a')

            row.append('{:.3f}'.format(np.median(patient.sample_mafs[sample_name])))

            row.append(patient.positives[sample_name])
            row.append(patient.negatives[sample_name])
            row.append(patient.unknowns[0][sample_name])
            row.append(patient.unknowns[1][sample_name])

            # write row to file
            csvwriter.writerow(row)


def print_genetic_distance_table(patient):
    """
    Print the genetic distance among the samples and the homogeneity index (Jaccard similarity coefficient)
    :param patient:
    :return:
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
