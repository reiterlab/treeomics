#!/usr/bin/python
"""Analyze genetic similarity across samples of a subject"""

import logging
import math
import utils.int_settings as def_sets
from utils.int_settings import NEG_UNKNOWN, POS_UNKNOWN

__author__ = 'jreiter'
__date__ = 'January, 2017'


# get logger for application
logger = logging.getLogger('treeomics')


def calculate_genetic_similarity(patient):
    """
    Calculate the genetic similarity among all samples in this patient based on
    significantly mutated or well-covered mutations (frequentist approach)
    (1) calculate the genetic distance between all pairs
    (2) calculate the Jaccard similarity coefficient between all pairs
    :param patient: data structure around sequencing data of a subject
    :return genetic distance matrix, Jaccard similarity coefficient matrix,
            Jaccard similarity coefficient matrix (excluding founders)
    """

    gen_dis = [[0 for _ in range(patient.n)] for _ in range(patient.n)]
    sim_coeff = [[0.0 for _ in range(patient.n)] for _ in range(patient.n)]
    sim_coeff_ex = [[0.0 for _ in range(patient.n)] for _ in range(patient.n)]  # excluding founders

    likely_founders = set(mut_idx for mut_idx in patient.data.keys()
                          if all(vaf > 0 or vaf < 0 for vaf in patient.data[mut_idx]))

    for s1_idx in range(len(patient.data[0])):
        for s2_idx in range(len(patient.data[0])):

            disagree = 0
            present_agree = 0
            no_known_variants = 0
            considered_vars = set()

            for mut_idx in range(len(patient.data)):

                if patient.data[mut_idx][s1_idx] > 0:  # mutation present
                    if patient.data[mut_idx][s2_idx] == 0:  # disagree
                        no_known_variants += 1
                        disagree += 1
                        considered_vars.add(mut_idx)

                    # mutation present in both
                    elif patient.data[mut_idx][s2_idx] > 0:  # agree
                        no_known_variants += 1
                        present_agree += 1
                        considered_vars.add(mut_idx)

                    elif patient.data[mut_idx][s2_idx] == NEG_UNKNOWN:
                        # no indication of a mutation at this position but low coverage
                        pass

                elif patient.data[mut_idx][s1_idx] == 0:  # mutation absent
                    if patient.data[mut_idx][s2_idx] > 0:
                        no_known_variants += 1
                        disagree += 1
                        considered_vars.add(mut_idx)

                    elif patient.data[mut_idx][s2_idx] == 0:  # agree
                        # mutation is absent in both samples and is therefore not relevant for these calculations
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

            gen_dis[s1_idx][s2_idx] = disagree
            sim_coeff[s1_idx][s2_idx] = 1.0 if no_known_variants == 0 \
                else (float(present_agree) / no_known_variants)

            # determine the number of likely founders among variants which are
            # known with certainty in the pair of samples
            no_founders = len(likely_founders.intersection(considered_vars))
            sim_coeff_ex[s1_idx][s2_idx] = 1.0 if no_known_variants - no_founders <= 0 \
                else ((float(present_agree) - no_founders) / (no_known_variants - no_founders))

    # Produce tables with the genetic distance between samples
    # print('Similarity index based on the fraction of shared mutations (including founders):')
    # print('Sample \t '+' \t '.join(patient.sample_names[sa_idx].replace('_', ' ') for sa_idx in range(patient.n))+'')
    # for s1_idx in range(patient.n):
    #     print('{} \t '.format(patient.sample_names[s1_idx].replace('_', ' ')) +
    #           ' \t '.join('{:.2f}'.format(patient.sim_coff[s1_idx][s2_idx]) for s2_idx in range(patient.n))+'')
    #
    # print('Similarity index based on the fraction of shared mutations (excluding founders):')
    # print('Sample \t '+' \t '.join(patient.sample_names[sa_idx].replace('_', ' ') for sa_idx in range(patient.n))+'')
    # for s1_idx in range(patient.n):
    #     print('{} \t '.format(patient.sample_names[s1_idx].replace('_', ' ')) +
    #           ' \t '.join('{:.2f}'.format(patient.sim_coff_ex[s1_idx][s2_idx]) for s2_idx in range(patient.n))+'')
    #
    # # Produce table with the genetic distance between samples
    # print('Genetic distance across the samples:')
    # print('Sample \t '+' \t '.join(patient.sample_names[sa_idx].replace('_', ' ') for sa_idx in range(patient.n)))
    # for s1_idx in range(patient.n):
    #     print('{} \t '.format(patient.sample_names[s1_idx].replace('_', ' ')) +
    #           ' \t '.join('{}'.format(patient.gen_dis[s1_idx][s2_idx]) for s2_idx in range(patient.n))+'')

    return gen_dis, sim_coeff, sim_coeff_ex


def calculate_bi_genetic_similarity(patient):
    """
    Calculate the genetic similarity among all samples in this patient based on the
    Bayesian inference model and a minimum confidence that a mutation is present or absent
    (1) calculate the genetic distance between all pairs
    (2) calculate the Jaccard similarity coefficient between all pairs
    :param patient: data structure around sequencing data of a subject
    """

    bi_sim_coeff = [[0.0 for _ in range(patient.n)] for _ in range(patient.n)]
    bi_gen_dist = [[0 for _ in range(patient.n)] for _ in range(patient.n)]

    # log probability to be classified with at least the given confidence threshold
    conf_clas_lpth = math.log(def_sets.CLA_CONFID_TH)

    for s1_idx in range(len(patient.data[0])):
        for s2_idx in range(s1_idx+1):

            if s1_idx == s2_idx:
                bi_sim_coeff[s1_idx][s2_idx] = 1.0
                continue

            no_present_vars = 0
            for mut_idx in range(len(patient.log_p01)):

                s1_p0, s1_p1 = patient.log_p01[mut_idx][s1_idx]
                s2_p0, s2_p1 = patient.log_p01[mut_idx][s2_idx]

                if s1_p1 > conf_clas_lpth or s2_p1 > conf_clas_lpth:  # variant present in at least one sample

                    # use only variants confidently classified as present or absent
                    if (s1_p1 > conf_clas_lpth or s1_p0 > conf_clas_lpth) and \
                            (s2_p1 > conf_clas_lpth or s2_p0 > conf_clas_lpth):

                        no_present_vars += 1
                        sim_prob = math.exp(s1_p1 + s2_p1) + math.exp(s1_p0 + s2_p0)
                        bi_sim_coeff[s1_idx][s2_idx] += sim_prob
                        if sim_prob < 0.5:
                            bi_gen_dist[s1_idx][s2_idx] += 1

            bi_sim_coeff[s1_idx][s2_idx] /= no_present_vars
            bi_sim_coeff[s2_idx][s1_idx] = bi_sim_coeff[s1_idx][s2_idx]
            bi_gen_dist[s2_idx][s1_idx] = bi_gen_dist[s1_idx][s2_idx]
            # logger.debug('Probabilistic similarity coefficient of {} and {}: {:.1%}'.format(
            #     patient.sample_names[s1_idx], patient.sample_names[s2_idx], bi_sim_coeff[s1_idx][s2_idx]))
            # logger.debug('Genetic distance of {} and {}: {:.0f}'.format(
            #     patient.sample_names[s1_idx], patient.sample_names[s2_idx], bi_gen_dist[s1_idx][s2_idx]))

    return bi_gen_dist, bi_sim_coeff
