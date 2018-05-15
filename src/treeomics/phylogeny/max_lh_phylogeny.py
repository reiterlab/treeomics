#!/usr/bin/python
"""Infers most reliable mutation pattern based on a bayesian inference model"""
import logging
import itertools
import math
from scipy.stats import binom
from collections import defaultdict
from itertools import combinations
import heapq
import numpy as np
import networkx as nx
import sys
from copy import deepcopy
import phylogeny.cplex_solver as cps
from phylogeny.phylogeny_utils import Phylogeny
from utils.statistics import get_log_p0
import utils.int_settings as def_sets

__author__ = 'Johannes REITER'
__date__ = 'July, 2015'

# get logger for application
logger = logging.getLogger('treeomics')


class MaxLHPhylogeny(Phylogeny):
    """
    Find likely (low reliability scores) technical and biological artifacts in the data
    Infer evolutionary tree from the resolved data
    """

    def __init__(self, patient, mps, loh_frequency=0.0):

        Phylogeny.__init__(self, patient, mps)

        # list of mutation pattern ids mapping to the corresponding mutation pattern
        self.idx_to_mp = None
        # dictionary from mutation patterns to the above column ids
        self.mp_col_ids = None
        # list of dictionaries [mut_idx, mp_col_idx]
        self.mp_weights = None
        # map from identified putative subclones to their original sample
        self.sc_sample_ids = None

        # most likely but also compatible mutation pattern for each variant
        self.max_lh_nodes = None
        self.max_lh_mutations = None
        # map from mut_idx to the weight of the highest ranked evolutionarily compatible MP
        self.max_lh_weights = None

        self.bootstrapping_values = None

        # weighted likelihood of the identified compatible mutation patterns based on the likelihood of the
        # individual solutions in the solution space
        self.weighted_node_lh = None

        # list of solutions as inferred by the MILP solution pool
        self.solutions = None

        self.mlh_tree = None

        self.max_no_mps = None

        # calculate the minimal number of variant reads k_min at the median coverage and purity such that p_1 > 50%
        k_mins = None
        called_ps = []
        missed_ps = []
        pres_lp = math.log(0.5)
        lh = 1.0
        for sa_idx, sample_name in enumerate(patient.sample_names):
            # calculate posterior according to prior, estimated purity and data
            for k in range(min(5000, int(np.nanmedian(patient.sample_coverages[sample_name])))):
                _, p1 = get_log_p0(np.nanmedian(patient.sample_coverages[sample_name]), k, self.patient.bi_error_rate,
                                   self.patient.bi_c0, cutoff_f=self.patient.get_cutoff_frequency(sample_name),
                                   pseudo_alpha=def_sets.PSEUDO_ALPHA, pseudo_beta=patient.betas[sample_name])
                if p1 > pres_lp:
                    k_mins = k
                    break

            if k_mins is None:
                if np.nanmedian(patient.sample_coverages[sample_name]) > 5000:
                    for k in range(5000, int(np.nanmedian(patient.sample_coverages[sample_name])),
                                   int(np.nanmedian(patient.sample_coverages[sample_name]) * 0.01)):

                        _, p1 = get_log_p0(np.nanmedian(patient.sample_coverages[sample_name]), k,
                                           self.patient.bi_error_rate,
                                           self.patient.bi_c0, cutoff_f=self.patient.get_cutoff_frequency(sample_name),
                                           pseudo_alpha=def_sets.PSEUDO_ALPHA, pseudo_beta=patient.betas[sample_name])
                        if p1 > pres_lp:
                            k_mins = k
                            break

                    else:
                        k_mins = int(np.nanmedian(patient.sample_coverages[sample_name]) * 0.05)
                        logger.warning('Treeomics struggles with the high median coverage {} and made some assumption.'
                                       .format(np.nanmedian(patient.sample_coverages[sample_name])))
                else:
                    logger.warning('Median coverage is very low. Perhaps exclude this sample: {}'.format(sample_name))
                    k_mins = 1

            logger.debug('{}: Minimum number of mutant reads such that presence probability is greater than 50%: {}.'
                         .format(sample_name, k_mins))
            called_ps.append(1.0 - binom.cdf(k_mins-1, int(np.nanmedian(patient.sample_coverages[sample_name])),
                                             self.patient.bi_error_rate))
            logger.debug('Probability to observe an incorrectly called variant: {:.3%}'.format(called_ps[-1]))

            if sample_name in patient.purities:
                missed_ps.append(binom.cdf(k_mins-1, int(np.nanmedian(patient.sample_coverages[sample_name])),
                                           self.patient.purities[sample_name] / 2.0))
            else:
                missed_ps.append(binom.cdf(k_mins-1, int(np.nanmedian(patient.sample_coverages[sample_name])),
                                           np.nanmedian(self.patient.sample_mafs[sample_name])))
            logger.debug('Probability to miss a clonal variant: {:.1e}'.format(missed_ps[-1]))

            # probability that all calls are correct
            # ensure that this probability is positive but also not too close to 1
            lh *= max((1.0 - max(called_ps[-1], 1.0 - def_sets.MAX_PRE_PROB) - missed_ps[-1] -
                      max(loh_frequency, 1.0 - def_sets.MAX_ABS_PROB)), 0.5)

        # probability that at least one call is wrong
        lh = 1 - lh
        subclone_th = 1.0
        self.min_score = -math.log(1.0 - lh) * subclone_th
        logger.debug('Likelihood of a pattern with at least one false-positive or false-negative' +
                     ' and an LOH probability along a linage of {}: {:.3e}'.format(loh_frequency, lh))

        logger.info('Minimum reliability score value to be considered as a potential subclone: {:.3e}'.format(
            self.min_score))

        lh_99 = 1.0 - math.pow(0.99, len(patient.sample_names))
        logger.debug('Reliability score of a pattern with 99% certainty in each call: {:.3e}'.format(
            -math.log(1.0 - lh_99)))

        lh_999 = 1.0 - math.pow(0.999, len(patient.sample_names))
        logger.debug('Reliability score of a pattern with 99.9% certainty in each call: {:.3e}'.format(
            -math.log(1.0 - lh_999)))

        lh_9999 = 1.0 - math.pow(0.9999, len(patient.sample_names))
        logger.debug('Reliability score of a pattern with 99.99% certainty in each call: {:.3e}'.format(
            -math.log(1.0 - lh_9999)))

    def infer_max_lh_tree(self, subclone_detection=False, pool_size=1, no_plotted_solutions=1, no_bootstrap_samples=0,
                          max_no_mps=None, time_limit=None, n_max_threads=0):
        """
        Infer maximum likelihood tree via calculation reliability scores for each
        possible mutation pattern from the likelihood that no variant has this pattern
        The inferred solution represents the reliable and evolutionary compatible mutation patterns
        The mutation pattern of each variant is given by the mp maximizing its likelihood
        :param subclone_detection: is subclone detection enabled?
        :param pool_size: number of best solutions explored by ILP solver to estimate confidence
        :param no_plotted_solutions: number of best solutions from the solution pool that will be plotted
        :param max_no_mps: only the given number of most likely (by joint likelihood) mutation patterns
            is explored per variant
        :param no_bootstrap_samples: number of samples with replacement for the bootstrapping
        :param time_limit: time limit for MILP solver in seconds
        :param n_max_threads: Sets the default maximal number of parallel threads that will be invoked by CPLEX
                              (0: default, let CPLEX decide; 1: single threaded; N: uses up to N threads)
        :return inferred evolutionary tree
        """

        if subclone_detection:
            self.patient.sc_names = deepcopy(self.patient.sample_names)
        else:
            self.patient.sc_names = self.patient.sample_names

        # only the given number of most likely (by joint likelihood) mutation patterns is explored per variant
        self.max_no_mps = max_no_mps

        # necessary to map from identified putative subclones to their original sample
        self.sc_sample_ids = dict()

        # compute various mutation patterns (nodes) and their reliability scores
        self.node_scores, self.idx_to_mp, self.mp_col_ids, self.mp_weights = infer_ml_graph_nodes(
            self.patient.log_p01, self.patient.sample_names, self.patient.mut_keys, gene_names=self.patient.gene_names,
            max_no_mps=max_no_mps)

        while True:
            # create conflict graph which forms the input to the ILP
            self.cf_graph = create_conflict_graph(self.node_scores)

            # translate the conflict graph into a minimum vertex cover problem
            # and solve this using integer linear programming
            self.solutions, self.weighted_node_lh = cps.solve_conflicting_phylogeny(
                self.cf_graph, len(self.patient.log_p01), pool_size, no_plotted_solutions,
                time_limit=time_limit, n_max_threads=n_max_threads)

            opt_sol = self.solutions[0]
            opt_sol.assign_variants(self, max_no_mps)
            opt_sol.find_artifacts(self)

            logger.info('Putative false-positives {}, putative false-negatives {}, put. false neg. unknowns {}'.format(
                sum(len(fps) for mut_idx, fps in opt_sol.false_positives.items()),
                sum(len(fns) for mut_idx, fns in opt_sol.false_negatives.items()),
                sum(len(fns) for mut_idx, fns in opt_sol.false_negative_unknowns.items())))
            if opt_sol.conflicting_mutations is not None and len(opt_sol.conflicting_mutations) > 0:
                logger.warning('{} evolutionarily incompatible variants occurred{}'.format(
                    len(opt_sol.conflicting_mutations),
                    ' in: '+', '.join(self.patient.gene_names[mut_idx] for mut_idx in opt_sol.conflicting_mutations) if
                    self.patient.gene_names is not None else '.'))

            # find parsimony-informative evolutionarily incompatible mps with high likelihood
            if subclone_detection:

                self.sc_sample_ids, found_subclones = self.find_subclones(opt_sol)
                if not found_subclones:
                    # not enough evidence for additional new subclones
                    break

            else:           # detection of subclones is disabled
                break

        if len(self.sc_sample_ids) > 0:
            logger.info('{} putative subclones have been detected: {}'.format(
                        len(self.patient.sc_names)-len(self.patient.sample_names),
                        ', '.join(self.patient.sc_names[sc_idx] for sc_idx in range(len(self.patient.sample_names),
                                                                                    len(self.patient.sc_names)))))
        if max_no_mps is not None:
            logger.info('{}/{} are compatible on the inferred evolutionary tree. '.format(
                len(opt_sol.max_lh_mutations), len(opt_sol.max_lh_mutations)+len(opt_sol.conflicting_mutations)) +
                '{} ({:.3%}) of the mutations are conflicting.'.format(
                len(opt_sol.conflicting_mutations), float(len(opt_sol.conflicting_mutations)) /
                    (len(opt_sol.max_lh_mutations)+len(opt_sol.conflicting_mutations))))

        # find founders, shared and unique mutations in updated mutation list
        opt_sol.infer_mutation_patterns(self)

        for sa_idx, sample_name in enumerate(self.patient.sample_names):
            logger.info('Inferred number of variants in sample {} in a perfect and persistent phylogeny: {}'.format(
                sample_name, len(opt_sol.variants[sa_idx])))

        # is bootstrapping enabled?
        if no_bootstrap_samples > 0:
            if subclone_detection:
                logger.error('Bootstrapping analysis does not support subclone detection! ')
                logger.error('Go to settings.py and set NO_BOOSTRAP_SAMPLES to 0 or SUBCLONE_DETECTION to False.')
                raise RuntimeError('Bootstrapping analysis does not support subclone detection! ')
            else:
                self.do_bootstrapping(no_bootstrap_samples, time_limit=time_limit)

        if self.bootstrapping_values is not None:
            confidence_values = self.bootstrapping_values
            logger.info('Confidence values for branching are given by bootstrapping.')

        elif self.weighted_node_lh is not None:
            confidence_values = self.weighted_node_lh
            logger.info('Confidence values for branching are given by exploring the whole solution space and '
                        'the weighted likelihoods of each solution.')
        else:
            confidence_values = None
            logger.info('Confidence values will not be provided.')

        # construct a phylogenetic tree from the maximum likelihood mutation patterns
        self.mlh_tree = self.infer_evolutionary_tree(opt_sol.shared_mlh_mps, opt_sol.mlh_founders,
                                                     opt_sol.mlh_unique_mutations, confidence=confidence_values)

        return self.mlh_tree

    def do_bootstrapping(self, no_samples, time_limit=None):
        """
        Validate the robustness of the identified most reliable mutation patterns through bootstrapping
        :param no_samples: Number of samples with replacement for the bootstrapping
        :param time_limit: time limit for MILP solver in seconds
        """

        # Most reliable mutation patterns need to be identified before
        if self.compatible_nodes is None:
            self.infer_max_lh_tree(subclone_detection=False, time_limit=time_limit)

        node_frequencies = cps.bootstrapping_solving(
            self.cf_graph, self.mp_weights, self.idx_to_mp, no_samples)

        self.bootstrapping_values = dict()

        # calculate the frequency with which compatible mutation patterns are reproduced
        # when bootstrapping is used
        for node in self.solutions[0].compatible_nodes:

            # consider only parsimony-informative MPs
            if len(node) <= 1 or len(node) == len(self.patient.sample_names):
                continue

            if node in node_frequencies:
                self.bootstrapping_values[node] = float(node_frequencies[node]) / no_samples
            else:
                self.bootstrapping_values[node] = 0.0

            logger.info('Bootstrapping value for node {}: {:.1%}'.format(node, self.bootstrapping_values[node]))

    def validate_node_robustness(self, no_replications, time_limit=None):
        """
        Validate the robustness of the identified most reliable mutation patterns
        through down-sampling
        :param no_replications: Number of replications for each used fraction of variants
        :param time_limit: time limit for MILP solver in seconds
        :return observed frequencies of the mutation patterns in each used fraction of variants
        """

        # Most reliable mutation patterns need to be identified before
        if self.compatible_nodes is None:
            self.infer_max_lh_tree(subclone_detection=False, time_limit=time_limit)

        node_frequencies = cps.solve_downsampled_nodes(
            self.cf_graph, self.mp_weights, self.idx_to_mp, no_replications)

        comp_node_frequencies = defaultdict(dict)

        # calculate the frequency with which compatible mutation patterns are reproduced
        # when only a subset of variants are used
        no_comp_pars_inf_mps = sum(1 for node in self.solutions[0].compatible_nodes
                                   if 1 < len(node) < len(self.patient.sample_names))
        for sample_fraction in sorted(node_frequencies.keys()):
            freq = 0
            for node in self.solutions[0].compatible_nodes:

                # consider only parsimony-informative MPs
                if len(node) <= 1 or len(node) == len(self.patient.sample_names):
                    continue

                if node in node_frequencies[sample_fraction]:
                    comp_node_frequencies[sample_fraction][node] = \
                        float(node_frequencies[sample_fraction][node]) / no_replications
                    freq += comp_node_frequencies[sample_fraction][node]
                else:
                    comp_node_frequencies[sample_fraction][node] = 0.0

            logger.info('Fraction of confirmed mutation patterns with {}% of the variants: {:.1%}'.format(
                sample_fraction, float(freq) / no_comp_pars_inf_mps))

        # produce latex table for paper
        # displayed_fractions = [10, 20, 50, 80, 90, 95]
        # print('\\textbf Mutation pattern  & {} \\\\'.format(' & '.join(
        #       '\\textbf {:.0f}$\%$'.format(fr) for fr in displayed_fractions)))
        # for node in sorted(self.compatible_nodes.keys(), key=lambda k: -self.cf_graph.node[k]['weight']):
        #
        #     print('({}) & {} \\\\'.format(','.join(str(n+1) for n in sorted(node)), ' & '.join(
        #         '{:.1f}$\%$'.format(comp_node_frequencies[fr][node]*100.0) for fr in displayed_fractions)))

        return comp_node_frequencies

    def find_subclones(self, optimal_solution):
        """

        :param optimal_solution:
        :return:
        """

        updated_nodes, self.sc_sample_ids = self.find_subclonal_mps(optimal_solution)

        if len(updated_nodes) == 0:
            logger.info('There are no more incompatible mutation patterns with a reliability score '
                        'of at least {:.1e}.'.format(self.min_score))
            return self.sc_sample_ids, False

        return self.sc_sample_ids, True

    def find_subclonal_mps(self, opt_sol):
        """
        Identify evolutionarily incompatible mutation patterns with a reliability score greater than min_score which
        where variants in some samples are subclonal
        :param opt_sol: optimal solution
        :return updated_nodes:
        """

        updated_nodes = dict()

        # find incompatible mutation pattern with the highest reliability score
        for mp in sorted(opt_sol.incompatible_nodes, key=lambda k: -self.node_scores[k]):

            if (self.node_scores[mp] < self.min_score or
                    len(self.patient.sc_names) > len(self.patient.sample_names) +
                    max(5, len(self.patient.sample_names)/2)):
                # no incompatible mutation patterns with high reliability score exist
                return updated_nodes, self.sc_sample_ids

            # determine the number of variants for which this mutation pattern would be
            # more likely than the before assigned evolutionarily compatible mp
            # no_sup_vars = sum(1 for mut_idx in range(len(self.mp_weights))
            #                   if self.mp_weights[mut_idx][self.mp_col_ids[mp]] > self.max_lh_weights[mut_idx])

            # # number of variants for which this mp would be the most likely one
            # no_sup_vars = sum(1 for mut_idx in range(len(self.mp_weights))
            #     if self.mp_weights[mut_idx][self.mp_col_ids[mp]] == max(self.mp_weights[mut_idx].values()))

            # number of variants for which the considered MP is more likely needs to be at least 2
            # and 1% of the total number of variants
            # if no_sup_vars < max(2, len(self.mp_weights) * 0.01):
            #     # mutation pattern is not the most likely for a significant fraction of variants
            #     return updated_nodes, self.sc_sample_ids

            # find highest ranked conflicting mutation pattern that is part of the current solution
            # step (a) in pseudo algorithm
            logger.debug('Highest ranked conflicting mutation pattern: {} and its neighbors'.format(mp))
            logger.debug('Its neighbors: {}'.format(self.cf_graph.neighbors(mp)))
            logger.debug('Conflicting nodes: {}'.format(opt_sol.incompatible_nodes))
            hcmp = max(set(self.cf_graph.neighbors(mp)).difference(opt_sol.incompatible_nodes),
                       key=lambda k: self.node_scores[k])
            logger.info('Highest ranked evolutionary incompatible mutation pattern of MP {} (w: {:.2e}): {} (w: {:.2e})'
                        .format(', '.join(self.patient.sc_names[sc] for sc in mp), self.cf_graph.node[mp]['weight'],
                                ', '.join(self.patient.sc_names[sc] for sc in hcmp),
                                self.cf_graph.node[hcmp]['weight']))

            # infer samples with putative mixed subclones
            msc_samples = list(mp.intersection(hcmp))
            if len(msc_samples) > 1:    # more than one subclone would be needed to generate, reconsider later
                logger.info('Multiple subclones are needed to explain the evolutionary incompatibility. '
                            'Continue search.')
                continue

            sc_sa = msc_samples[0]      # sample where subclone is created

            # if incompatibility is due to mutations in a generated subclone, ignore it
            if sc_sa >= len(self.patient.sample_names):
                logger.info('Multiple subclones within the same sample are needed to explain '
                            'the evolutionary incompatibility. Continue search.')
                continue

            # create new subclones if there is enough evidence
            created_scss = dict()
            # create new mutation pattern resulting from the new subclone
            new_mp = set(mp)

            # one could check if the VAFs of the subclonal mutation indicate subclones
            # no contradicting evidence for subclone found
            logger.info('Found evidence for subclones in sample {}. '.format(self.patient.sc_names[sc_sa]))

            # create new subclone in this sample and update mutation pattern
            created_scss[sc_sa] = len(self.patient.sc_names)
            sc_sa_idx = len(self.patient.sc_names)
            self.sc_sample_ids[sc_sa_idx] = sc_sa
            new_mp.remove(sc_sa)
            new_mp.add(sc_sa_idx)
            # did a subclone of this sample already get generated
            if '_SC' in self.patient.sc_names[sc_sa]:
                name = self.patient.sc_names[sc_sa][:self.patient.sc_names[sc_sa].find('_SC')]
            else:
                name = self.patient.sc_names[sc_sa]
            next_id = int(max(int(n.split('_SC')[1]) if len(n.split('_SC')) > 1 else 0
                              for n in self.patient.sc_names if n.split('_SC')[0] == name))+1
            self.patient.sc_names.append(name+'_SC'+str(next_id))

            # step (d): update mutation patterns in the patient
            new_mp = frozenset(new_mp)

            # update subclones in patient
            self.node_scores[new_mp] = self.node_scores[mp]
            del self.node_scores[mp]
            self.idx_to_mp[self.mp_col_ids[mp]] = new_mp
            self.mp_col_ids[new_mp] = self.mp_col_ids[mp]
            del self.mp_col_ids[mp]
            updated_nodes[mp] = new_mp

            # add subclone to nodes
            self.node_scores[frozenset([sc_sa_idx])] = self.node_scores[frozenset([sc_sa])]
            logger.info('Created new subclones {}'.format(
                        ', '.join(self.patient.sc_names[sc] for sc in new_mp.difference(mp))))

            # step (d) continued: check compatibility with other conflicting clones
            for related_mp in sorted(opt_sol.incompatible_nodes, key=lambda k: -self.cf_graph.node[k]['weight']):
                if related_mp == mp:
                    continue

                # is this conflicting node evolutionary compatible with the previously identified one
                # in other words, check for superset and subset relation
                if related_mp in self.cf_graph.neighbors(mp):
                    continue

                # some of the new subclonal samples have to be present in this conflicting node
                # otherwise there is no evidence that the new subclone is present in the conflicting node
                if not any(new_scsa in related_mp for new_scsa in created_scss):
                    continue

                # check if the MP (mutation pattern) would be evolutionary compatible if the subclones would be present
                # find highest ranked conflicting mutation pattern
                # hcmp = max(set(self.cf_graph.neighbors(related_mp)).difference(self.conflicting_nodes),
                #            key=lambda k: self.cf_graph.node[k]['weight'])
                # print('Highest ranked conflicting mutation pattern for {} (w: {:.1f}): {} (w: {:.1f})'.format(
                #      rel_cn, cf_graph.node[rel_cn]['weight'], hcmp, cf_graph.node[hcmp]['weight']))

                # if (related_mp.intersection(hcmp)).issubset(created_scss):
                # logger.debug('New subclone might also be present in this mutation pattern: {} (w: {:.1e})'.format(
                #              ', '.join(self.patient.sc_names[sc] for sc in related_mp), self.node_scores[related_mp]))
                updated_nodes = self._update_mutation_pattern(related_mp, created_scss, updated_nodes)
                # else:
                #     if self.node_scores[related_mp] > min_score/10.0:
                #         logger.debug('Generating the putative subclones in MP {} (w: {:.1e}) '.format(
                #             related_mp, self.node_scores[related_mp])
                #             + 'would not make it compatible to {} (w: {:.1e})'.format(hcmp, self.node_scores[hcmp]))

            # step (e): add new subclone to all compatible supersets of mp
            for anc in opt_sol.compatible_nodes:

                if anc.issuperset(mp) and anc != mp:
                    anc_sas = [sa for sa in anc]
                    anc_sas.append(sc_sa_idx)
                    new_anc = frozenset(anc_sas)

                    # update subclones in patient
                    self.node_scores[new_anc] = self.node_scores[anc]
                    del self.node_scores[anc]
                    self.idx_to_mp[self.mp_col_ids[anc]] = new_anc
                    self.mp_col_ids[new_anc] = self.mp_col_ids[anc]
                    del self.mp_col_ids[anc]
                    updated_nodes[anc] = new_anc

            break

        return updated_nodes, self.sc_sample_ids

    def _update_mutation_pattern(self, old_mp, created_scs, updated_nodes):
        """
        Update previously incompatible MPs with newly identified SCs
        :param old_mp: previous mutation pattern
        :param created_scs: newly created SCs
        :param updated_nodes:
        :return mapping from old to updated MPs
        """

        new_mp = set(old_mp)
        for new_sc in created_scs.keys():
            if new_sc in new_mp:
                new_mp.remove(new_sc)
                new_mp.add(created_scs[new_sc])

        new_mp = frozenset(new_mp)

        # update subclones in patient
        self.node_scores[new_mp] = self.node_scores[old_mp]
        del self.node_scores[old_mp]
        self.idx_to_mp[self.mp_col_ids[old_mp]] = new_mp
        self.mp_col_ids[new_mp] = self.mp_col_ids[old_mp]
        del self.mp_col_ids[old_mp]
        updated_nodes[old_mp] = new_mp
        # logger.debug('Replaced samples {} and updated mutation pattern (w: {:.1e}) to {}'.format(
        #              ', '.join(self.patient.sc_names[new_sc] for new_sc in created_scs.keys() if new_sc in old_mp),
        #              self.node_scores[new_mp], ', '.join(self.patient.sc_names[sc] for sc in new_mp)))

        return updated_nodes


def infer_ml_graph_nodes(log_p01, sample_names, mut_keys, gene_names=None, max_no_mps=None):
    """
    Infer maximum likelihood using bayesian inference for each possible mutation pattern
    :param log_p01: posterior: log probability that VAF = 0, log probability that VAF > 0
    :param sample_names:
    :param mut_keys: list with information about the variant
    :param gene_names: list with the names of the genes in which the variant occurred
    :param max_no_mps: maximal number of MPs per variant that are considered in the MILP; by default the full solution
                       space is considered and hence 2^(#samples) of MPs are generated
    :return dictionary of nodes and corresponding variants, pattern reliability scores, weights of patterns of variants
    """

    assert max_no_mps is None or max_no_mps > 0, 'At least one mutation pattern per variant has to be considered'

    n = len(sample_names)  # number of samples
    m = len(log_p01)       # number of variants

    # presence probability of a variant for calculating reliability score
    # is upper bounded because the same variant could have been independently acquired twice
    max_pre_llh = math.log(def_sets.MAX_PRE_PROB)
    not_max_pre_llh = 1.0 - max_pre_llh

    # absence probability of a variant for calculating reliability score
    # should be upper bounded because the variant could have been lost by LOH
    # for most sequencing depth this lower bound is irrelevant
    max_abs_llh = math.log(def_sets.MAX_ABS_PROB)
    not_max_abs_llh = 1.0 - max_abs_llh

    node_scores = dict()    # mutation patterns score summed over all variants

    # weight per inferred mutation pattern per variant given the p0's and p1's in each sample for a variant
    mp_weights = list()
    # mutation pattern to the corresponding column id in the weight matrix
    mp_col_ids = dict()
    # column id to corresponding mutation pattern
    idx_to_mp = list()

    trunk_mp = frozenset([sa_idx for sa_idx in range(n)])

    mp_idx = 0
    if max_no_mps is None:
        # generate all possible mutation patterns for <n> given samples and index them
        for no_pres_vars in range(0, n+1):     # number of present variants in the generated MPs (mutation patterns)
            # generate all mps with length no_pres_vars
            for mp in combinations(range(n), no_pres_vars):
                node = frozenset(mp)        # create mutation pattern
                idx_to_mp.append(node)
                mp_col_ids[node] = mp_idx

                mp_idx += 1
    else:
        # generate the <max_no_mps> most likely mutation patterns of each variant
        for mut_idx in range(m):
            mlog_pre_probs = list()
            mlog_abs_probs = list()
            for sa_idx in range(n):
                mlog_pre_probs.append(-min(log_p01[mut_idx][sa_idx][1], max_pre_llh))
                mlog_abs_probs.append(-min(log_p01[mut_idx][sa_idx][0], max_abs_llh))

            for i, (mlog_prob, ml_mp, not_flipped_sas) in enumerate(
                    _get_ml_mps(n, mlog_pre_probs, mlog_abs_probs, max_no_mps)):

                # check if this pattern already belongs to one of the most likely one for another variant
                if ml_mp not in mp_col_ids:
                    idx_to_mp.append(ml_mp)
                    mp_col_ids[ml_mp] = mp_idx
                    mp_idx += 1

                # logger.debug('{} {}th: {} {:.3f}'.format(gene_names[mut_idx], i + 1, ml_mp, math.exp(-mlog_prob)))

        logger.info('Only {} mutation patterns will be explored in total.'.format(mp_idx))

    # calculate the reliability scores for all as above as relevant determined mutation patterns
    for mut_idx in range(m):

        mp_weights.append(dict())

        for mp_idx, node in enumerate(idx_to_mp):

            log_ml = 0.0                # log maximum likelihood of the inferred pattern
            for sa_idx in node:                         # variant is present
                log_ml += min(log_p01[mut_idx][sa_idx][1], max_pre_llh)
            for sa_idx in trunk_mp.difference(node):    # variant is absent
                log_ml += min(log_p01[mut_idx][sa_idx][0], max_abs_llh)

            if log_ml == 0.0:   # numerical artifact, use approximation: ignore second order term
                for sa_idx in node:                         # variant is present
                    log_ml += math.exp(max(log_p01[mut_idx][sa_idx][0], not_max_pre_llh))     # sum probability
                for sa_idx in trunk_mp.difference(node):    # variant is absent
                    log_ml += math.exp(max(log_p01[mut_idx][sa_idx][1], not_max_abs_llh))

                log_ml = np.log1p(-log_ml)      # calculates log(1+argument)
                if gene_names is not None:
                    logger.debug('Approximated log probability of variant {} having pattern {} by {:.2e}.'.format(
                        gene_names[mut_idx], node, log_ml))
                else:
                    logger.debug('Approximated log probability of variant {} having pattern {} by {:.2e}.'.format(
                        mut_keys[mut_idx], node, log_ml))
                if log_ml == 0.0:
                    if len(node) == 0 or len(node) == 1 or len(node) == n:
                        logger.debug('Underflow warning. Set probability to minimal float value!')
                    else:
                        logger.warning('Underflow error. Set probability to minimal float value!')
                    log_ml = -200

                assert log_ml < 0.0, ('Underflow error while calculating the probability that the ' +
                                      'variant {} does not have pattern {}.'.format(
                                       mut_keys[mut_idx], ', '.join(sample_names[sa_idx] for sa_idx in node)))

            # if max_no_mps is not None:  # not the full solution space is explored
            #     # use heapq to keep track of the most likely <max_no_mps> of mutation patterns
            #     if len(heap) < max_no_mps:
            #         heapq.heappush(heap, (log_ml, mp_idx))
            #     elif log_ml > heap[0][0]:      # likelihood of currently considered MP is higher than smallest in heap
            #         # log likelihood of mp that is more likely for the currently considered variant
            #         heapq.heapreplace(heap, (log_ml, mp_idx))
            # else:                       # full solution space is explored, weight of every pattern is relevant
            # assign calculated log probability that this variant has this mutation pattern
            mp_weights[mut_idx][mp_idx] = log_ml

        # if max_no_mps is not None:          # most likely MPs for this variant if the solution space is limited
        #     # assign calculated log probability that this variant has this mutation pattern
        #     for log_ml, mp_idx in heap:
        #         mp_weights[mut_idx][mp_idx] = log_ml

        # run through all relevant MPs for this variant and sum their log likelihoods
        # to calculate the reliability scores
        for mp_idx, log_ml in mp_weights[mut_idx].items():
            node = idx_to_mp[mp_idx]
            # calculate the probability of a mp that no variant has this mutation pattern
            # product of (1 - the probability that a variant has this mp)
            if node in node_scores.keys():
                node_scores[node] -= math.log(-math.expm1(log_ml))
            else:
                node_scores[node] = -math.log(-math.expm1(log_ml))

    for mp_idx, node in enumerate(idx_to_mp):
        if node in node_scores and node_scores[node] == 0.0:
            if len(node) == 0 or len(node) == 1 or len(node) == len(sample_names):
                # logger.debug('Underflow warning for pattern {}. Set probability to minimal float value!'
                #       .format(', '.join(sample_names[sa_idx] for sa_idx in node)))
                pass
            else:
                # logger.warn('Underflow error for pattern {}. Set probability to minimal float value!'.format(
                #     ', '.join(sample_names[sa_idx] for sa_idx in node)))
                # raise RuntimeError(
                #     'Underflow error for pattern {}. Set probability to minimal float value!'.format(
                #         ', '.join(sample_names[sa_idx] for sa_idx in node)))
                pass

            node_scores[node] = sys.float_info.min

        # logger.debug('Variant {} has pattern {} with probability {:.1e}.'.format(
        #     gene_names[mut_idx], ', '.join(sample_names[sa_idx] for sa_idx in node), math.exp(log_ml)))

    # normalize reliability score by the number of processed variants (m)
    for node in node_scores.keys():
        node_scores[node] /= m
        if node_scores[node] == 0.0:
            node_scores[node] = sys.float_info.min

    # Show nodes with highest reliability score
    for node, score in itertools.islice(sorted(node_scores.items(), key=lambda k: -k[1]), 0, 30):
        logger.debug('Pattern {} has a normalized reliability score of {:.2e}.'.format(node, score))

    return node_scores, idx_to_mp, mp_col_ids, mp_weights


def create_conflict_graph(reliability_scores):
    """
    Create a graph where the nodes are given by the mutation patterns and
    the edges model the evolutionary conflicts among them
    :param reliability_scores: each node in the graph is weighted corresponding to the confidence
    in the sequencing data of the mutation modeled by the reliability scores
    :return conflict graph
    """

    # create empty conflict graph
    cf_graph = nx.Graph()
    incompatible_mp = defaultdict(set)

    for node, score in reliability_scores.items():
        if 1 < len(node):
            cf_graph.add_node(node, weight=score)

    # since each mutation occurs at most once running through this is
    # in O(m^2 * n) = O(|mutations|*|mutations|*|samples|) as
    # the number of clones is bounded by the number of distinct mutations
    for (node1, node2) in itertools.combinations(reliability_scores.keys(), 2):

        # variant need to appear at least in two samples and at most in n-1 samples
        # otherwise a conflict is not possible
        if not cf_graph.has_node(node1) or not cf_graph.has_node(node2):
            continue

        if len(node1.intersection(node2)) > 0:     # at least one sample where both characters are present (11)

            # check if some characters are present in one clone but not in the other and visa versa
            if len(node1.difference(node2)) > 0 and len(node2.difference(node1)) > 0:   # check for 01 and 10
                # => conflict exists among characters of mp 1 and mp 2
                # evolutionary incompatible mutation patterns

                incompatible_mp[node1].add(node2)
                incompatible_mp[node2].add(node1)

                # add edge between conflicting clones
                cf_graph.add_edge(node1, node2)

    logger.info('Created conflict graph with {} nodes of weight {:.2f} and {} evolutionary conflicts.'.format(
        cf_graph.order(), sum(data['weight'] for _, data in cf_graph.nodes_iter(data=True)), cf_graph.size()))

    return cf_graph


def _get_ml_mps(n, mlog_pre_probs, mlog_abs_probs, k):
    """
    Return k most likely mps of a variant where its presence
    probability in each sample is given in log_pre_probs and absence probability in log_abs_probs
    :param n: number of samples
    :param mlog_pre_probs: list of -log probabilities to be present in a sample
    :param mlog_abs_probs: list of -log probabilities to be absent in a sample
    :param k: number of most likely mutation patterns to be returned
    """
    heap = []
    explored_mps = set()

    # most likely pattern is given by the most likely classification in each sample
    mp = []
    mlog_prob = 0.0
    for sa_idx in range(n):
        if mlog_pre_probs[sa_idx] < mlog_abs_probs[sa_idx]:  # variant more likely to be present in sample sa_idx
            mlog_prob += mlog_pre_probs[sa_idx]
            mp.append(sa_idx)
        else:
            mlog_prob += mlog_abs_probs[sa_idx]

    heapq.heappush(heap, (mlog_prob, frozenset(mp), set([i for i in range(n)])))
    explored_mps.add(frozenset(mp))
    # logger.debug('Pushed ml pattern {}: {:.3f}'.format(mp, math.exp(-mlog_prob)))

    # return the k most likely mps
    for _ in range(k):
        yield _next_ml_mp(heap, mlog_pre_probs, mlog_abs_probs, explored_mps, n)


def _next_ml_mp(heap, mlog_pre_probs, mlog_abs_probs, explored_mps, n):
    """
    Get next less likely mutation pattern and further explore even less likely patterns
    :param heap:
    :param mlog_pre_probs: list of -log probabilities to be present in a sample
    :param mlog_abs_probs: list of -log probabilities to be absent in a sample
    :param explored_mps: set of mutation patterns that were already explored
    :param n: number of samples
    :return: next less likely mutation pattern
    """

    ml_p, ml_mp, not_flipped_sas = heapq.heappop(heap)
    # logger.debug('Popped {}: {:.3f}'.format(ml_mp, math.exp(-ml_p)))

    # investigate flipping individually all the not yet flipped samples
    # and put them on the heap (equivalent to breath-first-search)
    for sa_to_flip in not_flipped_sas:
        nfs = set(not_flipped_sas)
        nfs.remove(sa_to_flip)
        mlog_prob = 0.0
        mp = []
        for sa_idx in range(n):
            if sa_idx in nfs:
                if mlog_pre_probs[sa_idx] < mlog_abs_probs[sa_idx]:
                    mlog_prob += mlog_pre_probs[sa_idx]
                    mp.append(sa_idx)
                else:
                    mlog_prob += mlog_abs_probs[sa_idx]
            else:  # less likely classification is used in this sample
                if mlog_pre_probs[sa_idx] < mlog_abs_probs[sa_idx]:
                    mlog_prob += mlog_abs_probs[sa_idx]
                else:
                    mlog_prob += mlog_pre_probs[sa_idx]
                    mp.append(sa_idx)

        if frozenset(mp) not in explored_mps:
            heapq.heappush(heap, (mlog_prob, frozenset(mp), nfs))
            explored_mps.add(frozenset(mp))
            # logger.debug('Pushed {}: {:.3f}'.format(mp, mlog_prob))

    return ml_p, ml_mp, not_flipped_sas


def _subsets(mp):
    """
    Returns all subsets (descendant) mutation patterns of a given mutation pattern
    :param mp: ancestor mp
    :return: descendant mp
    """

    for l in range(len(mp)-1, 0, -1):
        for descendant_mp in itertools.combinations(mp, l):
            yield frozenset(descendant_mp)
