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

    def __init__(self, patient, mps):

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

        # confidence in each parsimony-informative branching
        # reliability score of node divided by its own score plus the sum of all evolutionarily incompatible patterns
        # (= sum of reliability scores of neighbors in the conflict graph)
        self.branch_confidence = None

        self.bootstrapping_values = None

        # false positives and false negatives compared to original classification
        # likely technical or biological artifacts in the data
        self.false_positives = None
        self.false_negatives = None
        self.false_negative_unknowns = None

        # updated clones after artifacts in the data have been removed
        self.shared_mlh_mps = None
        self.mlh_founders = None
        self.mlh_unique_mutations = None
        self.mlh_absent_mutations = None

        self.mlh_tree = None

        # calculate the minimal number of variant reads k_min at the median coverage and purity such that p_1 > 50%
        k_mins = []
        called_ps = []
        missed_ps = []
        pres_lp = math.log(0.5)
        lh = 1.0
        for sa_idx, sample_name in enumerate(patient.sample_names):
            for k in range(500):
                _, p1 = get_log_p0(np.median(patient.sample_coverages[sample_name]), k, self.patient.bi_error_rate,
                                   self.patient.bi_c0, pseudo_alpha=def_sets.PSEUDO_ALPHA,
                                   pseudo_beta=patient.betas[sample_name])
                if p1 > pres_lp:
                    k_mins.append(k)
                    break
            logger.debug('{}: Minimum number of mutant reads such that presence probability is greater than 50%: {}.'
                         .format(sample_name, k_mins[-1]))
            called_ps.append(1.0 - binom.cdf(k_mins[-1]-1, np.median(patient.sample_coverages[sample_name]),
                                             self.patient.bi_error_rate))
            logger.debug('Probability to observe an incorrectly called variant: {:.3%}'.format(called_ps[-1]))

            if sample_name in patient.estimated_purities:
                missed_ps.append(binom.cdf(k_mins[-1]-1, np.median(patient.sample_coverages[sample_name]),
                                           self.patient.estimated_purities[sample_name] / 2.0))
            else:
                missed_ps.append(binom.cdf(k_mins[-1]-1, np.median(patient.sample_coverages[sample_name]),
                                           np.median(self.patient.sample_mafs[sample_name])))
            logger.debug('Probability to miss a clonal variant: {:.1e}'.format(missed_ps[-1]))

            # probability that all calls are correct
            lh *= 1.0 - called_ps[-1] - missed_ps[-1]

        # probability that at least one call is wrong
        lh = 1 - lh
        self.min_score = -math.log(1.0 - lh)
        logger.debug('Likelihood of a pattern with at least one false-positive or false-negative: {:.3e}'.format(lh))
        logger.info('Minimum reliability score value to be considered as a potential subclone: {:.5f}'.format(
            self.min_score))

    def infer_max_lh_tree(self, subclone_detection=False, no_bootstrap_samples=0, max_no_mps=None):
        """
        Infer maximum likelihood tree via calculation reliability scores for each
        possible mutation pattern from the likelihood that no variant has this pattern
        The inferred solution represents the reliable and evolutionary compatible mutation patterns
        The mutation pattern of each variant is given by the mp maximizing its likelihood
        :param subclone_detection: is subclone detection enabled?
        :param max_no_mps: only the given maximal number of most likely (by joint likelihood) mutation patterns
            is explored per variant
        :param no_bootstrap_samples: number of samples with replacement for the bootstrapping
        :return inferred evolutionary tree
        """

        if subclone_detection:
            self.patient.sc_names = deepcopy(self.patient.sample_names)
        else:
            self.patient.sc_names = None

        # necessary to map from identified putative subclones to their original sample
        self.sc_sample_ids = dict()

        if max_no_mps is not None and max_no_mps < math.pow(2, self.patient.n):
            logger.warn('Some variants might be evolutionarily incompatible since the '
                        'solution space is only partially explored!')
            self.conflicting_mutations = set()

        # compute various mutation patterns (nodes) and their reliability scores
        self.node_scores, self.idx_to_mp, self.mp_col_ids, self.mp_weights = infer_ml_graph_nodes(
            self.patient.log_p01, self.patient.sample_names, self.patient.mut_keys, gene_names=self.patient.gene_names,
            max_no_mps=max_no_mps)

        while True:
            # create conflict graph which forms the input to the ILP
            self.cf_graph = create_conflict_graph(self.node_scores)

            # translate the conflict graph into a minimum vertex cover problem
            # and solve this using integer linear programming
            self.conflicting_nodes, self.compatible_nodes = cps.solve_conflicting_phylogeny(self.cf_graph)

            # # unused measure since it is not invariant to the number of samples
            # self.branch_confidence = dict()
            # # calculate confidence value for each branching
            # for node in self.compatible_nodes:
            #     # consider only parsimony-informative MPs
            #     if len(node) <= 1 or len(node) == len(self.patient.sample_names):
            #         continue
            #     self.branch_confidence[node] = \
            #         ((1.0 - math.exp(-self.node_scores[node])) *
            #          np.prod(np.array([math.exp(-self.node_scores[v]) for v in self.cf_graph.neighbors(node)])))
            #     self.branch_confidence[node] /= \
            #         (self.branch_confidence[node] + (math.exp(-self.node_scores[node]) *
            #          (1.0-np.prod(np.array([math.exp(-self.node_scores[v]) for v in self.cf_graph.neighbors(node)])))))
            #     # logger.info('Branching confidence of {}: {:.2e}'.format(node, self.branch_confidence[node]))

            # ##### assign each variant to the highest ranked evolutionarily compatible mutation pattern ########

            # maps from each evolutionarily compatible MP to the selected set of variants
            self.max_lh_nodes = defaultdict(set)
            # map from mut_idx to the highest ranked evolutionarily compatible MP
            self.max_lh_mutations = dict()
            # map from mut_idx to the weight of the highest ranked evolutionarily compatible MP
            self.max_lh_weights = dict()

            # dictionaries from mut_idx to putative artifacts
            self.false_positives = defaultdict(set)
            self.false_negatives = defaultdict(set)
            self.false_negative_unknowns = defaultdict(set)

            # find the highest ranked evolutionarily compatible mutation pattern for each variant
            for mut_idx in range(len(self.mp_weights)):

                # sort by descending log likelihood
                for sol_idx, (mp_col_idx, _) in enumerate(
                        sorted(self.mp_weights[mut_idx].items(), key=lambda k: -k[1]), 1):

                    if self.idx_to_mp[mp_col_idx] in self.compatible_nodes:
                        # found most likely pattern for this variant
                        self.max_lh_nodes[self.idx_to_mp[mp_col_idx]].add(mut_idx)
                        self.max_lh_mutations[mut_idx] = self.idx_to_mp[mp_col_idx]
                        self.max_lh_weights[mut_idx] = self.mp_weights[mut_idx][mp_col_idx]
                        # logger.debug('Chose {}th highest ranked pattern.'.format(sol_idx)
                        #     + 'Max LH pattern of variant in {} is {} with log likelihood {:.1e}.'.format(
                        #       self.patient.gene_names[mut_idx] if self.patient.gene_names is not None
                        #       else self.patient.mut_keys[mut_idx], self.idx_to_mp[mp_col_idx],
                        #       self.mp_weights[mut_idx][mp_col_idx]))

                        # determine false positives and false negatives compared to original classification
                        fps = set(self.patient.mutations[mut_idx].difference(self.idx_to_mp[mp_col_idx]))

                        # check if some of these false-positives are present in the newly created subclones
                        for sc_idx, sa_idx in self.sc_sample_ids.items():
                            if sc_idx in self.idx_to_mp[mp_col_idx] and sa_idx in fps:
                                fps.remove(sa_idx)
                        if len(fps) > 0:
                            self.false_positives[mut_idx] = fps

                        # distinguish between real false negatives and variants classified as unknown
                        for sc_idx in self.idx_to_mp[mp_col_idx].difference(self.patient.mutations[mut_idx]):

                            # map from identified putative subclones to their original sample
                            if sc_idx in self.sc_sample_ids.keys():
                                while sc_idx in self.sc_sample_ids.keys():
                                    sc_idx = self.sc_sample_ids[sc_idx]
                                if sc_idx in self.patient.mutations[mut_idx]:
                                    # mutation was already classified as present in original sample
                                    # => no false-negative
                                    continue
                                # mutation was not classified as present in original sample
                                # => must be a false-negative
                                if self.patient.data[mut_idx][sc_idx] < 0.0:      # unknown classified mutation
                                    self.false_negative_unknowns[mut_idx].add(sc_idx)
                                else:
                                    self.false_negatives[mut_idx].add(sc_idx)
                            else:
                                if self.patient.data[mut_idx][sc_idx] < 0.0:      # unknown classified mutation
                                    self.false_negative_unknowns[mut_idx].add(sc_idx)
                                else:
                                    self.false_negatives[mut_idx].add(sc_idx)

                        # search can be stopped since most likely pattern has been found for this variant
                        break

                # no evolutionarily compatible mutation pattern was among the <max_no_mps> most likely pattern
                # of this variant; variant will be incompatible to the inferred tree!
                else:

                    # check if indeed the solution space is only partially explored
                    assert max_no_mps is not None and max_no_mps < math.pow(2, self.patient.n), \
                        'Compatible mutation pattern must exist when the full solution space is explored!'

                    self.conflicting_mutations.add(mut_idx)

            logger.info('Putative false-positives {}, putative false-negatives {}, put. false neg. unknowns {}'.format(
                sum(len(fps) for mut_idx, fps in self.false_positives.items()),
                sum(len(fns) for mut_idx, fns in self.false_negatives.items()),
                sum(len(fns) for mut_idx, fns in self.false_negative_unknowns.items())))
            if self.conflicting_mutations is not None and len(self.conflicting_mutations) > 0:
                logger.warn('{} evolutionarily incompatible variants occurred{}'.format(
                    len(self.conflicting_mutations),
                    ' in: '+', '.join(self.patient.gene_names[mut_idx] for mut_idx in self.conflicting_mutations) if
                    self.patient.gene_names is not None else '.'))

            # find parsimony-informative evolutionarily incompatible mps with high likelihood
            if subclone_detection and max_no_mps is None:

                self.sc_sample_ids, found_subclones = self.find_subclones()
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
                len(self.max_lh_mutations), len(self.max_lh_mutations)+len(self.conflicting_mutations)) +
                '{} ({:.3%}) of the mutations are conflicting.'.format(
                len(self.conflicting_mutations), float(len(self.conflicting_mutations)) /
                    (len(self.max_lh_mutations)+len(self.conflicting_mutations))))

        # find founders and unique mutations in updated mutation list
        # build up a dictionary of resolved clones where
        # the likely sequencing errors have been updated
        self.mlh_founders = set()                       # set of founding mutations present in all samples
        self.mlh_unique_mutations = defaultdict(set)    # private mutations only appear in the leaves
        self.mlh_absent_mutations = set()               # mutations inferred to be absent in all samples
        self.shared_mlh_mps = defaultdict(set)          # parsimony-informative mps

        for mut_idx, samples in self.max_lh_mutations.items():
            if len(samples) == len(self.patient.sample_names) + len(self.sc_sample_ids):          # founder mut.
                self.mlh_founders.add(mut_idx)
            elif 1 < len(samples) < len(self.patient.sample_names) + len(self.sc_sample_ids):     # shared mut.
                self.shared_mlh_mps[frozenset(samples)].add(mut_idx)
            elif len(samples) == 1:                                                          # unique mut.
                for sa_idx in samples:
                    self.mlh_unique_mutations[sa_idx].add(mut_idx)
            else:
                self.mlh_absent_mutations.add(mut_idx)

        # compute the number of persistent and present mutations inferred
        # in the evolutionary trajectory of each sample
        for sa_idx in range(len(self.patient.sample_names)):
            for mut_idx in self.mlh_founders:
                self.patient.variants[sa_idx].append(mut_idx)
            for mut_idx in self.mlh_unique_mutations[sa_idx]:
                self.patient.variants[sa_idx].append(mut_idx)

        for mps, muts in self.shared_mlh_mps.items():
            for sa_idx in mps:
                for mut_idx in muts:
                    self.patient.variants[sa_idx].append(mut_idx)

        for sa_idx, sample_name in enumerate(self.patient.sample_names):
            logger.info('Inferred number of variants in sample {} in a perfect and persistent phylogeny: {}'.format(
                sample_name, len(self.patient.variants[sa_idx])))

        # is bootstrapping enabled?
        if no_bootstrap_samples > 0:
            if subclone_detection is not None and 0 < subclone_detection < 1:
                logger.error('Bootstrapping analysis does not support subclone detection! ')
                logger.error('Go to settings.py and set MIN_MP_LH to 1.')
            else:
                self.do_bootstrapping(no_bootstrap_samples)

        # construct a phylogenetic tree from the maximum likelihood mutation patterns
        self.mlh_tree = self.infer_evolutionary_tree(self.shared_mlh_mps, self.mlh_founders,
                                                     self.mlh_unique_mutations, confidence=self.bootstrapping_values)

        return self.mlh_tree

    def do_bootstrapping(self, no_samples):
        """
        Validate the robustness of the identified most reliable mutation patterns through bootstrapping
        :param no_samples: Number of samples with replacement for the bootstrapping
        """

        # Most reliable mutation patterns need to be identified before
        if self.compatible_nodes is None:
            self.infer_max_lh_tree(subclone_detection=False)

        node_frequencies = cps.bootstrapping_solving(
            self.cf_graph, self.mp_weights, self.idx_to_mp, no_samples)

        self.bootstrapping_values = dict()

        # calculate the frequency with which compatible mutation patterns are reproduced
        # when bootstrapping is used
        for node in self.compatible_nodes:

            # consider only parsimony-informative MPs
            if len(node) <= 1 or len(node) == len(self.patient.sample_names):
                continue

            if node in node_frequencies:
                self.bootstrapping_values[node] = float(node_frequencies[node]) / no_samples
            else:
                self.bootstrapping_values[node] = 0.0

            logger.info('Bootstrapping value for node {}: {:.1%}'.format(node, self.bootstrapping_values[node]))

    def validate_node_robustness(self, no_replications):
        """
        Validate the robustness of the identified most reliable mutation patterns
        through down-sampling
        :param no_replications: Number of replications for each used fraction of variants
        :return observed frequencies of the mutation patterns in each used fraction of variants
        """

        # Most reliable mutation patterns need to be identified before
        if self.compatible_nodes is None:
            self.infer_max_lh_tree(subclone_detection=False)

        node_frequencies = cps.solve_downsampled_nodes(
            self.cf_graph, self.mp_weights, self.idx_to_mp, no_replications)

        comp_node_frequencies = defaultdict(dict)

        # calculate the frequency with which compatible mutation patterns are reproduced
        # when only a subset of variants are used
        no_comp_pars_inf_mps = sum(1 for node in self.compatible_nodes
                                   if 1 < len(node) < len(self.patient.sample_names))
        for sample_fraction in sorted(node_frequencies.keys()):
            freq = 0
            for node in self.compatible_nodes:

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

    def find_subclones(self):

        updated_nodes, self.sc_sample_ids = self.find_subclonal_mps()

        if len(updated_nodes) == 0:
            logger.info('There are no more incompatible mutation patterns with a reliability score '
                        'of at least {:.1e}.'.format(self.min_score))
            return self.sc_sample_ids, False
        else:
            # anc = set()     # find the mutations present in the ancestor MP and assign them also to the new MP
            # for old_mp, new_mp in updated_nodes.items():
            #
            #     if old_mp.issuperset(anc):
            #         anc = old_mp
            #     elif (len(old_mp.intersection(anc)) > 0 and
            #             len(old_mp.difference(anc)) > 0 and len(anc.difference(old_mp)) > 0):
            #         logger.warn('Putative subclones are evolutionary incompatible '
            #                     'to each other: {} <> {}'.format(anc, old_mp))
            #
            # logger.debug('{} is the parental mutation pattern of the putative subclones {}.'.format(
            #     anc, updated_nodes.keys()))
            #
            # # the mutations present in the direct ancestor of the identified
            # # parent subclone should be present in all
            # parental_mps = set()
            # for sc in self.compatible_nodes:
            #     if sc.issuperset(anc):
            #         parental_mps.add(sc)
            #         logger.debug('Mutations present in {} should also be present in all putative subclones.'
            #                      .format(sc))
            #
            # # account for putative subclones in parental mutation patterns
            # for a in parental_mps:
            #     # mutations of parental clone should also be present in new putative subclones
            #     par_mp = a.union(updated_nodes[anc])
            #
            #     # update subclones in patient
            #     self.node_scores[par_mp] = self.node_scores[a]
            #     del self.node_scores[a]
            #     self.idx_to_mp[self.mp_col_ids[a]] = par_mp
            #     self.mp_col_ids[par_mp] = self.mp_col_ids[a]
            #     del self.mp_col_ids[a]
            #
            #     logger.info('Updated parental mutation pattern from {} to {}'.format(
            #         ', '.join(self.patient.sc_names[sc_idx] for sc_idx in a),
            #         ', '.join(self.patient.sc_names[sc_idx] for sc_idx in par_mp)))

            return self.sc_sample_ids, True

    def find_subclonal_mps(self):
        """
        Identify evolutionarily incompatible mutation patterns with a reliability score greater than min_score which
        where variants in some samples are subclonal
        :return updated_nodes:
        """

        updated_nodes = dict()

        # find incompatible mutation pattern with the highest reliability score
        for mp in sorted(self.conflicting_nodes, key=lambda k: -self.node_scores[k]):

            if self.node_scores[mp] < self.min_score:
                # no incompatible mutation patterns with high reliability score exist
                return updated_nodes, self.sc_sample_ids

            # determine the number of variants for which this mutation pattern would be
            # more likely than the before assigned evolutionarily compatible mp
            no_sup_vars = sum(1 for mut_idx in range(len(self.mp_weights))
                              if self.mp_weights[mut_idx][self.mp_col_ids[mp]] > self.max_lh_weights[mut_idx])

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
            hcmp = max(set(self.cf_graph.neighbors(mp)).difference(self.conflicting_nodes),
                       key=lambda k: self.node_scores[k])
            logger.info('Highest ranked evolutionary incompatible mutation pattern of MP {} (w: {:.2f}): {} (w: {:.2f})'
                        .format(', '.join(self.patient.sc_names[sc] for sc in mp), self.cf_graph.node[mp]['weight'],
                                ', '.join(self.patient.sc_names[sc] for sc in hcmp),
                                self.cf_graph.node[hcmp]['weight']))

            # infer samples with putative mixed subclones
            msc_samples = list(mp.intersection(hcmp))
            if len(msc_samples) > 1:    # more than one subclone would be needed to generate, reconsider later
                continue

            sc_sa = msc_samples[0]      # sample where subclone is created

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
            logger.info('Created new subclones {}'.format(
                        ', '.join(self.patient.sc_names[sc] for sc in new_mp.difference(mp))))

            # step (d) continued: check compatibility with other conflicting clones
            for related_mp in sorted(self.conflicting_nodes, key=lambda k: -self.cf_graph.node[k]['weight']):
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
            for anc in self.compatible_nodes:

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

    # generate all possible mutation patterns for <n> given samples and index them
    mp_idx = 0
    for no_pres_vars in range(0, n+1):     # number of present variants in the generated MPs (mutation patterns)
        # generate all mps with length no_pres_vars
        for mp in combinations(range(n), no_pres_vars):
            node = frozenset(mp)        # create mutation pattern
            idx_to_mp.append(node)
            mp_col_ids[node] = mp_idx

            mp_idx += 1

    # run through all mutations and generate either the <max_no_mps> of mutation patterns for each variant or
    # generate each possible mutation pattern for each variant
    for mut_idx in range(m):

        mp_weights.append(dict())

        if max_no_mps is not None:  # not the full solution space is explored
            # generate heap to keep track of the most likely <max_no_mps> of mutation patterns
            heap = list()   # element at heap[0] is always the minimal element hence lowest log likelihood

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
                        logger.warn('Underflow error. Set probability to minimal float value!')
                    log_ml = -200

                assert log_ml < 0.0, ('Underflow error while calculating the probability that the ' +
                                      'variant {} does not have pattern {}.'.format(
                                       mut_keys[mut_idx], ', '.join(sample_names[sa_idx] for sa_idx in node)))

            if max_no_mps is not None:  # not the full solution space is explored
                # use heapq to keep track of the most likely <max_no_mps> of mutation patterns
                if len(heap) < max_no_mps:
                    heapq.heappush(heap, (log_ml, mp_idx))
                elif log_ml > heap[0][0]:      # likelihood of currently considered MP is higher than smallest in heap
                    # log likelihood of mp that is more likely for the currently considered variant
                    heapq.heapreplace(heap, (log_ml, mp_idx))
            else:                       # full solution space is explored, weight of every pattern is relevant
                # assign calculated log probability that this variant has this mutation pattern
                mp_weights[mut_idx][mp_idx] = log_ml

        if max_no_mps is not None:          # most likely MPs for this variant if the solution space is limited
            # assign calculated log probability that this variant has this mutation pattern
            for log_ml, mp_idx in heap:
                mp_weights[mut_idx][mp_idx] = log_ml

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
    for node, score in itertools.islice(sorted(node_scores.items(), key=lambda k: -k[1]), 0, 50):
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
        cf_graph.add_node(node, weight=score)

    # since each mutation occurs at most once running through this is
    # in O(m^2 * n) = O(|mutations|*|mutations|*|samples|) as
    # the number of clones is bounded by the number of distinct mutations
    for (node1, node2) in itertools.combinations(reliability_scores.keys(), 2):

        # variant need to appear at least in two samples and at most in n-1 samples
        # otherwise a conflict is not possible
        if len(node1) < 2 or len(node2) < 2:
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


def _subsets(mp):
    """
    Returns all subsets (descendant) mutation patterns of a given mutation pattern
    :param mp: ancestor mp
    :return: descendant mp
    """

    for l in range(len(mp)-1, 0, -1):
        for descendant_mp in itertools.combinations(mp, l):
            yield frozenset(descendant_mp)
