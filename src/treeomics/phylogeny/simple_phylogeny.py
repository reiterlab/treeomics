#!/usr/bin/python
"""Infers most reliable mutation pattern based on standard binary classification"""
__author__ = 'Johannes REITER'
__date__ = 'April, 2014'

import logging
import math
import sys
import numpy as np
from collections import defaultdict
from itertools import islice
import phylogeny.cplex_solver as cps
from phylogeny.phylogeny_utils import Phylogeny, create_conflict_graph
# import conflict_solver as cf


# get logger for application
logger = logging.getLogger('treeomics')


class SimplePhylogeny(Phylogeny):

    def __init__(self, patient, mps):

        Phylogeny.__init__(self, patient, mps)

        # variables used for the minimum ignored mutation solution
        self.compatible_mutations = None
        self.conflicting_mutations = None

        self.nodes = None

        # perfect and persistent phylogeny tree according to the minimal number of ignored mutations
        self.compatible_tree = None

        # dictionary of the individual mutation pattern scores, dict key is given by the mut_idx
        self.mp_weights = None

    def find_max_compatible_tree(self):
        """
        Find the largest set of compatible characters and create from this largest
        set a phylogenetic tree.
        For variants which are not significantly mutated a minimum coverage can be required.
        If the minimum coverage is not met, the variant status is unknown and both patterns (presence and
        absence) will be explored.
        :return inferred evolutionary tree
        """

        # compute various mutation patterns (nodes) and their reliability scores
        self.nodes, self.node_scores, self.mp_weights = determine_graph_nodes(
            self.patient.log_p01, self.patient.sample_names, self.patient.mut_keys, self.patient.gene_names)

        # create conflict graph which forms the input to the ILP
        self.cf_graph = create_conflict_graph(self.nodes, weights=self.node_scores)

        # translate the conflict graph into a minimum vertex cover problem
        # and solve this using integer linear programming
        self.conflicting_nodes, _ = cps.solve_conflicting_phylogeny(self.cf_graph)

        self.compatible_nodes = dict()

        self.conflicting_mutations = set(mut for mp, muts in self.nodes.items() if len(mp) > 0 for mut in muts)
        self.compatible_mutations = [mut_idx for mut_idx in sorted(self.conflicting_mutations)]

        founders = set()  # find all founding mutations according to the inferred solution
        unique_mutations = defaultdict(set)  # find all private mutations for all samples

        # remove the minimum set of conflicting mutation patterns from the mutation lists
        # and compute the maximum compatible clones and mutation dictionaries
        for mp, muts in self.nodes.items():
            if mp not in self.conflicting_nodes and len(mp) > 0:
                self.compatible_nodes[mp] = muts
                if len(mp) == 1:      # unique (private) mutations
                    for sa_idx in mp:
                        unique_mutations[sa_idx] = muts
                elif len(mp) == self.patient.n:
                    founders = muts

        for compatible_clone, muts in self.compatible_nodes.items():
            for mut in muts:
                self.conflicting_mutations.remove(mut)
        for conflicting_mut in self.conflicting_mutations:
            self.compatible_mutations.remove(conflicting_mut)

        logger.info('{} out of {} are compatible on an evolutionary tree. '.format(
            len(self.compatible_mutations), len(self.compatible_mutations)+len(self.conflicting_mutations))
            + '{} ({:.3%}) of the mutations are conflicting.'.format(
            len(self.conflicting_mutations), float(len(self.conflicting_mutations)) /
                (len(self.compatible_mutations)+len(self.conflicting_mutations))))

        # construct a phylogenetic tree from the conflict-free set of clones
        self.compatible_tree = self.infer_evolutionary_tree(self.compatible_nodes, founders, unique_mutations)

        return self.compatible_tree

    def validate_node_robustness(self, no_replications):
        """
        Validate the robustness of the identified most reliable mutation patterns
        through down-sampling
        :param no_replications: Number of replications for each used fraction of variants
        :return observed frequencies of the mutation patterns in each used fraction of variants
        """

        # Most reliable mutation patterns need to be identified before
        if self.compatible_nodes is None:
            self.find_max_compatible_tree()

        node_frequencies = cps.solve_downsampled_binary_nodes(
            self.cf_graph, self.mp_weights, self.patient.shared_mutations,
            no_replications, len(self.patient.sample_names))

        comp_node_frequencies = defaultdict(dict)

        # infer the number of compatible parsimony-informative mutation patterns
        no_comp_nodes = sum(1 for node in self.compatible_nodes.keys()
                            if 1 < len(node) < len(self.patient.sample_names))
        for sample_fraction in sorted(node_frequencies.keys()):
            freq = 0
            for node in self.compatible_nodes.keys():

                # consider only parsimony-informative MPs
                if len(node) == 1 or len(node) == len(self.patient.sample_names):
                    continue

                if node in node_frequencies[sample_fraction]:
                    comp_node_frequencies[sample_fraction][node] = \
                        float(node_frequencies[sample_fraction][node]) / no_replications
                    freq += comp_node_frequencies[sample_fraction][node]
                else:
                    comp_node_frequencies[sample_fraction][node] = 0.0

            logger.info('Fraction of confirmed mutation patterns with {}% of the variants: {:.1%}'.format(
                sample_fraction, float(freq)/no_comp_nodes))

        # produce latex table for paper
        # displayed_fractions = [10, 20, 50, 80, 90, 95]
        # print('\\textbf Mutation pattern  & {} \\\\'.format(' & '.join(
        #       '\\textbf {:.0f}$\%$'.format(fr) for fr in displayed_fractions)))
        # for node in sorted(self.compatible_nodes.keys(), key=lambda k: -self.cf_graph.node[k]['weight']):
        #
        #     print('({}) & {} \\\\'.format(','.join(str(n+1) for n in sorted(node)), ' & '.join(
        #         '{:.1f}$\%$'.format(comp_node_frequencies[fr][node]*100.0) for fr in displayed_fractions)))

        return comp_node_frequencies


def determine_graph_nodes(log_p01, sample_names, mut_keys, gene_names=None):
    """
    Infer the mutation patterns (graph nodes) by the maximum likelihood using bayesian inference
    :param log_p01: posterior: log probability that VAF = 0
    :param sample_names: name of the samples
    :param mut_keys: unique key of variant
    :param gene_names: if available name of the gene where the variant occurred
    :return dictionary of nodes and corresponding variants, pattern reliability scores, weights of patterns of variants
    """

    # clones form the basis for the nodes in the graph
    nodes = defaultdict(set)
    node_scores = dict()    # mutation patterns score summed over all variants supporting this score
    mut_pattern_weights = dict()    # weight of the inferred mutation pattern given the p0 in each sample

    for mut_idx in range(len(log_p01)):
        node = set()
        log_ml = 0.0        # log maximum likelihood of the inferred pattern
        for sa_idx in range(len(sample_names)):
            if math.exp(log_p01[mut_idx][sa_idx][0]) < 0.5:   # variant is present
                node.add(sa_idx)
                log_ml += math.log(-np.expm1(log_p01[mut_idx][sa_idx][0]))
            else:                           # mutation is absent
                log_ml += log_p01[mut_idx][sa_idx][0]

        node = frozenset(node)
        if log_ml == 0.0:   # numerical artifact, use approximation: ignore second order term
            for sa_idx in range(len(sample_names)):
                if sa_idx in node:              # variant is present
                    log_ml += math.exp(log_p01[mut_idx][sa_idx][0])     # sum probability
                else:                           # variant is absent
                    log_ml += -np.expm1(log_p01[mut_idx][sa_idx][0])

            log_ml = np.log1p(-log_ml)      # calculates log(1+argument)
            logger.debug('Approximated log probability of variant {} having pattern {} by {:.2e}.'.format(
                gene_names[mut_idx] if gene_names is not None else mut_keys[mut_idx], node, log_ml))
            if log_ml == 0.0:
                if len(node) == 0 or len(node) == len(sample_names):
                    logger.debug('Underflow error. Set probability to minimal float value!')
                else:
                    logger.warn('Underflow error. Set probability to minimal float value!')
                log_ml = -sys.float_info.min
            assert log_ml < 0.0,\
                ('Underflow error while calculating the probability that the variant {} does not have pattern {}.'
                    .format(gene_names[mut_idx], ', '.join(sample_names[sa_idx] for sa_idx in node)))

        mut_pattern_weights[mut_idx] = log_ml

        nodes[node].add(mut_idx)
        # calculate the probability that no variant has pattern 'node'
        if node in node_scores.keys():
            node_scores[node] *= -np.expm1(log_ml)
        else:
            node_scores[node] = -np.expm1(log_ml)

        logger.debug('Variant {} {} has pattern {} with probability {:.1e}.'.format(
            mut_keys[mut_idx], '({})'.format(gene_names[mut_idx]) if gene_names is not None else '',
            ', '.join(sample_names[sa_idx] for sa_idx in node), math.exp(log_ml)))

    # calculate the final reliability score of a mutation pattern of the probability that no mutation has that pattern
    for node in node_scores.keys():
        if node_scores[node] == 0.0:
            if len(node) == 0 or len(node) == len(sample_names):
                logger.warn('Underflow error for pattern {}. Set probability to minimal float value!'.format(
                    ', '.join(sample_names[sa_idx] for sa_idx in node)))
            else:
                logger.debug('Underflow error for pattern {}. Set probability to minimal float value!'.format(
                    ', '.join(sample_names[sa_idx] for sa_idx in node)))

            node_scores[node] = sys.float_info.min
        node_scores[node] = -math.log(node_scores[node])

    # Show nodes with high reliability score
    for node, score in islice(sorted(node_scores.items(), key=lambda k: -k[1]), 0, 50):
        logger.debug('Pattern {} has a reliability score of {:.2f}.'.format(
            ', '.join(sample_names[sa_idx] for sa_idx in node), score))

    return nodes, node_scores, mut_pattern_weights
