#!/usr/bin/python
"""Infers most reliable mutation pattern based on a bayesian inference model"""
__author__ = 'Johannes REITER'
__date__ = 'July, 2015'

import logging
import itertools
import math
from collections import defaultdict
from itertools import combinations
import numpy as np
import networkx as nx
import sys
import phylogeny.cplex_solver as cps
from phylogeny.phylogeny import Phylogeny


# get logger for application
logger = logging.getLogger('treeomics')


class MaxLHPhylogeny(Phylogeny):
    """
    Find likely (low reliability scores) technical and biological artifacts in the data
    Infer evolutionary tree from the resolved data
    """

    def __init__(self, patient, mps):

        Phylogeny.__init__(self, patient, mps)

        # dictionary column ids in the mp_weights 2d array and the corresponding mutation pattern
        self.mp_col_ids = None
        # 2d np array [mut_idx, mp_col_idx]
        self.mp_weights = None

        # most likely but also compatible mutation pattern for each variant
        self.max_lh_nodes = None
        self.max_lh_mutations = None

        # false positives and false negatives compared to original classification
        self.false_positives = None
        self.false_negatives = None
        self.false_negative_unknowns = None

        # updated clones after artifacts in the data have been removed
        self.shared_mlh_mps = None
        self.mlh_founders = None
        self.mlh_unique_mutations = None
        self.mlh_absent_mutations = None

        # variables used for resolving technical and biological artifacts in the data
        self.no_incompatible_positions = None
        self.incompatible_positions = None

        self.mlh_tree = None

    def infer_max_lh_tree(self):
        """
        Infer maximum likelihood tree via calculation reliability scores for each
        possible mutation pattern from the likelihood that no variant has this pattern
        The inferred solution represents the reliable and evolutionary compatible mutation patterns
        The mutation pattern of each variant is given by the mp maximizing its likelihood
        :return inferred evolutionary tree
        """

        # compute various mutation patterns (nodes) and their reliability scores
        self.node_scores, self.mp_col_ids, self.mp_weights = infer_ml_graph_nodes(
            self.patient.log_p0, self.patient.sample_names, self.patient.gene_names)

        # create conflict graph which forms the input to the ILP
        self.cf_graph = create_conflict_graph_light(self.node_scores)

        # translate the conflict graph into a minimum vertex cover problem
        # and solve this using integer linear programming
        self.conflicting_nodes, self.compatible_nodes = cps.solve_conflicting_phylogeny(self.cf_graph)

        # find the most likely but also compatible mutation pattern for each variant
        self.max_lh_nodes = defaultdict(set)
        self.max_lh_mutations = dict()
        self.false_positives = defaultdict(set)
        self.false_negatives = defaultdict(set)
        self.false_negative_unknowns = defaultdict(set)

        for mut_idx in range(len(self.mp_weights)):

            for mp_col_idx in np.argsort(self.mp_weights[mut_idx])[::-1]:       # sort by descending log probability

                if self.mp_col_ids[mp_col_idx] in self.compatible_nodes:
                    # found most likely pattern for this variant
                    self.max_lh_nodes[self.mp_col_ids[mp_col_idx]].add(mut_idx)
                    self.max_lh_mutations[mut_idx] = self.mp_col_ids[mp_col_idx]
                    logger.info('Max LH pattern of variant in {} is {} with log likelihood {:.1e}.'.format(
                        self.patient.gene_names[mut_idx], self.mp_col_ids[mp_col_idx],
                        self.mp_weights[mut_idx][mp_col_idx]))

                    # determine false positives and false negatives compared to original classification
                    fps = self.patient.mutations[mut_idx].difference(self.mp_col_ids[mp_col_idx])
                    if len(fps) > 0:
                        self.false_positives[mut_idx] = fps

                    # distinguish between real false negatives and variants classified as unknown
                    for sa_idx in self.mp_col_ids[mp_col_idx].difference(self.patient.mutations[mut_idx]):
                        if self.patient.data[mut_idx][sa_idx] < 0.0:      # unknown classified mutation
                            self.false_negative_unknowns[mut_idx].add(sa_idx)
                        else:
                            self.false_negatives[mut_idx].add(sa_idx)

                    break

        logger.info('Putative false-positives {}, putative false-negatives {}, put. false neg. unknowns {}'.format(
            sum(len(fps) for mut_idx, fps in self.false_positives.items()),
            sum(len(fns) for mut_idx, fns in self.false_negatives.items()),
            sum(len(fns) for mut_idx, fns in self.false_negative_unknowns.items())))

        # find founders and unique mutations in updated mutation list
        # build up a dictionary of resolved clones where
        # the likely sequencing errors have been updated
        self.mlh_founders = set()                       # set of founding mutations present in all samples
        self.mlh_unique_mutations = defaultdict(set)    # private mutations only appear in the leaves
        self.mlh_absent_mutations = set()               # mutations inferred to be absent in all samples
        self.shared_mlh_mps = defaultdict(set)          # parsimony-informative mps

        for mut_idx, samples in self.max_lh_mutations.items():
            if len(samples) == len(self.patient.sample_names):          # founder mut.
                self.mlh_founders.add(mut_idx)
            elif 1 < len(samples) < len(self.patient.sample_names):     # shared mut.
                self.shared_mlh_mps[frozenset(samples)].add(mut_idx)
            elif len(samples) == 1:                                     # unique mut.
                for sa_idx in samples:
                    self.mlh_unique_mutations[sa_idx].add(mut_idx)
            else:
                self.mlh_absent_mutations.add(mut_idx)

        # construct a phylogenetic tree from the maximum likelihood mutation patterns
        self.mlh_tree = self.infer_evolutionary_tree(self.shared_mlh_mps, self.mlh_founders,
                                                     self.mlh_unique_mutations)

        return self.mlh_tree


def infer_ml_graph_nodes(log_p0, sample_names, gene_names):
    """
    Infer maximum likelihood using bayesian inference for each possible mutation pattern
    :param log_p0: posterior: log probability that VAF = 0
    :param sample_names:
    :param gene_names:
    :return dictionary of nodes and corresponding variants, pattern reliability scores, weights of patterns of variants
    """

    n = len(sample_names)  # number of samples

    node_scores = dict()    # mutation patterns score summed over all variants

    # weight per inferred mutation pattern per variant given the p0's in each sample for a variant
    mp_weights = np.empty([len(log_p0), math.pow(2, n)], dtype=np.float64)
    # mutation pattern to the corresponding column id in the weight matrix
    mp_col_ids = dict()

    trunk_mp = frozenset([sa_idx for sa_idx in range(n)])
    col_idx = 0
    for no_pres_vars in range(n+1):     # number of present variants
        for mp in combinations(range(n), no_pres_vars):  # generate all mps with length no_pres_vars

            node = frozenset(mp)        # create mutation pattern
            mp_col_ids[col_idx] = node

            for mut_idx in range(len(log_p0)):
                log_ml = 0.0                # log maximum likelihood of the inferred pattern
                for sa_idx in node:                         # variant is present
                    log_ml += math.log(-np.expm1(log_p0[mut_idx][sa_idx]))
                for sa_idx in trunk_mp.difference(node):    # variant is absent
                    log_ml += log_p0[mut_idx][sa_idx]

                if log_ml == 0.0:   # numerical artifact, use approximation: ignore second order term
                    for sa_idx in node:                         # variant is present
                        log_ml += math.exp(log_p0[mut_idx][sa_idx])     # sum probability
                    for sa_idx in trunk_mp.difference(node):    # variant is absent
                        log_ml += -np.expm1(log_p0[mut_idx][sa_idx])

                    log_ml = np.log1p(-log_ml)      # calculates log(1+argument)
                    logger.debug('Approximated log probability of variant {} having pattern {} by {:.2e}.'.format(
                        gene_names[mut_idx], node, log_ml))
                    if log_ml == 0.0:
                        if len(node) == 0 or len(node) == n:
                            logger.debug('Underflow error. Set probability to minimal float value!')
                        else:
                            logger.warn('Underflow error. Set probability to minimal float value!')
                        log_ml = -sys.float_info.min
                    assert log_ml < 0.0, ('Underflow error while calculating the probability that the '
                                          + 'variant {} does not have pattern {}.'.format(
                                            gene_names[mut_idx], ', '.join(sample_names[sa_idx] for sa_idx in node)))

                # assign calculated log probability that this variant has this mutation pattern
                mp_weights[mut_idx][col_idx] = log_ml

                # calculate the probability of a mp that no variant has this mutation pattern
                # product of (1 - the probability that a variant has this mp)
                if node in node_scores.keys():
                    node_scores[node] *= -np.expm1(log_ml)
                else:
                    node_scores[node] = -np.expm1(log_ml)

                # logger.debug('Variant {} has pattern {} with probability {:.1e}.'.format(
                #     gene_names[mut_idx], ', '.join(sample_names[sa_idx] for sa_idx in node), math.exp(log_ml)))

            col_idx += 1      # next mutation pattern

    # calculate the final reliability score of a mutation pattern from the probability that no mutation has that pattern
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

    for node, score in sorted(node_scores.items(), key=lambda k: -k[1]):
        logger.info('Pattern {} has a reliability score of {:.2e}.'.format(node, score))

    return node_scores, mp_col_ids, mp_weights


def create_conflict_graph_light(reliability_scores):
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
