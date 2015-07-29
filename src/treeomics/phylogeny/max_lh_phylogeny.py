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
from copy import deepcopy
import phylogeny.cplex_solver as cps
from phylogeny.phylogeny_utils import Phylogeny


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
        self.col_ids_mp = None
        # dictionary from mutation patterns to the above column ids
        self.mp_col_ids = None
        # 2d np array [mut_idx, mp_col_idx]
        self.mp_weights = None
        # map from identified putative subclones to their original sample
        self.sc_sample_ids = None

        # most likely but also compatible mutation pattern for each variant
        self.max_lh_nodes = None
        self.max_lh_mutations = None
        # map from mut_idx to the weight of the highest ranked evolutionarily compatible MP
        self.max_lh_weights = None

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

    def infer_max_lh_tree(self, min_mp_lh=None):
        """
        Infer maximum likelihood tree via calculation reliability scores for each
        possible mutation pattern from the likelihood that no variant has this pattern
        The inferred solution represents the reliable and evolutionary compatible mutation patterns
        The mutation pattern of each variant is given by the mp maximizing its likelihood
        :param min_mp_lh: minimum likelihood that at least one variant has an incompatible mp
            such that this mp is considered as a subclone
        :return inferred evolutionary tree
        """

        # convert minimum likelihood into minimum reliability score
        if min_mp_lh is not None:
            min_score = -math.log(1.0 - min_mp_lh)
            self.patient.sc_names = deepcopy(self.patient.sample_names)
        else:
            min_score = None

        # necessary to map from identified putative subclones to their original sample
        self.sc_sample_ids = dict()

        # compute various mutation patterns (nodes) and their reliability scores
        self.node_scores, self.col_ids_mp, self.mp_col_ids, self.mp_weights = infer_ml_graph_nodes(
            self.patient.log_p01, self.patient.sample_names, self.patient.mut_keys, gene_names=self.patient.gene_names)

        while True:
            # create conflict graph which forms the input to the ILP
            self.cf_graph = create_conflict_graph_light(self.node_scores)

            # translate the conflict graph into a minimum vertex cover problem
            # and solve this using integer linear programming
            self.conflicting_nodes, self.compatible_nodes = cps.solve_conflicting_phylogeny(self.cf_graph)

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

                for mp_col_idx in np.argsort(self.mp_weights[mut_idx])[::-1]:       # sort by descending log probability

                    if self.col_ids_mp[mp_col_idx] in self.compatible_nodes:
                        # found most likely pattern for this variant
                        self.max_lh_nodes[self.col_ids_mp[mp_col_idx]].add(mut_idx)
                        self.max_lh_mutations[mut_idx] = self.col_ids_mp[mp_col_idx]
                        self.max_lh_weights[mut_idx] = self.mp_weights[mut_idx][mp_col_idx]
                        logger.info('Max LH pattern of variant in {} is {} with log likelihood {:.1e}.'.format(
                            self.patient.gene_names[mut_idx] if self.patient.gene_names is not None
                            else self.patient.mut_keys[mut_idx], self.col_ids_mp[mp_col_idx],
                            self.mp_weights[mut_idx][mp_col_idx]))

                        # determine false positives and false negatives compared to original classification
                        fps = set(self.patient.mutations[mut_idx].difference(self.col_ids_mp[mp_col_idx]))

                        # check if some of these false-positives are present in the newly created subclones
                        for sc_idx, sa_idx in self.sc_sample_ids.items():
                            if sc_idx in self.col_ids_mp[mp_col_idx] and sa_idx in fps:
                                fps.remove(sa_idx)
                        if len(fps) > 0:
                            self.false_positives[mut_idx] = fps

                        # distinguish between real false negatives and variants classified as unknown
                        for sa_idx in self.col_ids_mp[mp_col_idx].difference(self.patient.mutations[mut_idx]):

                            # map from identified putative subclones to their original sample
                            if sa_idx in self.sc_sample_ids:
                                sa_idx = self.sc_sample_ids[sa_idx]
                                if sa_idx in self.patient.mutations[mut_idx]:
                                    # mutation was already classified as present in original sample
                                    # => no false-negative
                                    continue
                                # mutation was not classified as present in original sample
                                # => must be a false-negative

                            if self.patient.data[mut_idx][sa_idx] < 0.0:      # unknown classified mutation
                                self.false_negative_unknowns[mut_idx].add(sa_idx)
                            else:
                                self.false_negatives[mut_idx].add(sa_idx)

                        break

            logger.info('Putative false-positives {}, putative false-negatives {}, put. false neg. unknowns {}'.format(
                sum(len(fps) for mut_idx, fps in self.false_positives.items()),
                sum(len(fns) for mut_idx, fns in self.false_negatives.items()),
                sum(len(fns) for mut_idx, fns in self.false_negative_unknowns.items())))

            # find parsimony-informative evolutionarily incompatible mps with high likelihood
            if min_mp_lh is not None:

                self.sc_sample_ids, found_subclones = self.find_subclones(min_score)
                if not found_subclones:
                    # not enough evidence for additional new subclones
                    break

            else:           # no detection of subclones
                break

        if len(self.sc_sample_ids) > 0:
            logger.info('{} putative subclones have been detected: {}'.format(
                        len(self.patient.sc_names)-len(self.patient.sample_names),
                        ', '.join(self.patient.sc_names[sc_idx] for sc_idx in range(len(self.patient.sample_names),
                                                                                    len(self.patient.sc_names)))))

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

        # construct a phylogenetic tree from the maximum likelihood mutation patterns
        self.mlh_tree = self.infer_evolutionary_tree(self.shared_mlh_mps, self.mlh_founders,
                                                     self.mlh_unique_mutations)

        return self.mlh_tree

    def validate_node_robustness(self, no_replications):
        """
        Validate the robustness of the identified most reliable mutation patterns
        through down-sampling
        :param no_replications: Number of replications for each used fraction of variants
        :return observed frequencies of the mutation patterns in each used fraction of variants
        """

        # Most reliable mutation patterns need to be identified before
        if self.compatible_nodes is None:
            self.infer_max_lh_tree(min_mp_lh=1.0)

        node_frequencies = cps.solve_downsampled_nodes(
            self.cf_graph, self.mp_weights, self.col_ids_mp, no_replications)

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

    def find_subclones(self, min_score):

        updated_nodes, self.sc_sample_ids = self.find_subclonal_mps(min_score)

        if len(updated_nodes) == 0:
            logger.info('There are no more incompatible mutation patterns with a reliability score '
                        'of at least {:.1e}.'.format(min_score))
            return self.sc_sample_ids, False
        else:
            anc = set()     # find the mutations present in the ancestor MP and assign them also to the new MP
            for old_mp, new_mp in updated_nodes.items():

                if old_mp.issuperset(anc):
                    anc = old_mp
                elif (len(old_mp.intersection(anc)) > 0 and
                        len(old_mp.difference(anc)) > 0 and len(anc.difference(old_mp)) > 0):
                    logger.warn('Putative subclones are evolutionary incompatible '
                                'to each other: {} <> {}'.format(anc, old_mp))

            logger.debug('{} is the parental mutation pattern of the putative subclones {}.'.format(
                anc, updated_nodes.keys()))

            # the mutations present in the direct ancestor of the identified
            # parent subclone should be present in all
            parental_mps = set()
            for sc in self.compatible_nodes:
                if sc.issuperset(anc):
                    parental_mps.add(sc)
                    logger.debug('Mutations present in {} should also be present in all putative subclones.'
                                 .format(sc))

            # account for putative subclones in parental mutation patterns
            for a in parental_mps:
                # mutations of parental clone should also be present in new putative subclones
                par_mp = a.union(updated_nodes[anc])

                # update subclones in patient
                self.node_scores[par_mp] = self.node_scores[a]
                del self.node_scores[a]
                self.col_ids_mp[self.mp_col_ids[a]] = par_mp
                self.mp_col_ids[par_mp] = self.mp_col_ids[a]
                del self.mp_col_ids[a]

                logger.info('Updated parental mutation pattern from {} to {}'.format(
                    ', '.join(self.patient.sc_names[sc_idx] for sc_idx in a),
                    ', '.join(self.patient.sc_names[sc_idx] for sc_idx in par_mp)))

            return self.sc_sample_ids, True

    def find_subclonal_mps(self, min_score, max_sc_maf=0.6):
        """
        Identify evolutionarily incompatible mutation patterns with a reliability score greater than min_score which
        where variants in some samples are subclonal
        :param min_score: minimum reliability score of a mutation pattern with putative subclones
        :param max_sc_maf: sum of median MAFs of conflicting MPs must be below this threshold
        :return updated_nodes:
        """

        updated_nodes = dict()

        # find incompatible mutation pattern with the highest reliability score
        for mp in sorted(self.conflicting_nodes, key=lambda k: -self.node_scores[k]):

            if self.node_scores[mp] < min_score:
                # no incompatible mutation patterns with high reliability score exist
                return updated_nodes, self.sc_sample_ids

            # find highest ranked conflicting mutation pattern
            hcmp = max(set(self.cf_graph.neighbors(mp)).difference(self.conflicting_nodes),
                       key=lambda k: self.node_scores[k])
            logger.info('Highest ranked evolutionary incompatible mutation pattern of MP {} (w: {:.1f}): {} (w: {:.1f})'
                        .format(', '.join(self.patient.sc_names[sc] for sc in mp), self.cf_graph.node[mp]['weight'],
                                ', '.join(self.patient.sc_names[sc] for sc in hcmp),
                                self.cf_graph.node[hcmp]['weight']))

            # infer samples with putative mixed subclones
            msc_samples = mp.intersection(hcmp)

            # check if all subsets (descendants) of these putative subclone are subclonal and hence have low MAFs
            mafs_desc = defaultdict(list)

            descendant_mps = dict()
            for desc_mp in _subsets(msc_samples):
                # is this pattern a subset (descendant)?
                if desc_mp.issubset(msc_samples):  # TODO: delete
                    descendant_mps[desc_mp] = None
                else:
                    raise RuntimeError('Should never happen!!!')

            # for desc_mp in descendant_mps:
            #     # infer variants for which the new pattern would be the highest ranked one
            #     # note that smaller mean more likely (log space)
            #     d_muts = [mut_idx for mut_idx in range(len(self.mp_weights))
            #               if (self.mp_weights[mut_idx][self.mp_col_ids[desc_mp]] < self.max_lh_weights[mut_idx])
            #               and all(self.mp_weights[mut_idx][self.mp_col_ids[desc_mp]]
            #                       <= self.mp_weights[mut_idx][self.mp_col_ids[d_mp]] for d_mp in descendant_mps)]
            #     descendant_mps[desc_mp] = d_muts
            #     for s in desc_mp:
            #         for mut_idx in d_muts:
            #             mafs_desc[s].append(
            #                 self.patient.data[mut_idx][s] if self.patient.data[mut_idx][s] > 0.0 else 0)

            # check if the mutations of these patterns are subclonal in the common samples
            mafs_mp = defaultdict(list)
            mafs_hcmp = defaultdict(list)

            # create new subclones if there is enough evidence
            created_scs = dict()
            # create new mutation pattern resulting from the new subclone
            new_mp = set(mp)
            for s in msc_samples:
                mp_muts = [mut_idx for mut_idx in range(len(self.mp_weights))
                           if (self.mp_weights[mut_idx][self.mp_col_ids[mp]] < self.max_lh_weights[mut_idx])]
                for mut_idx in mp_muts:
                    mafs_mp[s].append(self.patient.data[mut_idx][s] if self.patient.data[mut_idx][s] > 0 else 0)

                hcmp_muts = [mut_idx for mut_idx in range(len(self.mp_weights))
                             if (self.mp_weights[mut_idx][self.mp_col_ids[hcmp]] < self.max_lh_weights[mut_idx])]
                for mut_idx in hcmp_muts:
                    mafs_hcmp[s].append(self.patient.data[mut_idx][s] if self.patient.data[mut_idx][s] > 0 else 0)

                # check if overlapping mutations are indeed subclonal supported by low MAFs
                if np.median(mafs_mp[s]) + np.median(mafs_hcmp[s]) > max_sc_maf:
                    logger.info('Found not enough evidence for subclones in sample {}. '.format(
                        self.patient.sample_names[s]))
                    logger.info('Median MAFs of subclones too high: {:.1%} and {:.1%}'.format(
                        np.mean(mafs_mp[s]), np.mean(mafs_hcmp[s])))
                    continue
                elif s in mafs_desc.keys():
                    if np.median(mafs_mp[s]) + np.median(mafs_desc[s]) > max_sc_maf:
                        if len(mafs_desc[s]) < len(self.patient.data)/50.0:
                            logger.info('In sample {} a few mutations have suspiciously high MAF (median: {:.1%}) '
                                        'for being subclonal: {}'.format(self.patient.sample_names[s],
                                                                         np.median(mafs_desc[s])))
                        else:
                            logger.info('Found not enough evidence for subclones in sample {}. '.format(
                                self.patient.sample_names[s]))
                            logger.info('Median MAFs of descending subclones too high: {:.1%} and {:.1%}'.format(
                                np.mean(mafs_mp[s]), np.mean(mafs_desc[s])))
                            continue
                    else:
                        if len(mafs_desc[s]) > 0:
                            logger.info('Low MAFs of mutations in descending sample {} of putative subclone '.format(
                                self.patient.sample_names[s])
                                + 'provide additional evidence for its existance: median {:.1%}'.format(
                                np.median(mafs_desc[s])))

                # no contradicting evidence for subclone found
                logger.info('Found some evidence for subclones in sample {}. '.format(self.patient.sc_names[s]))
                logger.info('Median MAFs of subclones in the samples: {:.1%} and {:.1%}'.format(
                    np.mean(mafs_mp[s]), np.mean(mafs_hcmp[s])))

                # create new subclone in this sample and update mutation pattern
                created_scs[s] = len(self.patient.sc_names)
                self.sc_sample_ids[len(self.patient.sc_names)] = s
                new_mp.remove(s)
                new_mp.add(len(self.patient.sc_names))
                # did a subclone of this sample already get generated
                if '_SC' in self.patient.sc_names[s]:
                    name = self.patient.sc_names[s][:self.patient.sc_names[s].find('_SC')]
                else:
                    name = self.patient.sc_names[s]
                next_id = int(max(n.split('_SC')[1] if len(n.split('_SC')) > 1 else 0
                                  for n in self.patient.sc_names if n.split('_SC')[0] == name))+1
                self.patient.sc_names.append(name+'_SC'+str(next_id))

            # update mutation patterns in the patient
            if len(new_mp.difference(mp)) > 0:
                new_mp = frozenset(new_mp)

                # update subclones in patient
                self.node_scores[new_mp] = self.node_scores[mp]
                del self.node_scores[mp]
                self.col_ids_mp[self.mp_col_ids[mp]] = new_mp
                self.mp_col_ids[new_mp] = self.mp_col_ids[mp]
                del self.mp_col_ids[mp]
                updated_nodes[mp] = new_mp
                logger.info('Created new subclones {}'.format(
                            ', '.join(self.patient.sc_names[sc] for sc in new_mp.difference(mp))))
                break

        # check compatibility with other conflicting clones
        for related_mp in sorted(self.conflicting_nodes, key=lambda k: -self.cf_graph.node[k]['weight']):
            if related_mp == mp:
                continue

            # is this conflicting node evolutionary compatible with the previously identified one
            if related_mp in self.cf_graph.neighbors(mp):
                continue

            # some of the new subclones have to be present in this conflicting node
            if not any(new_sc in related_mp for new_sc in created_scs):
                continue

            # check if the MP (mutation pattern) would be evolutionary compatible if the subclones would be present
            # find highest ranked conflicting mutation pattern
            hcmp = max(set(self.cf_graph.neighbors(related_mp)).difference(self.conflicting_nodes),
                       key=lambda k: self.cf_graph.node[k]['weight'])
            # print('Highest ranked conflicting mutation pattern for {} (w: {:.1f}): {} (w: {:.1f})'.format(
            #      rel_cn, cf_graph.node[rel_cn]['weight'], hcmp, cf_graph.node[hcmp]['weight']))

            if (related_mp.intersection(hcmp)).issubset(created_scs):
                logger.info('New subclone might also be present in this mutation pattern: {} (w: {:.1e})'.format(
                            ', '.join(self.patient.sc_names[sc] for sc in related_mp), self.node_scores[related_mp]))
                updated_nodes = self._update_mutation_pattern(related_mp, created_scs, updated_nodes)
            else:
                if self.node_scores[related_mp] > min_score/10.0:
                    logger.debug('Generating the putative subclones in MP {} (w: {:.1e}) '.format(
                        related_mp, self.node_scores[related_mp])
                        + 'would not make it compatible to {} (w: {:.1e})'.format(hcmp, self.node_scores[hcmp]))

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
        self.col_ids_mp[self.mp_col_ids[old_mp]] = new_mp
        self.mp_col_ids[new_mp] = self.mp_col_ids[old_mp]
        del self.mp_col_ids[old_mp]
        updated_nodes[old_mp] = new_mp
        logger.info('Replaced samples {} and updated mutation pattern (w: {:.1e}) to {}'.format(
                    ', '.join(self.patient.sc_names[new_sc] for new_sc in created_scs.keys() if new_sc in old_mp),
                    self.node_scores[new_mp], ', '.join(self.patient.sc_names[sc] for sc in new_mp)))

        return updated_nodes


def infer_ml_graph_nodes(log_p01, sample_names, mut_keys, gene_names=None):
    """
    Infer maximum likelihood using bayesian inference for each possible mutation pattern
    :param log_p01: posterior: log probability that VAF = 0, log probability that VAF > 0
    :param sample_names:
    :param mut_keys: list with information about the variant
    :param gene_names: list with the names of the genes in which the variant occurred
    :return dictionary of nodes and corresponding variants, pattern reliability scores, weights of patterns of variants
    """

    n = len(sample_names)  # number of samples

    node_scores = dict()    # mutation patterns score summed over all variants

    # weight per inferred mutation pattern per variant given the p0's in each sample for a variant
    mp_weights = np.empty([len(log_p01), math.pow(2, n)], dtype=np.float64)
    # mutation pattern to the corresponding column id in the weight matrix
    mp_col_ids = dict()
    # column id to corresponding mutation pattern
    col_ids_mp = dict()

    trunk_mp = frozenset([sa_idx for sa_idx in range(n)])
    col_idx = 0
    for no_pres_vars in range(n+1):     # number of present variants
        for mp in combinations(range(n), no_pres_vars):  # generate all mps with length no_pres_vars

            node = frozenset(mp)        # create mutation pattern
            col_ids_mp[col_idx] = node
            mp_col_ids[node] = col_idx

            for mut_idx in range(len(log_p01)):
                log_ml = 0.0                # log maximum likelihood of the inferred pattern
                for sa_idx in node:                         # variant is present
                    log_ml += log_p01[mut_idx][sa_idx][1]
                for sa_idx in trunk_mp.difference(node):    # variant is absent
                    log_ml += log_p01[mut_idx][sa_idx][0]

                if log_ml == 0.0:   # numerical artifact, use approximation: ignore second order term
                    for sa_idx in node:                         # variant is present
                        log_ml += math.exp(log_p01[mut_idx][sa_idx][0])     # sum probability
                    for sa_idx in trunk_mp.difference(node):    # variant is absent
                        log_ml += math.exp(log_p01[mut_idx][sa_idx][1])

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
                        log_ml = -1.0e-150
                    assert log_ml < 0.0, ('Underflow error while calculating the probability that the '
                                          + 'variant {} does not have pattern {}.'.format(
                                            mut_keys[mut_idx], ', '.join(sample_names[sa_idx] for sa_idx in node)))

                # assign calculated log probability that this variant has this mutation pattern
                mp_weights[mut_idx][col_idx] = log_ml

                # calculate the probability of a mp that no variant has this mutation pattern
                # product of (1 - the probability that a variant has this mp)
                if node in node_scores.keys():
                    node_scores[node] -= math.log(-math.expm1(log_ml))
                else:
                    node_scores[node] = -math.log(-math.expm1(log_ml))

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
            raise RuntimeError('Underflow error for pattern {}. Set probability to minimal float value!'.format(
                ', '.join(sample_names[sa_idx] for sa_idx in node)))
            # should never happen, delete this check

    # Show nodes with high reliability score
    for node, score in itertools.islice(sorted(node_scores.items(), key=lambda k: -k[1]), 0, 50):
        logger.info('Pattern {} has a reliability score of {:.2e}.'.format(node, score))

    return node_scores, col_ids_mp, mp_col_ids, mp_weights


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


def _subsets(mp):
    """
    Returns all subsets (descendant) mutation patterns of a given mutation pattern
    :param mp: ancestor mp
    :return: descendant mp
    """

    for l in range(len(mp)-1, 0, -1):
        for descendant_mp in itertools.combinations(mp, l):
            yield frozenset(descendant_mp)
