"""Data structure for compatible solutions of the MILP"""
import logging
import math
import numpy as np
from collections import defaultdict
import settings
import utils.int_settings as def_sets


__author__ = 'Johannes REITER'
__date__ = 'March 6, 2017'


# get logger for application
logger = logging.getLogger('treeomics')


class Solution:

    # scale all values to avoid numerical issues with CPLEX
    # calculated in CPLEX solver
    SCALING_FACTOR = 1

    def __init__(self, rank, obj_func_vals, obj_val, solution_values, cf_graph, no_muts):

        self.rank = rank            # solution ranking
        self.obj_val = obj_val      # value of objective function

        # calculate log likelihood
        self._calculate_sol_log_likelihood(solution_values, obj_func_vals, no_muts)
        # logger.debug('Solution of rank {} (obj-val: {:.2e}) with log-likelihood {:.3e} ({:.3e}).'.format(
        #     rank, obj_val / Solution.SCALING_FACTOR, self.llh, math.exp(self.llh)))

        self._translate_solution(solution_values, cf_graph, rank=rank)

        # maps from each evolutionarily compatible MP to the selected set of variants
        self.max_lh_nodes = None
        # map from mut_idx to the highest ranked evolutionarily compatible MP
        self.max_lh_mutations = None
        # map from mut_idx to the weight of the highest ranked evolutionarily compatible MP
        self.max_lh_weights = None

        self.mlh_founders = None            # set of founding mutations present in all samples
        self.mlh_unique_mutations = None    # private mutations only appear in the leaves
        self.mlh_absent_mutations = None    # mutations inferred to be absent in all samples
        self.shared_mlh_mps = None          # parsimony-informative mps

        # mutations present in each sample inferred by maximum likelihood tree
        self.variants = None

        # only relevant if not the whole solution space is explored
        self.conflicting_mutations = None

        # false positives and false negatives compared to original classification
        # likely technical or biological artifacts in the data
        self.false_positives = None
        self.false_negatives = None
        self.false_negative_unknowns = None

    def _translate_solution(self, sol_values, cf_graph, rank=None):
        # set of conflicting mutation patterns, set of compatible mutation patterns
        # translate solution of the ILP back to the phylogeny problem
        # removing the minimum vertex cover (conflicting mutation patterns) gives the maximum compatible set of mps
        self.compatible_nodes = set()
        self.incompatible_nodes = set()
        conflicting_mutations_weight = 0

        for col_idx, (node, data) in enumerate(cf_graph.nodes(data=True)):

            # if round(sol.get_values(str(node)), 5) == 0:
            if round(sol_values[col_idx], 5) == 0:              # compatible
                self.compatible_nodes.add(node)
            else:                                               # incompatible
                self.incompatible_nodes.add(node)
                conflicting_mutations_weight += data['weight']

        if rank is not None and rank <= 30:
            logger.debug('Solution values: '
                         + ', '.join('{}: {:.1f} (w:{:.2e})'.format(node, sol_values[col_idx], data['weight'])
                                     for col_idx, (node, data) in enumerate(cf_graph.nodes(data=True))))

        assert round(conflicting_mutations_weight, 4) == round(self.obj_val / Solution.SCALING_FACTOR, 4), \
            "As long as weights are given by the number of mutations: {} == {}".format(
                conflicting_mutations_weight, self.obj_val / Solution.SCALING_FACTOR)

        if not settings.SHOW_BI_ABSENT_MUTS:
            self.compatible_nodes.add(frozenset([]))

        if rank is not None and rank <= 30:
            logger.debug('Identified {} compatible mutation patterns: {}'.format(
                len(self.compatible_nodes), self.compatible_nodes))

    def _calculate_sol_log_likelihood(self, sol_values, obj_func_vals, no_muts):
        """
        Calcualte log likelihood of solution
        :param sol_values: list of almost binary variables of a mutation pattern is present or absent in the solution
        :param obj_func_vals: list of reliability scores for each mutation pattern
        :param no_muts: number of considered mutations
        :return: log likelihood of solution
        """

        self.llh = 0.0  # log likelihood
        # llh_comp_mps = 0.0

        for ilp_col_idx, val in enumerate(sol_values):
            # o = obj_vals[ilp_col_idx]
            # p_p = -np.expm1(-obj_vals[ilp_col_idx] * no_muts)  # probability to exhibit at least one mutation
            # p_a = math.exp(-obj_vals[ilp_col_idx] * no_muts)   # probability to not exhibit any mutation

            # mutation pattern was not selected and hence is evolutionary incompatible with the nodes in the solution
            if round(val, 5) != 0:
                # unnormalize reliability score by multiplying with no_muts
                # probability that no mutation has this pattern
                self.llh += -obj_func_vals[ilp_col_idx] * no_muts

            # mutation pattern was selected and hence is evolutionary compatible
            # else:
            #     # numpy.expm1(x[, out]): exp(x) - 1 =>> -(exp(x) - 1) =>> 1 - exp(x)
            #     # unnormalize reliability score by multiplying with no_muts
            #     # probability that some mutation has this pattern
            #     # self.llh += math.log(-np.expm1(-obj_func_vals[ilp_col_idx] * no_muts))
            #     llh_comp_mps += math.log(-np.expm1(-obj_func_vals[ilp_col_idx] * no_muts))

        # logger.debug('Compatible mutation patterns likelihood: {:.3e} (log: {:.3e})'.format(
        #     math.exp(llh_comp_mps), llh_comp_mps))

    def get_unscaled_obj_value(self):
        # scaling factor is used that python does not run into numerical precision errors with CPLEX
        return self.obj_val / Solution.SCALING_FACTOR

    def assign_variants(self, phylogeny, max_no_mps):

        if max_no_mps is not None and max_no_mps < math.pow(2, phylogeny.patient.n):
            logger.warning('Some variants might be evolutionarily incompatible since the '
                           'solution space is only partially explored!')
            self.conflicting_mutations = set()

        # ##### assign each variant to the highest ranked evolutionarily compatible mutation pattern ########
        # merge all identified compatible nodes with the parsimony-uninformative nodes (founders, unique)
        sol_nodes = set(list(self.compatible_nodes))

        for node in phylogeny.node_scores.keys():
            if len(node) == 1 or len(node) == len(phylogeny.patient.sample_names):
                sol_nodes.add(node)

        # maps from each evolutionarily compatible MP to the selected set of variants
        self.max_lh_nodes = defaultdict(set)
        # map from mut_idx to the highest ranked evolutionarily compatible MP
        self.max_lh_mutations = dict()
        # map from mut_idx to the weight of the highest ranked evolutionarily compatible MP
        self.max_lh_weights = dict()

        # find the highest ranked evolutionarily compatible mutation pattern for each variant
        for mut_idx in range(len(phylogeny.mp_weights)):

            # sort by descending log likelihood
            for sol_idx, (mp_col_idx, _) in enumerate(
                    sorted(phylogeny.mp_weights[mut_idx].items(), key=lambda k: -k[1]), 1):

                if phylogeny.idx_to_mp[mp_col_idx] in sol_nodes:
                    # found most likely pattern for this variant
                    self.max_lh_nodes[phylogeny.idx_to_mp[mp_col_idx]].add(mut_idx)
                    self.max_lh_mutations[mut_idx] = phylogeny.idx_to_mp[mp_col_idx]
                    self.max_lh_weights[mut_idx] = phylogeny.mp_weights[mut_idx][mp_col_idx]
                    # logger.debug('Chose {}th highest ranked pattern.'.format(sol_idx)
                    #     + 'Max LH pattern of variant in {} is {} with log likelihood {:.1e}.'.format(
                    #       phylogeny.patient.gene_names[mut_idx] if phylogeny.patient.gene_names is not None
                    #       else phylogeny.patient.mut_keys[mut_idx], phylogeny.idx_to_mp[mp_col_idx],
                    #       phylogeny.mp_weights[mut_idx][mp_col_idx]))

                    # search can be stopped since most likely pattern has been found for this variant
                    break

            # no evolutionarily compatible mutation pattern was among the <max_no_mps> most likely pattern
            # of this variant; variant will be incompatible to the inferred tree!
            else:

                # check if indeed the solution space is only partially explored
                assert max_no_mps is not None and max_no_mps < math.pow(2, phylogeny.patient.n), \
                    'Compatible mutation pattern must exist when the full solution space is explored!'

                self.conflicting_mutations.add(mut_idx)

    def infer_mutation_patterns(self, phylogeny):
        """
        Find founders, shared and unique mutations in updated mutation list
        Build up a dictionary of resolved subclones where the likely sequencing errors have been updated
        """

        self.mlh_founders = set()  # set of founding mutations present in all samples
        self.mlh_unique_mutations = defaultdict(set)  # private mutations only appear in the leaves
        self.mlh_absent_mutations = set()  # mutations inferred to be absent in all samples
        self.shared_mlh_mps = defaultdict(set)  # parsimony-informative mps

        for mut_idx, samples in self.max_lh_mutations.items():
            if len(samples) == len(phylogeny.patient.sample_names) + len(phylogeny.sc_sample_ids):  # founder mut.
                self.mlh_founders.add(mut_idx)
            elif 1 < len(samples) < len(phylogeny.patient.sample_names) + len(phylogeny.sc_sample_ids):  # shared mut.
                self.shared_mlh_mps[frozenset(samples)].add(mut_idx)
            elif len(samples) == 1:  # unique mut.
                for sa_idx in samples:
                    self.mlh_unique_mutations[sa_idx].add(mut_idx)
            else:
                self.mlh_absent_mutations.add(mut_idx)

        # mutations present in each sample inferred by maximum likelihood tree
        self.variants = defaultdict(list)

        # compute the number of persistent and present mutations inferred
        # in the evolutionary trajectory of each sample
        for sa_idx in range(len(phylogeny.patient.sample_names)):
            for mut_idx in self.mlh_founders:
                self.variants[sa_idx].append(mut_idx)
            for mut_idx in self.mlh_unique_mutations[sa_idx]:
                self.variants[sa_idx].append(mut_idx)

        for mps, muts in self.shared_mlh_mps.items():
            for sa_idx in mps:
                for mut_idx in muts:
                    self.variants[sa_idx].append(mut_idx)

    def find_artifacts(self, phylogeny):

        # dictionaries from mut_idx to putative artifacts
        self.false_positives = defaultdict(set)
        self.false_negatives = defaultdict(set)
        self.false_negative_unknowns = defaultdict(set)

        # log probability to be classified with at least the given confidence threshold
        conf_clas_lpth = math.log(def_sets.CLA_CONFID_TH)
        lpth_half = math.log(0.5)

        # find the highest ranked evolutionarily compatible mutation pattern for each variant
        for mut_idx in range(len(phylogeny.mp_weights)):

            # artifact calculation based on BI model from the p-value based model in new Treeomics version
            fps = set()
            # distinguish between real false negatives and variants classified as unknown
            for sa_idx, (p0, p1) in enumerate(phylogeny.patient.log_p01[mut_idx]):
                if p1 > conf_clas_lpth and sa_idx not in self.max_lh_mutations[mut_idx]:
                    fps.add(sa_idx)

                if sa_idx in phylogeny.sc_sample_ids.keys():
                    while sc_idx in phylogeny.sc_sample_ids.keys():
                        sc_idx = phylogeny.sc_sample_ids[sc_idx]
                    # if sc_idx in self.patient.mutations[mut_idx]:
                    if phylogeny.patient.log_p01[mut_idx][sc_idx][1] > conf_clas_lpth:
                        # mutation was already classified as present in original sample
                        # => no false-negative
                        continue
                else:
                    sc_idx = sa_idx

                if p0 > conf_clas_lpth and sc_idx in self.max_lh_mutations[mut_idx]:
                    self.false_negatives[mut_idx].add(sc_idx)
                elif p0 <= lpth_half and p1 <= conf_clas_lpth and sc_idx in self.max_lh_mutations[mut_idx]:
                    self.false_negative_unknowns[mut_idx].add(sc_idx)

            # check if some of these false-positives are present in the newly created subclones
            for sc_idx, sa_idx in phylogeny.sc_sample_ids.items():
                if sc_idx in self.max_lh_mutations[mut_idx] and sa_idx in fps:
                    fps.remove(sa_idx)

            if len(fps) > 0:
                self.false_positives[mut_idx] = fps
