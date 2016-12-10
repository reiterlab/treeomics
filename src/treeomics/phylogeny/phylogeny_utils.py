#!/usr/bin/python
"""Provides abstract class for all subclasses (different methods to infer compatible patterns)"""
__author__ = 'Johannes REITER'
__date__ = 'April, 2014'

import logging
import itertools
import os
from copy import deepcopy
from collections import defaultdict
import networkx as nx
from networkx.readwrite import json_graph
import json


# get logger for application
logger = logging.getLogger('treeomics')

TREE_ROOT = 'germline'      # root node of the evolutionary tree


class Phylogeny(object):
    """
    Managing the derivation of perfect (and persistent) trees
    representing the clonal evolution of cancer
    """
    def __init__(self, patient, mps):

        # reference to the patient data
        self.patient = patient

        # dictionary of mutation patterns with the corresponding set of shared mutations
        self.mps = mps

        self._tree_root = None

        self.node_scores = None
        self.compatible_nodes = None
        self.conflicting_nodes = None

        # only relevant if not the whole solution space is explored
        self.conflicting_mutations = None

        # mutation pattern conflict graph
        self.cf_graph = None

        # path to generated PNG of inferred tree for HTML report
        self.tree_plot = None

    def infer_evolutionary_tree(self, updated_mps, founders, unique_mutations, confidence=None):
        """
        Construct a phylogenetic tree from the maximum set of compatible mutations
        :param updated_mps: evolutionary-conflict-free mutation patterns from which the tree is inferred
        :param founders: set of mutations which are present in all samples
        :param unique_mutations: dictionary of mutations which are only present in a sample
        :param confidence: confidence value in parsimony-informative branch
        :return inferred tree
        """

        # create new persistent tree and the root node
        tree = nx.DiGraph()
        tree.add_node(TREE_ROOT, name=TREE_ROOT, muts=set())

        founding_mp = frozenset([sa_idx for sa_idx in range(
            len(self.patient.sample_names) if self.patient.sc_names is None else len(self.patient.sc_names))])
        tree.add_node(founding_mp, name='', muts=founders)

        # add edge from the germline to the founding clone
        tree.add_edge(TREE_ROOT, founding_mp, muts=tree.node[founding_mp]['muts'])

        # add compatible mutation patterns to the evolutionary tree
        for mp, muts in updated_mps.items():
            if mp == founding_mp:
                for mut in muts:
                    tree.edge[TREE_ROOT][mp]['muts'].add(mut)
            elif len(mp) > 1:
                _add_evolutionary_node(tree, founding_mp, mp, '', muts, confidence=confidence)
            else:
                # leaves (private mutations) are considered separately
                pass

        # Insert the leaves (samples) into the tree with all its private mutations
        for sc_idx in range(
                len(self.patient.sample_names) if self.patient.sc_names is None else len(self.patient.sc_names)):
            # add leaves to the evolutionary tree
            _add_evolutionary_node(tree, founding_mp, frozenset([sc_idx]),
                                   self.patient.sample_names[sc_idx] if self.patient.sc_names is None
                                   else self.patient.sc_names[sc_idx], unique_mutations[sc_idx])

        # _add_subclone_mutations(tree, TREE_ROOT, set())
        # prepare tree for meaningful output in the figures
        self._add_evolutionary_information(tree)

        logger.info('Inferred evolutionary tree from the compatible mutation patterns.')

        return tree

    def _add_evolutionary_information(self, tree):
        """
        Traverse tree and assign meaningful names to the nodes
        """

        # use first in first out queue to traverse the tree in level order: start with the root
        fifo_queue = list()
        fifo_queue.append((TREE_ROOT, set()))
        subclone_idx = 0

        while len(fifo_queue):

            # take first node of the first in first out queue and process it
            cur_node, acquired_mutations = fifo_queue.pop(0)
            tree.node[cur_node]['muts'] = acquired_mutations

            # add meaningful names to internal nodes
            if len(tree.successors(cur_node)):      # internal node

                if subclone_idx == 0:
                    tree.node[cur_node]['name'] = 'Germline'+' '+self.patient.name
                else:
                    # subclones are numbered in level order
                    tree.node[cur_node]['name'] = 'SC '+str(subclone_idx)
                subclone_idx += 1

                # add children of the current node to the queue
                for child in tree.successors(cur_node):
                    total_mutations = acquired_mutations.copy()
                    fifo_queue.append((child, total_mutations.union(tree[cur_node][child]['muts'])))

            # else:       # this node is a leave node (=sample)

        logger.debug("Added evolutionary information to the derived tree.")

    @staticmethod
    def save_json_tree(filepath, tree):
        """
        Transform inferred phylogeny to JSON object and save it to a file
        :param filepath: path to the output file
        :param tree: reconstructed phylogeny given as networkx DiGraph object
        """

        # create simplified JSON tree
        out_ids = dict()
        out_tree = nx.DiGraph()
        for out_id, node in enumerate(tree.nodes_iter()):
            out_ids[node] = out_id
            out_tree.add_node(out_id, name=tree.node[node]['name'])

        for u, v in tree.edges_iter():
            out_tree.add_edge(out_ids[u], out_ids[v], value=len(tree.edge[u][v]['muts']))

        # create json output from reconstructed phylogeny
        json_data = json_graph.node_link_data(out_tree)

        # json object to output file
        with open(filepath, 'w') as json_file:
            json.dump(json_data, json_file, indent=4)
            logger.info('Create JSON file from reconstructed phylogeny: {}.'.format(filepath))

    @staticmethod
    def write_html_file(outfilepath, json_filepath):
        with open(outfilepath, 'w') as html_file, open(_get_html_template()) as temp_file:
            for line in temp_file:
                html_file.write(line.replace('SETFILENAME', json_filepath))


def create_conflict_graph(nodes, weights=None):
    """
    Create a graph where the nodes are given by the mutation patterns and
    the edges model the evolutionary conflicts among them
    :param nodes: dictionary of frozensets describing the samples where a set of mutations is present
    :param weights: each node in the graph is weighted corresponding to the confidence
    in the sequencing data of the mutation modeled by the reliability scores
    :return conflict graph
    """

    # create empty conflict graph
    cf_graph = nx.Graph()
    incompatible_mp = defaultdict(set)

    # since each mutation occurs at most once running through this is
    # in O(m^2 * n) = O(|mutations|*|mutations|*|samples|) as
    # the number of clones is bounded by the number of distinct mutations
    for (node1, node2) in itertools.combinations(nodes.keys(), 2):

        # characters need to appear at least in two samples
        # otherwise a conflict is not possible
        if len(node1) < 2 or len(node2) < 2:
            continue

        if len(node1.intersection(node2)) > 0:     # at least one sample where both characters are present (11)

            # check if some characters are present in one clone but not in the other and visa versa
            if len(node1.difference(node2)) > 0 and len(node2.difference(node1)) > 0:   # check for 01 and 10
                # => conflict exists among characters of clone 1 and clone 2

                incompatible_mp[node1].add(node2)
                incompatible_mp[node2].add(node1)

                if node1 not in cf_graph:
                    if weights is not None:
                        cf_graph.add_node(node1, weight=weights[node1], muts=deepcopy(nodes[node1]))
                    else:
                        cf_graph.add_node(node1, weight=len(nodes[node1]), muts=deepcopy(nodes[node1]))

                if node2 not in cf_graph:
                    if weights is not None:
                        cf_graph.add_node(node2, weight=weights[node2], muts=deepcopy(nodes[node2]))
                    else:
                        cf_graph.add_node(node2, weight=len(nodes[node2]), muts=deepcopy(nodes[node2]))

                # add edge between conflicting clones
                cf_graph.add_edge(node1, node2)

    logger.info('Created conflict graph with {} nodes of weight {:.2f} and {} evolutionary conflicts.'.format(
        cf_graph.order(), sum(data['weight'] for _, data in cf_graph.nodes_iter(data=True)), cf_graph.size()))

    return cf_graph


def compute_graph_nodes(mps, sample_names, mut_names, present_p_values, absent_p_values,
                        coverage=None, min_absent_cov=0):
    """
    Compute the nodes (given by the pattern of the clones) in the conflict graph
    If the minimum required coverage at non-significantly mutated mutations is 0,
    the number of clones is identical to the number of nodes (i.e., no unknown positions)
    :param mps: mutation patterns
    :param sample_names
    :param mut_names: mutation keys
    :param present_p_values: p-values for the confidence that a variant is present
    :param absent_p_values: p-values for the confidence that a variant is absent
    :param coverage: dictionary of mutation keys and dictionary of sample names for the coverage
    :param min_absent_cov: minimum required coverage at non-significantly mutated position otherwise unknown
    :return dictionary of nodes and corresponding variants, node weights, dictionary of clones which are unknown
    """

    # clones form the basis for the nodes in the graph
    nodes = deepcopy(mps)
    node_weights = dict()
    mut_pattern_scores = dict()

    unknown_positions = defaultdict(set)
    # dictionary of mutations which form different mutation patterns due to unknown positions
    unknown_muts = defaultdict(dict)
    for node in nodes.keys():

        # calculate weight of each mutation pattern by going through all variants having the same pattern
        node_weights[node] = 0

        for mut_idx in nodes[node]:
            # t = []
            weight = 1.0
            for sa_idx, sample_name in enumerate(sample_names):
                if sa_idx in node:          # variants are present in these samples

                    # check if data for p-value calculation was provided
                    if mut_names[mut_idx] in present_p_values[sample_name]:
                        weight *= 1.0 - present_p_values[sample_name][mut_names[mut_idx]]
                    else:
                        weight *= 1.0
                    # t.append(present_p_values[sample_name][mut_names[mut_idx]])
                else:                       # variants are absent in these samples

                    # not enough coverage at this position to sure about the absence
                    # check if the presence, reduces the number of conflicts
                    if coverage is not None and mut_names[mut_idx] in coverage and \
                            sample_name in coverage[mut_names[mut_idx]] and \
                            0 <= coverage[mut_names[mut_idx]][sample_name] < min_absent_cov:

                        unknown_positions[(node, mut_idx)].add(sa_idx)
                        logger.debug('Variant {} has low coverage of {} and absent-p-value {:.3f}'.format(
                            mut_names[mut_idx], coverage[mut_names[mut_idx]][sample_name],
                            min(absent_p_values[sample_name][mut_names[mut_idx]], 0.5)))

                    # p-value can be at most 0.5 (variant is present or absent with 50% chance each)
                    # come positions have coverage 0; presence or absence has to be random
                    # check if data for p-value calculation was provided
                    if mut_names[mut_idx] in absent_p_values[sample_name]:
                        weight *= 1.0 - min(absent_p_values[sample_name][mut_names[mut_idx]], 0.5)
                    # else:
                    #     weight *= 0.9       # coverage data is not available!
                    # t.append(absent_p_values[sample_name][mut_names[mut_idx]])

            # logger.debug('Variant {} has weight {:.3e} from p-values {}.'.format(
            #     mut_names[mut_idx], weight, ', '.join('{:.2e}'.format(a) for a in t)))
            mut_pattern_scores[mut_idx] = weight
            node_weights[node] += weight
            if (node, mut_idx) in unknown_positions:
                unknown_muts[mut_idx][node] = weight

    # run through all unknown positions and generate nodes with all possible patterns
    for node, mut_idx in unknown_positions.keys():

        # generates all combinations (present-absent) of all unknown positions in a given clone
        for length in range(1, len(unknown_positions[(node, mut_idx)])+1):
            for unknown_samples in itertools.combinations(unknown_positions[(node, mut_idx)], length):

                # form new mutation pattern
                new_node = frozenset(node.union(unknown_samples))
                weight = 1.0
                for sa_idx, sample_name in enumerate(sample_names):
                    if sa_idx in node:          # variants are present in these samples
                        if mut_names[mut_idx] in present_p_values[sample_name]:
                            weight *= 1.0 - present_p_values[sample_name][mut_names[mut_idx]]
                        # t.append(present_p_values[sample_name][mut_names[mut_idx]])
                    elif sa_idx in unknown_samples:
                        # p-value can be at most 0.5 (variant is present or absent with 50% chance each)
                        if mut_names[mut_idx] in absent_p_values[sample_name]:
                            weight *= min(absent_p_values[sample_name][mut_names[mut_idx]], 0.5)
                        else:
                            # will never be used, since unknowns can only appear if there is coverage data
                            # and if there is coverage data, then there are also p-values
                            weight *= 0.5
                    else:                       # variants are absent in these samples
                        # p-value can be at most 0.5 (variant is present or absent with 50% chance each)
                        if mut_names[mut_idx] in absent_p_values[sample_name]:
                            weight *= 1.0 - min(absent_p_values[sample_name][mut_names[mut_idx]], 0.5)
                        # else:
                        #     weight *= 0.9       # coverage data is not available!

                # add mutation to the newly formed pattern
                nodes[new_node].add(mut_idx)
                # add weight to the newly formed pattern
                if new_node in node_weights:
                    node_weights[new_node] += weight
                else:
                    node_weights[new_node] = weight
                unknown_muts[mut_idx][new_node] = weight

            if len(unknown_positions[(node, mut_idx)]) > 3 and length == 3:
                logger.warn('Variant {} has low coverage in too many samples ({}) such '.format(
                    mut_names[mut_idx], len(unknown_positions[(node, mut_idx)]))
                    + 'that all possible combinations can be explored.')
                break

    return nodes, node_weights, unknown_muts, mut_pattern_scores


def _add_evolutionary_node(tree, parent, mp, node_name, mutations, confidence=None):
    """

    Find the right place in the tree to insert the new node (clone)
    There is exactly one such place. Depth-first-search is used.
    :param tree:
    :param parent:
    :param mp:
    :param node_name:
    :param mutations: set of mutation uniquely present in this node
    :param confidence: confidence in branching
    :return:
    """

    # logger.debug('Parent {} - clone {} - mutations {}'.format(parent, clone, mutations))
    # find a node where the given clone is a subset of the node's samples
    # but none of its children is a superset of the given clone => insertion place is found
    for child in tree.successors(parent):
        if mp.issubset(child):
            if _add_evolutionary_node(tree, child, mp, node_name, mutations, confidence=confidence):
                return True

    # clone is not a subset of any existing children
    # hence, a new node has to be created on this level
    tree.add_node(mp, name=node_name, muts=mutations.union(tree.node[parent]['muts']))
    # add confidence value of this branching
    if confidence is not None and mp in confidence:
        tree.node[mp]['conf'] = confidence[mp]

    # check if any of the existing children of the parent are subsets of the new node
    new_grandchildren = set()
    for child in tree.successors(parent):
        if child.issubset(mp):
            new_grandchildren.add(child)

    for grandchild in new_grandchildren:
        helper_mutations = tree[parent][grandchild]['muts']
        tree.remove_edge(parent, grandchild)
        tree.add_edge(mp, grandchild, muts=helper_mutations)
        tree.node[grandchild]['muts'] = tree.node[grandchild]['muts'].union(mutations)
        # logger.debug('{} became a child of {}.'.format(grandchild, clone))

    # insert edge to the newly added clone (node)
    tree.add_edge(parent, mp, muts=mutations)
    # logger.debug('Successfully inserted {} with {} mutations as child of {}.'.format(clone,
    #              len(tree.node[clone]['muts']), parent))

    return True


def _get_html_template():

    filename = "tree_template.html"
    input_dir = "input"
    # derive correct input directory
    if not os.getcwd().endswith('treeomics'):
        template_path = os.path.join(input_dir, filename)
    else:
        template_path = os.path.join('..', input_dir, filename)

    logger.debug('Path to HTML template: '.format(template_path))

    return template_path
