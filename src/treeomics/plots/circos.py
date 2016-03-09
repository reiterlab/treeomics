#!/usr/bin/python
"""Helper functions to generate input files for circos plots"""
__author__ = 'Johannes REITER'

import logging
import os
from itertools import chain
from plots.plots_utils import _format_gene_name
from phylogeny.simple_phylogeny import SimplePhylogeny
from phylogeny.max_lh_phylogeny import MaxLHPhylogeny

# get logger for application
logger = logging.getLogger('treeomics')


def create_raw_data_file(raw_data_filename, mutations, mut_pos, data=None, sample_names=None):
    """
    Create a space separated file with full list of mutations in the following format:
    chr start end value options
    :param raw_data_filename: output filename
    :param mutations: dictionary with all mutations mapping to the set of samples where it is present
    :param mut_pos: array of tuples with the mutation position (chr, start_pos, end_pos)
    """

    logger.debug('Creating raw mutation data file: {}'.format(raw_data_filename))

    # create directory for data files
    if not os.path.exists(os.path.dirname(raw_data_filename)):
        os.makedirs(os.path.dirname(raw_data_filename))

    # open output file
    with open(raw_data_filename, 'w') as data_file:

        data_file.write('# chr start end value [options]\n')

        # run through all mutations
        for mut_idx in range(len(mutations)):

            if data is None:        # no data for unknown mutations
                # write a line with the mutation to the file for each sample where it is present
                for sa_idx in mutations[mut_idx]:
                    data_file.write('hs{} {} {} {} shared={}\n'.format(
                        mut_pos[mut_idx][0], mut_pos[mut_idx][1], mut_pos[mut_idx][2], sa_idx+1,
                        len(mutations[mut_idx])))
            elif len(mutations[mut_idx]) > 0:       # add information about unknown positions
                for sa_idx in range(len(data[mut_idx])):
                    if data[mut_idx][sa_idx] > 0:
                        data_file.write('hs{} {} {} {} shared={},unknown=0\n'.format(
                            mut_pos[mut_idx][0], mut_pos[mut_idx][1], mut_pos[mut_idx][2], sa_idx+1,
                            len(mutations[mut_idx])))
                    elif data[mut_idx][sa_idx] < 0:
                        data_file.write('hs{} {} {} {} shared={},unknown=1\n'.format(
                            mut_pos[mut_idx][0], mut_pos[mut_idx][1], mut_pos[mut_idx][2], sa_idx+1,
                            len(mutations[mut_idx])))

        if sample_names is not None:    # add supplementary mapping from value to sample name
            data_file.write('# Mapping from value to sample name: {} \n'
                            .format('; '.join(
                                'val{}: {}'.format(sa_idx+1, sa_name) for sa_idx, sa_name in enumerate(sample_names))))

    logger.info('Circos raw mutation data file created: {}'.format(os.path.abspath(
        raw_data_filename)))


def create_mutation_labels_file(labels_filename, mutations, gene_names, mut_pos, driver_pathways):
    """
    Create a space separated file with a full list of mutations in the following format:
    chr start end name options
    :param labels_filename: output filename
    :param mutations: dictionary with all mutations mapping to the set of samples where it is present
    :param gene_names: name/label to each mutation
    :param mut_pos: array of tuples with the mutation position (chr, start_pos, end_pos)
    :param driver_pathways: is this mutation a known driver
    """

    logger.debug('Creating circos mutation labels file: {}'.format(labels_filename))

    # create directory for data files
    if not os.path.exists(os.path.dirname(labels_filename)):
        os.makedirs(os.path.dirname(labels_filename))

    # open output file
    with open(labels_filename, 'w') as data_file:

        data_file.write('# chr start end name [options]\n')

        # run through all mutations and write their position and their label to the file
        for mut_idx in range(len(mutations)):

            # show only shared labels or drivers
            if len(mutations[mut_idx]) > 0 or mut_idx in driver_pathways:
                data_file.write('hs{} {} {} {} shared={},driver={}\n'.format(
                    mut_pos[mut_idx][0], mut_pos[mut_idx][1], mut_pos[mut_idx][2],
                    _format_gene_name(gene_names[mut_idx]), len(mutations[mut_idx]),
                    (str(1) if mut_idx in driver_pathways else str(0))))

    logger.info('Circos mutation label file created: {}'.format(os.path.abspath(
        labels_filename)))


def create_mutation_links_file(links_filename, phylogeny, mut_pos):
    """
    Create a space separated file with the following information:
    chr1 start1 end1 chr2 start2 end2 [options]
    e.g.: hs5 50 150 hs4 500000 500100 color=blue,thickness=5p
    :param links_filename: output filename
    :param phylogeny: data structure around the phylogenetic tree
    :param mut_pos: array of tuples with the mutation position: (chr, start_pos, end_pos)
    """

    logger.debug('Creating circos links file: {}'.format(links_filename))

    # open file
    with open(links_filename, 'w') as links_file:

        links_file.write('# chr1 start1 end1 chr2 start2 end2 [options]\n')

        # run through all edges of conflicting mutation patterns
        for source, sink in phylogeny.cf_graph.edges_iter():

            if source not in phylogeny.patient.mps or sink not in phylogeny.patient.mps:
                # only generate the conflict links among the present and absent mutation patterns
                # do not consider clones determined by unknown positions which are given by the nodes
                # in the actual conflict graph
                continue

            # create links between all mutations present in the two mutation patterns
            for source_mut in phylogeny.patient.mps[source]:

                for sink_mut in phylogeny.patient.mps[sink]:
                    # write link to file in the following format
                    # hs1 450 550 hs2 800100 800200

                    links_file.write('hs{} {} {} '.format(mut_pos[source_mut][0], mut_pos[source_mut][1],
                                                          mut_pos[source_mut][2]))
                    links_file.write('hs{} {} {}\n'.format(mut_pos[sink_mut][0], mut_pos[sink_mut][1],
                                                           mut_pos[sink_mut][2]))

        logger.info('Circos links file created: {}'.format(os.path.abspath(links_filename)))


def create_conflict_graph_files(cfg_nodes_filename, cfg_mutnode_labels_filename, cfg_mutnode_data_filename,
                                cfg_links_filename, phylogeny, gene_names, driver_pathways, data=None,
                                min_node_weight=None, max_no_mps=50):
    """
    Create all input data files for circular conflict graph plots with circos
    :param cfg_nodes_filename: path to output file circos nodes
    :param cfg_mutnode_labels_filename: path to output file circos node labels
    :param cfg_mutnode_data_filename: path to output file circos node data
    :param cfg_links_filename: path to output file containing all links (conflicts) among the nodes
    :param phylogeny: instance of the class phylogeny
    :param gene_names:
    :param driver_pathways:
    :param data:
    :param min_node_weight: minimal reliability score of a mutation pattern to be displayed
    :param max_no_mps: apply min_node_weight if there are more than this number of MPs in the data
    """

    # creates three data files: (i) basic nodes (clones), (ii) labels for mutations per node,
    # (iii) mutation pattern per node
    cfg_nodes = _create_cfg_nodes_files(cfg_nodes_filename, cfg_mutnode_labels_filename, cfg_mutnode_data_filename,
                                        phylogeny, gene_names, driver_pathways, data=data,
                                        min_node_weight=min_node_weight, max_no_mps=max_no_mps)

    _create_cfg_links_file(cfg_links_filename, phylogeny, cfg_nodes)


def _create_cfg_nodes_files(cfg_nodes_filename, cfg_mutnode_labels_filename, cfg_mutnode_data_filename, phylogeny,
                            gene_names, driver_pathways, data=None, min_node_weight=None, max_no_mps=50):
    """
    Create a space separated file with all nodes (clones) in the reduced conflict graph
    The format is as given here:
    # chr - ID LABEL START END COLOR
    chr - 1 (1,2,3) 0 1 lpred
    :param cfg_nodes_filename: output filename for the nodes
    :param cfg_mutnode_labels_filename: output filename for the labels of the mutation present in the node
    :param cfg_mutnode_data_filename: output filename for mutation pattern data in each clone
    :param phylogeny: data structure around the phylogenetic tree
    :param gene_names: gene names of the position of the mutation
    :param driver_pathways: is this mutation a known driver
    :param data:
    :param min_node_weight: minimal reliability score of a mutation pattern to be displayed
    :param max_no_mps: apply min_node_weight if there are more than this number of MPs in the data
    :return: cfg_nodes: mapping from clones (nodes) to node id in the circos files
    """

    logger.debug('Creating circos conflict graph nodes files: {}, {}, {}'.format(cfg_nodes_filename,
                 cfg_mutnode_labels_filename, cfg_mutnode_data_filename))

    # dictionary with mapping from clones (nodes) to node id in the circos files
    cfg_nodes = dict()

    # open file
    with open(cfg_nodes_filename, 'w') as cfg_nodes_file, \
            open(cfg_mutnode_labels_filename, 'w') as cfg_labels_file, \
            open(cfg_mutnode_data_filename, 'w') as cfg_data_file:

        cfg_nodes_file.write('# chr - ID LABEL(=samples) START END COLOR\n')
        cfg_labels_file.write('# NODE START END GENE-NAME [options]\n')
        cfg_data_file.write('# chr - ID START END VALUE [options]\n')

        # run through all nodes of the reduced conflicting conflict graph
        # and generate the circus input with the node size according to the
        # the number of mutations in the clone
        for node_idx, node in enumerate(sorted(phylogeny.nodes.keys(),
                                               key=lambda k: (-len(k), -phylogeny.node_scores[k])), 1):

            cfg_nodes[node] = node_idx
            # write node id to file
            cfg_nodes_file.write('chr - n'+str(node_idx)+' ')

            # if there many nodes then don't show the ones with very small weights

            if len(phylogeny.nodes.keys()) > max_no_mps and \
                    phylogeny.node_scores[node] < min_node_weight:
                cfg_nodes_file.write('<{:.2f} 0 1 grey\n'.format(min_node_weight))
                cfg_labels_file.write('n{} 0 1 {}-{} driver=0,conflicting=0,chosen=0,unknown=1\n'.format(
                    node_idx, node_idx+1, len(phylogeny.cf_graph.nodes())))
                del cfg_nodes[node]
                break

            # write node (clone) id to file which is given by the samples where
            # the mutation is present
            cfg_nodes_file.write('{:.2f} '.format(phylogeny.node_scores[node]))

            # write start and end position of the node to the file
            # end position is given by the number of mutations present in this clone
            cfg_nodes_file.write('0 {} '.format(len(phylogeny.nodes[node])))
            # write node color to the file which is given by the fact if the clone is
            # removed by the solution or kept
            if node in phylogeny.compatible_nodes:
                cfg_nodes_file.write('blue')       # prefix p for pure color, e.g. pblue
            elif node in phylogeny.conflicting_nodes:
                cfg_nodes_file.write('red')         # prefix l or d for light or dark, e.g. lpblue
            else:
                cfg_nodes_file.write('lgrey')
                logger.warn('Mutation pattern {} is not conflicting but only contains unknown variants.'.format(node))

            cfg_nodes_file.write('\n')

            for pos, mut_idx in enumerate(sorted(phylogeny.nodes[node],
                                                 key=lambda k: phylogeny.patient.gene_names[k]
                                                 if len(phylogeny.patient.gene_names) else 0)):
                # create label entry in the according label file if available
                if len(gene_names):
                    cfg_labels_file.write('n{} {} {} {} driver={},conflicting={},chosen={}'.format(
                        node_idx, pos, pos+1, _format_gene_name(gene_names[mut_idx]),
                        (str(1) if mut_idx in driver_pathways else str(0)),
                        1 if mut_idx in phylogeny.conflicting_mutations else 0,
                        1 if node in phylogeny.compatible_nodes and mut_idx in phylogeny.nodes[node] else 0))
                    if node in phylogeny.nodes and mut_idx in phylogeny.nodes[node]:
                        cfg_labels_file.write(',unknown=0\n')
                    else:
                        cfg_labels_file.write(',unknown=1\n')

                # create mutation pattern entry in the according data file
                for sa_idx in node:
                    cfg_data_file.write('n{} {} {} {} shared={}'.format(
                        node_idx, pos, pos+1, sa_idx+1, len(node)))
                    if data is None:
                        cfg_data_file.write('\n')
                    elif data[mut_idx][sa_idx] < 0:
                        cfg_data_file.write(',unknown=1\n')
                    else:
                        cfg_data_file.write(',unknown=0\n')

        logger.info('Circos conflict graph nodes file created: {}'.format(os.path.abspath(cfg_nodes_filename)))

    return cfg_nodes


def create_mlh_graph_files(res_nodes_filename, res_mutnode_labels_filename, res_mutnode_data_filename,
                           data, phylogeny, gene_names, driver_pathways):
    """
    Create a space separated files with all mutation and their corresponding resolved mutation pattern
    The format is as given here:
    # chr - ID LABEL START END COLOR
    chr - 1 (1,2,3) 0 1 lpgreen
    :param res_nodes_filename: output filename for the nodes
    :param res_mutnode_labels_filename: output filename for the labels of the mutation present in the node
    :param res_mutnode_data_filename: output filename for mutation pattern data in each clone
    :param data: holds the detailed mutation data
    :param phylogeny: data structure around the phylogenetic tree
    :param gene_names: gene names of the position of the mutation
    :param driver_pathways: is this mutation a known driver
    :return: res_nodes: mapping from clones (nodes) to node id in the circos files
    """

    logger.debug('Creating circos conflict graph nodes files: {}, {}, {}'.format(res_nodes_filename,
                 res_mutnode_labels_filename, res_mutnode_data_filename))

    # dictionary with mapping from clones (nodes) to node id in the circos files
    res_nodes = dict()

    # open file
    with open(res_nodes_filename, 'w') as res_nodes_file, \
            open(res_mutnode_labels_filename, 'w') as res_labels_file, \
            open(res_mutnode_data_filename, 'w') as res_data_file:

        res_nodes_file.write('# chr - ID LABEL(=samples) START END COLOR\n')
        res_labels_file.write('# NODE START END GENE-NAME [options]\n')
        res_data_file.write('# chr - ID START END VALUE [options]\n')

        founders = dict()
        founding_mp = frozenset(sa for sa in range(phylogeny.patient.n))
        founders[founding_mp] = phylogeny.mlh_founders
        unique_mutations = dict()
        for sa_idx, muts in phylogeny.mlh_unique_mutations.items():
            if len(muts) > 0:
                unique_mutations[frozenset([sa_idx])] = muts

        for node_idx, (mp, muts) in enumerate(
                chain(founders.items(), sorted(phylogeny.shared_mlh_mps.items(), key=lambda k: -len(k[0])),
                      sorted(unique_mutations.items(), key=lambda k: k[0])), 1):

            # write node id to file
            res_nodes_file.write('chr - n'+str(node_idx)+' ')
            # write node (clone) id to file which is given by the samples where
            # the mutation is present
            res_nodes_file.write(','.join(str(sa_idx) for sa_idx in mp)+' ')
            # write start and end position of the node to the file
            # end position is given by the number of mutations present in this clone
            res_nodes_file.write('0 '+str(len(muts))+' ')
            # write node color to the file which is blue since all conflicts are resolved
            res_nodes_file.write('blue')
            res_nodes_file.write('\n')

            # order mutations according to their names within a clone
            for pos, mut_idx in enumerate(
                    sorted(muts, key=lambda k:
                           phylogeny.patient.gene_names[k] if phylogeny.patient.gene_names is not None else 0)):
                # create label entry in the according label file
                res_labels_file.write('n{} {} {} {} driver={} \n'.format(
                    node_idx, pos, pos+1,
                    _format_gene_name(gene_names[mut_idx]) if gene_names is not None else str(mut_idx),
                    str(1) if mut_idx in driver_pathways else str(0)))

                # create mutation pattern entry in the according data file
                for sa_idx in chain(mp, phylogeny.false_positives[mut_idx], phylogeny.false_negatives[mut_idx],
                                    phylogeny.false_negative_unknowns[mut_idx]):
                    res_data_file.write('n'+str(node_idx)+' '+str(pos)+' '+str(pos+1)
                                        + ' '+str(sa_idx+1) + ' '+'shared='+str(len(mp)))
                    # check if this data was given or resolved
                    if mut_idx in phylogeny.false_negatives.keys() \
                            and sa_idx in phylogeny.false_negatives[mut_idx]:

                        # this mutation has been added
                        res_data_file.write(',resolved=1')

                    elif mut_idx in phylogeny.false_positives.keys() \
                            and sa_idx in phylogeny.false_positives[mut_idx]:

                        # mutation has been removed
                        res_data_file.write(',resolved=-1')
                    elif mut_idx in phylogeny.false_negative_unknowns.keys() \
                            and sa_idx in phylogeny.false_negative_unknowns[mut_idx]:

                        # this mutation has been added
                        res_data_file.write(',resolved=1')
                    else:                       # mutation data has not been changed
                        res_data_file.write(',resolved=0')

                    # map from identified putative subclones to their original sample
                    while sa_idx in phylogeny.sc_sample_ids:
                        sa_idx = phylogeny.sc_sample_ids[sa_idx]
                    if data[mut_idx][sa_idx] < 0:
                        res_data_file.write(',unknown=1')
                    else:
                        res_data_file.write(',unknown=0')

                    res_data_file.write('\n')

        logger.info('Circos resolved graph nodes file created: {}'.format(os.path.abspath(res_nodes_filename)))

    return res_nodes


def create_mp_graph_files(mp_nodes_filename, mp_mutnode_data_filename, mp_links_filename,
                          phylogeny, mps, mp_weights, min_node_weight=0.25, max_no_mps=50, pars_infor=False):
    """
    Create the circos data files for the mutation pattern overview graph
    :param mp_nodes_filename: output filename for the nodes
    :param mp_mutnode_data_filename: output filename for mutation pattern data in each clone
    :param mp_links_filename: circos links output filename
    :param phylogeny: data structure around the phylogenetic tree
    :param mps: mutation patterns
    :param mp_weights: mutation pattern reliability scores
    :param min_node_weight: minimal reliability score of a mutation pattern to be displayed
                            if there are more than max_no_mps of MPs
    :param max_no_mps: apply min_node_weight if there are more than this number of MPs
    :param pars_infor: show only parsimony informative mutation patterns
    """

    # creates three data files: (i) basic nodes (clones), (ii) labels for mutations per node,
    # (iii) mutation pattern per node
    mp_nodes = _create_mp_nodes_files(mp_nodes_filename, mp_mutnode_data_filename, phylogeny, mps, mp_weights,
                                      min_node_weight=min_node_weight, max_no_mps=max_no_mps, pars_infor=pars_infor)

    _create_cfg_links_file(mp_links_filename, phylogeny, mp_nodes)


def _create_mp_nodes_files(mp_nodes_filename, mp_mutnode_data_filename, phylogeny, mps, mp_weights,
                           min_node_weight=0.25, max_no_mps=50, pars_infor=False):
    """
    Create a space separated file with all nodes (clones) in the mutation pattern overview graph
    The format is as given here:
    # chr - ID LABEL START END COLOR
    chr - 1 (1,2,3) 0 1 lpred
    :param mp_nodes_filename: output filename for the nodes
    :param mp_mutnode_data_filename: output filename for mutation pattern data in each clone
    :param phylogeny: data structure around the phylogenetic tree
    :param mps: mutation patterns
    :param mp_weights: mutation pattern reliability scores
    :param min_node_weight: minimal reliability score of a mutation pattern to be displayed
                            if there are more than max_no_mps of MPs
    :param max_no_mps: apply min_node_weight if there are more than this number of MPs
    :param pars_infor: show only parsimony informative mutation patterns
    :return: mp_nodes: mapping from mutation patterns (nodes) to node id in the circos files
    """

    logger.debug('Creating circos mutation pattern overview graph nodes files: {}, {}'.format(
        mp_nodes_filename, mp_mutnode_data_filename))

    # dictionary with mapping from clones (nodes) to node id in the circos files
    mp_nodes = dict()

    # open file
    with open(mp_nodes_filename, 'w') as mp_nodes_file, \
            open(mp_mutnode_data_filename, 'w') as mp_data_file:

        mp_nodes_file.write('# chr - ID LABEL(=samples) START END COLOR\n')
        mp_data_file.write('# chr - ID START END VALUE [options]\n')

        # run through all nodes of the reduced conflicting conflict graph
        # and generate the circus input with the node size according to the
        # the number of mutations in the clone
        for node_idx, node in enumerate(sorted(mps, key=lambda k: -mp_weights[k]), 1):

            if pars_infor:
                if len(node) < 2:
                    continue

            mp_nodes[node] = node_idx
            # write node id to file
            mp_nodes_file.write('chr - n'+str(node_idx)+' ')

            # show only a maximal number of mutation pattern
            if node_idx > max_no_mps:
                mp_nodes_file.write(('<{:.2f}'.format(mp_weights[node]) if mp_weights[node] < 10
                                     else '<{:.1f} '.format(mp_weights[node]) if mp_weights[node] < 100
                                     else '<{:.0f} '.format(mp_weights[node]))+' 0 1 grey\n')
                del mp_nodes[node]
                break
            # show only the nodes with a minimum reliability score
            elif mp_weights[node] < min_node_weight:
                mp_nodes_file.write(('<{:.2f}'.format(min_node_weight) if mp_weights[node] < 10
                                     else '<{:.1f} '.format(mp_weights[node]) if mp_weights[node] < 100
                                     else '<{:.0f} '.format(mp_weights[node]))+' 0 1 grey\n')
                del mp_nodes[node]
                break

            # write mutation pattern reliability score (weight) to file
            # mp_nodes_file.write('{}-{:.2f} '.format(len(phylogeny.nodes[node]), mp_weights[node]))
            mp_nodes_file.write('{:.2f} '.format(mp_weights[node]) if mp_weights[node] < 10
                                else '{:.1f} '.format(mp_weights[node]) if mp_weights[node] < 100
                                else '{:.0f} '.format(mp_weights[node]))

            # write start and end position of the node to the file
            # end position is given by the number of mutations present in this clone
            mp_nodes_file.write('0 {} '.format(1))
            # write node color to the file which is given by the fact if the clone is
            # removed by the solution or kept
            if ((isinstance(phylogeny, SimplePhylogeny) and node in phylogeny.compatible_nodes) or
                    (isinstance(phylogeny, MaxLHPhylogeny) and node in phylogeny.compatible_nodes) or len(node) == 0):
                mp_nodes_file.write('blue')       # prefix p for pure color, e.g. pblue
            elif node in phylogeny.conflicting_nodes:
                mp_nodes_file.write('red')         # prefix l or d for light or dark, e.g. lpblue
            else:
                mp_nodes_file.write('lgrey')
                logger.warn('Mutation pattern {} is neither conflicting nor compatible! It has to be either one.'
                            .format(node))

            mp_nodes_file.write('\n')

            # create mutation pattern entry in the according data file
            for sa_idx in node:
                mp_data_file.write('n{} {} {} {} shared={}\n'.format(
                    node_idx, 0, 1, sa_idx+1, len(node)))

        logger.info('Circos conflict graph nodes file created: {}'.format(os.path.abspath(mp_nodes_filename)))

    return mp_nodes


def _create_cfg_links_file(cfg_links_filename, phylogeny, cfg_nodes):
    """
    Create a space separated file with the following information:
    n1 start1 end1 n2 start2 end2 [options]
    e.g.: n1 0 5 n4 0 1 color=blue,thickness=5p
    :param cfg_links_filename: output filename
    :param phylogeny: data structure around the phylogenetic tree
    :param cfg_nodes: mapping from clones (nodes) to the node id in the circos file
    :return:
    """

    logger.debug('Creating circos conflict graph links file: {}'.format(cfg_links_filename))

    # open file
    with open(cfg_links_filename, 'w') as links_file:

        links_file.write('# n1 start1 end1 n2 start2 end2 [options]\n')

        # run through all edges of conflicting clones
        for source, sink in phylogeny.cf_graph.edges_iter():

            if source not in cfg_nodes or sink not in cfg_nodes:
                # logger.debug('Edge endnode is not shown in the circos file and hence the link is also not depicted.')
                continue

            # write link to file
            links_file.write('n{} 0 {} '.format(
                cfg_nodes[source],
                len(phylogeny.cf_graph.node[source]['muts']) if 'muts' in phylogeny.cf_graph.node[source] else '1'))
            links_file.write('n{} 0 {}\n'.format(
                cfg_nodes[sink],
                len(phylogeny.cf_graph.node[sink]['muts'] if 'muts' in phylogeny.cf_graph.node[sink] else '1')))

        logger.info('Circos conflict graph links file created: {}'.format(os.path.abspath(cfg_links_filename)))
