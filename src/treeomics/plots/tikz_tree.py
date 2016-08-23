#!/usr/bin/python
"""Function to generate a tree figure with latex and tikz"""
import logging
import copy
import networkx as nx
import numpy as np
from itertools import islice
import utils.latex_output as latex_output
from phylogeny.phylogeny_utils import TREE_ROOT

__author__ = 'jreiter'

# get logger for application
logger = logging.getLogger('treeomics')


def create_figure_file(tree, tree_root_key, filename, patient, phylogeny, figure_caption, drivers=set(),
                       germline_distance=2.0, standalone=False):
    """
    Takes a networkx tree and creates a tikz source file such that a pdf can be generated with latex
    :param tree:
    :param tree_root_key:
    :param filename: path to output file
    :param patient: data structure around patient
    :param phylogeny:
    :param figure_caption: caption for tikz figure
    :param drivers: optional set of known driver gene names highlighted on each edge
    :return path to the generated latex/tikz tree
    """

    try:
        with open(filename, 'w') as latex_file:

            # generate latex file headers
            latex_output.write_tikz_header(latex_file, germline_distance=germline_distance, standalone=standalone)
            latex_output.write_figure_header(latex_file, standalone)

            # generate tikz tree and write it to the opened file
            _write_tikz_tree(tree, tree_root_key, latex_file, 0, patient, phylogeny,
                             gene_names=patient.gene_names, mut_keys=patient.mut_keys, drivers=drivers)

            # generate latex file footers
            latex_output.write_figure_footer(latex_file, figure_caption, standalone)

            latex_file.write('\\end{document}\n')
            logger.info('Generated tree in tikz latex format: {}'.format(filename))
            return filename

    except OSError:
        logger.exception("Unexpected error occurred during tikz file generation: ")
        return None


def _write_tikz_tree(tree, cur_node, latex_file, level, patient, pg,
                     gene_names=None, mut_keys=None, drivers=set()):
    """
    Run recursively through the tree and write the tree in tikz format to the opened file
    :param tree:
    :param cur_node:
    :param latex_file:
    :param level:
    :param patient:
    :param pg: phylogeny
    :param gene_names:
    :param mut_keys:
    :param drivers:
    :return:
    """

    # maximal number of depicted driver gene names per branch
    MAX_DR_NAS = 6

    tree.node[cur_node]['level'] = level
    # calculate indentation dependent on the level of the node
    pre = ''.join('\t' for _ in range(level))

    if len(tree.neighbors(cur_node)) > 0 or \
            ('name' in tree.node[cur_node] and (cur_node == TREE_ROOT or tree.node[cur_node]['name'] == TREE_ROOT)):

        if tree.node[cur_node]['name'] == TREE_ROOT or cur_node == TREE_ROOT:
            # Remove underscores in patient names as tikz can not process underscores
            tree_name = str(tree.node[cur_node]['name']).replace(' ', '~')
            latex_file.write('\\Tree [.'+tree_name.replace('_', '~')+' \n')
        # internal subclone
        else:
            node_name = str(tree.node[cur_node]['name']).replace(' ', '~')
            if node_name.startswith('Pam'):
                node_name = node_name[5:]
            latex_file.write(pre+'[.\\small{'+node_name+'} \n')

        # add edges to its children and information about the acquired mutations
        # for child in tree.successors(parent):     had to change to be able to handle undirected graph
        for child in sorted(tree.neighbors(cur_node), key=lambda k: tree.node[k]['name']):

            # check if not has been covered already
            if 'level' in tree.node[child] and tree.node[child]['level'] < level:
                continue

            # is any of the acquired mutations in a driver gene?
            if gene_names is not None:
                acquired_drivers = [gene_names[m] for m in tree[cur_node][child]['muts'] if gene_names[m] in drivers]
            else:
                acquired_drivers = list()

            latex_file.write(
                pre+'\t\\edge node[above, NavyBlue]{' +
                (('\\footnotesize{+'+str(len(tree[cur_node][child]['muts'])) +
                 (' ('+(','.join(d for d in sorted(islice(acquired_drivers, 0, MAX_DR_NAS)))) +
                 (',...+{})'.format(len(acquired_drivers)-MAX_DR_NAS) if len(acquired_drivers) > MAX_DR_NAS
                  else ')') if len(acquired_drivers) > 0 else '') + '}')
                 if 'muts' in tree[cur_node][child] and len(tree[cur_node][child]['muts']) > 0 else '') + '}' +
                (' node[below, red!70]{\\footnotesize{'+str(-len(tree[cur_node][child]['dels']))+'}}'
                 if 'dels' in tree[cur_node][child] and len(tree[cur_node][child]['dels']) > 0 else '') +
                (' node[below, black!60]{\\footnotesize{'+'{:.0f}\\%'.format(
                    tree.node[child]['conf']*100.0)+'}}' if 'conf' in tree.node[child] else '') + ';\n ')

            # mutations which have been acquired on the edge from the parent to this child
            if 'muts' in tree[cur_node][child] and len(tree[cur_node][child]['muts']) > 0:
                if gene_names is not None:
                    latex_file.write(pre+'\t% Acquired mutations ({}): '.format(len(tree[cur_node][child]['muts'])) +
                                     ', '.join(sorted(gene_names[m] for m in tree[cur_node][child]['muts']))+'\n')
                if mut_keys is not None:
                    latex_file.write(pre+'\t% Acquired mutations ({}): '.format(len(tree[cur_node][child]['muts'])) +
                                     ', '.join(sorted(mut_keys[m] for m in tree[cur_node][child]['muts']))+'\n')
                latex_file.write(pre+'\t% Acquired mutations ({}): '.format(len(tree[cur_node][child]['muts'])) +
                                 ','.join(sorted((str(m) for m in tree[cur_node][child]['muts']), key=int))+'\n')
                latex_file.write(pre + '\t% VAF of acquired mutations: mean {:.2%}; median {:.2%}'.format(
                    np.mean([patient.vafs[m, pg.sc_sample_ids[sa_idx] if sa_idx in pg.sc_sample_ids else sa_idx]
                             for m in tree[cur_node][child]['muts'] for sa_idx in child]),
                    np.median([patient.vafs[m, pg.sc_sample_ids[sa_idx] if sa_idx in pg.sc_sample_ids else sa_idx]
                               for m in tree[cur_node][child]['muts'] for sa_idx in child])) +
                                 '\n')

            _write_tikz_tree(tree, child, latex_file, level+1, patient, pg, gene_names=gene_names,
                             mut_keys=mut_keys, drivers=drivers)

        latex_file.write(pre+']\n')

    else:                                   # node is a leaf

        sa_idx = (patient.sample_names.index(tree.node[cur_node]['name']) if patient.sc_names is None
                  else patient.sc_names.index(tree.node[cur_node]['name']))
        reported_muts = len(patient.samples[sa_idx])

        # Calculate number of acquired mutations in this sample
        # and compare it to the reported number
        root = None
        for node in tree.nodes_iter():
            if node == TREE_ROOT or tree.node[node]['name'] == TREE_ROOT:
                root = node
        path = nx.shortest_path(tree, source=root, target=cur_node)

        acquired_muts = 0
        lost_muts = 0
        for i in range(len(path[:-1])):
            if 'muts' in tree[path[i]][path[i+1]]:
                acquired_muts += len(tree[path[i]][path[i+1]]['muts'])
            if 'dels' in tree[path[i]][path[i+1]]:
                lost_muts += len(tree[path[i]][path[i+1]]['dels'])

        # Remove underscores in sample names as tikz can not process underscores
        node_name = str(tree.node[cur_node]['name']).replace('_', '~')
        if node_name.startswith('Pam'):
            node_name = node_name[5:]

        if acquired_muts > 0:
            sample_text = '\\node[black,draw,text width=1.09cm,inner sep=2pt,align=center]{' + node_name
            # sample_text += '\\\\ \\small{'+'{}'.format(acquired_muts-lost_muts)+'}'
        else:
            sample_text = '\\node[black,draw,text width=1.09cm,inner sep=2pt,align=center]{' + node_name
            # sample_text += '\\\\ \\small{'+'{}'.format(reported_muts)+'}'

        latex_file.write(pre+sample_text+'}; '+' \n')
        if acquired_muts > 0:
            latex_file.write(pre+'% Present mutations: {}; Reported mutations: {} \n'.format(
                acquired_muts-lost_muts, reported_muts))


def name_internal_nodes(tree, node, level, root_path, prename, idx):
    """
    Recursively run through the tree and annotate the nodes with levels and meaningful names
    Internal nodes get assigned a subclone identifier
    :param tree:
    :param node:
    :param level:
    :param root_path:
    :param prename: prefix to the subclone identifier
    :param idx:
    :return: subclone identifier
    """

    tree.node[node]['level'] = level
    path = copy.deepcopy(root_path)
    path.append(node)
    tree.node[node]['root_path'] = path

    for child in tree.neighbors(node):

        if 'level' not in tree.node[child]:
            if len(tree.neighbors(child)) > 1:
                tree.node[child]['name'] = prename+str(idx)
                idx += 1

            idx = name_internal_nodes(tree, child, level+1, path, prename, idx)

    return idx


def get_dendrogram_clusters(tree, tree_root_key, no_samples, include_germline=True):
    """
    Returns the clusters (matrix with n-1 rows and 4 columns) of a given tree for the scipy dendrogram function
    Each row records which two clusters were combined as the neighbor joining was performed.
    For instance, if the first row is [5., 8., 3.8, 2.], then clusters 5 and 8 were combined
    because they had a distance of 3.8.
    The 2 means that there are two original observations (samples) in the newly formed cluster.
    If the second row is [3., 9., 5.6, 2.], then clusters 3 and 9 were combined
    because they had a distance of 5.6 and there are two original observations in the new cluster etc.
    :param tree: networkx tree (label n is the germline sample)
    :return matrix with n-1 rows and 4 columns
    """

    clusters = []
    tree.node[tree_root_key]['name'] == TREE_ROOT
    _generate_dendrogram_clusters(tree, clusters, tree_root_key, 0, no_samples, include_germline=include_germline)

    return clusters


def _generate_dendrogram_clusters(tree, clusters, node, level, no_samples, include_germline=True):
    """
    Traverse the tree and generate the dendrogram clustering according to the tree topology
    :return triple of cluster id, distance (height), number of samples within the cluster
    """

    if (len(tree.neighbors(node)) > 1 or
            ('name' in tree.node[node] and tree.node[node]['name'] == TREE_ROOT)):     # internal node

        cl_ids = []
        distances = []
        observations = []
        for child in tree.neighbors(node):

            # if the node is an ancestor continue
            if 'level' in tree.node[child] and tree.node[child]['level'] < level:
                continue

            cl_id, d, o = _generate_dendrogram_clusters(tree, clusters, child, level+1,
                                                        no_samples, include_germline=include_germline)
            cl_ids.append(cl_id)
            distances.append(d)
            observations.append(o)

        if tree.node[node]['name'] != TREE_ROOT or include_germline:
            clusters.append([cl_ids[0], cl_ids[1], max(distances)+1, sum(observations)])

        return len(clusters)-1+no_samples, max(distances)+1.0, sum(observations)

    else:                                       # leave node
        return node, 0, 1           # return id of sample, height 0, and number of samples 1
