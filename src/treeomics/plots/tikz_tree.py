#!/usr/bin/python
"""Function to generate a tree figure with latex and tikz"""
import logging
import copy
import csv
import math
import networkx as nx
import numpy as np
from itertools import islice
import utils.latex_output as latex_output
from phylogeny.phylogeny_utils import TREE_ROOT
from patient import get_variant_details
import os

__author__ = 'jreiter'

# get logger for application
logger = logging.getLogger('treeomics')

# maximal number of depicted driver gene names per branch
MAX_DR_NAS = 6


def create_figure_file(tree, tree_root_key, filename, patient, phylogeny, figure_caption, driver_vars=None,
                       germline_distance=2.0, standalone=False, variant_filepath=None):
    """
    Takes a networkx tree and creates a tikz source file such that a pdf can be generated with latex
    :param tree:
    :param tree_root_key:
    :param filename: path to output file
    :param patient: data structure around patient
    :param phylogeny:
    :param figure_caption: caption for tikz figure
    :param driver_vars: defaultdict with mutation IDs and instance of driver class to be highlighted on each edge
    :param germline_distance: distance from the germ line (root node) to the first branch in the tikz figure
    :param standalone: Latex/tikz produces pdf only with the figure
    :param variant_filepath: path to output file with information about variants and where they were acquired
    :return path to the generated latex/tikz tree
    """

    try:
        with open(filename, 'w') as latex_file:

            # generate latex file headers
            latex_output.write_tikz_header(latex_file, germline_distance=germline_distance, standalone=standalone)
            latex_output.write_figure_header(latex_file, standalone)

            if variant_filepath is not None:
                var_file = open(variant_filepath, 'w')
                logger.info('Write variant analysis file: {}'.format(variant_filepath))
                var_writer = csv.writer(var_file)

                # writer header of file with variants
                header = ['Chromosome', 'StartPosition', 'EndPosition', 'RefAllele', 'AltAllele', 'GeneSymbol',
                          'MutType', 'Driver', 'CGC_region', 'Phylogeny', 'BranchName']
                header += ['p1_'+sample_name for sample_name in patient.sample_names]
                header.append('BISharingStatus')

                var_writer.writerow(header)
            else:
                var_writer = None
                var_file = None

            # generate tikz tree and write it to the opened file
            _write_tikz_tree(tree, tree_root_key, latex_file, 0, patient, phylogeny, variant_writer=var_writer,
                             gene_names=patient.gene_names, mut_keys=patient.mut_keys, drivers=driver_vars)

            if var_file is not None:
                var_file.close()

            # generate latex file footers
            latex_output.write_figure_footer(latex_file, figure_caption, standalone)

            latex_file.write('\\end{document}\n')
            logger.info('Generated tree in tikz latex format: {}'.format(filename))

            # study mutation signatures by investigating the signatures on the individual branches
            # mat_file, extension = os.path.splitext(filename)
            # mat_file += '_matrix.txt'
            # logger.info("Writing matrix file!")
            # with open(mat_file, 'w') as matrix_file:
            #     _write_assigned_mutation_matrix(tree, tree_root_key, matrix_file, 0, patient, drivers=drivers,
            #                                     gene_names=patient.gene_names, mut_keys=patient.mut_keys)

            return filename

    except OSError:
        logger.exception("Unexpected error occurred during tikz file generation: ")
        return None


def _write_tikz_tree(tree, cur_node, latex_file, level, patient, pg,
                     gene_names=None, mut_keys=None, drivers=None, variant_writer=None):
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
    :param variant_writer: CSV writer to store relevant information about a variant and its presence to a file
    :return:
    """

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
            if gene_names is not None and drivers is not None:
                acquired_drivers = [gene_names[m] for m in tree[cur_node][child]['muts'] if m in drivers]
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
                    np.nanmean([patient.vafs[m, pg.sc_sample_ids[sa_idx] if sa_idx in pg.sc_sample_ids else sa_idx]
                               for m in tree[cur_node][child]['muts'] for sa_idx in child]),
                    np.nanmedian([patient.vafs[m, pg.sc_sample_ids[sa_idx] if sa_idx in pg.sc_sample_ids else sa_idx]
                                 for m in tree[cur_node][child]['muts'] for sa_idx in child])) +
                                 '\n')

                # write acquired variants to variant analysis file
                if variant_writer is not None:
                    _write_muts(variant_writer, pg, patient, tree[cur_node][child]['muts'], child, driver_vars=drivers)

            _write_tikz_tree(tree, child, latex_file, level+1, patient, pg, gene_names=gene_names,
                             mut_keys=mut_keys, drivers=drivers, variant_writer=variant_writer)

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


def _write_assigned_mutation_matrix(tree, cur_node, matrix_file, level, patient, gene_names=None, mut_keys=None,
                                    drivers=set()):
    """
    Run recursively through the tree and write mutations per branch to a file to study mutation signatures
    :author Jeff Gerold
    """

    tree.node[cur_node]['level'] = level
    # calculate indentation dependent on the level of the node
    pre = ''.join('\t' for _ in range(level))

    if len(tree.neighbors(cur_node)) > 0 or \
            ('name' in tree.node[cur_node] and (cur_node == TREE_ROOT or tree.node[cur_node]['name'] == TREE_ROOT)):

        if tree.node[cur_node]['name'] == TREE_ROOT or cur_node == TREE_ROOT:
            # Remove underscores in patient names as tikz can not process underscores
            tree_name = str(tree.node[cur_node]['name']).replace(' ', '~')

        # internal subclone
        else:
            node_name = str(tree.node[cur_node]['name']).replace(' ', '~')
            if node_name.startswith('Pam'):
                node_name = node_name[5:]

        # add edges to its children and information about the acquired mutations
        # for child in tree.successors(parent):     had to change to be able to handle undirected graph
        for child in sorted(tree.neighbors(cur_node), key=lambda k: tree.node[k]['name']):

            # check if not has been covered already
            if 'level' in tree.node[child] and tree.node[child]['level'] < level:
                continue

            # mutations which have been acquired on the edge from the parent to this child
            if 'muts' in tree[cur_node][child] and len(tree[cur_node][child]['muts']) > 0:
                mut_key_list = [patient.mut_keys[m] for m in tree[cur_node][child]['muts']]
                for x in mut_key_list:
                    if len(x) >= 3:
                        reference = x[-3]
                        alt = x[-1]
                        if cur_node == 'germline':
                            parent_str = cur_node
                        else:
                            parent_str = ' '.join(sorted(patient.sample_names[int(s)] for s in cur_node))
                        child_str = ' '.join(sorted(patient.sample_names[int(s)] for s in child))
                        outstr = parent_str + '\t' + child_str + '\t' + reference + '\t' + alt + '\t' + x + '\n'
                        matrix_file.write(outstr)

            _write_assigned_mutation_matrix(tree, child, matrix_file, level+1, patient, gene_names=gene_names,
                                            mut_keys=mut_keys, drivers=drivers)

    else:
        pass  # node is a leaf


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


def _write_muts(var_writer, pg, patient, muts, sample_ids, driver_vars=None):
    """
    Write variants that were acquired on the given branch to a CSV file
    :param var_writer:  CSV writer
    :param pg: phylogeny
    :param patient: instance of class Patient
    :param muts: list of integer ids of the mutations that were acquired on this branch
    :param sample_ids: set of sample ids of this branch
    :param driver_vars: dictionary with mut IDs as keys
    """

    # CSV file has the following columns
    # 'Chromosome', 'StartPosition', 'EndPosition', 'RefAllele', 'AltAllele', 'GeneSymbol',
    # 'MutType', 'Driver', 'CGC_region', 'Phylogeny', 'BranchName']

    for mut_idx in muts:

        row = list()
        mut_key = patient.mut_keys[mut_idx]
        chrom, start_pos, end_pos, ref, alt = get_variant_details(mut_key)
        row.append(chrom)
        row.append(start_pos)
        row.append(end_pos)
        row.append(ref)
        row.append(alt)

        if patient.gene_names is not None:
            row.append(patient.gene_names[mut_idx])
        else:
            row.append('')

        if patient.mut_types is not None:
            row.append(patient.mut_types[mut_idx])
        else:
            row.append('')

        if driver_vars is not None:
            if mut_idx in driver_vars.keys():
                row.append(str(True))
                row.append(str(driver_vars[mut_idx].cgc_driver))
            else:
                row.append(str(False))
                row.append('')
        else:
            row.append(str(False))
            row.append('')

        phylogeny_annotations = list()
        # determine branch status
        if len(sample_ids) == 1:
            sa_name = patient.sample_names[list(sample_ids)[0]]
            if 'TM' in sa_name:
                phylogeny_annotations.append('TrMetPrivate')
            elif 'M' in sa_name and 'TM' not in sa_name:
                phylogeny_annotations.append('UntrMetPrivate')
            elif 'PT' in sa_name or 'Primary' in sa_name:
                phylogeny_annotations.append('PTPrivate')
            else:
                phylogeny_annotations.append('Private')

        elif len(sample_ids) == 0:
            phylogeny_annotations.append('Absent')
            raise RuntimeError('Mutation {} is not assigned to any samples!'.format(mut_key))

        else:
            if all(sa_idx in sample_ids for sa_idx in range(len(patient.sample_names))):
                phylogeny_annotations.append('Trunk')
            else:
                if (all(sa_idx in sample_ids for sa_idx, sa_name in enumerate(patient.sample_names)
                        if 'M' in sa_name and 'TM' not in sa_name)
                        and any('M' in sa_name and 'TM' not in sa_name for sa_name in patient.sample_names)):
                    phylogeny_annotations.append('UntrMetTrunk')

                elif (all(sa_idx in sample_ids for sa_idx, sa_name in enumerate(patient.sample_names) if 'M' in sa_name)
                      and any('M' in sa_name for sa_name in patient.sample_names)):
                    phylogeny_annotations.append('MetTrunk')

                elif any('M' in patient.sample_names[sa_idx] and 'TM' not in patient.sample_names[sa_idx]
                         for sa_idx in sample_ids):
                    phylogeny_annotations.append('UntrMetShared')

                elif any('M' in patient.sample_names[sa_idx] for sa_idx in sample_ids):
                    phylogeny_annotations.append('MetShared')

                if (all(sa_idx in sample_ids for sa_idx, sa_name in enumerate(patient.sample_names)
                        if 'PT' in sa_name or 'Primary' in sa_name)
                        and any('PT' in sa_name or 'Primary' in sa_name for sa_name in patient.sample_names)):
                    phylogeny_annotations.append('PTTrunk')

                if len(phylogeny_annotations) == 0:
                    phylogeny_annotations.append('Shared')

        row.append('|'.join(pa for pa in phylogeny_annotations))

        branch_name = '__'.join(sorted(
            patient.sample_names[pg.sc_sample_ids[sa_idx] if sa_idx in pg.sc_sample_ids else sa_idx]
            for sa_idx in sample_ids))
        row.append(branch_name)

        row += ['{:.9g}'.format(math.exp(patient.log_p01[mut_idx][sa_idx][1])) for sa_idx in range(
            len(patient.sample_names))]

        row.append(patient.sharing_status[mut_idx])

        var_writer.writerow(row)
