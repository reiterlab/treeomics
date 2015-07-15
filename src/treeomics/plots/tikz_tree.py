__author__ = 'jreiter'

import logging
import copy
import networkx as nx
import utils.latex_output as latex_output
from phylogeny.phylogeny import TREE_ROOT
from utils.int_settings import NEG_UNKNOWN, POS_UNKNOWN

# get logger for application
logger = logging.getLogger('treeomics')


def create_figure_file(tree, tree_root_key, filename, patient, figure_caption, standalone=False):
    """
    Takes a networkx tree and creates a tikz source file such that a pdf can be generated with latex
    """

    try:
        with open(filename, 'w') as latex_file:

            # generate latex file headers
            latex_output.write_tikz_header(latex_file, standalone)
            latex_output.write_figure_header(latex_file, standalone)

            # generate tikz tree and write it to the opened file
            _write_tikz_tree(tree, tree_root_key, latex_file, 0, patient, patient.gene_names)

            # generate latex file footers
            latex_output.write_figure_footer(latex_file, figure_caption, standalone)

            latex_file.write('\\end{document}\n')
            logger.info('Generated tree in tikz latex format: {}'.format(filename))

    except OSError:
        logger.exception("Unexpected error occurred during tikz file generation: ")


def _write_tikz_tree(tree, cur_node, latex_file, level, patient, gene_names):
    """
    Run recursively through the tree and write the tree in tikz format to the opened file
    """

    tree.node[cur_node]['level'] = level
    # calculate indentation dependent on the level of the node
    pre = ''.join('\t' for _ in range(level))

    if len(tree.neighbors(cur_node)) > 0 or \
            ('name' in tree.node[cur_node] and (cur_node == TREE_ROOT or tree.node[cur_node]['name'] == TREE_ROOT)):

        if tree.node[cur_node]['name'] == TREE_ROOT or cur_node == TREE_ROOT:
            # Remove underscores in patient names as tikz can not process underscores
            tree_name = str(tree.node[cur_node]['name']).replace(' ', '~')
            latex_file.write('\\Tree [.'+tree_name.replace('_', '')+' \n')
        # internal subclone
        else:
            node_name = str(tree.node[cur_node]['name']).replace(' ', '~')
            if node_name.startswith('Pam'):
                node_name = node_name[5:]
            latex_file.write(pre+'[.\\small{'+node_name+'} \n')

        # add edges to its children and information about the acquired mutations
        # for child in tree.successors(parent):     had to change to be able to handle undirected graph
        for child in tree.neighbors(cur_node):

            # check if not has been covered already
            if 'level' in tree.node[child] and tree.node[child]['level'] < level:
                continue

            latex_file.write(pre+'\t\\edge node[above, blue!70]{'
                             + ('\\footnotesize{+'+str(len(tree[cur_node][child]['muts']))+'}'
                                if 'muts' in tree[cur_node][child] and len(tree[cur_node][child]['muts']) > 0 else '')
                             + '}'
                             + (' node[below, red!70]{\\footnotesize{'+str(-len(tree[cur_node][child]['dels']))+'}}'
                                if 'dels' in tree[cur_node][child] and len(tree[cur_node][child]['dels']) > 0 else '')
                             + ';\n ')

            # mutations which have been acquired on the edge from the parent to this child
            if 'muts' in tree[cur_node][child] and len(tree[cur_node][child]['muts']) > 0:
                if len(gene_names):
                    latex_file.write(pre+'\t% Acquired mutations ({}): '.format(len(tree[cur_node][child]['muts']))
                                     + ', '.join(sorted(gene_names[m] for m in tree[cur_node][child]['muts']))+'\n')
                latex_file.write(pre+'\t% Acquired mutations ({}): '.format(len(tree[cur_node][child]['muts']))
                                 + ','.join(sorted((str(m) for m in tree[cur_node][child]['muts']), key=int))+'\n')

            # mutations which have been lost on the edge from the parent to this child
            if 'dels' in tree[cur_node][child] and len(tree[cur_node][child]['dels']) > 0:
                if len(gene_names):
                    latex_file.write(pre+'\t% Lost mutations ({}): '.format(len(tree[cur_node][child]['dels']))
                                     + ', '.join(sorted(gene_names[m] for m in tree[cur_node][child]['dels']))+'\n')
                latex_file.write(pre+'\t% Lost mutations ({}): '.format(len(tree[cur_node][child]['dels']))
                                 + ','.join(sorted((str(m) for m in tree[cur_node][child]['dels']), key=int))+'\n')

            _write_tikz_tree(tree, child, latex_file, level+1, patient, gene_names)

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
        node_name = str(tree.node[cur_node]['name']).replace('_', '')
        if node_name.startswith('Pam'):
            node_name = node_name[5:]

        if acquired_muts > 0:
            sample_text = '\\node[black,draw,text width=1.09cm,inner sep=2pt,align=center]{' + node_name + '\\\\'
            sample_text += '\\small{'+'{}'.format(acquired_muts-lost_muts)+'}'
        else:
            sample_text = '\\node[black,draw,text width=1.09cm,inner sep=2pt,align=center]{' + node_name + '\\\\'
            sample_text += '\\small{'+'{}'.format(reported_muts)+'}'

        latex_file.write(pre+sample_text+'}; '+' \n')
        if acquired_muts > 0:
            latex_file.write(pre+'% Present mutations: {}; Reported mutations: {} \n'.format(
                acquired_muts-lost_muts, reported_muts))


def compute_parsimony_score(tree, patient):
    """
    Compute the parsimony score for a given tree topology and sequencing data of a patient
    :param tree: tree topology describing the evolution of the samples
    :param patient: sequencing for multiple samples
    :return: tuple of the parsimony score and the number of required deletions in a perfect but not persistent phylogeny
    """

    deletions = 0
    score = 0

    for v1, v2 in tree.edges_iter():
        tree[v1][v2]['muts'] = set()
        tree[v1][v2]['dels'] = set()

    # update sample nodes with unknowns
    for sa_idx in range(patient.n):
        tree.node[sa_idx]['pos_unknowns'] = set()
        tree.node[sa_idx]['neg_unknowns'] = set()

        tree.node[sa_idx]['pos_u_pos'] = set()
        tree.node[sa_idx]['pos_u_neg'] = set()
        tree.node[sa_idx]['neg_u_pos'] = set()
        tree.node[sa_idx]['neg_u_neg'] = set()

    # assign unknown mutations to their sample nodes
    for mut in patient.mutations.keys():
        for sa_idx in range(patient.n):
            if patient.data[mut][sa_idx] == POS_UNKNOWN:
                tree.node[sa_idx]['pos_unknowns'].add(mut)
            elif patient.data[mut][sa_idx] == NEG_UNKNOWN:
                tree.node[sa_idx]['neg_unknowns'].add(mut)

    # founder mutations appear on the edge from the germline to the founding clone
    founding_node = [v for v in tree.neighbors(patient.n)][0]
    for mut in patient.founders:
        tree[patient.n][founding_node]['muts'].add(mut)
    # mutations present in only one sample appear on the edge to this one sample
    for mut in patient.shared_muts[1]:
        sa_idx = [sa for sa in patient.mutations[mut]][0]
        parent = [v for v in tree.neighbors(sa_idx)][0]
        tree[parent][sa_idx]['muts'].add(mut)

    # run through all clones (set of mutations which share the samples - mutation pattern)
    # and calculate their appearance and their deletions
    for clone, muts in patient.clones.items():

        anc = _get_common_ancestor(tree, clone)
        parent = [father for father in tree.neighbors(anc)
                  if tree.node[father]['level'] < tree.node[anc]['level']]
        assert len(parent) == 1

        for mut in muts:        # mutations appear at the latest common ancestor
            tree[parent[0]][anc]['muts'].add(mut)

        # run recursively from the latest common ancestor to the leaves and calculate the deletions
        if 1 < len(clone) < patient.n:
            sc, dels, stats = _compute_variant_score(tree, anc, parent[0], clone, muts)

            score += sc
            deletions += sum(dels)

            # update events in all descendants to variant is present if it was unknown
            descendants = [child for child in tree.neighbors(anc)
                           if tree.node[child]['level'] > tree.node[anc]['level']]

            while len(descendants):
                tmp_descendants = []
                for desc in descendants:
                    for mut in muts:
                        if 'pos_unknowns' in tree.node[desc] and mut in tree.node[desc]['pos_unknowns']:
                            tree.node[desc]['pos_unknowns'].remove(mut)
                            tree.node[desc]['pos_u_pos'].add(mut)
                            # score does not change

                        elif 'neg_unknowns' in tree.node[desc] and mut in tree.node[desc]['neg_unknowns']:
                            tree.node[desc]['neg_unknowns'].remove(mut)
                            tree.node[desc]['neg_u_pos'].add(mut)
                            # score has already been changed in the leave

                    for child in tree.neighbors(desc):
                        if tree.node[child]['level'] > tree.node[desc]['level']:
                            tmp_descendants.append(child)

                descendants = tmp_descendants

    logger.debug('Tree needs {} individual deletions. Score: {}'.format(deletions, score))
    tree.graph['total_dels'] = deletions

    return score, deletions


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


def _compute_variant_score(tree, node, parent, clone, muts):
    """
    Recursively compute where in the tree mutations appear and get deleted
    :param tree: tree topology
    :param node: node in the tree
    :param parent: parent of the given node
    :param clone: set of samples where these mutations are present
    :param muts: set of mutations present in the given clone
    :return: tuple of parsimony score, number of deletions, and mutation status
    """

    children = [child for child in tree.neighbors(node) if tree.node[child]['level'] > tree.node[node]['level']]

    score = 0
    # node represents a leave
    if len(children) == 0:      # leave node corresponding to a sample
        if node in clone:       # mutations are present in this sample
            return 0, [0 for _ in range(len(muts))], [1 for _ in range(len(muts))]
        else:
            # return tuple of score, list with deletions, list of status
            deletions = []
            statuses = []

            for mut in muts:
                if mut in tree.node[node]['pos_unknowns']:        # mutation is believed to be present
                    deletions.append(0)
                    statuses.append(-1)
                    # score does not change

                elif mut in tree.node[node]['neg_unknowns']:      # mutation is believed to be absent
                    deletions.append(0)
                    statuses.append(-1)
                    score += EPSILON
                else:                                             # mutation is absent
                    tree[parent][node]['dels'].add(mut)
                    deletions.append(1)
                    statuses.append(0)
                    score += 1

            return score, deletions, statuses

    # internal node
    else:
        # process internal node
        deletions = [0 for _ in range(len(muts))]      # start with no deletions
        statuses = [-1 for _ in range(len(muts))]      # start with status unknown

        # calculate recursively events in the children
        for child in children:
            sc, d, st = _compute_variant_score(tree, child, node, clone, muts)

            score += sc
            for i in range(len(muts)):
                deletions[i] += d[i]
                if st[i] == 1:       # mutation is present downstream
                    statuses[i] = 1
                elif st[i] == 0 and statuses[i] == -1:   # mutation is absent and was unknown in the siblings
                    statuses[i] = 0

        for i, mut in enumerate(muts):
            # if mutation is absent in some and unknown in the remaining descendants,
            # => only one deletion is required
            if statuses[i] == 0:
                score -= deletions[i] + 1
                deletions[i] = 1
                tree[parent][node]['dels'].add(mut)

                # update events in all descendants to variant is absent
                for desc in children:
                    score = _upate_mutation_status(tree, desc, node, mut, score)

        return score, deletions, statuses


def _upate_mutation_status(tree, node, parent, mut, score):

    if mut in tree[parent][node]['dels']:
        tree[parent][node]['dels'].remove(mut)

    elif 'pos_unknowns' in tree.node[node] and mut in tree.node[node]['pos_unknowns']:
        tree.node[node]['pos_unknowns'].remove(mut)
        tree.node[node]['pos_u_neg'].add(mut)
        score += EPSILON

    elif 'neg_unknowns' in tree.node[node] and mut in tree.node[node]['neg_unknowns']:
        tree.node[node]['neg_unknowns'].remove(mut)
        tree.node[node]['neg_u_neg'].add(mut)
        # score has already been changed in the leave

    for child in tree.neighbors(node):
        if tree.node[child]['level'] > tree.node[node]['level']:
            _upate_mutation_status(tree, child, node, mut, score)

    return score


def _get_common_ancestor(tree, nodes):
    """
    Returns the latest common ancestor of a given list of nodes
    :param tree: graph (tree topology annoted with levels)
    :param nodes: list of nodes
    :return: latest common ancestor
    """

    if len(nodes) == 1:
        for node in nodes:
            return node

    root_paths = []
    shortest_path = tree.order()
    for node in nodes:
        root_paths.append(tree.node[node]['root_path'])
        if len(tree.node[node]['root_path']) < shortest_path:
            shortest_path = len(tree.node[node]['root_path'])

    anc = -1
    for level in range(shortest_path):
        anc = root_paths[0][level]
        for path in root_paths:

            # check if this node on the path from the root to the leave is different
            # if yes, last node was the latest common ancestor
            if path[level] != anc:
                logger.debug('Latest common ancestor of nodes {} is {}.'.format(
                    ', '.join(str(n) for n in nodes), path[level - 1]))
                return path[level - 1]

    return anc


def get_newick_representation(tree):
    """
    Returns the newick representation of a given tree
    """

    return _generate_newick_representation(tree, TREE_ROOT)


def _generate_newick_representation(tree, node):
    """
    Traverse the tree and generate a newick representation
    """

    if len(tree.successors(node)) > 0:          # internal node
        output = '('
        for child_idx, child in enumerate(tree.successors(node)):
            output += (',' if child_idx > 0 else '') + _generate_newick_representation(tree, child)

        output += ')'
        if len(tree.predecessors(node)):
            node_name = ''.join(str(leave) for leave in sorted(node))
        else:       # if it is the root (=germline)
            node_name = node

        output += node_name
        if len(tree.predecessors(node)):
            output += ':'+str(sum(1 for _ in tree[tree.predecessors(node)[0]][node]['muts']))

        return output

    else:                                       # leave node
        return ''.join(str(leave) for leave in sorted(node))+':'\
               + str(sum(1 for _ in tree[tree.predecessors(node)[0]][node]['muts']))


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