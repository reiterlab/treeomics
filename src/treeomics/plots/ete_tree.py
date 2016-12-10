#!/usr/bin/python
"""Function to generate a tree figure with ete3 and Qt4"""
import logging
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, TextFace, CircleFace
import networkx as nx
import numpy as np
from itertools import islice
from phylogeny.phylogeny_utils import TREE_ROOT

__author__ = 'jreiter'

# helpful documentation:
# http://etetoolkit.org/docs/latest/reference/reference_treeview.html
# http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html

# get logger for application
logger = logging.getLogger('treeomics')

# maximal number of depicted driver gene names per branch
MAX_DR_NAS = 6


def create_tree(tree, tree_root_key, filepath, patient, phylogeny, drivers=set()):
    """
    Takes a networkx tree and creates an ete3 tree and creates a PDF, SVG or PNG according to the given filename at
    the given path
    :param tree: inferred cancer phylogeny encoded as networkx tree
    :param tree_root_key: key of root in networkx tree
    :param filepath: path to output file (excluding file extension)
    :param patient: data structure around patient
    :param phylogeny: inferred cancer phylogeny
    :param drivers: optional set of known driver gene names highlighted on each edge
    :return path to the generated ete3 PNG file
    """

    # generate ete3 tree
    rt = Tree()   # generate root
    rt.name = 'Germline {}'.format(patient.name)
    _generate_ete_tree(tree, tree_root_key, rt, 0, patient, phylogeny,
                       gene_names=patient.gene_names, drivers=drivers)

    ts = TreeStyle()
    ts.layout_fn = ete_layout
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.extra_branch_line_color = 'Black'
    # XX pixels per branch length unit
    ts.scale = 5 if len(patient.mut_keys) < 500 else (1 if len(patient.mut_keys) < 500 else 0.1)
    ts.branch_vertical_margin = 15  # XX pixels between adjacent branches

    # rt.show(tree_style=ts)

    ete_tree_plot = filepath+'.png'
    rt.render(ete_tree_plot, w=183, units="mm", tree_style=ts, dpi=300)
    rt.render(filepath+'.pdf', w=183, units="mm", tree_style=ts, dpi=300)

    return ete_tree_plot


def _generate_ete_tree(tree, cur_node, ete_cur_node, level, patient, pg,
                       gene_names=None, drivers=set()):
    """
    Run recursively through the tree and write the tree in tikz format to the opened file
    :param tree: inferred cancer phylogeny encoded as networkx tree
    :param cur_node: currently processed node in the networkx tree
    :param ete_cur_node: currently being generated node in the ete tree
    :param level: level in the networkx tree starting with 0 at the root (normal sample)
    :param patient: data structure around the input data
    :param pg: instance of the phylogeny
    :param gene_names: list with the gene names associated with the given list of mutations
    :param drivers: names of genes associated with a likely selective advantage
    :return:
    """

    tree.node[cur_node]['level'] = level

    if len(tree.neighbors(cur_node)) > 0 or \
            ('name' in tree.node[cur_node] and (cur_node == TREE_ROOT or tree.node[cur_node]['name'] == TREE_ROOT)):

        # generate children and information about the acquired mutations
        # for child in tree.successors(parent):     had to change to be able to handle undirected graph
        for child in sorted(tree.neighbors(cur_node), key=lambda k: tree.node[k]['name']):

            # check if node has been covered already
            if 'level' in tree.node[child] and tree.node[child]['level'] < level:
                continue

            node_name = str(tree.node[child]['name']).replace('_', ' ')
            if node_name.startswith('Pam'):
                node_name = node_name[5:]
            new_n = ete_cur_node.add_child(name=node_name, dist=len(tree[cur_node][child]['muts']))
            if 'conf' in tree.node[child]:          # set confidence or support value for this branch
                # new_n.support = tree.node[child]['conf']
                str_conf = ('{:.1%}'.format(tree.node[child]['conf']) if tree.node[child]['conf'] < 0.999
                            else '>99.9%')
            else:
                # new_n.support = float('nan')
                str_conf = ''
            new_n.add_features(confidence=str_conf)

            # is any of the acquired mutations in a driver gene?
            if gene_names is not None:
                acquired_drivers = [gene_names[m] for m in tree[cur_node][child]['muts'] if gene_names[m] in drivers]
                formatted_drs = ((','.join(d for d in sorted(islice(acquired_drivers, 0, MAX_DR_NAS)))) +
                                 (',...+{})'.format(len(acquired_drivers)-MAX_DR_NAS)
                                  if len(acquired_drivers) > MAX_DR_NAS else '') if len(acquired_drivers) > 0 else '')
                new_n.add_features(drivers=formatted_drs+'  ')

            new_n.add_features(total_muts=len(tree[cur_node][child]['muts']))
            mean_vaf = np.nanmean([patient.vafs[m, pg.sc_sample_ids[sa_idx] if sa_idx in pg.sc_sample_ids
                                   else sa_idx] for m in tree[cur_node][child]['muts'] for sa_idx in child])
            new_n.add_features(mean_vaf=mean_vaf)
            median_vaf = np.nanmedian([patient.vafs[m, pg.sc_sample_ids[sa_idx] if sa_idx in pg.sc_sample_ids
                                       else sa_idx] for m in tree[cur_node][child]['muts'] for sa_idx in child])
            new_n.add_features(median_vaf=median_vaf)

            _generate_ete_tree(tree, child, new_n, level + 1, patient, pg, gene_names=gene_names,
                               drivers=drivers)

    else:                                   # node is a leaf
        # mutations inferred to be present in that sample by the statistical model ignoring cancer evolution
        sa_idx = (patient.sample_names.index(tree.node[cur_node]['name']) if patient.sc_names is None
                  else patient.sc_names.index(tree.node[cur_node]['name']))
        classified_muts = len(patient.samples[sa_idx])
        ete_cur_node.add_features(classified_muts=classified_muts)

        # Calculate number of acquired mutations in this sample
        # and compare it to the classified number
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
        ete_cur_node.add_features(acquired_muts=acquired_muts)


def ete_layout(node):
    """
    Formatting of tree nodes while tree is rendered
    :param node: ete node
    """

    nstyle = NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 0

    nstyle["vt_line_width"] = 1     # line width of vertical branches
    nstyle["hz_line_width"] = 1     # line width of horizontal branches

    node.set_style(nstyle)

    if 'median_vaf' in node.features and not np.isnan(node.median_vaf) and logger.isEnabledFor(logging.DEBUG):
        # # Creates a sphere face whose size is proportional to the given VAF or CCF of the acquired muts
        # cf = CircleFace(radius=node.median_vaf*50, color="RoyalBlue", style="sphere")
        # # Let's make the sphere transparent
        # cf.opacity = 0.6
        # # And place as a float face over the tree
        # faces.add_face_to_node(cf, node, column=0, position="float-behind")

        vaf_face = TextFace('({:.0%}) '.format(node.median_vaf), fsize=9, fgcolor='SlateGray')
        vaf_face.opacity = 0.8  # from 0 to 1
        faces.add_face_to_node(vaf_face, node, column=2, position="branch-top")

    if node.is_root():
        # tf_root = TextFace(node.name, fsize=14, fgcolor="Black")
        # tf_root.margin_bottom = 25
        # faces.add_face_to_node(tf_root, node, column=0, position='branch-right')
        return

    bt_face_muts = AttrFace("dist", fsize=11, fgcolor="Black", text_prefix=' ', text_suffix=' ', formatter='%d')
    bt_face_muts.margin_bottom = 2
    bt_face_muts.margin_top = 2
    faces.add_face_to_node(bt_face_muts, node, column=1, position="branch-top")

    if 'drivers' in node.features and len(node.drivers) > 0:
        drivers_face = TextFace(node.drivers, fsize=13, fstyle='italic', fgcolor='Orange')
        drivers_face.opacity = 0.8  # from 0 to 1
        faces.add_face_to_node(drivers_face, node, column=0, position="branch-top")

    # If node is a leaf, add the nodes name and a its scientific name
    if node.is_leaf():

        # for colors see: http://etetoolkit.org/docs/latest/reference/reference_treeview.html#color-names
        leaf_face = AttrFace("name", fsize=14, fgcolor="Black", text_prefix=' ', text_suffix=' ')
        if node.name.startswith('PT'):          # primary tumor sample
            leaf_face.fgcolor = 'SteelBlue'
        elif node.name.startswith('LiM'):       # liver met
            leaf_face.fgcolor = 'SaddleBrown'
        elif node.name.startswith('LuM'):       # lung met
            leaf_face.fgcolor = 'DarkOliveGreen'
        elif node.name.startswith('NoM'):       # lymph node met
            leaf_face.fgcolor = 'Magenta'
        elif node.name.startswith('BrM'):       # brain met
            leaf_face.fgcolor = 'Crimson'

        leaf_face.border.type = 0
        leaf_face.border.width = 1
        leaf_face.margin_bottom = 1
        leaf_face.margin_top = 1
        faces.add_face_to_node(leaf_face, node, column=1)  # , position="aligned"

    else:     # inner node
        # if not np.isnan(node.support):
            # bt_face_conf = AttrFace("support", fsize=12, fgcolor="DimGrey", text_prefix=' ',
            #                         text_suffix=' ', formatter='%0.2f')
        if node.confidence is not None and node.confidence != '':
            tf_conf = TextFace(node.confidence, fsize=12, fgcolor='DimGrey')
            tf_conf.margin_left = 3
            tf_conf.margin_right = 3
            faces.add_face_to_node(tf_conf, node, column=1, position="branch-bottom")
