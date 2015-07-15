#!/usr/bin/python
"""Helper functions to generate latex source code for a tikz tree plot"""
__author__ = 'jreiter'

import logging
from phylogeny.phylogeny import TREE_ROOT

# get logger for application
logger = logging.getLogger('treeomics')


def write_tikz_header(latex_file, standalone=False):
    """
    Write the latex and tikz header to the given and already opened file
    :param latex_file: already opened file writer instance
    :param standalone: plot cut on the border of the figure instead of regular paper format
    """

    latex_file.write('\\documentclass'+('[border=1pt]{standalone}' if standalone else '[11pt]{article}')+'\n')
    latex_file.write('\\usepackage[usenames,dvipsnames]{xcolor}\n')
    latex_file.write('\\usepackage{tikz,array,comment}\n')
    latex_file.write('\\usepackage{tikz-qtree}\n')
    latex_file.write('\\usepackage{sansmath}\n')
    # if standalone:
    #     latex_file.write('\\usepackage{varwidth}\n')
    if not standalone:
        # tikz output in landscape
        latex_file.write('\\usepackage[a4paper,landscape]{geometry}\n')
        latex_file.write('\\usepackage[cm]{fullpage}\n')
    latex_file.write('\\usetikzlibrary{arrows,shapes,automata,patterns,trees,decorations}\n')
    latex_file.write('\n')
    latex_file.write('\\renewcommand{\\familydefault}{\\sfdefault}\n')
    latex_file.write('\\sffamily\\sansmath\n')
    latex_file.write('\n')
    latex_file.write('\\begin{document}\n\n')

    latex_file.write('% Set the overall layout of the tree\n')
    latex_file.write('\\tikzset{level 1/.style={level distance=2cm}}\n')
    latex_file.write('\\tikzset{level 2+/.style={level distance=1.6cm}}\n')
    latex_file.write('\\tikzset{edge from parent/.append style={->,line width=1.2,'
                     + 'shorten >=2.5pt,shorten <=2.5pt}}\n')
    latex_file.write('% Define styles for leafs\n')
    latex_file.write('\\tikzset{every leaf node/.style={draw,text width=1.09cm,inner sep=2pt,align=center}}\n\n')

    if not standalone:
        latex_file.write('\\begin{flushright}\n')
        latex_file.write('\\today\n')
        latex_file.write('\\end{flushright}\n\n')


def write_figure_header(latex_file, standalone=False):
    """
    Write the latex figure header to the given and already opened file
    :param latex_file: already opened file writer instance
    :param standalone: plot cut on the border of the figure instead of regular paper format
    """

    # if standalone:
    #     latex_file.write('\\begin{varwidth}{20cm}\n')
    if not standalone:
        latex_file.write('\\begin{figure}[h]\n')
        latex_file.write('\\centering\n')

    # write tikz figure to file
    latex_file.write('\\begin{tikzpicture}[sloped,font=\\sansmath\\sffamily\\normalsize]\n\n')


def write_figure_footer(latex_file, figure_caption, standalone=False):
    """
    Write figure footer to the given file with the given figure caption.
    :param latex_file: already opened file writer instance
    :param figure_caption: figure caption string
    :param standalone: plot cut on the border of the figure instead of regular paper format
    """

    latex_file.write('\n\\end{tikzpicture}\n')

    if standalone:
        # latex_file.write('\\end{varwidth}\n')
        latex_file.write('% \\caption{'+figure_caption+'}\n')
    else:
        latex_file.write('\\caption{'+figure_caption+'}\n')
        latex_file.write('\\end{figure}\n\n')


def add_ignored_mut_info(filename, phylogeny, tree):
    """
    Add additional information to the figure file about ignored mutations and variation events on the edges
    :param filename: name of the file to which the information should be appended
    :param phylogeny: data around the phylogeny
    :param tree: inferred evolutionary tree
    """

    try:
        with open(filename, 'a') as latex_file:

            latex_file.write('\\begin{comment} \n')

            gene_names = phylogeny.patient.gene_names
            driver_pathways = phylogeny.patient.driver_pathways

            # print the list of mutations which need to be excluded from the derivation
            # such that no deletions are necessary
            if len(gene_names):
                latex_file.write('\\noindent \n')
                latex_file.write('\\textbf{Ignored mutations ('+str(len(phylogeny.conflicting_mutations))+'): } \n')
                latex_file.write(', '.join(sorted('\\textcolor{orange}{'+gene_names[mut]+'} ('+driver_pathways[mut]
                                           + ')' for mut in phylogeny.conflicting_mutations if mut in driver_pathways)
                                           + sorted(gene_names[mut] for mut in phylogeny.conflicting_mutations
                                                    if mut not in driver_pathways)))

            # add a list with all edges and all their variation events
            if len(gene_names):
                latex_file.write('\n\n\\noindent \n')
                latex_file.write('\\textbf{Variation events on each edge:} \n')
                latex_file.write('\\footnotesize \n')
                latex_file.write('\\begin{itemize} \n')

                # process edges in through a first in first out queue
                fifo_queue = list()

                # add edges to queue in level order
                for child in tree.successors(TREE_ROOT):
                    fifo_queue.append((TREE_ROOT, child))

                while len(fifo_queue):

                    # take first node of the first in first out queue and process it
                    parent, child = fifo_queue.pop(0)

                    latex_file.write('\t\\item \\textbf{'+str(tree.node[parent]['name'])
                                     + ' to '+str(tree.node[child]['name'])+'}: ')
                    if len(tree[parent][child]['muts']) > 0:
                        latex_file.write(
                            '\\\\ \nAcquired mutations: ' + ', '.join(
                                sorted(('\\textcolor{orange}{'+gene_names[m]+'} (' + driver_pathways[m] + ')')
                                       for m in tree[parent][child]['muts'] if m in driver_pathways)
                                + sorted(gene_names[m] for m in tree[parent][child]['muts']
                                         if m not in driver_pathways)))
                        latex_file.write('\n\n')

                    if len(tree.successors(child)):
                        # add edges to queue in level order
                        for grandchild in tree.successors(child):
                            fifo_queue.append((child, grandchild))

                latex_file.write('\\end{itemize} \n')
                latex_file.write('\\end{comment} \n\n')

            # add a table with detailed mutation number translation
            # latex_file.write('\n\n\n% Mutation number translation table:\n')
            # for mut_idx, mut_name in enumerate(mut_names):
            #    latex_file.write('% {}: {}\n'.format(mut_idx, mut_name))

    except OSError:
        logger.exception("Unexpected error occurred during tikz file generation: ")


def add_artifact_info(file_path, phylogeny):
    """
    Add additional information to the figure file about the resolved mutation positions on the edges
    :param file_path: path to the file
    :param phylogeny: data structure around the phylogenetic tree
    """

    try:

        with open(file_path, 'a') as latex_file:

            latex_file.write('\n')
            latex_file.write('\\begin{comment} \n')

            pat = phylogeny.patient

            # print putative false positives
            for mut_idx, samples in sorted(phylogeny.false_positives.items(),
                                           key=lambda x: pat.gene_names[x[0]]):
                for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                    # putative false-positive
                    latex_file.write('Putative false-positive {} ({}) in sample {}'.format(
                        pat.gene_names[mut_idx], pat.mut_keys[mut_idx],
                        pat.sample_names[sa_idx]))
                    latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                        pat.phred_coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                        pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))

            # print putative false negatives
            for mut_idx, samples in sorted(phylogeny.false_negatives.items(),
                                           key=lambda x: pat.gene_names[x[0]]):
                for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                    latex_file.write('Putative false-negative {} ({}) in sample {}'.format(
                        pat.gene_names[mut_idx], pat.mut_keys[mut_idx],
                        pat.sample_names[sa_idx]))
                    latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                        pat.phred_coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                        pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))

            # print putative false negatives with too low coverage (unknowns)
            for mut_idx, samples in sorted(phylogeny.false_negative_unknowns.items(),
                                           key=lambda x: pat.gene_names[x[0]]):
                for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                    latex_file.write('Putative positive unknown {} ({}) in sample {}'.format(
                        pat.gene_names[mut_idx], pat.mut_keys[mut_idx],
                        pat.sample_names[sa_idx]))
                    latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                        pat.phred_coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                        pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))

            latex_file.write('\\end{comment} \n\n')

    except OSError:
        logger.exception("Unexpected error occurred during tikz file generation: ")
