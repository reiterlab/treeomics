#!/usr/bin/python
"""Helper functions to generate latex source code for a tikz tree plot"""
__author__ = 'jreiter'

import logging
import itertools
import numpy as np
from phylogeny.phylogeny_utils import TREE_ROOT

# get logger for application
logger = logging.getLogger('treeomics')


def write_tikz_header(latex_file, germline_distance=2.0, standalone=False):
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
    latex_file.write('\\tikzset{level 1/.style={level distance='+'{:.1f}'.format(germline_distance)+'cm}}\n')
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


def add_branch_mut_info(filename, phylogeny, tree):
    """
    Add additional information to the figure file about ignored mutations and variation events on the edges
    :param filename: name of the file to which the information should be appended
    :param phylogeny: data around the phylogeny
    :param tree: inferred evolutionary tree
    """

    try:
        with open(filename, 'a') as latex_file:

            latex_file.write('\\begin{comment} \n')

            mut_keys = phylogeny.patient.mut_keys
            gene_names = phylogeny.patient.gene_names
            driver_pathways = phylogeny.patient.driver_pathways

            # print the list of mutations which need to be excluded from the derivation
            # such that all variants allow a perfect and persistent phylogeny
            # latex_file.write('\\noindent \n')
            # latex_file.write('\\textbf{Ignored mutations ('+str(len(phylogeny.conflicting_mutations))+'): } \n')
            # if gene_names is not None:
            #     latex_file.write(', '.join(sorted('\\textcolor{orange}{'+gene_names[mut]+'} ('+driver_pathways[mut]
            #                                + ')' for mut in phylogeny.conflicting_mutations if mut in driver_pathways)
            #                                + sorted(gene_names[mut] for mut in phylogeny.conflicting_mutations
            #                                         if mut not in driver_pathways)))
            #
            # latex_file.write(', '.join(sorted('\\textcolor{orange}{'+mut_keys[mut]+'} ('+driver_pathways[mut]
            #                            + ')' for mut in phylogeny.conflicting_mutations if mut in driver_pathways)
            #                            + sorted(mut_keys[mut] for mut in phylogeny.conflicting_mutations
            #                                     if mut not in driver_pathways)))

            # add a list with all edges and all their variation events
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
                            sorted(('\\textcolor{orange}{'+mut_keys[m]
                                    + ('({})'.format(gene_names[m]) if gene_names is not None else '')
                                    + '} (' + driver_pathways[m] + ')')
                                   for m in tree[parent][child]['muts'] if m in driver_pathways)
                            + sorted(mut_keys[m]+(' ({})'.format(gene_names[m]) if gene_names is not None else '')
                                     for m in tree[parent][child]['muts'] if m not in driver_pathways)))
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

            if pat.gene_names is not None:
                if phylogeny.conflicting_mutations is not None and len(phylogeny.conflicting_mutations) > 0:
                    # print evolutionarily incompatible (unclassified artifacts) due to limited solution space
                    latex_file.write('Evolutionarily incompatible variants due to limited search space: {} \n\n'.format(
                        ', '.join('{} ({})'.format(pat.gene_names[mut_idx], pat.mut_keys[mut_idx])
                                  for mut_idx in sorted(phylogeny.conflicting_mutations,
                                                        key=lambda x: pat.gene_names[x].lower()))))

                    for mut_idx in sorted(phylogeny.conflicting_mutations, key=lambda x: pat.gene_names[x].lower()):
                        latex_file.write('Evolutionarily incompatible {} ({}): '.format(
                                         pat.gene_names[mut_idx], pat.mut_keys[mut_idx]))
                        for sa_idx in sorted(range(len(pat.sample_names)), key=lambda x: pat.sample_names[x]):
                            latex_file.write('{} (reads {}/{}); '.format(
                                pat.sample_names[sa_idx],
                                pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                                pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))
                        latex_file.write('\n')

                # print putative false positives
                for mut_idx, samples in sorted(phylogeny.false_positives.items(),
                                               key=lambda x: pat.gene_names[x[0]].lower()):
                    for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                        # putative false-positive
                        latex_file.write('Putative false-positive {} ({}) in sample {}'.format(
                            pat.gene_names[mut_idx], pat.mut_keys[mut_idx],
                            pat.sample_names[sa_idx]))
                        latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))

                fp_cov = [pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]
                          for mut_idx, samples in phylogeny.false_positives.items() for sa_idx in samples]
                fp_var = [pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]
                          for mut_idx, samples in phylogeny.false_positives.items() for sa_idx in samples]
                fp_vaf = [(pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]] /
                          pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]])
                          if pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]] > 0 else 0.0
                          for mut_idx, samples in phylogeny.false_positives.items() for sa_idx in samples]
                latex_file.write('False-positives: coverage: mean {}, median {}; vars: mean {}, median {}\n'.format(
                    np.mean(fp_cov), np.median(fp_cov), np.mean(fp_var), np.median(fp_var)))
                latex_file.write('False-positives: VAF: mean {:.2%}, median {:.2%}\n\n'.format(
                    np.mean(fp_vaf), np.median(fp_vaf)))

                # print putative false negatives
                for mut_idx, samples in sorted(phylogeny.false_negatives.items(),
                                               key=lambda x: pat.gene_names[x[0]].lower()):
                    for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                        latex_file.write('Putative powered false-negative {} ({}) in sample {}'.format(
                            pat.gene_names[mut_idx], pat.mut_keys[mut_idx],
                            pat.sample_names[sa_idx]))
                        latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))
                pfn_cov = [pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]
                           for mut_idx, samples in phylogeny.false_negatives.items() for sa_idx in samples]
                latex_file.write('Powered false-negatives: Coverage: mean {}, median {}\n\n'.format(
                    np.mean(pfn_cov), np.median(pfn_cov)))

                # print putative false negatives with too low coverage (unknowns)
                for mut_idx, samples in sorted(phylogeny.false_negative_unknowns.items(),
                                               key=lambda x: pat.gene_names[x[0]].lower()):
                    for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                        latex_file.write('Putative under-powered false-negative {} ({}) in sample {}'.format(
                            pat.gene_names[mut_idx], pat.mut_keys[mut_idx],
                            pat.sample_names[sa_idx]))
                        latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))
                upfn_cov = [pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]
                            for mut_idx, samples in phylogeny.false_negative_unknowns.items() for sa_idx in samples]
                latex_file.write('Under-powered false-negatives: Coverage: mean {}, median {}\n\n'.format(
                    np.mean(upfn_cov), np.median(upfn_cov)))

            else:
                if phylogeny.conflicting_mutations is not None and len(phylogeny.conflicting_mutations) > 0:
                    # print evolutionarily incompatible (unclassified artifacts) due to limited solution space
                    for mut_idx in sorted(phylogeny.conflicting_mutations, key=lambda x: pat.gene_names[x].lower()):
                        latex_file.write('Evolutionarily incompatible {}: '.format(pat.mut_keys[mut_idx]))
                        for sa_idx in sorted(range(len(pat.sample_names)), key=lambda x: pat.sample_names[x]):
                            latex_file.write('{} (reads {}/{}); '.format(
                                pat.sample_names[sa_idx],
                                pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                                pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))
                        latex_file.write('\n')
                # print putative false positives
                for mut_idx, samples in sorted(phylogeny.false_positives.items(),
                                               key=lambda x: pat.mut_keys[x[0]]):
                    for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                        # putative false-positive
                        latex_file.write('Putative false-positive {} in sample {}'.format(
                            pat.mut_keys[mut_idx], pat.sample_names[sa_idx]))
                        latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))

                # print putative false negatives
                for mut_idx, samples in sorted(phylogeny.false_negatives.items(),
                                               key=lambda x: pat.mut_keys[x[0]]):
                    for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                        latex_file.write('Putative false-negative {} in sample {}'.format(
                            pat.mut_keys[mut_idx], pat.sample_names[sa_idx]))
                        latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))

                # print putative false negatives with too low coverage (unknowns)
                for mut_idx, samples in sorted(phylogeny.false_negative_unknowns.items(),
                                               key=lambda x: pat.mut_keys[x[0]]):
                    for sa_idx in sorted(samples, key=lambda x: pat.sample_names[x]):
                        latex_file.write('Putative positive unknown {} in sample {}'.format(
                            pat.mut_keys[mut_idx], pat.sample_names[sa_idx]))
                        latex_file.write(' (Cov: {}, var-reads: {}).\n'.format(
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]))

            latex_file.write('\\end{comment} \n\n')

            latex_file.write('\\begin{comment} \n')
            # produce latex table for paper
            latex_file.write(' & '.join('\\textbf{'+col_name.replace('_','~')+'}' for col_name in itertools.chain(
                ['Sample'], pat.sample_names)) + '\\\\\n')
            latex_file.write('\\hline \n')

            ufns = [0 for _ in range(pat.n)]
            for mut_idx, samples in phylogeny.false_negative_unknowns.items():
                for sa_idx in samples:
                    ufns[sa_idx] += 1
            latex_file.write('\\textbf{Under-powered false-negatives} & '
                             + ' & '.join('$'+str(f)+'$' for f in ufns) + ' \\\\\n')

            fns = [0 for _ in range(pat.n)]
            for mut_idx, samples in phylogeny.false_negatives.items():
                for sa_idx in samples:
                    fns[sa_idx] += 1
            latex_file.write('\\textbf{Powered false-negatives} & '
                             + ' & '.join('$'+str(f)+'$' for f in fns) + ' \\\\\n')

            fps = [0 for _ in range(pat.n)]
            for mut_idx, samples in phylogeny.false_positives.items():
                for sa_idx in samples:
                    fps[sa_idx] += 1
            latex_file.write('\\textbf{False-positives} & '
                             + ' & '.join('$'+str(f)+'$' for f in fps) + ' \\\\\n')

            latex_file.write('\\end{comment} \n\n')

    except OSError:
        logger.exception("Unexpected error occurred during tikz file generation: ")
