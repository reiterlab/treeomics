"""Infer evolutionary trees by using various methods"""
import logging
import os
import numpy as np
from subprocess import call
import sys
from phylogeny.simple_phylogeny import SimplePhylogeny
from phylogeny.max_lh_phylogeny import MaxLHPhylogeny
from utils.mutation_matrix import write_mutation_matrix
from utils.data_tables import write_mutation_patterns
import plots.tikz_tree as tikz
import utils.latex_output as latex

__author__ = 'jreiter'
__date__ = 'July 11, 2015'


# get logger for application
logger = logging.getLogger('treeomics')


def create_max_lh_tree(patient, tree_filepath=None, mm_filepath=None, mp_filepath=None, subclone_detection=False,
                       loh_frequency=0.0, driver_vars=set(), max_no_mps=None, time_limit=None, plots=True,
                       pool_size=0, no_bootstrap_samples=0, variant_filepath=None, n_max_threads=0):
    """
    Create an evolutionary tree based on the maximum likelihood mutation patterns of each variant
    :param patient: data structure around the patient
    :param tree_filepath: tree is written to the given file
    :param mm_filepath: path to mutation matrix output file
    :param mp_filepath: path to mutation pattern output file
    :param subclone_detection: is subclone detection enabled?
    :param loh_frequency: probability that a SNV along a lineage is lost due loss of heterozygosity
    :param driver_vars: defaultdict with mutation IDs and instance of driver class to be highlighted on edges
    :param max_no_mps: only the given maximal number of most likely (by joint likelihood) mutation patterns
            is explored per variant; limits the solution space
    :param time_limit: time limit for MILP solver in seconds
    :param plots: generate pdf from tex file
    :param pool_size: number of best solutions explored by ILP solver to estimate confidence
    :param n_max_threads: Sets the default maximal number of parallel threads that will be invoked by CPLEX
                          (0: default, let CPLEX decide; 1: single threaded; N: uses up to N threads)
    :param no_bootstrap_samples: number of samples with replacement for the bootstrapping
    :param variant_filepath: path to output file with information about variants and where they were acquired
    :return: evolutionary tree as graph
    """

    mlh_pg = MaxLHPhylogeny(patient, patient.mps, loh_frequency=loh_frequency)

    mlh_tree = mlh_pg.infer_max_lh_tree(subclone_detection=subclone_detection, max_no_mps=max_no_mps,
                                        pool_size=pool_size, time_limit=time_limit, n_max_threads=n_max_threads,
                                        no_bootstrap_samples=no_bootstrap_samples)

    if mlh_tree is not None:

        # if more solutions were inferred, plot their trees as well
        for sol in mlh_pg.solutions[1:]:
            sol.assign_variants(mlh_pg, max_no_mps)
            sol.find_artifacts(mlh_pg)
            sol.infer_mutation_patterns(mlh_pg)
            # construct a phylogenetic tree from the maximum likelihood mutation patterns
            tree = mlh_pg.infer_evolutionary_tree(sol.shared_mlh_mps, sol.mlh_founders, sol.mlh_unique_mutations,
                                                  confidence=None if mlh_pg.weighted_node_lh is None else
                                                  mlh_pg.weighted_node_lh)

            if tree_filepath is not None:
                _create_tree_plots(sol, mlh_pg, tree, plots, tree_filepath + '_{}'.format(sol.rank), driver_vars)

        # create tree for most likely solution
        if tree_filepath is not None:
            _create_tree_plots(mlh_pg.solutions[0], mlh_pg, mlh_tree, plots,
                               tree_filepath + '_1', driver_vars, variant_filepath=variant_filepath)

        # create mutation matrix for benchmarking
        if mm_filepath is not None:
            write_mutation_matrix(mlh_pg, mm_filepath)

        if mp_filepath is not None:
            write_mutation_patterns(mlh_pg, mp_filepath)

    else:
        logger.warning('Conflicts could not be resolved. No evolutionary tree has been created.')

    return mlh_pg


def _create_tree_plots(solution, mlh_pg, mlh_tree, plots, tree_filepath, driver_vars, variant_filepath=None):
    """
    Produce evolutionary tree plots for the given solution
    :param solution:
    :param mlh_pg:
    :param mlh_tree:
    :param plots:
    :param tree_filepath:
    :param driver_vars: defaultdict with mutation IDs and instance of driver class
    :param variant_filepath:
    """

    # ignore mutations which are not in any sample which passed the filtering
    present_mutations = mlh_pg.patient.present_mutations

    no_fps = sum(len(fps) for mut_idx, fps in solution.false_positives.items())
    no_fns = sum(len(fns) for mut_idx, fns in solution.false_negatives.items())
    classification_info = ('Putative false-positives {}, put. false-negatives {}, put. false neg.-unknowns {}. '
                           .format(no_fps, no_fns, sum(
                            len(fns) for mut_idx, fns in solution.false_negative_unknowns.items())))

    if solution.conflicting_mutations is not None:
        compatibility_info = ('{} ({:.1%})'.format(
            len(solution.conflicting_mutations), float(len(solution.conflicting_mutations)) /
            (len(solution.max_lh_mutations) + len(solution.conflicting_mutations))) +
            ' variants were evolutionarily incompatible due to the limited search space.')
    else:
        compatibility_info = ''

    logger.info(classification_info)
    caption = ('Phylogenetic tree illustrating the evolutionary history of the cancer. ' +
               'The derivation of an evolutionarily-compatible maximum likelihood tree identified ' +
               '{} putative false-positives or false-negatives (out of {}; {:.1%}). '.format(
                   no_fps + no_fns, len(mlh_pg.patient.sample_names) * len(present_mutations),
                   float(no_fps + no_fns) / (len(mlh_pg.patient.sample_names) * len(present_mutations))) +
               classification_info + compatibility_info)
    # calculate median of number of persistent and present mutations inferred
    # in the evolutionary trajectory
    if len(solution.variants.values()) > 0:
        median_no_muts = np.median([len(muts) for muts in solution.variants.values()])
    else:
        median_no_muts = float('nan')

    if plots and tree_filepath is not None:
        try:
            # create ETE tree
            from plots.ete_tree import create_tree
            mlh_pg.tree_plot = create_tree(mlh_tree, tikz.TREE_ROOT, tree_filepath, mlh_pg.patient, mlh_pg,
                                           drivers=driver_vars)
            logger.info('Generated ete tree as PNG: {}'.format(mlh_pg.tree_plot))

        except ImportError as ie:
            logger.warn('ImportError! ete3 is not installed! {}'.format(ie))

        # create Latex/TIkZ tree
        tikz_tree = tikz.create_figure_file(
            mlh_tree, tikz.TREE_ROOT, tree_filepath + '_full.tex', mlh_pg.patient, mlh_pg, caption,
            driver_vars=driver_vars, germline_distance=10.0 * max(1.0, len(solution.mlh_founders) / median_no_muts),
            standalone=True, variant_filepath=variant_filepath)
        # add information about the ignored mutations and the position of the acquired mutations
        latex.add_branch_mut_info(tree_filepath, mlh_pg, mlh_tree)

        # add information about the resolved mutation positions
        # which are likely sequencing errors
        latex.add_artifact_info(tree_filepath, mlh_pg)

        tikz_path, tikz_file = os.path.split(tikz_tree)
        logger.debug('Tikzpath: {} {}'.format(tikz_path, tikz_file))
        # increase buffer size on mac: 'buf_size=1000000 pdflatex {}'.format(tikz_file)
        # increase buffer size on windows: 'pdflatex --buf-size=1000000 {}'.format(tikz_file)
        # check which operating system is used
        if sys.platform == 'darwin':
            pdflatex_cmd = 'buf_size=1000000 pdflatex {}'.format(tikz_file)
        elif sys.platform == 'win32':
            pdflatex_cmd = 'pdflatex --buf-size=1000000 {}'.format(tikz_file)
        else:
            pdflatex_cmd = 'pdflatex {}'.format(tikz_file)

        fnull = open(os.devnull, 'w')
        return_code = call(pdflatex_cmd, shell=True, cwd=tikz_path, stdout=fnull)

        if return_code == 0:
            pdf_tree = tikz_tree.replace('.tex', '.pdf')
            logger.info('Successfully called pdflatex to create pdf of the evolutionary tree at {}'.format(
                pdf_tree))
        else:
            logger.error('PDF of the evolutionary tree was not created. Is Latex/tikz installed?')


def infer_max_compatible_tree(filepath, patient, drivers=set(), time_limit=None):
    """
    Create an evolutionary tree where most conflicting mutations have been ignored
    due to the ambiguous binary classification of variants being present/absent
    :param filepath: tree is writen to the given file
    :param patient: data structure around the patient
    :param drivers: known driver genes in considered cancer type
    :param time_limit: time limit for MILP solver in seconds
    :return: evolutionary tree as graph
    """

    phylogeny = SimplePhylogeny(patient, patient.mps)

    # infer the tree which as to ignore the least number of mutation
    # to derive a conflict-free tree
    simple_tree = phylogeny.find_max_compatible_tree(time_limit=time_limit)

    # number of mutations inferred to be present in at least one sample
    no_present_muts = len(phylogeny.compatible_mutations)+len(phylogeny.conflicting_mutations)
    caption = ('Phylogenetic tree illustrating the clonal evolution of cancer. ' +
               'The derivation of an evolutionary conflict-free tree required the exclusion ' +
               'of {} out of {} ({:.1%}) mutations.'.format(
                 len(phylogeny.conflicting_mutations), no_present_muts,
                 float(len(phylogeny.conflicting_mutations)) / no_present_muts))

    # create tikz figure
    tikz.create_figure_file(simple_tree, tikz.TREE_ROOT, filepath,
                            patient, phylogeny, caption, driver_vars=drivers, standalone=True)
    # add information about the ignored mutations and the position of the acquired mutations
    latex.add_branch_mut_info(filepath, phylogeny, simple_tree)

    # ensure that the derived tree has the correct number of mutations on all leaves
    if logging.DEBUG == logger.getEffectiveLevel():
        for sa_idx, sa_name in enumerate(patient.sample_names):
            logger.debug("Compatible mutations present in sample {}: {}, {}".format(sa_name,
                         sum(1 for mut in patient.samples[sa_idx] if mut in phylogeny.compatible_mutations),
                         ', '.join(str(mut) for mut in patient.samples[sa_idx]
                                   if mut in phylogeny.compatible_mutations)))

        #     assert (len(pers_tree.node[frozenset([sa_idx])]['muts'])
        #             >= sum(1 for mut in patient.samples[sa_idx] if mut in phylogeny.compatible_mutations)), \
        #         'Mutations are incorrect for sample {}: {} != {}'.format(sa_idx,
        #         len(pers_tree.node[frozenset([sa_idx])]['muts']),
        #         sum(1 for mut in patient.samples[sa_idx] if mut in phylogeny.compatible_mutations))
        #
        # assert (sum(len(pers_tree[v1][v2]['muts']) for (v1, v2) in pers_tree.edges_iter())
        #         == len(phylogeny.compatible_mutations)), \
        #     'Total number of acquired mutations equals the number of compatible mutations: {} == {}'.format(
        #         sum(len(pers_tree[v1][v2]['muts']) for (v1, v2) in pers_tree.edges_iter()),
        #         len(phylogeny.compatible_mutations))

    return phylogeny
