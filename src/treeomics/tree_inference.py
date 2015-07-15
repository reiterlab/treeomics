"""Infer evolutionary trees by using various methods"""
__author__ = 'jreiter'
__date__ = 'July 11, 2015'

import logging
from phylogeny.simple_phylogeny import SimplePhylogeny
from phylogeny.max_lh_phylogeny import MaxLHPhylogeny
#from subclonal_phylogeny import SubclonalPhylogeny
from utils.vaf_data import get_present_mutations
import plots.tikz_tree as tikz
import utils.latex_output as latex


# get logger for application
logger = logging.getLogger('treeomics')


def infer_max_compatible_tree(filepath, patient):
    """
    Create an evolutionary tree where most conflicting mutations have been ignored
    due to the ambiguous binary classification of variants being present/absent
    :param filepath: tree is writen to the given file
    :param patient: data structure around the patient
    :return: evolutionary tree as graph
    """

    phylogeny = SimplePhylogeny(patient, patient.clones)

    # infer the tree which as to ignore the least number of mutation
    # to derive a conflict-free tree
    simple_tree = phylogeny.find_max_compatible_tree()
    caption = ('Phylogenetic tree illustrating the clonal evolution of cancer. '
               + 'The derivation of an evolutionary conflict-free tree required the exclusion '
               + 'of {} out of {} ({:.1%}) mutations.'.format(
                 len(phylogeny.conflicting_mutations), len(get_present_mutations(patient.data, include_unknowns=False)),
                 float(len(phylogeny.conflicting_mutations)) /
                 len(get_present_mutations(patient.data, include_unknowns=False))))

    # TODO: show tree with mutation patterns determined by bayesian inference

    # create tikz figure
    tikz.create_figure_file(simple_tree, tikz.TREE_ROOT, filepath,
                            patient, caption, True)
    # add information about the ignored mutations and the position of the acquired mutations
    latex.add_ignored_mut_info(filepath, phylogeny, simple_tree)

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


def create_subclonal_tree(filepath, patient, min_sc_score, min_absent_cov=0):
    """
    :param filepath: tree is writen to the given file
    :param patient: data structure around a patient
    :param min_sc_score: minimum reliability score of  an incomp. mutation pattern that a subclone is considered
    :param min_absent_cov: minimum coverage for absent variant otherwise status is unknown
    """
    phylogeny = SubclonalPhylogeny(patient, patient.clones)

    # infer the tree which as to ignore the least number of mutation
    # to derive a conflict-free tree
    pers_tree = phylogeny.find_subclonal_tree(min_sc_score, min_absent_cov=min_absent_cov)
    caption = ('Phylogenetic tree illustrating the subclonal evolution of cancer. '
               + 'The derivation of an evolutionary conflict-free tree required the exclusion '
               + 'of {} out of {} ({:.1%}) mutations.'
               .format(len(phylogeny.conflicting_mutations),
                       len(get_present_mutations(patient.data, include_unknowns=False)),
                       float(len(phylogeny.conflicting_mutations)) / len(
                           get_present_mutations(patient.data, include_unknowns=False))))

    # create tikz figure
    tikz.create_figure_file(pers_tree, tikz.TREE_ROOT, filepath,
                            patient, caption, True)
    # add information about the ignored mutations and the position of the acquired mutations
    latex.add_ignored_mut_info(filepath, phylogeny, pers_tree)

    return phylogeny


def create_max_lh_tree(file_path, patient):
    """
    Create an evolutionary tree based on the maximum likelihood mutation patterns of each variant
    :param file_path: tree is writen to the given file
    :param patient: data structure around the patient
    :return: evolutionary tree as graph
    """

    mlh_pg = MaxLHPhylogeny(patient, patient.clones)

    mlh_tree = mlh_pg.infer_max_lh_tree()

    if mlh_tree is not None:

        # ignore mutations which are not in any sample which passed the filtering
        present_mutations = get_present_mutations(patient.data, include_unknowns=False)

        no_fps = sum(len(fps) for mut_idx, fps in mlh_pg.false_positives.items())
        no_fns = sum(len(fns) for mut_idx, fns in mlh_pg.false_negatives.items())
        classification_info = ('Putative false-positives {}, put. false-negatives {}, put. false neg.-unknowns {}.'
                               .format(no_fps, no_fns,
                                       sum(len(fns) for mut_idx, fns in mlh_pg.false_negative_unknowns.items())))
        logger.info(classification_info)

        caption = ('Phylogenetic tree illustrating the clonal evolution of cancer. '
                   + 'The derivation of an evolutionarily-compatible maximum likelihood tree identified '
                   + '{} suspicious positions (out of {}; {:.1%}). '.format(
                     no_fps+no_fns, len(patient.sample_names) * len(present_mutations),
                     float(no_fps+no_fns) / (len(patient.sample_names) * len(present_mutations)))
                   + classification_info)

        tikz.create_figure_file(mlh_tree, tikz.TREE_ROOT, file_path,
                                patient, caption, True)

        # add information about the resolved mutation positions
        # which are likely sequencing errors
        latex.add_artifact_info(file_path, mlh_pg)
    else:
        logger.warn('Conflicts could not be resolved. No evolutionary tree has been created.')

    return mlh_pg
