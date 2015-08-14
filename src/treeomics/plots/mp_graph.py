#!/usr/bin/python
"""Automatically generate circos configuration files for a circos mutation pattern overview graph"""
__author__ = 'Johannes REITER'
__date__ = 'March 6, 2015'


import logging
import plots.circos as circos
import os
from subprocess import call

# get logger for application
logger = logging.getLogger('treeomics')


def create_mp_graph(fn_pattern, phylogeny, mps, mp_weights, output_directory='',
                    min_node_weight=None, max_no_mps=50, pars_informative=False):
    """
    Create mutation pattern overview plot showing only the different patterns and not the individual variants
    :param fn_pattern: common name template for output files
    :param phylogeny: data structure around the phylogenetic tree
    :param mps: mutation patterns
    :param mp_weights: mutation pattern reliability scores
    :param output_directory: path to the output directory for the generated files
    :param min_node_weight: minimal reliability score of a mutation pattern to be displayed
    :param max_no_mps: apply min_node_weight if there are more than this number of MPs in the data
    :param pars_informative: show only parsimony informative mutation patterns
    """

    if min_node_weight is None:
        min_node_weight = max(0.5, len(phylogeny.patient.present_mutations)/200.0)
        logger.info('Minimal mutation pattern weight for MP overview graph set to {}'.format(min_node_weight))
    else:
        logger.info('Minimal mutation pattern weight for MP overview graph: {}'.format(min_node_weight))

    # create directory for data files
    data_dir = os.path.join(output_directory, 'data', '')
    try:
        if not os.path.exists(os.path.dirname(data_dir)):
            os.makedirs(os.path.dirname(data_dir))
        logger.debug('Data directory: {}'.format(os.path.abspath(data_dir)))
    except OSError:
        logger.exception("Could not create mutation pattern overview graph directory {} ".format(data_dir))

    circos.create_mp_graph_files(
        os.path.join(data_dir, 'mp_nodes_'+fn_pattern+'.txt'),
        os.path.join(data_dir, 'mp_mutnode_data_'+fn_pattern+'.txt'),
        os.path.join(data_dir, 'mp_links_'+fn_pattern+'.txt'), phylogeny, mps, mp_weights,
        min_node_weight=min_node_weight, max_no_mps=max_no_mps, pars_infor=pars_informative)

    _create_mp_gr_confs(phylogeny, output_directory, fn_pattern)

    mp_filename = 'fig_mp_{}.png'.format(phylogeny.patient.name)
    mp_cmd = 'circos -outputfile {}'.format(mp_filename)
    return_code = call(mp_cmd, shell=True, cwd=output_directory)

    if return_code == 0:
        logger.info('Successfully called circos to create mutation pattern overview graph at {}'.format(mp_filename))
        return mp_filename
    else:
        logger.error('Mutation pattern overview graph was not created. Is circos installed?')
        return None


def _create_mp_gr_confs(phylogeny, output_directory, fn_pattern):
    """
    Create configuration files for the mutation pattern overview graph
    :param phylogeny: data structure around the phylogenetic tree
    :param output_directory: path to the output directory for the generated files
    """

    # create directory for configuration files
    conf_dir = os.path.join(output_directory, 'etc', '')
    try:
        if not os.path.exists(os.path.dirname(conf_dir)):
            os.makedirs(os.path.dirname(conf_dir))
        logger.debug('Data directory: {}'.format(os.path.abspath(conf_dir)))
    except OSError:
        logger.exception("Could not create mutation pattern overview graph directory {} ".format(conf_dir))

    logger.debug('Creating circos mutation pattern overview graph configuration files at {}'.format(conf_dir))

    _create_mp_gr_conf_cir(phylogeny, os.path.join(conf_dir, 'circos.conf'), fn_pattern)
    _create_mp_gr_conf_ideo(os.path.join(conf_dir, 'ideogram.conf'))
    _create_mp_gr_conf_ticks(os.path.join(conf_dir, 'ticks.conf'))


def _create_mp_gr_conf_cir(phylogeny, filepath, fn_pattern):
    """
    Create circos.conf
    :param phylogeny: data structure around the phylogenetic tree
    :param filepath: path to the output file
    :param fn_pattern: common name template for output files
    """

    # open file
    with open(filepath, 'w') as f:

        f.write('karyotype = data/mp_nodes_{}.txt'.format(fn_pattern)+'\n')
        f.write('chromosomes_units = 1'+'\n')
        f.write('\n')

        f.write('<links>'+'\n')
        f.write('<link>'+'\n')
        f.write('file = data/mp_links_{}.txt'.format(fn_pattern)+'\n')
        f.write('radius        = 0.45r'+'\n')
        f.write('bezier_radius = 0r'+'\n')
        f.write('color         = dred'+'\n')
        f.write('thickness     = 3'+'\n')
        f.write('</link>'+'\n')
        f.write('</links>'+'\n')
        f.write('\n')

        f.write('<plots>'+'\n')
        f.write('# configuration for scatter plot'+'\n')
        f.write('<plot>'+'\n')
        f.write('type = scatter'+'\n')
        f.write('file = data/mp_mutnode_data_{}.txt'.format(fn_pattern)+'\n')
        f.write('r1 = 0.95r'+'\n')
        f.write('r0 = 0.5r'+'\n')
        f.write('max = {}'.format(len(phylogeny.patient.sample_names))+'\n')
        f.write('min = 1.0'+'\n')
        f.write('\n')

        f.write('glyph = circle'+'\n')
        f.write('glyph_size = 15'+'\n')
        f.write('color = black'+'\n')
        f.write('stroke_color = black'+'\n')
        f.write('stroke_thickness = 3'+'\n')
        f.write('\n')

        f.write('<backgrounds>'+'\n')
        f.write('<background>'+'\n')
        f.write('color = vvlgrey'+'\n')
        f.write('</background>'+'\n')
        f.write('</backgrounds>'+'\n')
        f.write('\n')

        f.write('<axes>'+'\n')
        f.write('<axis>'+'\n')
        f.write('color     = grey'+'\n')
        f.write('thickness = 1'+'\n')
        f.write('spacing   = 1'+'\n')
        f.write('</axis>'+'\n')
        f.write('</axes>'+'\n')
        f.write('\n')

        f.write('</plot>'+'\n')
        f.write('</plots>'+'\n')
        f.write('\n\n')

        f.write('<<include ideogram.conf>>'+'\n')
        f.write('<ideogram>'+'\n')
        f.write('<spacing>'+'\n')
        f.write('# spacing between ideograms'+'\n')
        f.write('default = {}r'.format(max(0.01, 1.0/len(phylogeny.mps)))+'\n')
        f.write('</spacing>'+'\n')
        f.write('\n')

        f.write('# ideogram position, thickness and fill'+'\n')
        f.write('radius = 0.89r'+'\n')
        f.write('thickness = 30p'+'\n')
        f.write('fill = yes'+'\n')
        f.write('show_label = yes'+'\n')
        f.write('label_with_tag = yes'+'\n')
        f.write('label_font = condensed'+'\n')
        f.write('label_radius = dims(ideogram,radius) + 0.08r'+'\n')
        f.write('label_center = yes'+'\n')
        f.write('label_size = 28p'+'\n')
        f.write('label_color = dgrey'+'\n')
        f.write('label_parallel = yes'+'\n')
        f.write('label_format = eval(sprintf("%s",var(label)))'+'\n')
        f.write('</ideogram>'+'\n')
        f.write('\n')

        f.write('<image>'+'\n')
        f.write('<<include etc/image.conf>>'+'\n')
        f.write('file* = mp_graph.png'+'\n')
        f.write('radius* = {}p'.format(500 if len(phylogeny.patient.sample_names) < 12 else
                                       len(phylogeny.patient.sample_names)*50)+'\n')
        f.write('</image>'+'\n')
        f.write('<<include etc/colors_fonts_patterns.conf>>'+'\n')
        f.write('<<include etc/housekeeping.conf>>'+'\n')
        f.write('\n')


def _create_mp_gr_conf_ideo(filepath):
    """
    Create ideogram.conf
    :param filepath: path to the output file
    """

    # open file
    with open(filepath, 'w') as f:

        f.write('<ideogram>'+'\n')
        f.write('<spacing>'+'\n')
        f.write('default = 0.004r'+'\n')
        f.write('</spacing>'+'\n')
        f.write('\n')

        f.write('# Ideogram position, fill and outline'+'\n')
        f.write('radius = 0.90r'+'\n')
        f.write('thickness = 20p'+'\n')
        f.write('fill = yes'+'\n')
        f.write('stroke_color = dgrey'+'\n')
        f.write('stroke_thickness = 2p'+'\n')
        f.write('\n')

        f.write('# Minimum definition for ideogram labels.'+'\n')
        f.write('show_label = yes'+'\n')
        f.write('# see etc/fonts.conf for list of font names'+'\n')
        f.write('label_font = default '+'\n')
        f.write('label_radius = dims(image,radius) - 60p'+'\n')
        f.write('label_size = 30'+'\n')
        f.write('label_parallel = yes'+'\n')
        f.write('</ideogram>'+'\n')
        f.write('\n')


def _create_mp_gr_conf_ticks(filepath):
    """
    Create ticks.conf
    :param filepath: path to the output file
    """

    # open file
    with open(filepath, 'w') as f:

        f.write('show_ticks = no'+'\n')
        f.write('show_tick_labels = no'+'\n')
