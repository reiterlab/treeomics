#!/usr/bin/python
"""Main file to run treeomics"""
__author__ = 'jreiter'
__date__ = 'March 31, 2014'

import logging
import sys
import os
import argparse
from patient import Patient
import settings
import utils.int_settings as int_sets
import tree_inference as ti
from utils.html_report import HTMLReport
import plots.plots_utils as plts
import plots.mp_graph as mp_graph
import plots.circos as circos
import utils.analysis as analysis


# create logger for application
logger = logging.getLogger('treeomics')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('treeomics.log')
fh.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s %(levelname)s %(pathname)s %(lineno)s: %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


def init_output(patient_name=None):
    """
    Set up output directory.
    :param patient_name: if patient name is given generate separate directory for each subject
    :return: path to output directory
    """

    # derive correct output directory
    if os.getcwd().endswith('treeomics'):
        # application has been started from this directory
        # set output directory to the parent
        output_directory = os.path.join('..', settings.OUTPUT_FOLDER)
    else:
        output_directory = os.path.join(settings.OUTPUT_FOLDER)

    if patient_name is not None:
        output_directory = os.path.join(output_directory, patient_name, '')
    else:
        output_directory = os.path.join(output_directory, '')

    # create directory for data files
    try:
        if not os.path.exists(os.path.dirname(output_directory)):
            os.makedirs(os.path.dirname(output_directory))
        logger.info('Output directory: {}'.format(os.path.abspath(output_directory)))
    except OSError:
        logger.exception("Could not create output folder {} ".format(settings.OUTPUT_FOLDER))

    return output_directory


def get_patients_name(name):
    """
    Guess/extract patient name from input file
    :param name: input filename
    :return: patient name
    """
    # extract patient's name from filename or path
    basename = os.path.basename(name)
    patient_name, _ = os.path.splitext(basename)

    return patient_name


def get_output_fn_template(name, total_no_samples, fpr, fdr, min_absent_coverage, min_sa_coverage,
                           min_sa_maf, bi_e=None, bi_c0=None):
    """
    Generate common name template for output files
    :param bi_e:  sequencing error for bayesian inference
    :param bi_c0: prior mixture parameter of delta function and uniform distribution for bayesian inference
    :return: output filename pattern
    """

    pattern = ('{}_{}'.format(name, total_no_samples)
               + ('_st={}'.format(min_sa_coverage) if min_sa_coverage > 0 else '')
               + ('_mf={}'.format(min_sa_maf) if min_sa_maf > 0 else '')
               + '_fpr={:.1e}_fdr={:.2f}'.format(fpr, fdr)
               + ('_at={}'.format(min_absent_coverage))
               + ('_e={}'.format(bi_e) if bi_e is not None else '')
               + ('_c0={}'.format(bi_c0) if bi_c0 is not None else ''))

    # replace points with commas because latex cannot handle points in file names (interprets it as file type)
    pattern = pattern.replace('.', '_')

    return pattern


def usage():
    """
    Give the user feedback on how to call the tool
    Terminates the tool afterwards
    """
    logger.warn("Usage: python treeomics.zip -r <mut-reads table> -s <coverage table> | -v vcf_file | "
                "-d vcf_file_directory  [-n normal_sample_name] \n")
    logger.warn("Example: python treeomics.zip -r Pam03_mutant_reads.txt -s Pam03_phredcoverage.txt ")
    sys.exit(2)


# -------------------------------- main function ---------------------------------
def main():

    parser = argparse.ArgumentParser(description='Infers the evolution of cancer.')

    parser.add_argument("-m", "--mode", help="running mode: 1...fast (one mutation pattern per variant), "
                                             "2...complete (explore full solution space).",
                        type=int, default=2)

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--vcf_file", help="path to the VCF file", type=str)
    group.add_argument("-d", "--directory", help="directory with multiple VCF files", type=str)

    parser.add_argument("-n", "--normal", help="name of normal sample (excluded from analysis)", type=str)

    parser.add_argument("-r", "--mut_reads", help="table with the number of reads with a mutation", type=str)
    parser.add_argument("-s", "--mut_cov", help="table with read phred coverage at the mutated positions", type=str)
    parser.add_argument("-t", "--mut_dis_cov",
                        help="table with read distinct phred coverage at the mutated positions", type=str)

    # parameter for unknown variant thresholds
    parser.add_argument("-c", "--min_median_coverage",
                        help="minimum median phred coverage of a sample",
                        type=int, default=settings.SAMPLE_COVERAGE_THRESHOLD)
    parser.add_argument("-f", "--min_median_maf",
                        help="minimum median mutant allele frequency of a sample",
                        type=float, default=settings.MAF_THRESHOLD)
    parser.add_argument("-p", "--false_positive_rate",
                        help="false positive rate for the statistical test",
                        type=float, default=settings.FPR)
    parser.add_argument("-i", "--false_discovery_rate",
                        help="false discovery rate for the statistical test",
                        type=float, default=settings.FDR)
    parser.add_argument("-y", "--min_absent_coverage",
                        help="minimum coverage for a true negative (negative threshold)",
                        type=int, default=settings.MIN_ABSENT_COVERAGE)

    parser.add_argument('-o', '--down', help='Replications for robustness analysis through down-sampling',
                        type=int, default=0)

    parser.add_argument("-u", "--min_sc_score",
                        help="minimum reliability score of a mutation pattern with putative subclones",
                        type=float, default=settings.MIN_MP_LH)

    parser.add_argument("-e", "--error_rate", help="data error rate for bayesian inference",
                        type=float, default=settings.BI_E)

    args = parser.parse_args()
    plots_report = True    # for debugging set to False
    plots_paper = False

    if args.normal:
        normal_sample_name = args.normal
        logger.info('Exclude normal sample with name: {}'.format(normal_sample_name))
    else:
        normal_sample_name = ''

    if args.min_median_coverage > 0:
        logger.info('Minimum sample median coverage (otherwise discarded): {}'.format(
            args.min_median_coverage))
    if args.min_median_maf > 0:
        logger.info('Minimum sample median mutant allele frequency (otherwise discarded): {}'.format(
            args.min_median_maf))

    fpr = args.false_positive_rate
    logger.info('False positive rate for the statistical test: {}.'.format(fpr))
    fdr = args.false_discovery_rate
    logger.info('False discovery rate for the statistical test: {}.'.format(fdr))
    min_absent_cov = args.min_absent_coverage
    logger.info('Minimum coverage for an absent variant: {} (otherwise unknown)'.format(min_absent_cov))
    logger.info('Replications for robustness analysis through down-sampling: {}'.format(args.down))

    # ##########################################################################################################
    # ############################################### LOAD DATA ################################################
    # ##########################################################################################################

    # take mutant read and coverage tables to calculate positives, negatives, and unknowns
    if args.mut_reads and args.mut_cov:

        read_table = args.mut_reads
        cov_table = args.mut_cov

        if args.mut_dis_cov:
            dis_cov_table = args.mut_dis_cov
        else:
            logger.warn('No distinct phred coverage data was provided!')
            dis_cov_table = args.mut_cov
            logger.warn('Using phred coverage data for negative/unknown classification.')

        patient_name = get_patients_name(read_table)
        if patient_name.find('_') != -1:
            patient_name = patient_name[:patient_name.find('_')]
        logger.debug('Patient name: {}'.format(patient_name))
        patient = Patient(patient_name, min_absent_cov=min_absent_cov, error_rate=args.error_rate, c0=settings.BI_C0)
        read_no_samples = patient.read_raw_data(
            read_table, cov_table, dis_cov_table, fpr, fdr, min_absent_cov, args.min_median_coverage,
            args.min_median_maf, excluded_columns={normal_sample_name})      # excluded (=normal) samples
        # 'LiM_2' Pam01
        # 'LiM_7', 'LiM_8', 'PT_18' Pam02
        # 'LiM_5', 'LuM_2', 'LuM_3', 'PT_12' Pam03
        # 'PeM_6', 'PT_26', 'PT_27'   Pam04

    elif args.vcf_file:      # take path to the input VCF file
        vcf_file = args.vcf_file
        if not os.path.isfile(vcf_file):
            logger.error("Provided VCF file {} does not exist.".format(vcf_file))
            usage()

        patient_name = get_patients_name(vcf_file)
        patient = Patient(patient_name, error_rate=args.error_rate, c0=settings.BI_C0)
        read_no_samples = patient.read_vcf_file(vcf_file, args.min_median_coverage, args.min_median_maf,
                                                fpr, fdr, args.min_absent_coverage,
                                                normal_sample_name=normal_sample_name)

    elif args.directory:      # take path to the directory with all VCF files
        vcf_directory = args.directory

        if not os.path.isdir(vcf_directory):
            logger.error("Directory named {} does not exist ({}).".format(
                args.directory, os.path.abspath(vcf_directory)))
            usage()

        patient_name = get_patients_name(
            vcf_directory[:-1] if vcf_directory.endswith('/') else vcf_directory)
        patient = Patient(patient_name, error_rate=args.error_rate, c0=settings.BI_C0)
        read_no_samples = patient.read_vcf_directory(vcf_directory, args.min_median_coverage, args.min_median_maf,
                                                     fpr, fdr, min_absent_cov, normal_sample_name)

    else:
        raise RuntimeError('No input files were provided!')

    output_directory = init_output(patient_name=patient_name)
    # create output filename pattern
    fn_pattern = get_output_fn_template(patient.name, read_no_samples, fpr, fdr,
                                        min_absent_cov, args.min_median_coverage, args.min_median_maf,
                                        bi_e=patient.bi_error_rate, bi_c0=patient.bi_c0)
    fn_pattern += '_s' if args.mode == 1 else ''

    # do basic analysis on provided input data
    patient.analyze_data()

    if plots_report:   # deactivate plot generation for debugging

        # generate mutation table plot
        # show only mutations which are present in at least one sample
        if patient.gene_names is not None:
            col_labels = patient.gene_names
        else:
            col_labels = patient.mut_keys

        if len(col_labels) < int_sets.MAX_MUTS_TABLE_PLOT and plots_report:
            mut_table_name = 'mut_table_'+fn_pattern
            plts.hinton(patient.data, os.path.join(output_directory, mut_table_name),
                        row_labels=patient.sample_names, column_labels=col_labels,
                        displayed_mutations=patient.present_mutations)
        else:
            mut_table_name = None
            logger.warn('There are too many variants to create a mutation table plot: {}'.format(len(col_labels)))
    else:
        mut_table_name = None
        col_labels = None

    if plots_paper:     # deactivate plot generation for regular analysis
        # generate scatter plot about raw sequencing data (mutant reads vs. coverage)
        plts.reads_plot(os.path.join(output_directory, 'fig_reads_'+patient.name+'.pdf'), patient)

        # generate box plots with MAFs per sample
        plts.boxplot(os.path.join(output_directory, 'fig_mafs_'+patient.name+'.pdf'), patient)

    # create raw data analysis file
    analysis.create_data_analysis_file(patient, os.path.join(output_directory, 'data_'+fn_pattern+'.txt'))
    # utils.analysis.print_genetic_distance_table(patient)

    # create HTML analysis report
    html_report = HTMLReport(os.path.join(output_directory, 'report_'+fn_pattern+'.html'), patient_name)
    html_report.start_report()
    html_report.add_sequencing_information(
        patient, mut_table_path=mut_table_name+'.png' if mut_table_name is not None else None)
    # html_report.add_similarity_information(patient)       # not relevant for this report

    # ############################################################################################
    # infer evolutionary compatible mutation patterns and subsequently evolutionary trees based on
    # different principles divided into three modes
    # ############################################################################################
    if args.mode == 1 or args.mode == 2:

        phylogeny = None
        comp_node_frequencies = None

        if args.mode == 1:   # find evolutionary incompatible mutation patterns based on standard binary classification

            phylogeny = ti.infer_max_compatible_tree(os.path.join(output_directory, 'btree_'+fn_pattern+'.tex'),
                                                     patient)

            if plots_report:
                # create mutation pattern overview plot
                # show only the different patterns and not the individual variants
                # (convenient for large numbers of variants)
                mp_graph_name = mp_graph.create_mp_graph(
                    fn_pattern, phylogeny, phylogeny.nodes, phylogeny.node_scores,
                    output_directory=output_directory, min_node_weight=settings.MIN_MP_SCORE,
                    max_no_mps=settings.MAX_NO_MPS)

                if mp_graph_name is not None:
                    html_report.add_mp_overview_graph(patient, phylogeny, mp_graph_name)

                # create plot only if there is enough space for all the incompatible mutations
                if len(phylogeny.conflicting_mutations) < int_sets.MAX_MUTS_TABLE_PLOT:
                    # illustrative mutation table plot of incompatible mutation patterns
                    x_length, y_length = plts.create_incompatible_mp_table(
                        patient, os.path.join(output_directory, settings.incomp_mps_plot_prefix+fn_pattern),
                        phylogeny, row_labels=patient.sample_names, column_labels=col_labels)

                    # add information about evolutionarily incompatible mutation patterns to the HTML report
                    html_report.add_inc_mp_information(
                        phylogeny, incomp_mps_plot_filepath=settings.incomp_mps_plot_prefix+fn_pattern+'.png',
                        plot_width=x_length*7)
                else:
                    # too many evolutionarily incompatible mutation to create mutation table plot
                    # add information about evolutionarily incompatible mutation patterns to the HTML report
                    html_report.add_inc_mp_information(phylogeny)

            # do mutation pattern robustness analysis through down-sampling
            if args.down > 0:
                if min_absent_cov > 0:
                    logger.warn('Down-sampling analysis can only run if the minimal absent coverage is 0.')
                else:
                    comp_node_frequencies = phylogeny.validate_node_robustness(args.down)

                    # generate mutation pattern robustness plot
                    plts.robustness_plot(
                        os.path.join(output_directory, 'robustness_{}_{}.pdf'.format(args.down, fn_pattern)),
                        comp_node_frequencies)

        elif args.mode == 2:     # find likely sequencing artifacts based on a bayesian inference model

            # determine mutation patterns based on standard binary classification to generate an overview graph

            phylogeny = ti.create_max_lh_tree(os.path.join(output_directory, 'mlhtree_'+fn_pattern+'.tex'), patient,
                                              min_mp_lh=settings.MIN_MP_LH if 0 < settings.MIN_MP_LH < 1 else None)

            # determine mutation patterns based on standard binary classification to generate an overview graph
            if plots_report:
                # create mutation pattern overview plot
                # show only the different patterns and not the individual variants
                # (convenient for large numbers of variants)
                mp_graph_name = mp_graph.create_mp_graph(
                    fn_pattern, phylogeny, phylogeny.node_scores.keys(), phylogeny.node_scores,
                    output_directory=output_directory, min_node_weight=settings.MIN_MP_SCORE,
                    max_no_mps=settings.MAX_NO_MPS)

                if mp_graph_name is not None:
                    html_report.add_mp_overview_graph(patient, phylogeny, mp_graph_name)

                # create plot only if there is enough space for all the incompatible mutations
                if (len(phylogeny.false_positives) + len(phylogeny.false_negatives)
                        + len(phylogeny.false_negative_unknowns) < int_sets.MAX_MUTS_TABLE_PLOT):

                    # illustrative mutation table plot of incompatible mutation patterns and their putative artifacts
                    x_length, y_length = plts.create_incompatible_mp_table(
                        patient, os.path.join(output_directory, settings.artifacts_plot_prefix+fn_pattern),
                        phylogeny, row_labels=patient.sample_names, column_labels=col_labels)
                    # add information about putative false-positives and false-negatives to the HTML report
                    html_report.add_artifacts_information(
                        phylogeny, artifacts_plot_filepath=settings.artifacts_plot_prefix+fn_pattern+'.png',
                        plot_width=x_length*7)
                else:
                    html_report.add_artifacts_information(phylogeny)

            # do mutation pattern robustness analysis through down-sampling
            if args.down > 0:
                if 0 < settings.MIN_MP_LH < 1:
                    logger.error('Down-sampling analysis does not support subclone detection! ')
                    logger.error('Please set MIN_MP_LH in settings.py to 1.')
                else:
                    comp_node_frequencies = phylogeny.validate_node_robustness(args.down)

                    # generate mutation pattern robustness plot
                    plts.robustness_plot(
                        os.path.join(output_directory, 'robustness_{}_{}.pdf'.format(args.down, fn_pattern)),
                        comp_node_frequencies)

        # generate analysis file to provide an overview about the derived results
        analysis.create_analysis_file(patient, args.min_median_coverage,
                                      os.path.join(output_directory, 'analysis_'+fn_pattern+'.txt'),
                                      phylogeny, comp_node_frequencies, args.down)

        # create input data file for circos conflict graph plots
        if plots_paper:
            circos.create_raw_data_file(os.path.join(output_directory, 'fig_data_'+fn_pattern+'_mutdata.txt'),
                                        patient.mutations, patient.mut_positions, data=patient.data,
                                        sample_names=patient.sample_names)

            # create labels file if gene names are available
            # typically not available in VCF files
            if len(patient.gene_names):
                circos.create_mutation_labels_file(
                    os.path.join(output_directory, 'fig_data_'+fn_pattern+'_mutlabels.txt'),
                    patient.mutations, patient.gene_names, patient.mut_positions, patient.driver_pathways)

            # if there are less than 10000 edges in the conflict graph
            if args.mode == 1:
                if phylogeny.cf_graph.size() < 50000:
                    circos.create_mutation_links_file(
                        os.path.join(output_directory, 'fig_data_'+fn_pattern+'_conflicts.txt'),
                        phylogeny, patient.mut_positions)
                else:
                    logger.warn('Circos mutation conflicts data has not been created as there '
                                'are {} edges in the graph.'.format(phylogeny.cf_graph.size()))

                # create input data files for circular plots with circos: conflict graph
                circos.create_conflict_graph_files(
                    os.path.join(output_directory, 'cfg_nodes_'+fn_pattern+'.txt'),
                    os.path.join(output_directory, 'cfg_mutnode_labels_'+fn_pattern+'.txt'),
                    os.path.join(output_directory, 'cfg_mutnode_data_'+fn_pattern+'.txt'),
                    os.path.join(output_directory, 'cfg_links_'+fn_pattern+'.txt'),
                    phylogeny, patient.gene_names, patient.driver_pathways, data=patient.data,
                    min_node_weight=settings.MIN_MP_SCORE, max_no_mps=settings.MAX_NO_MPS)

            if args.mode == 2:
                # create input data files for circular plots with circos: conflict graph
                circos.create_mlh_graph_files(
                    os.path.join(output_directory, 'mlh_nodes_'+fn_pattern+'.txt'),
                    os.path.join(output_directory, 'mlh_mutnode_labels_'+fn_pattern+'.txt'),
                    os.path.join(output_directory, 'mlh_mutnode_data_'+fn_pattern+'.txt'),
                    patient.data, phylogeny, patient.gene_names, patient.driver_pathways)

    else:
        raise RuntimeError("No mode was provided (e.g. -m 1) to infer the phylogenetic tree.")

    # finalize HTML report
    html_report.end_report(fpr, fdr, min_absent_cov, args.min_median_coverage, args.min_median_maf)
    logger.info('Treeomics finished evolutionary analysis.')

if __name__ == '__main__':

    main()
