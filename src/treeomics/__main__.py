#!/usr/bin/python
import logging
import sys
import os
import argparse
from collections import Counter
from utils.driver import get_drivers, can_be_driver
from utils.mutation_effects import get_top_effect_name
from patient import Patient
import settings
import utils.int_settings as def_sets
import tree_inference as ti
from utils.html_report import HTMLReport
from utils.analysis import analyze_data
import plots.plots_utils as plts
import plots.mp_graph as mp_graph
import plots.circos as circos
import utils.analysis as analysis


"""Main file to run Treeomics"""
__author__ = 'Johannes REITER, IST Austria - Harvard University'
__date__ = 'March 31, 2014'


# create logger for application
logger = logging.getLogger('treeomics')
logger.setLevel(logging.INFO)
# create file handler which logs even debug messages
fh = logging.FileHandler('treeomics.log')
fh.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
# formatter = logging.Formatter('%(asctime)s %(levelname)s %(pathname)s %(lineno)s: %(message)s')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(filename)s %(lineno)s: %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


def init_output(patient_name=None, output_dir=None):
    """
    Set up output directory.
    :param patient_name: if patient name is given generate separate directory for each subject
    :param output_dir: name of the output directory
    :return: path to output directory
    """

    # derive correct output directory
    if output_dir is None:
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
    else:
        output_directory = os.path.join(output_dir, '')

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


def get_output_fn_template(name, total_no_samples, subclone_detection=False, fpr=None, fdr=None,
                           min_absent_coverage=None, min_sa_coverage=None, min_sa_vaf=None,
                           bi_e=None, bi_c0=None, max_absent_vaf=None, no_boot=None, mode=None, max_no_mps=None):
    """
    Generate common name template for output files
    :param name: subject name or id
    :param total_no_samples: number of considered sequencing samples
    :param subclone_detection: is subclone detection enabled?
    :param fpr: false-positive rate
    :param fdr: false-discovery rate
    :param min_absent_coverage: minimum coverage such that variant can be classified as absent in the classical model
    :param min_sa_coverage: minimum median coverage per sample
    :param min_sa_vaf: minimum median variant allele frequency per sample
    :param bi_e: sequencing error for bayesian inference
    :param bi_c0: prior mixture parameter of delta function and uniform distribution for bayesian inference
    :param max_absent_vaf: maximal absent VAF before considering estimated purity
    :param no_boot: number of bootstrapping samples
    :param mode: mode of Treeomics
    :param max_no_mps: maximal number of considered most likely distinct mutation patterns per variant
    :return: output filename pattern
    """

    pattern = ('{}_{}{}'.format(name, total_no_samples, '_SC' if subclone_detection else '') +
               ('_st={}'.format(min_sa_coverage) if min_sa_coverage > 0 else '') +
               ('_mf={}'.format(min_sa_vaf) if min_sa_vaf > 0 else '') +
               ('_fpr={:.1e}_fdr={:.2f}'.format(fpr, fdr) if fpr is not None or fdr is not None else '') +
               ('_at={}'.format(min_absent_coverage) if min_absent_coverage is not None else '') +
               ('_vaf={}'.format(settings.MIN_VAF) if settings.MIN_VAF is not None and settings.MIN_VAF > 0 else '') +
               ('_var={}'.format(settings.MIN_VAR_READS) if settings.MIN_VAR_READS is not None and
                settings.MIN_VAR_READS > 0 else '') +
               ('_e={}'.format(bi_e) if bi_e is not None else '') +
               ('_c0={}'.format(bi_c0) if bi_c0 is not None else '') +
               ('_af={}'.format(max_absent_vaf) if max_absent_vaf is not None else '') +
               ('_mps={}'.format(max_no_mps) if max_no_mps is not None else '') +
               ('_b={}'.format(no_boot) if no_boot is not None and no_boot > 0 else '') +
               ('_s' if mode is not None and mode == 2 else ''))

    # replace points with commas because latex cannot handle points in file names (interprets it as file type)
    pattern = pattern.replace('.', '_')

    return pattern


def usage():
    """
    Give the user feedback on how to call the tool
    Terminates the tool afterwards
    """
    logger.warn("Usage: python treeomics.zip -r <mut-reads table> -s <coverage table> | -v vcf_file | "
                "-d vcf_file_directory  [-n <normal sample name>] [-e <sequencing error rate] "
                "[-z <prior absent probability>] [-c <minimum sample median coverage>] [] \n")
    logger.warn("Example: python treeomics.zip -r input/Makohon2015/Pam03_mutant_reads.txt "
                "-s input/Makohon2015/Pam03_phredcoverage.txt ")
    sys.exit(2)


# -------------------------------- main function ---------------------------------
def main():

    parser = argparse.ArgumentParser(description='Infers the evolution of cancer.')

    parser.add_argument("-m", "--mode", help="running mode: 1...run Treeomics and explore full solution space, "
                                             "2...fast (one mutation pattern per variant).",
                        type=int, default=1)

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--csv_file", help="path to the CSV file", type=str)
    group.add_argument("-v", "--vcf_file", help="path to the VCF file", type=str)
    group.add_argument("-d", "--directory", help="directory with multiple VCF files", type=str)

    parser.add_argument("-n", "--normal", help="name of normal sample (excluded from analysis)", type=str, default=None)

    parser.add_argument("-r", "--mut_reads", help="path table with the number of reads with a mutation", type=str)
    parser.add_argument("-s", "--coverage", help="path to table with read coverage at the mutated positions", type=str)

    # specify output directory
    parser.add_argument("-o", "--output", help="output directory", type=str, default=settings.OUTPUT_FOLDER)

    # read in parameter value
    parser.add_argument("-e", "--error_rate", help="sequencing error rate for bayesian inference",
                        type=float, default=settings.BI_E)
    parser.add_argument("-z", "--prob_zero", help="prior probability of being absent",
                        type=float, default=settings.BI_C0)
    parser.add_argument("-a", "--max_absent_vaf", help="maximal absent vaf before considering purity",
                        type=float, default=settings.MAX_ABSENT_VAF)

    parser.add_argument("-c", "--min_median_coverage",
                        help="minimum median coverage of a sample",
                        type=int, default=settings.SAMPLE_COVERAGE_THRESHOLD)
    parser.add_argument("-f", "--min_median_vaf",
                        help="minimum median mutant allele frequency of a sample",
                        type=float, default=settings.MAF_THRESHOLD)

    parser.add_argument("-g", "--ref_genome",
                        help="to which reference genome was the sequencing data aligned",
                        type=str, default=settings.REF_GENOME)

    # DEPRECATED FROM VERSION 1.7.0 ONWARD
    parser.add_argument("-p", "--false_positive_rate",
                        help="false positive rate for the statistical test",
                        type=float, default=settings.FPR)
    parser.add_argument("-i", "--false_discovery_rate",
                        help="false discovery rate for the statistical test",
                        type=float, default=settings.FDR)
    parser.add_argument("-y", "--min_absent_coverage",
                        help="minimum coverage for a true negative (negative threshold)",
                        type=int, default=settings.MIN_ABSENT_COVERAGE)

    # -v, --verbose Verbose mode to provide more information while running treeomics

    plots_parser = parser.add_mutually_exclusive_group(required=False)
    plots_parser.add_argument('--plots', dest='plots', action='store_true', help="Is plot generation enabled?")
    plots_parser.add_argument('--no_plots', dest='plots', action='store_false', help="Is plot detection disabled?")
    parser.set_defaults(plots=True)

    parser.add_argument('-b', '--boot', help='Number of bootstrapping samples', type=int,
                        default=settings.NO_BOOTSTRAP_SAMPLES)
    parser.add_argument('--pool_size', help='number of best solutions explored by ILP solver', type=int,
                        default=settings.POOL_SIZE)

    # limit search space exploration to decrease the run time
    parser.add_argument("-t", "--time_limit",
                        help="maximum running time for CPLEX to solve the MILP",
                        type=int, default=settings.TIME_LIMIT)
    parser.add_argument("-l", "--max_no_mps", help="limit the solution space size by the maximal number of " +
                                                   "explored mutation patterns per variant",
                        type=int, default=settings.MAX_NO_MPS)

    feature_parser = parser.add_mutually_exclusive_group(required=False)
    feature_parser.add_argument('-u', '--subclone_detection', dest='subclone_detection', action='store_true',
                                help="Is subclone detection enabled?")
    feature_parser.add_argument('--no_subclone_detection', dest='subclone_detection', action='store_false',
                                help="Is subclone detection disabled?")
    parser.set_defaults(subclone_detection=settings.SUBCLONE_DETECTION)

    args = parser.parse_args()

    # disable plot generation if the script is called from another script (benchmarking)
    # set log level to info
    if not os.getcwd().endswith('treeomics') and not os.getcwd().endswith('src'):
        logger.setLevel(logging.INFO)
        fh.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)

    if args.subclone_detection:
        logger.info('Subclone detection is enabled.')
    else:
        logger.info('Subclone detection is disabled.')

    if args.plots:
        logger.info('Plot generation is enabled.')
    else:
        logger.info('Plot generation is disabled.')

    plots_report = args.plots    # for debugging set to False
    plots_paper = logger.isEnabledFor(logging.DEBUG)
    # plots_report = False  # for debugging set to False
    # plots_paper = False

    if args.normal:
        normal_sample_name = args.normal
        logger.info('Exclude normal sample with name: {}'.format(normal_sample_name))
    else:
        normal_sample_name = None

    if args.min_median_coverage > 0:
        logger.info('Minimum sample median coverage (otherwise discarded): {}'.format(
            args.min_median_coverage))
    if args.min_median_vaf > 0:
        logger.info('Minimum sample median mutant allele frequency (otherwise discarded): {}'.format(
            args.min_median_vaf))

    if args.boot is not None:
        if args.boot < 0:
            raise AttributeError('Number of bootstrapping samples can not be negative!')

    if args.time_limit is not None:
        if args.time_limit <= 0:
            raise AttributeError('Time limit for the MILP solver needs to be positive!')
        logger.info('MILP solver running time is limited to {} seconds. Obtained solution may not be optimal.'.format(
            args.time_limit))

    if args.max_no_mps is not None:
        if args.max_no_mps <= 0:
            raise AttributeError('Solution space can only be limited to a positive number of mutation patterns!')
        logger.info('Solution space is limited to the {} most likely mutation patterns per variant.'.format(
            args.max_no_mps))

    if args.max_absent_vaf < args.error_rate:
        raise AttributeError('The maximal absent VAF has to be larger than the error rate in the Bayesian model!')

    if args.pool_size < 1:
        raise AttributeError('Solution pool size of ILP solver has to be at least 1!')

    # DEPRECATED FROM VERSION 1.7.0 ONWARD
    fpr = args.false_positive_rate
    logger.debug('False positive rate for the statistical test: {}.'.format(fpr))
    fdr = args.false_discovery_rate
    logger.debug('False discovery rate for the statistical test: {}.'.format(fdr))
    min_absent_cov = args.min_absent_coverage
    logger.debug('Minimum coverage for an absent variant: {} (otherwise unknown)'.format(min_absent_cov))

    ref_genome = args.ref_genome.lower()
    logger.info('Sequencing data was aligned to reference genome: {}'.format(ref_genome))

    if args.boot > 0:
        logger.info('Number of samples for bootstrapping analysis: {}'.format(args.boot))

    if args.subclone_detection and args.max_no_mps is not None:
        logger.error('Subclone and partial solution space search are not supported to be performed at the same time! ')
        usage()

    # ##########################################################################################################
    # ############################################### LOAD DATA ################################################
    # ##########################################################################################################

    # take mutant read and coverage tables to calculate positives, negatives, and unknowns
    if args.mut_reads and args.coverage:

        patient_name = get_patients_name(args.mut_reads)
        if patient_name.find('_') != -1:
            patient_name = patient_name[:patient_name.find('_')]
        logger.debug('Patient name: {}'.format(patient_name))
        patient = Patient(error_rate=args.error_rate, c0=args.prob_zero, max_absent_vaf=args.max_absent_vaf,
                          pat_name=patient_name, min_absent_cov=min_absent_cov, reference_genome=ref_genome)
        read_no_samples = patient.process_raw_data(
            fpr, fdr, min_absent_cov, args.min_median_coverage, args.min_median_vaf, var_table=args.mut_reads,
            cov_table=args.coverage, normal_sample=normal_sample_name)

    elif args.csv_file:

        patient_name = get_patients_name(args.csv_file)
        if patient_name.find('_') != -1:
            patient_name = patient_name[:patient_name.find('_')]
        logger.debug('Patient name: {}'.format(patient_name))
        patient = Patient(error_rate=args.error_rate, c0=args.prob_zero, max_absent_vaf=args.max_absent_vaf,
                          pat_name=patient_name, min_absent_cov=min_absent_cov, reference_genome=ref_genome)
        read_no_samples = patient.process_raw_data(
            fpr, fdr, min_absent_cov, args.min_median_coverage, args.min_median_vaf, csv_file=args.csv_file,
            normal_sample=normal_sample_name)

    elif args.vcf_file:      # take path to the input VCF file
        vcf_file = args.vcf_file
        if not os.path.isfile(vcf_file):
            logger.error("Provided VCF file {} does not exist.".format(vcf_file))
            usage()

        patient_name = get_patients_name(vcf_file)
        if patient_name.find('_') != -1:
            patient_name = patient_name[:patient_name.find('_')]
        patient = Patient(error_rate=args.error_rate, c0=args.prob_zero, max_absent_vaf=args.max_absent_vaf,
                          pat_name=patient_name, reference_genome=ref_genome)
        read_no_samples = patient.read_vcf_file(vcf_file, fpr, fdr, min_sa_cov=args.min_median_coverage,
                                                min_sa_maf=args.min_median_vaf, min_absent_cov=args.min_absent_coverage,
                                                normal_sample_name=normal_sample_name)

    elif args.directory:      # take path to the directory with all VCF files
        vcf_directory = args.directory

        if not os.path.isdir(vcf_directory):
            logger.error("Directory named {} does not exist ({}).".format(
                args.directory, os.path.abspath(vcf_directory)))
            usage()

        patient_name = get_patients_name(
            vcf_directory[:-1] if vcf_directory.endswith('/') else vcf_directory)
        patient = Patient(error_rate=args.error_rate, c0=args.prob_zero, max_absent_vaf=args.max_absent_vaf,
                          pat_name=patient_name, reference_genome=ref_genome)
        read_no_samples = patient.read_vcf_directory(vcf_directory, args.min_median_coverage, args.min_median_vaf,
                                                     fpr, fdr, min_absent_cov, normal_sample_name)

    else:
        raise RuntimeError('No input files were provided!')

    cgc_drivers, user_drivers = get_drivers(settings.CGC_PATH, settings.DRIVER_PATH, ref_genome)

    # check if any of the variants are a putative driver
    likely_drivers = dict()
    unlikely_driver_mut_effects = Counter()
    if patient.gene_names is not None or patient.ensembl_data is not None:

        for mut_idx, mut_pos in enumerate(patient.mut_positions):

            if patient.gene_names[mut_idx] in user_drivers.keys():

                # is there a reference genome to predict the effect of the mutation?
                if patient.ensembl_data is not None:
                    mut_effect = get_top_effect_name(patient.vc_variants[mut_idx])
                    if not can_be_driver(mut_effect):
                        unlikely_driver_mut_effects[mut_effect] += 1
                        logger.info('{} variant found in likely driver {}.'.format(
                            mut_effect, patient.gene_names[mut_idx]))
                        continue
                else:
                    mut_effect = None

                # are there known positions for this driver gene?
                if cgc_drivers is not None and patient.gene_names[mut_idx] in cgc_drivers:
                    dri_pos = cgc_drivers[patient.gene_names[mut_idx]].genomic_location
                else:
                    dri_pos = None

                if dri_pos is None:
                    # no positions provided => assume it's a driver
                    likely_drivers[patient.gene_names[mut_idx]] = \
                        (mut_effect, user_drivers[patient.gene_names[mut_idx]].sources, False)

                # check if variant is at the same chromosome and is within given region
                elif (dri_pos[0] == mut_pos[0] and
                      (dri_pos[1] <= int(mut_pos[1]) <= dri_pos[2] or dri_pos[1] <= int(mut_pos[2]) <= dri_pos[2])):
                    # variant is among CGC region
                    likely_drivers[patient.gene_names[mut_idx]] = \
                        (mut_effect, user_drivers[patient.gene_names[mut_idx]].sources, True)

    output_directory = init_output(patient_name=patient_name,
                                   output_dir=args.output if args.output is not settings.OUTPUT_FOLDER else None)
    # create output filename pattern
    fn_pattern = get_output_fn_template(
        patient.name, read_no_samples, mode=args.mode, min_sa_coverage=args.min_median_coverage,
        min_sa_vaf=args.min_median_vaf, bi_e=patient.bi_error_rate, bi_c0=patient.bi_c0,
        max_absent_vaf=patient.max_absent_vaf)

    # do basic analysis on provided input data
    analyze_data(patient, post_table_filepath=os.path.join(output_directory,
                 get_output_fn_template(patient.name, read_no_samples, min_sa_coverage=args.min_median_coverage,
                                        min_sa_vaf=args.min_median_vaf, bi_e=patient.bi_error_rate,
                                        bi_c0=patient.bi_c0, max_absent_vaf=patient.max_absent_vaf) +
                                        '_posterior.txt'))

    if plots_report:   # deactivate plot generation for debugging and benchmarking

        # generate mutation table plot
        # show only mutations which are present in at least one sample
        if patient.gene_names is not None:
            col_labels = patient.gene_names
        else:
            col_labels = patient.mut_keys

        if len(col_labels) < def_sets.MAX_MUTS_TABLE_PLOT and plots_report:
            mut_table_name = \
                (get_output_fn_template(
                    patient.name, read_no_samples, min_sa_coverage=args.min_median_coverage,
                    min_sa_vaf=args.min_median_vaf, bi_e=patient.bi_error_rate, bi_c0=patient.bi_c0,
                    max_absent_vaf=patient.max_absent_vaf)+'_bayesian_data_table')

            plts.bayesian_hinton(patient.log_p01, output_directory, mut_table_name,
                                 row_labels=patient.sample_names, column_labels=col_labels,
                                 displayed_mutations=patient.present_mutations, drivers=likely_drivers)
        else:
            logger.warning('Too many reported variants for a detailed mutation table plot: {}'.format(len(col_labels)))
            mut_table_name = None
    else:
        mut_table_name = None
        col_labels = None

    if plots_paper:     # deactivate plot generation for regular analysis
        # generate violin coverage distribution plot
        plts.coverage_plot(os.path.join(output_directory, patient.name+'_coverage_distr.pdf'), patient)

        # generate box plots with MAFs per sample
        # plts.boxplot(os.path.join(output_directory, 'fig_mafs_'+patient.name+'.pdf'), patient)
        # generate violin VAF distribution plot
        plts.vaf_distribution_plot(os.path.join(output_directory, patient.name+'_vaf_distr.pdf'), patient)

    # create raw data analysis file
    analysis.create_data_analysis_file(patient, os.path.join(output_directory, fn_pattern+'_data.txt'))
    # utils.analysis.print_genetic_distance_table(patient)

    # create output filename pattern
    fn_pattern = get_output_fn_template(
        patient.name, read_no_samples, mode=args.mode, min_sa_coverage=args.min_median_coverage,
        min_sa_vaf=args.min_median_vaf, bi_e=patient.bi_error_rate, bi_c0=patient.bi_c0,
        max_absent_vaf=patient.max_absent_vaf, max_no_mps=args.max_no_mps)

    # create HTML analysis report
    html_report = HTMLReport(os.path.join(output_directory, fn_pattern+'_report.html'), patient_name)
    html_report.start_report()
    html_report.add_sequencing_information(
        patient, mut_table_path=mut_table_name+'.png' if mut_table_name is not None else None,
        likely_drivers=likely_drivers, unlikely_driver_mut_effects=unlikely_driver_mut_effects)

    html_report.add_similarity_information(patient)

    # ############################################################################################
    # infer evolutionary compatible mutation patterns and subsequently evolutionary trees based on
    # different principles divided into three modes
    # ############################################################################################
    if args.mode == 1 or args.mode == 2:

        phylogeny = None
        # comp_node_frequencies = None

        # ### RUN TREEOMICS ###
        if args.mode == 1:     # find likely sequencing artifacts based on a bayesian inference model

            # generate filename for tree
            fn_tree = get_output_fn_template(
                patient.name, read_no_samples, subclone_detection=args.subclone_detection,
                min_sa_coverage=args.min_median_coverage, min_sa_vaf=args.min_median_vaf, no_boot=args.boot,
                max_no_mps=args.max_no_mps, bi_e=patient.bi_error_rate, bi_c0=patient.bi_c0,
                max_absent_vaf=patient.max_absent_vaf, mode=args.mode)
            fn_matrix = get_output_fn_template(
                patient.name, read_no_samples, subclone_detection=args.subclone_detection,
                min_sa_coverage=args.min_median_coverage, min_sa_vaf=args.min_median_vaf, max_no_mps=args.max_no_mps,
                bi_e=patient.bi_error_rate, bi_c0=patient.bi_c0, max_absent_vaf=patient.max_absent_vaf, mode=args.mode)
            # infer maximum likelihood tree
            phylogeny = ti.create_max_lh_tree(
                patient, tree_filepath=os.path.join(output_directory, fn_tree+'_mlhtree'),
                mm_filepath=os.path.join(output_directory, fn_matrix+'_treeomics_mm.csv'),
                mp_filepath=os.path.join(output_directory, fn_matrix+'_treeomics_mps.tsv'),
                subclone_detection=args.subclone_detection, loh_frequency=settings.LOH_FREQUENCY,
                drivers=likely_drivers.keys(), pool_size=args.pool_size, no_bootstrap_samples=args.boot,
                max_no_mps=args.max_no_mps, time_limit=args.time_limit, plots=plots_report)

            # previously used for benchmarking
            # if plots_paper:     # generate Java Script D3 trees
            #     json_file = 'mlhtree_'+fn_tree+'.json'
            #     Phylogeny.save_json_tree(os.path.join(output_directory, json_file), phylogeny.mlh_tree)
            #     Phylogeny.write_html_file(os.path.join(output_directory, 'mlhtree_'+fn_tree+'.html'), json_file)

            # determine mutation patterns and generate an overview graph
            if plots_report:
                if phylogeny.tree_plot is not None:     # add generated tree plot to HTML report
                    html_report.add_tree_plot(patient, phylogeny)

                # create mutation pattern overview plot
                # show only the different patterns and not the individual variants
                # (convenient for large numbers of variants)
                fn_mp_plot = get_output_fn_template(
                    patient.name, read_no_samples, min_sa_coverage=args.min_median_coverage, max_no_mps=args.max_no_mps,
                    min_sa_vaf=args.min_median_vaf, bi_e=patient.bi_error_rate, bi_c0=patient.bi_c0,
                    max_absent_vaf=patient.max_absent_vaf, mode=args.mode)

                if args.subclone_detection:
                    pg = ti.create_max_lh_tree(
                        patient, tree_filepath=None, mm_filepath=None, mp_filepath=None, subclone_detection=False,
                        loh_frequency=settings.LOH_FREQUENCY, drivers=likely_drivers.keys(), pool_size=args.pool_size,
                        no_bootstrap_samples=0, max_no_mps=args.max_no_mps, time_limit=args.time_limit, plots=False)
                else:
                    pg = phylogeny

                mp_graph_name = mp_graph.create_mp_graph(
                    fn_mp_plot, pg, pg.node_scores.keys(), pg.node_scores,
                    output_directory=output_directory, min_node_weight=settings.MIN_MP_SCORE,
                    circos_max_no_mps=settings.CIRCOS_MAX_NO_MPS)

                if mp_graph_name is not None:
                    html_report.add_conflict_graph(patient, mp_graph_name, phylogeny=pg)

                opt_sol = phylogeny.solutions[0]  # optimal solution
                # create plot only if there is enough space for all the incompatible mutations
                if (0 < len(opt_sol.false_positives) + len(opt_sol.false_negatives) +
                        len(opt_sol.false_negative_unknowns) < def_sets.MAX_MUTS_TABLE_PLOT):

                    # illustrative mutation table plot of incompatible mutation patterns and their putative artifacts
                    x_length, y_length = plts.create_incompatible_mp_table(
                        patient, os.path.join(output_directory, fn_pattern+settings.artifacts_plot_suffix),
                        phylogeny, row_labels=patient.sample_names, column_labels=col_labels)
                    if x_length > 0 and y_length > 0:
                        # add information about putative false-positives and false-negatives to the HTML report
                        html_report.add_artifacts_information(
                            phylogeny, artifacts_plot_filepath=fn_pattern+settings.artifacts_plot_suffix+'.png',
                            plot_width=x_length*7)
                else:
                    html_report.add_artifacts_information(phylogeny)

        # find evolutionary incompatible mutation patterns based on standard binary classification
        elif args.mode == 2:

            phylogeny = ti.infer_max_compatible_tree(os.path.join(output_directory, fn_pattern+'_btree.tex'),
                                                     patient, drivers=likely_drivers.keys())

            if plots_report:
                # create mutation pattern overview plot
                # show only the different patterns and not the individual variants
                # (convenient for large numbers of variants)
                mp_graph_name = mp_graph.create_mp_graph(
                    fn_pattern, phylogeny, phylogeny.nodes, phylogeny.node_scores,
                    output_directory=output_directory, min_node_weight=settings.MIN_MP_SCORE,
                    circos_max_no_mps=settings.CIRCOS_MAX_NO_MPS)

                if mp_graph_name is not None:
                    html_report.add_conflict_graph(patient, mp_graph_name)

                # create plot only if there is enough space for all the incompatible mutations
                if len(phylogeny.conflicting_mutations) < def_sets.MAX_MUTS_TABLE_PLOT:
                    # illustrative mutation table plot of incompatible mutation patterns
                    x_length, y_length = plts.create_incompatible_mp_table(
                        patient, os.path.join(output_directory, fn_pattern+settings.incomp_mps_plot_suffix),
                        phylogeny, row_labels=patient.sample_names, column_labels=col_labels)

                    if x_length > 0 and y_length > 0:
                        # add information about evolutionarily incompatible mutation patterns to the HTML report
                        html_report.add_inc_mp_information(
                            phylogeny, incomp_mps_plot_filepath=fn_pattern+settings.incomp_mps_plot_suffix + '.png',
                            plot_width=x_length*7)
                else:
                    # too many evolutionarily incompatible mutation to create mutation table plot
                    # add information about evolutionarily incompatible mutation patterns to the HTML report
                    html_report.add_inc_mp_information(phylogeny)

        # generate analysis file to provide an overview about the derived results
        analysis.create_analysis_file(patient, args.min_median_coverage,
                                      os.path.join(output_directory, fn_pattern+'_analysis.txt'), phylogeny)

        # create input data file for circos conflict graph plots
        if plots_report:
            circos.create_raw_data_file(os.path.join(output_directory, 'fig_data_'+fn_pattern+'_mutdata.txt'),
                                        patient.mutations, patient.mut_positions, data=patient.data,
                                        sample_names=patient.sample_names)

            # create labels file if gene names are available
            # typically not available in VCF files
            if patient.gene_names is not None and len(patient.gene_names):
                circos.create_mutation_labels_file(
                    os.path.join(output_directory, 'fig_data_'+fn_pattern+'_mutlabels.txt'),
                    patient.mutations, patient.gene_names, patient.mut_positions, patient.driver_pathways)

            if args.mode == 1 and len(patient.data) < 500:
                # create input data files for circular plots with circos: conflict graph
                circos.create_mlh_graph_files(
                    os.path.join(output_directory, 'mlh_nodes_'+fn_pattern+'.txt'),
                    os.path.join(output_directory, 'mlh_mutnode_labels_'+fn_pattern+'.txt'),
                    os.path.join(output_directory, 'mlh_mutnode_data_'+fn_pattern+'.txt'),
                    patient.data, phylogeny, patient.gene_names, patient.driver_pathways)

            # if there are less than 10000 edges in the conflict graph
            elif args.mode == 2:
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
                    min_node_weight=settings.MIN_MP_SCORE, max_no_mps=settings.CIRCOS_MAX_NO_MPS)

    else:
        raise RuntimeError("No mode was provided (e.g. -m 1) to infer the phylogenetic tree.")

    # finalize HTML report
    html_report.end_report(patient.bi_error_rate, patient.bi_c0, patient.max_absent_vaf, settings.LOH_FREQUENCY,
                           fpr, fdr, min_absent_cov, args.min_median_coverage, args.min_median_vaf,
                           max_no_mps=args.max_no_mps)
    logger.info('Treeomics finished evolutionary analysis.')

if __name__ == '__main__':

    main()
