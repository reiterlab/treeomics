#!/usr/bin/python
"""Data structure around sequencing data of a subject"""
import logging
from collections import defaultdict, Counter
import re
import heapq
import numpy as np
import math
import settings
import utils.int_settings as def_sets
from utils.int_settings import NEG_UNKNOWN, POS_UNKNOWN
from utils.vcf_parser import read_vcf_files, read_vcf_file
from utils.filtering import is_intronic, is_incompletetranscript, is_intergenic
from utils.data_tables import read_mutation_table, read_csv_file
from utils.statistics import calculate_present_pvalue, find_significant_mutations
from utils.vaf_data import calculate_p_values
from utils.statistics import get_log_p0

try:    # check if varcode and pyensembl is available (necessary for Windows)
    from varcode import Variant as VCVariant  # https://github.com/hammerlab/varcode
    from pyensembl import EnsemblRelease
    from utils.mutation_effects import get_top_effect_name

except ImportError:
    # mutation effect prediction will not be performed since VarCode is not available
    get_top_effect_name = None

__author__ = 'jreiter'
__date__ = 'April, 2014'


# get logger for application
logger = logging.getLogger('treeomics')
re_sample = re.compile(r'(?:Pam[0-9]{2})?PT[0-9]{1,2}[B]?')


class Patient(object):
    """ Patient: sample data processing
    """

    def __init__(self, error_rate, c0, max_absent_vaf, pat_name='Patient', min_absent_cov=0, reference_genome='grch37',
                 purities=None):
        """
        Initialize patient
        :param error_rate:
        :param c0:
        :param max_absent_vaf:
        :param pat_name:
        :param min_absent_cov:
        :param reference_genome: reference genome name, see https://github.com/hammerlab/pyensembl
        :param purities: dictionary of externally estimated purities for each included sample
        """

        self.name = pat_name

        # minimum coverage for an absent variant
        self.min_absent_cov = min_absent_cov

        # observed Variant Allele Frequency
        self.vafs = None
        # allele frequency data in each numbered sample (frequentist)
        self.data = defaultdict(list)

        # list of tuple posterior (log probability that VAF = 0, log probability that VAF > 0) for each variant
        self.log_p01 = defaultdict(list)
        self.bi_error_rate = error_rate     # sequencing error rate for bayesian inference
        self.bi_c0 = c0    # prior mixture parameter of delta function and uniform distribution for bayesian inference
        self.max_absent_vaf = max_absent_vaf  # maximal absent VAF before considering estimated purity
        self.betas = None
        logger.info('Bayesian inference model: error rate e {}, prior weight c0 {}, max absent vaf {}.'.format(
            self.bi_error_rate, self.bi_c0, self.max_absent_vaf))

        # raw sequencing data of the variants
        self.mut_reads = None
        self.coverage = None
        # median coverages per sample and median MAFs of all confirmed present mutations per sample
        self.sample_coverages = None
        self.sample_mafs = None
        # either externally estimated purities or estimated purities based on VAF of shared variants
        self.purities = purities

        self.positives = None
        self.unknowns = None
        self.negatives = None

        self.present_p_values = None

        # sample contain which numbered mutations
        self.samples = defaultdict(set)
        self.discarded_samples = 0

        # mutation is present in which numbered samples
        self.mutations = defaultdict(set)
        # list of mutations which are present in some of the samples
        self.present_mutations = None

        # number of present and absent mutations according to the Bayesian inference model per sample
        self.no_present_vars = None
        self.no_absent_vars = None

        # variant statistics: 1...passed filtering, -1 did not reach significant level, -2...intronic,
        #                     -3...intergenic, -4...incomplete transcript annotation, -5...common variant in normals
        #                     -6...present in normal samples
        self.variant_stats = None

        self.n = 0      # read samples
        # holds the names for the numbered samples
        self.sample_names = list()
        # holds the names of the numbered mutations
        self.mut_keys = []
        # holds the gene names of the numbered mutations
        self.gene_names = None
        # holds the type/class of the numbered mutation, e.g. missense, nonsense, etc.
        self.mut_types = None
        # holds the data about the mutation: chromosome, start position, and end position
        self.mut_positions = []
        # dictionary with the name of the driver pathway if the mutation is a driver gene mutation
        self.driver_pathways = dict()

        # list of varcode variant class instances
        self.vc_variants = None

        try:
            from varcode import Variant as VCVariant  # https://github.com/hammerlab/varcode
            from pyensembl import EnsemblRelease

            # release 77 uses human reference genome GRCh38/hg38
            # release 75 uses human reference genome GRCh37/hg19
            # release 54 uses human reference genome GRCh36/hg18
            if reference_genome is None:
                self.ensembl_data = None
            else:
                logger.info('VarCode is used for mutation effect prediction.')

            if reference_genome == 'grch36':
                self.ensembl_data = EnsemblRelease(54)

            elif reference_genome == 'grch37':
                self.ensembl_data = EnsemblRelease(75)

            elif reference_genome == 'grch38':
                self.ensembl_data = EnsemblRelease(77)

            else:
                self.ensembl_data = None
                logger.error('Provided reference genome {} is not supported. '.format(reference_genome)
                             + 'Please use GRCh36 (hg18), GRCh37 (hg19), or GRCh38 (hg38).')

        except ImportError as ie:
            self.ensembl_data = None
            logger.warning('ImportError! VarCode is not installed! No mutation effect prediction. {}'.format(ie))

        # mutations present in all samples
        self.founders = set()
        # mutations ordered via the numbered of shared samples
        self.shared_muts = defaultdict(set)
        # sharing status according to Bayesian inference model (list indexed by mut IDs)
        self.sharing_status = None

        # holds a dictionary with a frozenset of samples mapping to a set of mutations
        # present in exactly the same set
        self.mps = None         # mutation patterns

        self.subclones = None
        self.sc_names = None
        self.updated_clones = None

        self.common_muts = []
        self.add_muts = []

        # based on frequentist approach
        self.gen_dis = None         # genetic distance between any pair of samples
        self.sim_coeff = None        # Jaccard similarity coefficient between any pair of samples

        # based on Bayesian inference model
        self.bi_gen_dis = None         # genetic distance between any pair of samples
        self.bi_sim_coeff = None        # Jaccard similarity coefficient between any pair of samples

        # same as above just in the form of a pandas dataframe
        self.df_bi_gen_dist = None      # genetic distance between any pair of samples
        self.df_bi_jsc = None           # Jaccard similarity coefficient between any pair of samples

    def process_raw_data(self, false_positive_rate, false_discovery_rate, min_absent_cov, min_sa_cov, min_sa_maf,
                         mut_reads_normal_th, vaf_normal_th,
                         var_table=None, cov_table=None, csv_file=None, normal_sample=None, excluded_columns=set(),
                         considered_samples=None, wes_filtering=False, artifacts=None):
        """
        Read raw sequencing data from tsv files
        :param false_positive_rate: false positive read of the used sequencing technology
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_absent_cov: minimum coverage at a non-significantly mutated position for absence classification
        :param min_sa_cov: minimum median coverage per sample
        :param min_sa_maf: minimum median mutant allele frequency per sample
        :param mut_reads_normal_th: exclude variants reaching this number of mutant reads with a given VAF in the normal
        :param vaf_normal_th: exclude variants reaching this VAF with a given number of mutant reads in the normal
        :param var_table: path to file with mutant reads
        :param cov_table: path to file with phred coverage
        :param csv_file: path to CSV file with sequencing data of all samples
        :param normal_sample: name of normal sample
        :param excluded_columns: matched normal sample or other samples to be excluded from the analysis
        :param considered_samples: if not None then only samples included in this set will be considered
        :param wes_filtering: remove intronic and intergenic variants due to frequent sequencing artifacts in whole
                              exome sequencing data
        :param artifacts: likely sequencing artifacts that are excluded from the analysis
        """

        if var_table is not None and cov_table is not None:
            # read sequencing data from tsv files
            self.mut_reads, gene_names, norm_var = read_mutation_table(
                var_table, normal_sample=normal_sample, excluded_columns=excluded_columns,
                considered_samples=considered_samples)
            self.coverage, _, norm_cov = read_mutation_table(cov_table, normal_sample=normal_sample,
                                                             excluded_columns=excluded_columns)

        elif csv_file is not None:
            # read sequencing data from csv file
            self.coverage, self.mut_reads, gene_names, norm_cov, norm_var = read_csv_file(
                csv_file, normal_sample=normal_sample, excluded_columns=excluded_columns,
                considered_samples=considered_samples)

        else:
            raise AttributeError('Either TSV or CSV files need to be provided!')

        assert len(self.mut_reads) == len(self.coverage), \
            'Number of reported variants is different in the data files.'

        if artifacts is None:
            # no potential dictionary of sequencing artifacts was given
            artifacts = dict()

        # - - - classifier for positives, negatives, positive unknowns and negative unknown - - -
        # calculate median coverage and median MAF of all confirmed present mutations in each sample
        self.sample_coverages = defaultdict(list)
        self.sample_mafs = defaultdict(list)
        vc_vars = dict()

        putative_sequencing_artifacts = 0
        low_vaf_artifacts = 0
        # #################### APPLY FILTERS ###################
        self.variant_stats = Counter()
        for mut_key in list(self.mut_reads.keys()):

            # check for minimum VAF in at least on of the samples
            if all(float(self.mut_reads[mut_key][sample_name]) / self.coverage[mut_key][sample_name] <
                   settings.MIN_VAF for sample_name in self.mut_reads[mut_key].keys()
                   if self.mut_reads[mut_key][sample_name] >= max(1, settings.MIN_VAR_READS)):

                vafs = [float(self.mut_reads[mut_key][sample_name]) / self.coverage[mut_key][sample_name]
                        for sample_name in self.mut_reads[mut_key].keys()
                        if self.mut_reads[mut_key][sample_name] >= max(1, settings.MIN_VAR_READS)]
                if len(vafs) == 0:
                    logger.debug('Excluded variant {} ({}) as it has in no sample at least {} variant reads.'
                                 .format(gene_names[mut_key], mut_key, max(1, settings.MIN_VAR_READS)))
                else:
                    logger.debug('Excluded variant {} ({}) present with highest VAF of {:.1%}.'
                                 .format(gene_names[mut_key], mut_key, max(vafs)))
                low_vaf_artifacts += 1

                # exclude these variants
                del self.mut_reads[mut_key]
                del self.coverage[mut_key]
                del gene_names[mut_key]
                self.variant_stats[-1] += 1
                continue

            # check if variant is present in the normal sample
            if normal_sample is not None:
                if (norm_var[mut_key] >= mut_reads_normal_th
                        and float(norm_var[mut_key]) / norm_cov[mut_key] >= vaf_normal_th):
                    logger.warning('Excluded variant {} ({}) present at a VAF of {:.1%} ({}/{}) in the normal sample.'
                                   .format(gene_names[mut_key], mut_key, float(norm_var[mut_key]) / norm_cov[mut_key],
                                           norm_var[mut_key], norm_cov[mut_key]))
                    putative_sequencing_artifacts += 1
                    # exclude these variants
                    del self.mut_reads[mut_key]
                    del self.coverage[mut_key]
                    del gene_names[mut_key]
                    self.variant_stats[-6] += 1
                    continue

            # remove intronic and intergenic variants due to frequent artifacts in whole exome sequencing data
            if self.ensembl_data is not None:
                chrom, start_pos, end_pos, ref, alt = get_variant_details(mut_key)

                # varcode variant, see https://github.com/hammerlab/varcode
                variant = VCVariant(contig=chrom, start=start_pos, ref=ref, alt=alt, ensembl=self.ensembl_data)

                if wes_filtering:
                    # remove intronic variants
                    if is_intronic(variant):
                        logger.debug('Excluded intronic variant {} ({}).'.format(gene_names[mut_key], mut_key))
                        # exclude this variant
                        del self.mut_reads[mut_key]
                        del self.coverage[mut_key]
                        del gene_names[mut_key]
                        self.variant_stats[-2] += 1
                        continue

                    # remove intergenic variants
                    elif is_intergenic(variant):
                        logger.debug('Excluded intergenic variant {} ({}).'.format(gene_names[mut_key], mut_key))
                        # exclude this variant
                        del self.mut_reads[mut_key]
                        del self.coverage[mut_key]
                        del gene_names[mut_key]
                        self.variant_stats[-3] += 1
                        continue

                    # remove variants with incomplete transcript annotation (likely introns)
                    elif is_incompletetranscript(variant):
                        logger.debug('Excluded variant with incomplete transcript annotation {} ({}).'.format(
                            gene_names[mut_key], mut_key))
                        # exclude this variant
                        del self.mut_reads[mut_key]
                        del self.coverage[mut_key]
                        del gene_names[mut_key]
                        self.variant_stats[-4] += 1
                        continue

                if mut_key in artifacts.keys():
                    logger.debug('Excluded variant {} ({}) from analysis as it was found in the common variant list.'
                                 .format(mut_key, artifacts[mut_key][1]))
                    # exclude this variant
                    del self.mut_reads[mut_key]
                    del self.coverage[mut_key]
                    del gene_names[mut_key]
                    self.variant_stats[-5] += 1
                    continue

                vc_vars[mut_key] = variant

            # variant passed all filters
            self.variant_stats[1] += 1
            for sample_name in self.mut_reads[mut_key].keys():

                if self.coverage[mut_key][sample_name] >= 0:
                    self.sample_coverages[sample_name].append(self.coverage[mut_key][sample_name])

                if self.mut_reads[mut_key][sample_name] > 2:
                    maf = float(self.mut_reads[mut_key][sample_name]) / self.coverage[mut_key][sample_name]
                    if maf > 0.01:      # ensure it's not due to sequencing errors
                        self.sample_mafs[sample_name].append(maf)

        if putative_sequencing_artifacts > 0:
            logger.warning('{} variants were detected as significantly present in the normal sample.'.format(
                putative_sequencing_artifacts))
        if low_vaf_artifacts > 0:
            logger.warning('{} variants did not reach a VAF of {:.1%} and at least {} var reads in any of the samples.'
                           .format(low_vaf_artifacts, settings.MIN_VAF, max(1, settings.MIN_VAR_READS)))

        self._print_filter_statistics()
        if wes_filtering and self.ensembl_data is None:
            logger.warning('Filtering of intronic and intergenic variants could not be performed because VarCode'
                           'is not available!')

        logger.info('{} variants passed the filtering.'.format(len(gene_names)))

        # calculate p-values for presence and absence for all given variants
        self.present_p_values = calculate_p_values(self.mut_reads, self.coverage, false_positive_rate)

        # remove low quality samples
        self.discarded_samples = self._filter_samples(min_sa_cov, min_sa_maf)
        # samples having passed the filtering
        self.n = len(self.sample_names)
        if self.n == 0:
            raise RuntimeError('No sample passed the filtering.')
        logger.info('{} samples passed filtering. {} samples have been discarded ({}).'.format(
            self.n, len(self.discarded_samples), ', '.join(self.discarded_samples)))

        # merge all p-values of a patient into one directory
        merged_p_values = dict()
        for sample_name in self.sample_names:
            for mut_key, p_value in self.present_p_values[sample_name].items():
                merged_p_values[(sample_name, mut_key)] = p_value

        # find significantly mutated genes using the Benjamini Hochberg procedure
        # returns a set of sample_idx - mutation key tuples
        sig_muts = find_significant_mutations(merged_p_values, false_discovery_rate)

        # merge mutant reads, phred coverage and distinct phred coverage
        # set MAF and determine positives (>0), negatives (0),
        # and unknowns (-1: likely positive; -2: likely negative)
        self.positives = Counter()
        self.unknowns = [Counter() for _ in range(2)]
        self.negatives = Counter()

        self.mut_keys = []
        self.gene_names = []
        if self.ensembl_data is not None:
            self.vc_variants = []
            self.mut_types = []

        # posterior log probability if no data was reported
        non_log_p0 = math.log(def_sets.NO_DATA_P0)
        non_log_p1 = math.log(1.0 - def_sets.NO_DATA_P0)

        self._calculate_hyperparameters()

        self.vafs = np.zeros((len(gene_names), len(self.sample_names)))

        # ##################################################################################
        # - - - - - - - - CLASSIFY MUTATIONS with BAYESIAN INFERENCE MODEL - - - - - - - - -
        # ##################################################################################
        for mut_key, gene_name in sorted(gene_names.items(), key=lambda k: k[1].lower()):

            chrom, start_pos, end_pos, ref, alt = get_variant_details(mut_key)

            if self.ensembl_data is not None:
                # varcode variants were previously generated to check for possible filtering
                self.vc_variants.append(vc_vars[mut_key])
                self.mut_types.append(get_top_effect_name(vc_vars[mut_key]))

            # add mutation name
            self.mut_keys.append(mut_key)
            # add gene name
            self.gene_names.append(gene_names[mut_key])
            # position data of this variant
            self.mut_positions.append((chrom, start_pos, end_pos))

            # - - - - - - - - CLASSIFY MUTATIONS - - - - - - - - -
            for sa_id, sample_name in enumerate(self.sample_names):

                # add VAF
                if self.coverage[mut_key][sample_name] > 0:
                    self.vafs[len(self.mut_keys)-1, sa_id] = (float(self.mut_reads[mut_key][sample_name]) /
                                                              self.coverage[mut_key][sample_name])
                else:
                    self.vafs[len(self.mut_keys) - 1, sa_id] = 0.0

                # calculate posterior: log probability that VAF = 0
                if self.coverage[mut_key][sample_name] < 0:   # no sequencing data in this sample
                    self.log_p01[len(self.mut_keys)-1].append([non_log_p0, non_log_p1])

                else:                          # calculate posterior according to prior, estimated purity and data

                    # calculate posterior: log probability that VAF = 0, log probability that VAF > 0
                    p0, p1 = get_log_p0(self.coverage[mut_key][sample_name], self.mut_reads[mut_key][sample_name],
                                        self.bi_error_rate, self.bi_c0, cutoff_f=self.get_cutoff_frequency(sample_name),
                                        pseudo_alpha=def_sets.PSEUDO_ALPHA, pseudo_beta=self.betas[sample_name])
                    self.log_p01[len(self.mut_keys)-1].append([p0, p1])
                # logger.debug('p0: {;.2e}, k: {}, n: {}.'.format(
                #     self.log_p0[len(self.mut_keys)-1], self.mut_reads[mut_key][sample_name],
                #     self.phred_coverage[mut_key][sample_name]))

                # conventional binary present/absent classification
                # not used in inference model, just for artifact calculations

                # logger.debug('Confidence p-value values {:.2e}, {:.2e} for {} mut-reads at {}x coverage'.format(
                #     self.present_p_values[sample_name][mut_key], self.absent_p_values[sample_name][mut_key],
                #     self.mut_reads[mut_key][sample_name], self.phred_coverage[mut_key][sample_name]))

                # was the null hypothesis rejected (declared as significant)?
                if (sample_name, mut_key) in sig_muts:

                    maf = float(self.mut_reads[mut_key][sample_name]) / self.coverage[mut_key][sample_name]

                    self.data[len(self.mut_keys)-1].append(maf)
                    self.positives[sample_name] += 1

                # is there enough coverage supporting a conclusion
                elif min_absent_cov == 0 or self.coverage[mut_key][sample_name] >= min_absent_cov \
                        or self.coverage[mut_key][sample_name] == -1:     # coverage has not been reported

                    self.data[len(self.mut_keys)-1].append(0)
                    self.negatives[sample_name] += 1

                # not enough coverage at this position => unknown
                else:
                    if self.mut_reads[mut_key][sample_name] > 0:        # unknown present
                        self.data[len(self.mut_keys)-1].append(POS_UNKNOWN)
                        self.unknowns[0][sample_name] += 1
                    else:                                               # unknown absent
                        self.data[len(self.mut_keys)-1].append(NEG_UNKNOWN)
                        self.unknowns[1][sample_name] += 1

        for sample_name in self.sample_names:
            logger.debug('Sample {} conventional classifications: '.format(sample_name) +
                         '{} positives; {} negatives; {} unknowns;'.format(
                         self.positives[sample_name], self.negatives[sample_name],
                         self.unknowns[0][sample_name]+self.unknowns[1][sample_name]))

        # for sample_name in self.sample_names:
        #     logger.debug('Sample {} classifications: '.format(sample_name)
        #                  + '{} positives; {} negatives; {} positive unknowns, {} negative unknowns;'.format(
        #                  self.positives[sample_name], self.negatives[sample_name],
        #                  self.unknowns[0][sample_name], self.unknowns[1][sample_name]))

        logger.info("{} samples passed the filtering and have been processed. ".format(self.n))
        logger.info("Completed reading of allele frequencies at {} mutated positions. ".format(len(self.mut_keys)))

        return len(self.discarded_samples) + self.n

    def read_maf_data(self, maf_filename):
        """
        Read allele frequencies for all found mutations in the samples
        of the given MAF file
        :param maf_filename: path to the MAF file
        """

        with open(maf_filename, 'r') as fData:

            for line in fData:
                if line[0] == '#':
                    continue

                # header
                elif line.startswith("Position"):
                    self.sample_names = re_sample.findall(line)
                    logger.debug("Read sample ids: {}".format(self.sample_names))
                    col_chr = col_data_start = col_data_end = col_function = col_gene = col_pathway = -1
                    # Find columns with the sample data
                    for col_idx, entry in enumerate(line.split('\t')):
                        if entry == self.sample_names[0]:
                            col_data_start = col_idx
                        elif entry == self.sample_names[-1]:
                            col_data_end = col_idx
                        elif entry == 'Chromosome':
                            col_chr = col_idx
                        elif entry == 'Start':
                            col_start_pos = col_idx
                        elif entry == 'End':
                            col_end_pos = col_idx
                        elif entry == 'Gene':
                            col_gene = col_idx
                        elif entry == 'Driver_Pathway' or entry == 'Driver':
                            col_pathway = col_idx
                        elif entry == 'Function':
                            col_function = col_idx
                            self.mut_types = []

                    # check if the columns with needed information were found
                    if col_function == -1:
                        logger.info('Data file {} does not contain information '.format(maf_filename) +
                                    'about the function of the mutated region. ')
                    if col_chr == -1 or col_start_pos == -1 or col_end_pos == -1:
                        logger.info('Data file {} does not contain information '.format(maf_filename) +
                                    'about the exact mutation position on the chromosome. ')

                # ignore line
                elif line.startswith("#") or line == '\n':
                    continue
                # process line with mutation allele frequencies
                else:

                    afrs = (line[:-1].split("\t"))[:col_data_end+1]

                    # check if the mutation is present in one of the samples
                    # if not, ignore this mutation
                    if any(float(freq) for freq in afrs[col_data_start:]):

                        self.mut_keys.append(afrs[0])      # add detailed mutation information
                        if col_gene != -1:
                            # add gene name of the mutation
                            self.gene_names.append(afrs[col_gene].replace(',', ', '))

                        if col_function != -1:      # check if function information is in the data file
                            self.mut_types.append(afrs[col_function])
                        else:
                            self.mut_types.append('')

                        # start position data of this mutation
                        if col_chr != -1 and col_start_pos != -1 and col_end_pos != -1:
                            self.mut_positions.append((afrs[col_chr], afrs[col_start_pos], afrs[col_end_pos]))
                        else:
                            self.mut_positions.append((0, 0, 0))

                        if col_pathway != -1 and afrs[col_pathway] != '.':                  # driver pathway information
                            self.driver_pathways[len(self.mut_keys)-1] = afrs[col_pathway]

                        # process MAF of the current mutation for each sample
                        for sa_idx in range(col_data_start, len(afrs)):

                            if float(afrs[sa_idx]) > 1:
                                logger.warning('Found allele frequency of {} indicating different file format. '.format(
                                    float(afrs[sa_idx])))
                                if sa_idx + 1 == len(afrs):
                                    logger.warning('Column {} is not considered as a sample column from now on'.format(
                                        self.sample_names[sa_idx-col_data_start]))
                                    col_data_end -= 1
                                    del self.sample_names[sa_idx-col_data_start]
                                else:
                                    raise ValueError('Cannot read input file {}'.format(maf_filename))
                            else:
                                self.data[len(self.mut_keys)-1].append(float(afrs[sa_idx]))

                    # otherwise move to the next mutation
                    else:
                        logger.debug('Mutation {} is not present in any of the samples and hence ignored.'
                                     .format(afrs[0]))

            self.n = len(self.sample_names)
            logger.info("{} samples have been processed. ".format(self.n))

            logger.info("Completed reading of allele frequencies at {} mutated positions. ".format(len(self.mut_keys)))

    def read_vcf_directory(self, vcf_directory, min_sa_cov, min_sa_maf, false_positive_rate, false_discovery_rate,
                           min_absent_cov, normal_sample_name=None, excluded_samples=None, considered_samples=None,
                           wes_filtering=False, artifacts=None):
        """
        Read allele frequencies for all variants in the samples in the files of the given directory
        :param vcf_directory: directory with VCF files
        :param min_sa_cov: minimum median coverage per sample
        :param min_sa_maf: minimum median mutant allele frequency per sample
        :param false_positive_rate: false positive rate of the used sequencing technology
        :param normal_sample_name: do not consider given normal sample
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_absent_cov: minimum coverage at a non-significantly mutated position for absence classification
        :param excluded_samples: set of sample names to exclude
        :param considered_samples: if not None then only samples included in this set will be considered
        :param wes_filtering: remove intronic and intergenic variants due to frequent sequencing artifacts in whole
                              exome sequencing data
        :param artifacts: likely sequencing artifacts that are excluded from the analysis
        :return: number of samples which were processed (independent of filtering)
        """
        excl_samples = set()
        if normal_sample_name is not None:
            excl_samples.add(normal_sample_name)
        if excluded_samples is not None:
            for sa in excluded_samples:
                excl_samples.add(sa)
        samples = read_vcf_files(vcf_directory, excluded_samples=excl_samples, considered_samples=considered_samples)

        processed_samples = self._process_samples(samples, min_sa_cov, min_sa_maf, false_positive_rate,
                                                  false_discovery_rate, min_absent_cov, wes_filtering=wes_filtering,
                                                  artifacts=artifacts)

        return processed_samples

    def read_vcf_file(self, vcf_file, false_positive_rate, false_discovery_rate, min_sa_cov=0, min_sa_maf=0.0,
                      min_absent_cov=0, normal_sample_name=None, excluded_samples=None, considered_samples=None,
                      wes_filtering=False, artifacts=None):
        """
        Read allele frequencies for all variants in the samples in the given VCF file
        :param vcf_file: path to vcf file
        :param false_positive_rate: false positive rate of the used sequencing technology
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_sa_cov: minimum median coverage per sample (default 0), if below sample discarded
        :param min_sa_maf: minimum median mutant allele frequency per sample (default 0.0), if below sample discarded
        :param min_absent_cov: min coverage at a non-significantly mutated position for absence classification (def 0)
        :param normal_sample_name: name of the normal sample to discard
        :param excluded_samples: set of sample names to exclude
        :param considered_samples: if not None then only samples included in this set will be considered
        :param wes_filtering: remove intronic and intergenic variants due to frequent sequencing artifacts in whole
                              exome sequencing data
        :param artifacts: likely sequencing artifacts that are excluded from the analysis
        :return: number of processed samples
        """

        excl_samples = set()
        if normal_sample_name is not None:
            excl_samples.add(normal_sample_name)
        if excluded_samples is not None:
            for sa in excluded_samples:
                excl_samples.add(sa)
        samples = read_vcf_file(vcf_file, excluded_samples=excl_samples, considered_samples=considered_samples)
        processed_samples = self._process_samples(samples, min_sa_cov, min_sa_maf, false_positive_rate,
                                                  false_discovery_rate, min_absent_cov, wes_filtering=wes_filtering,
                                                  artifacts=artifacts)

        return processed_samples

    def _process_samples(self, samples, min_sa_cov, min_sa_maf, fpr, false_discovery_rate, min_absent_cov,
                         wes_filtering=False, artifacts=None):
        """
        Determine which variants are shared among the samples
        Remove low quality samples not meeting the median coverage or median VAF
        Find significantly mutated variants while controlling for the false discovery rate
        :param samples: list of relevant samples
        :param min_sa_cov: minimum median coverage per sample
        :param min_sa_maf: minimum median mutant allele frequency per sample
        :param fpr: false positive rate of the used sequencing technology
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_absent_cov: minimum coverage at a non-significantly mutated position for absence classification
        :param wes_filtering: remove intronic and intergenic variants due to frequent sequencing artifacts in whole
                              exome sequencing data
        :return: processed samples
        """

        # raw sequencing data
        self.mut_reads = defaultdict(dict)
        self.coverage = defaultdict(dict)

        # calculate median coverage and median MAF of all confirmed present mutations in each sample
        self.sample_coverages = defaultdict(list)
        self.sample_mafs = defaultdict(list)

        # Calculate p-values for confidence in presence and absence
        self.present_p_values = defaultdict(dict)

        heap = []
        tmp_vars = []

        self.mut_keys = []
        self.gene_names = []
        if self.ensembl_data is not None:
            self.vc_variants = []
            self.mut_types = []

        self.variant_stats = Counter()

        if artifacts is None:
            # no potential dictionary of sequencing artifacts was given
            artifacts = dict()

        # take first variant from all samples
        for sample_name, sample in sorted(samples.items(), key=lambda x: x[0]):
            heapq.heappush(heap, (heapq.heappop(sample.variants), sample_name))
            self.sample_names.append(sample.name)

        logger.info('Start processing of {} samples.'.format(len(self.sample_names)))

        while len(heap) > 0:

            variant, sample_name = heapq.heappop(heap)
            # check if heap with variants of this sample is not yet empty
            if samples[sample_name].variants:
                heapq.heappush(heap, (heapq.heappop(samples[sample_name].variants), sample_name))

            if len(tmp_vars) > 0:
                if tmp_vars[0][0].CHROM == variant.CHROM \
                        and tmp_vars[0][0].POS == variant.POS \
                        and tmp_vars[0][0].REF == variant.REF \
                        and tmp_vars[0][0].ALT[0] == variant.ALT[0]:

                    # add new identical variant
                    tmp_vars.append((variant, sample_name))

                    if len(heap) == 0:
                        # process identical variants and remove them from the temporal list
                        status = self._add_variant(tmp_vars, fpr, artifacts, wes_filtering=wes_filtering)
                        self.variant_stats[status] += 1

                else:       # variant at different position
                    # process old identical variants and remove them from the temporal list
                    status = self._add_variant(tmp_vars, fpr, artifacts, wes_filtering=wes_filtering)
                    self.variant_stats[status] += 1

                    # add the new variant to the list to process
                    tmp_vars.append((variant, sample_name))

            else:
                # add new variant
                tmp_vars.append((variant, sample_name))

        logger.info('{} samples have been processed. '.format(len(samples.keys())))
        logger.info('Completed reading of allele frequencies at {} mutated positions. '.format(len(self.mut_keys)))

        self._print_filter_statistics()

        if wes_filtering and self.ensembl_data is None:
            logger.warning('Filtering of intronic and intergenic variants could not be performed because VarCode'
                           'is not available!')

        self.discarded_samples = self._filter_samples(min_sa_cov, min_sa_maf)
        for sample_name in self.discarded_samples:
            samples.pop(sample_name)

        # samples having passed the filtering
        self.n = len(self.sample_names)
        if self.n == 0:
            raise RuntimeError('No sample passed the filtering.')
        logger.info('{} samples passed filtering. {} samples have been discarded ({}).'.format(
            self.n, len(self.discarded_samples), ', '.join(self.discarded_samples)))

        # merge all p-values of a patient into one directory
        merged_p_values = dict()
        for sample_name in self.sample_names:
            for mut_key, p_value in self.present_p_values[sample_name].items():
                merged_p_values[(sample_name, mut_key)] = p_value

        # find significantly mutated genes using the Benjamini Hochberg procedure
        # returns a set of sample_idx - mutation key tuples
        sig_muts = find_significant_mutations(merged_p_values, false_discovery_rate)

        # merge mutant reads, phred coverage and distinct phred coverage
        # set MAF and determine positives (>0), negatives (0),
        # and unknowns (-1: likely positive; -2: likely negative)
        self.positives = Counter()
        self.unknowns = [Counter() for _ in range(2)]
        self.negatives = Counter()

        # posterior log probability if no data was reported
        non_log_p0 = math.log(def_sets.NO_DATA_P0)
        non_log_p1 = math.log(1.0 - def_sets.NO_DATA_P0)

        # analyze VAF distribution per sample and calculate a sample-specific prior
        # based on an estimated purity of shared mutations
        self._calculate_hyperparameters()

        self.vafs = np.zeros((len(self.mut_keys), len(self.sample_names)))

        # ##################################################################################
        # - - - - - - - - CLASSIFY MUTATIONS with BAYESION INFERENCE MODEL - - - - - - - - -
        # ##################################################################################
        for mut_id, mut_key in enumerate(self.mut_keys):

            for sa_id, sample_name in enumerate(self.sample_names):

                # add VAF
                if self.coverage[mut_key][sample_name] > 0:
                    self.vafs[mut_id, sa_id] = (float(self.mut_reads[mut_key][sample_name]) /
                                                self.coverage[mut_key][sample_name])
                else:
                    self.vafs[mut_id, sa_id] = 0.0

                # calculate posterior: log probability that VAF = 0
                if self.coverage[mut_key][sample_name] < 0 or math.isnan(self.coverage[mut_key][sample_name]):
                    # no sequencing data in this sample
                    self.log_p01[mut_id].append([non_log_p0, non_log_p1])

                else:                        # calculate posterior according to prior, estimated purity and data

                    p0, p1 = get_log_p0(self.coverage[mut_key][sample_name], self.mut_reads[mut_key][sample_name],
                                        self.bi_error_rate, self.bi_c0, cutoff_f=self.get_cutoff_frequency(sample_name),
                                        pseudo_alpha=def_sets.PSEUDO_ALPHA, pseudo_beta=self.betas[sample_name])
                    self.log_p01[mut_id].append([p0, p1])

                # logger.debug('p0: {:.2e}, k: {}, n: {}.'.format(
                #     self.log_p0[mut_id][len(self.log_p0[mut_id])-1], self.mut_reads[mut_key][sample_name],
                #     self.phred_coverage[mut_key][sample_name]))

                # conventional binary present/absent classification
                # not used in inference model, just for artifact calculations
                # logger.debug('Confidence p-value values {:.2e}, {:.2e} for {} mut-reads at {}x coverage'.format(
                #     self.present_p_values[sample_name][mut_key], self.absent_p_values[sample_name][mut_key],
                #     self.mut_reads[mut_key][sample_name], self.phred_coverage[mut_key][sample_name]))

                # was the null hypothesis rejected (declared as significant)?
                if (sample_name, mut_key) in sig_muts:

                    maf = float(self.mut_reads[mut_key][sample_name]) / self.coverage[mut_key][sample_name]

                    self.data[mut_id].append(maf)
                    self.positives[sample_name] += 1

                # is there enough coverage supporting a conclusion
                elif min_absent_cov == 0 or sample_name not in self.coverage[mut_key] \
                        or self.coverage[mut_key][sample_name] >= min_absent_cov:

                    self.data[mut_id].append(0)
                    self.negatives[sample_name] += 1

                # not enough coverage at this position => unknown
                else:
                    if self.mut_reads[mut_key][sample_name] > 0:        # unknown present
                        self.data[mut_id].append(POS_UNKNOWN)
                        self.unknowns[0][sample_name] += 1
                    else:                                               # unknown absent
                        self.data[mut_id].append(NEG_UNKNOWN)
                        self.unknowns[1][sample_name] += 1

        for sample_name in self.sample_names:
            logger.debug('Sample {} classifications: '.format(sample_name) +
                         '{} positives; {} negatives; {} unknowns;'.format(
                         self.positives[sample_name], self.negatives[sample_name],
                         self.unknowns[0][sample_name]+self.unknowns[1][sample_name]))

        for sample_name in self.sample_names:
            logger.debug('Sample {} classifications: '.format(sample_name) +
                         '{} positives; {} negatives; {} positive unknowns, {} negative unknowns;'.format(
                         self.positives[sample_name], self.negatives[sample_name],
                         self.unknowns[0][sample_name], self.unknowns[1][sample_name]))

        # calculate median coverage and median allele frequency across all samples
        coverages = []
        for mut_key in self.mut_reads.keys():
            for sample_name in self.mut_reads[mut_key].keys():
                if self.coverage[mut_key][sample_name] >= 0:
                    coverages.append(self.coverage[mut_key][sample_name])

        logger.info('Median coverage in patient {}: {} (mean: {:.1f})'.format(self.name, np.median(coverages),
                                                                              np.mean(coverages)))

        return len(self.discarded_samples) + self.n

    def _add_variant(self, tmp_vars, fpr, artifacts, wes_filtering=False):
        """
        Process identical variants and remove them from the given list
        :param tmp_vars: list of identical variants in different samples
        :param fpr: false positive rate of the used sequencing technology
        :param wes_filtering: remove intronic and intergenic variants due to frequent sequencing artifacts in whole
                              exome sequencing data
        :return 1 if variant passed all filters, -1 if variant never reached a significant level,
                -2 intronic, -3 intergenic, -4 incomplete transcript annotation
        """

        # process identical variants
        mut_key = '{}__{}__{}>{}'.format(tmp_vars[0][0].CHROM, tmp_vars[0][0].POS,
                                         tmp_vars[0][0].REF, tmp_vars[0][0].ALT[0])

        self.gene_names.append(tmp_vars[0][0].GENE_NAME if tmp_vars[0][0].GENE_NAME is not None else 'unknown')

        # generate varcode Variant class instance (useful to infer gene name and mutation effect)
        if self.ensembl_data is not None:

            # only for SNVs and small InDels reference allele make sense
            if tmp_vars[0][0].REF != 'N/A':
                # varcode variant, see https://github.com/hammerlab/varcode
                variant = VCVariant(contig=tmp_vars[0][0].CHROM, start=tmp_vars[0][0].POS, ref=tmp_vars[0][0].REF,
                                    alt=tmp_vars[0][0].ALT[0], ensembl=self.ensembl_data)

                # remove intronic and intergenic variants due to frequent artifacts in whole exome sequencing data
                if wes_filtering:
                    # remove intronic variants
                    if is_intronic(variant):
                        logger.debug('Excluded intronic variant in {} ({}).'.format(self.gene_names[-1], mut_key))
                        # exclude this variant
                        del self.gene_names[-1]
                        del tmp_vars[:]
                        return -2

                    # remove intergenic variants
                    elif is_intergenic(variant):
                        logger.debug('Excluded intergenic variant in {} ({}).'.format(self.gene_names[-1], mut_key))
                        # exclude this variant
                        del self.gene_names[-1]
                        del tmp_vars[:]
                        return -3

                    # remove variants with incomplete transcript annotation (likely introns)
                    elif is_incompletetranscript(variant):
                        logger.debug('Excluded variant with incomplete transcript annotation in {} ({}).'.format(
                            self.gene_names[-1], mut_key))
                        # exclude this variant
                        del self.gene_names[-1]
                        del tmp_vars[:]
                        return -4

                if mut_key in artifacts.keys():
                    logger.debug('Excluded variant {} ({}) from analysis as it was found in the provided artifact list.'
                                 .format(mut_key, artifacts[mut_key][1]))
                    # exclude this variant
                    del self.gene_names[-1]
                    del tmp_vars[:]
                    return -5

                # infer gene name if not given
                if self.gene_names[-1] == 'unknown':
                    # query missing data from Ensembl data
                    if len(variant.gene_names) == 1:
                        self.gene_names[-1] = variant.gene_names[0]
                    elif len(variant.gene_names) > 1:
                        self.gene_names[-1] = ','.join(g for g in variant.gene_names)
                        # else:
                        # No gene name could be inferred

                self.vc_variants.append(variant)
                self.mut_types.append(get_top_effect_name(variant))

            # for SVs and CNVs no mutation types can be inferred using VARCODE
            else:
                self.vc_variants.append(None)
                self.mut_types.append(None)

        self.mut_keys.append(mut_key)      # add detailed mutation information

        # start position data of this mutation
        self.mut_positions.append((tmp_vars[0][0].CHROM, tmp_vars[0][0].POS,
                                   str(int(tmp_vars[0][0].POS)+len(tmp_vars[0][0].REF))))

        tmp_vars.sort(key=lambda x: x[1], reverse=True)

        # calculate p-values and store coverages and MAFs of each sample
        var, var_sample_name = tmp_vars.pop()
        for sample_name in self.sample_names:

            if var_sample_name == sample_name:

                # set raw sequencing data
                self.mut_reads[mut_key][sample_name] = var.AD[1]
                self.coverage[mut_key][sample_name] = var.DP

                self.sample_coverages[sample_name].append(var.DP)
                # only consider variants with at least three supporting reads
                # for median VAF calculation
                if var.AD[1] > 2:
                    self.sample_mafs[sample_name].append(var.BAF)

                if var.DP > 0:
                    self.present_p_values[sample_name][mut_key] = calculate_present_pvalue(var.AD[1], var.DP, fpr)

                # get variant in the next called sample
                if len(tmp_vars) > 0:
                    var, var_sample_name = tmp_vars.pop()
                else:
                    var_sample_name = -1

            else:
                # sequencing data information was not provided in this sample for this variant
                self.mut_reads[mut_key][sample_name] = -1       # -1 = unknown which is different from 0
                self.coverage[mut_key][sample_name] = -1   # -1 = unknown which is different from 0

        # check for minimum VAF in at least on of the samples
        if all(float(self.mut_reads[mut_key][sample_name]) / self.coverage[mut_key][sample_name] <
                settings.MIN_VAF for sample_name in self.mut_reads[mut_key].keys()
               if self.mut_reads[mut_key][sample_name] >= max(1, settings.MIN_VAR_READS)):

            logger.debug('Variant {}{} did not pass filtering.'.format(
                mut_key, ' ({})'.format(self.gene_names[-1]) if self.gene_names is not None else ''))

            # exclude these variants
            del self.mut_reads[mut_key]
            del self.coverage[mut_key]
            del self.mut_keys[-1]
            del self.mut_positions[-1]
            del self.gene_names[-1]
            if self.ensembl_data is not None:
                del self.vc_variants[-1]
                del self.mut_types[-1]

            return -1

        else:
            # logger.debug('Variant {}{} passed filtering.'.format(
            #     mut_key, ' ({})'.format(self.gene_names[-1]) if self.gene_names is not None else ''))
            return 1

    def _filter_samples(self, min_sa_cov, min_sa_maf):

        # filter low median coverage and low median MAF samples
        discarded_samples = []
        self.sample_names = []

        for sample_name in sorted(self.sample_coverages.keys(),
                                  key=lambda item: (item.split('_')[0], int(item.split('_')[1]))
                                  if len(item.split('_')) > 1 and item.split('_')[1].isdigit() else (item, item)):

            # discard a sample if its median coverage is lower than the threshold
            if np.nanmedian(self.sample_coverages[sample_name]) < min_sa_cov:

                logger.warning('Sample {} is discarded! Median coverage of {} is below the threshold {}.'.format(
                    sample_name, np.nanmedian(self.sample_coverages[sample_name]),
                    min_sa_cov))

                self.sample_coverages.pop(sample_name)
                if sample_name in self.sample_mafs.keys():
                    self.sample_mafs.pop(sample_name)

                self.present_p_values.pop(sample_name)

                discarded_samples.append(sample_name)

            # discard a sample if its median MAF is lower than the threshold
            elif np.nanmedian(self.sample_mafs[sample_name]) < min_sa_maf:

                logger.warning('Sample {} is discarded!'.format(sample_name) +
                               ' Median mutant allele frequency (MAF) of {:.3f} is below the threshold {}.'.format(
                               np.nanmedian(self.sample_mafs[sample_name]), min_sa_maf))

                self.sample_coverages.pop(sample_name)
                self.sample_mafs.pop(sample_name)

                self.present_p_values.pop(sample_name)

                discarded_samples.append(sample_name)

            # sample passed filtering
            else:
                self.sample_names.append(sample_name)
                logger.info('Sample {}: median coverage {:.1f}, median VAF {:.3f}.'.format(
                    sample_name, np.nanmedian(self.sample_coverages[sample_name]) if
                    len(self.sample_coverages[sample_name]) > 0 else float('nan'),
                    np.nanmedian(self.sample_mafs[sample_name]) if len(self.sample_mafs[sample_name]) > 0
                    else float('nan')))
                # logger.info('Median distinct phred coverage in sample {}: {}'.format(
                #     sample_name, np.median(self.sample_dis_phred_coverages[sample_name])))

        return discarded_samples

    def _print_filter_statistics(self):
        """
        Print statistics of the applied filters
        """
        logger.info('{} variants passed filtering.'.format(self.variant_stats[1]))
        if self.variant_stats[-1] > 0:
            logger.info('{} variants did not reach a significant level.'.format(self.variant_stats[-1]))
        if self.variant_stats[-2] > 0:
            logger.info('{} variants were intronic and excluded.'.format(self.variant_stats[-2]))
        if self.variant_stats[-3] > 0:
            logger.info('{} variants were intergenic and excluded.'.format(self.variant_stats[-3]))
        if self.variant_stats[-4] > 0:
            logger.info('{} variants were excluded due to incomplete transcript annotation (likely introns).'.format(
                self.variant_stats[-4]))
        if self.variant_stats[-5] > 0:
            logger.info('{} variants were in provided artifact list and excluded.'.format(self.variant_stats[-5]))

    def get_cutoff_frequency(self, sample_name):
        """
        Calculate cutoff frequency (max absent frequency given an estimated purity)
        Cutoff frequency lower bound is 1%
        :param sample_name:
        :return:
        """

        # what if purity estimation is not available? use mean of other samples, otherwise use 80%
        if sample_name in self.purities.keys():
            cutoff_f = max(self.max_absent_vaf * self.purities[sample_name], 0.01,
                           self.bi_error_rate)

        # no purity was estimated: use 80% for each
        elif np.isnan(np.median([p for p in self.purities.values()])):
            cutoff_f = max(self.max_absent_vaf * 0.8, 0.01, self.bi_error_rate)

        # purity was estimated for some samples: use that one
        else:
            cutoff_f = max(self.max_absent_vaf * np.median([p for p in self.purities.values()]),
                           0.01, self.bi_error_rate)

        return cutoff_f

    def _calculate_hyperparameters(self):
        """
        Compute hyperparameters for the prior in the Bayesian inference model based on the estimated purities
        of each sample
        """

        # get purities
        if self.purities is None:
            self._estimate_purities()

        # set hyperparameter beta according to the estimated purities
        self.betas = dict()
        for sample_name in self.sample_names:
            if sample_name in self.purities:
                self.betas[sample_name] = 1.0 / self.purities[sample_name]
                logger.debug('Beta for prior in sample {}: {:.1f}'.format(sample_name, self.betas[sample_name]))
            else:
                self.betas[sample_name] = def_sets.PSEUDO_BETA
                logger.warning('Purity could not be estimated. Used default beta for prior in sample {}: {:.1f}'.format(
                               sample_name, self.betas[sample_name]))

    def _estimate_purities(self):
        """
        Estimate purity from shared (none-private) variants present in multiple samples and
        take their median VAF in a given sample; assumes diploid cancer cells
        """

        pres_lp = math.log(0.5)
        clonal_vafs = defaultdict(lambda: defaultdict(list))
        for mut_key in self.mut_reads.keys():
            present_samples = list()        # samples where variant is likely to be present
            for sample_name in self.mut_reads[mut_key].keys():
                if self.mut_reads[mut_key][sample_name] > 2:

                    _, p1 = get_log_p0(self.coverage[mut_key][sample_name], self.mut_reads[mut_key][sample_name],
                                       self.bi_error_rate, self.bi_c0, cutoff_f=0.05,
                                       pseudo_alpha=def_sets.PSEUDO_ALPHA, pseudo_beta=def_sets.PSEUDO_BETA)

                    if p1 > pres_lp:        # probability to be present is greater than 50%
                        present_samples.append(sample_name)

            if len(present_samples) >= max(2.0, len(self.sample_names)/3.0):    # variant is present in multiple samples
                # hence, not a private mutation and therefore helpful to estimate purity
                for sample_name in present_samples:
                    clonal_vafs[sample_name][len(present_samples)].append(
                        float(self.mut_reads[mut_key][sample_name]) / self.coverage[mut_key][sample_name])

        self.purities = dict()
        for sample_name in self.sample_names:
            if sample_name in clonal_vafs.keys() and \
                            sum(len(clonal_vafs[sample_name][k]) for k in clonal_vafs[sample_name].keys()) > 5:

                shared_vafs = list()
                for no_shared, vafs in sorted(clonal_vafs[sample_name].items(), key=lambda k: -k[0]):
                    shared_vafs += vafs
                    if len(shared_vafs) > max(min(10, len(self.mut_reads)*0.3), len(self.mut_reads)*0.02):
                        break
                self.purities[sample_name] = 2 * np.median(shared_vafs)
                logger.info('Identified {} shared variants in sample {}. Estimated purity: {:.1%}.'.format(
                    len(shared_vafs), sample_name, self.purities[sample_name]))
                if self.purities[sample_name] > 0.995:
                    logger.warning('Sample {} has unusual high estimated purity of {:.1%}.'.format(
                        sample_name, self.purities[sample_name]))
                    self.purities[sample_name] = 0.99

            elif 0.05 < np.median(self.sample_mafs[sample_name]) < 0.5:
                self.purities[sample_name] = 2 * np.median(self.sample_mafs[sample_name])
                logger.warning('Insufficient shared variants in '
                               'sample {} to reliably estimate purity. Used median VAF for estimation: {:.1%}'.format(
                                sample_name, self.purities[sample_name]))
            else:
                logger.warning('Only {} private mutations with median VAF of {:.1%} identified in sample {}. '.format(
                    len(self.sample_mafs[sample_name]), np.median(self.sample_mafs[sample_name]),
                    sample_name) + 'Unable to estimate purity!')

    def calculate_no_present_vars(self):
        """
        Calculate number of present and absent variants according to the Bayesian inference model for each sample
        """

        self.no_present_vars = Counter()
        self.no_absent_vars = Counter()

        pres_lp = math.log(0.5)  # log probability of 50%

        for ps in self.log_p01.values():
            for sa_idx, (_, log_p1) in enumerate(ps):
                if log_p1 > pres_lp:
                    self.no_present_vars[sa_idx] += 1
                else:
                    self.no_absent_vars[sa_idx] += 1


def get_variant_details(mut_key):
    """
    Extract the details of the variant from the mutation key
    :param mut_key:
    :return: tuple of chromosome, start position, end position, reference allele, alternative allel
    """
    # determine exact position of variant
    var = mut_key.split('__')
    if var[0].lower().startswith('chr'):
        if max(var[0].find('p'), var[0].find('q')) > -1:
            chrom = var[0][3:max(var[0].find('p'), var[0].find('q'))]
        else:
            chrom = var[0][3:]
    else:
        chrom = var[0]  # chromosome
    start_pos = int(var[1])
    if '>' in var[2]:
        ref, alt = var[2].split('>')

    # insertion/duplication
    elif 'dup' in var[2] or 'ins' in var[2]:
        ref = '-'
        alt = var[2][var[2].find('dup') + 3:]

    # deletion
    elif 'del' in var[2]:
        alt = '-'
        ref = var[2][var[2].find('del') + 3:]
    else:
        ref = alt = var[2]
        logger.warning('Change format not supported: {}'.format(var[2]))

    end_pos = start_pos + (len(ref) if ref != '-' else 0) - 1

    return chrom, start_pos, end_pos, ref, alt
