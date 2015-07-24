#!/usr/bin/python
"""Data structure around sequencing data of a subject"""
__author__ = 'jreiter'
__date__ = 'April, 2014'

import logging
from collections import defaultdict, Counter
from itertools import islice
import re
import heapq
import numpy as np
import math
import settings
import utils.int_settings as int_sets
from utils.int_settings import NEG_UNKNOWN, POS_UNKNOWN
from utils.vcf_parser import read_vcf_files
from utils.data_tables import read_mutation_table
from utils.statistics import calculate_present_pvalue, find_significant_mutations
from utils.vaf_data import get_shared_mutations, get_present_mutations, calculate_p_values
from utils.statistics import get_log_p0

# get logger for application
logger = logging.getLogger('treeomics')
re_sample = re.compile(r'(?:Pam[0-9]{2})?PT[0-9]{1,2}[B]?')


class Patient(object):
    """ Patient: sample data processing
    """

    def __init__(self, pat_name='Patient', min_absent_cov=0, error_rate=0.01, c0=0.5):
        """
        Constructor
        """

        self.name = pat_name

        # minimum coverage for an absent variant
        self.min_absent_cov = min_absent_cov

        # allele frequency data in each numbered sample
        self.data = defaultdict(list)

        # tuple posterior: log probability that VAF = 0, lob probability that VAF > 0
        self.log_p01 = defaultdict(list)
        self.bi_error_rate = error_rate     # sequencing error rate for bayesian inference
        self.bi_c0 = c0    # prior mixture parameter of delta function and uniform distribution for bayesian inference
        logger.info('Bayesian inference model uses an error rate of {} and p0 prior weight of {}.'.format(
            self.bi_error_rate, self.bi_c0))

        # raw sequencing data of the variants
        self.mut_reads = None
        self.phred_coverage = None
        self.dis_phred_coverage = None
        # median coverages per sample and median MAFs of all confirmed present mutations per sample
        self.sample_dis_phred_coverages = None
        self.sample_phred_coverages = None
        self.sample_mafs = None

        self.positives = None
        self.unknowns = None
        self.negatives = None

        self.present_p_values = None
        self.absent_p_values = None

        # sample contain which numbered mutations
        self.samples = defaultdict(set)
        self.discarded_samples = 0

        # mutation is present in which numbered samples
        self.mutations = defaultdict(set)
        # list of mutations which are present in some of the samples (excluding unknowns)
        self.present_mutations = None
        # list of shared mutations, excluding unknowns
        self.shared_mutations = None

        self.n = 0      # read samples
        # holds the names for the numbered samples
        self.sample_names = None
        # holds the names of the numbered mutations
        self.mut_keys = None
        # holds the gene names of the numbered mutations
        self.gene_names = None
        # holds the function of the numbered mutation, e.g. intronic, extronic, etc.
        self.mut_functions = []
        # holds the data about the mutation: chromosome, start position, and end position
        self.mut_positions = []
        # dictionary with the name of the driver pathway if the mutation is a driver gene mutation
        self.driver_pathways = dict()

        # mutations present in all samples
        self.founders = set()
        # mutations ordered via the numbered of shared samples
        self.shared_muts = defaultdict(set)

        # holds a dictionary with a frozenset of samples mapping to a set of mutations
        # present in exactly the same set
        self.mps = None         # mutation patterns

        self.subclones = None
        self.sc_names = None
        self.updated_clones = None

        self.common_muts = []
        self.add_muts = []

        self.gen_dis = None         # genetic distance between any pair of samples
        self.sim_coff = None        # Jaccard similarity coefficient between any pair of samples

    def read_raw_data(self, read_table, cov_table, dis_cov_table, false_positive_rate, false_discovery_rate,
                      min_absent_cov, min_sa_cov, min_sa_maf, excluded_columns=set()):
        """
        Read raw sequencing data from tsv files
        :param read_table: path to file with mutant reads
        :param cov_table: path to file with phred coverage
        :param dis_cov_table: path to file with distinct phred coverage
        :param false_positive_rate: false positive read of the used sequencing technology
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_absent_cov: minimum coverage at a non-significantly mutated position for absence classification
        :param min_sa_cov: minimum median coverage per sample
        :param min_sa_maf: minimum median mutant allele frequency per sample
        :param excluded_columns: matched normal sample or other samples to be excluded from the analysis
        """

        # read sequencing data from tsv files
        self.mut_reads, gene_names = read_mutation_table(read_table, excluded_columns)
        self.phred_coverage, _ = read_mutation_table(cov_table, excluded_columns)
        if dis_cov_table is not None:
            self.dis_phred_coverage, _ = read_mutation_table(dis_cov_table, excluded_columns)
        else:
            self.dis_phred_coverage = self.phred_coverage

        assert len(self.mut_reads.keys()) == len(self.phred_coverage), \
            'Number of reported variants is different in the data files.'

        # - - - classifier for positives, negatives, positive unknowns and negative unknown - - -
        # calculate median coverage and median MAF of all confirmed present mutations in each sample
        self.sample_dis_phred_coverages = defaultdict(list)
        self.sample_phred_coverages = defaultdict(list)
        self.sample_mafs = defaultdict(list)

        for mut_key in self.mut_reads.keys():
            for sample_name in self.mut_reads[mut_key].keys():

                if self.phred_coverage[mut_key][sample_name] >= 0:
                    self.sample_phred_coverages[sample_name].append(self.phred_coverage[mut_key][sample_name])

                if self.mut_reads[mut_key][sample_name] > 0:
                    maf = float(self.mut_reads[mut_key][sample_name]) / self.phred_coverage[mut_key][sample_name]
                    self.sample_mafs[sample_name].append(maf)

                if self.dis_phred_coverage[mut_key][sample_name] >= 0:
                    self.sample_dis_phred_coverages[sample_name].append(self.dis_phred_coverage[mut_key][sample_name])

        # calculate p-values for presence and absence for all given variants
        self.present_p_values = calculate_p_values(self.mut_reads, self.phred_coverage, false_positive_rate)

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

        for mut_key, gene_name in gene_names.items():

            # add mutation name
            self.mut_keys.append(mut_key)
            # add gene name
            self.gene_names.append(gene_names[mut_key])

            # determine exact position of variant
            var = mut_key.split('__')
            if var[0].lower().startswith('chr'):
                if max(var[0].find('p'), var[0].find('q')) > -1:
                    chrom = var[0][3:max(var[0].find('p'), var[0].find('q'))]
                else:
                    chrom = var[0][3:]
            else:
                chrom = var[0]              # chromosome
            start_pos = int(var[1])
            ref, alt = var[2].split('>')
            end_pos = start_pos + len(ref) - 1

            # position data of this variant
            self.mut_positions.append((chrom, start_pos, end_pos))

            # - - - - - - - - CLASSIFY MUTATIONS - - - - - - - - -
            for sample_name in self.sample_names:

                # logger.debug('Confidence p-value values {:.2e}, {:.2e} for {} mut-reads at {}x coverage'.format(
                #     self.present_p_values[sample_name][mut_key], self.absent_p_values[sample_name][mut_key],
                #     self.mut_reads[mut_key][sample_name], self.phred_coverage[mut_key][sample_name]))

                # was the null hypothesis rejected (declared as significant)?
                if (sample_name, mut_key) in sig_muts:

                    maf = float(self.mut_reads[mut_key][sample_name]) / self.phred_coverage[mut_key][sample_name]

                    self.data[len(self.mut_keys)-1].append(maf)
                    self.positives[sample_name] += 1

                # is there enough coverage supporting a conclusion
                elif min_absent_cov == 0 or self.phred_coverage[mut_key][sample_name] >= min_absent_cov \
                        or self.phred_coverage[mut_key][sample_name] == -1:     # coverage has not been reported

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

                # calculate posterior: log probability that VAF = 0, log probability that VAF > 0
                p0, p1 = get_log_p0(self.phred_coverage[mut_key][sample_name], self.mut_reads[mut_key][sample_name],
                                    self.bi_error_rate, self.bi_c0)
                self.log_p01[len(self.mut_keys)-1].append([p0, p1])
                # logger.debug('p0: {;.2e}, k: {}, n: {}.'.format(
                #     self.log_p0[len(self.mut_keys)-1], self.mut_reads[mut_key][sample_name],
                #     self.phred_coverage[mut_key][sample_name]))

        for sample_name in self.sample_names:
            logger.info('Sample {} classifications: '.format(sample_name)
                        + '{} positives; {} negatives; {} unknowns;'.format(
                        self.positives[sample_name], self.negatives[sample_name],
                        self.unknowns[0][sample_name]+self.unknowns[1][sample_name]))

        for sample_name in self.sample_names:
            logger.debug('Sample {} classifications: '.format(sample_name)
                         + '{} positives; {} negatives; {} positive unknowns, {} negative unknowns;'.format(
                         self.positives[sample_name], self.negatives[sample_name],
                         self.unknowns[0][sample_name], self.unknowns[1][sample_name]))

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
                        elif entry == 'Driver_Pathway':
                            col_pathway = col_idx
                        elif entry == 'Function':
                            col_function = col_idx

                    # check if the columns with needed information were found
                    if col_function == -1:
                        logger.info('Data file {} does not contain information '.format(maf_filename)
                                    + 'about the function of the mutated region. ')
                    if col_chr == -1 or col_start_pos == -1 or col_end_pos == -1:
                        logger.info('Data file {} does not contain information '.format(maf_filename)
                                    + 'about the exact mutation position on the chromosome. ')

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
                            self.mut_functions.append(afrs[col_function])
                        else:
                            self.mut_functions.append('')

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
                                logger.warn('Found allele frequency of {} indicating a different file format. '.format(
                                    float(afrs[sa_idx])))
                                if sa_idx + 1 == len(afrs):
                                    logger.warn('Column {} is not considered as a sample column from now on'.format(
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

    def read_vcf_directory(self, vcf_directory, min_sa_cov, min_sa_maf, false_positive_rate,
                           false_discovery_rate, min_absent_cov, normal_sample_name=None):
        """
        Read allele frequencies for all variants in the samples in the files of the given directory
        :param vcf_directory: directory with VCF files
        :param min_sa_cov: minimum median coverage per sample
        :param min_sa_maf: minimum median mutant allele frequency per sample
        :param false_positive_rate: false positive rate of the used sequencing technology
        :param normal_sample_name: do not consider given normal sample
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_absent_cov: minimum coverage at a non-significantly mutated position for absence classification
        :return: number of samples which were processed (independent of filtering)
        """

        samples = read_vcf_files(vcf_directory, excluded_samples=[normal_sample_name])

        processed_samples = self._process_samples(samples, min_sa_cov, min_sa_maf, false_positive_rate,
                                                  false_discovery_rate, min_absent_cov)

        return processed_samples

    def read_vcf_file(self, vcf_file, min_sa_cov, min_sa_maf, false_positive_rate,
                      false_discovery_rate, min_absent_cov, normal_sample_name=None):
        """
        Read allele frequencies for all variants in the samples in the given VCF file
        :param vcf_file: path to vcf file
        :param min_sa_cov: minimum median coverage per sample
        :param min_sa_maf: minimum median mutant allele frequency per sample
        :param false_positive_rate: false positive rate of the used sequencing technology
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_absent_cov: minimum coverage at a non-significantly mutated position for absence classification
        :param normal_sample_name:
        :return: number of processed samples
        """

        samples = read_vcf_files(vcf_file, excluded_samples=[normal_sample_name])
        processed_samples = self._process_samples(samples, min_sa_cov, min_sa_maf, false_positive_rate,
                                                  false_discovery_rate, min_absent_cov)

        return processed_samples

    def _process_samples(self, samples, min_sa_cov, min_sa_maf, fpr,
                         false_discovery_rate, min_absent_cov):
        """
        Determine which variants are shared among the samples
        Remove low quality samples not meeting the median coverage or MAF
        Find significantly mutated variants while controlling for the false discovery rate
        :param samples: list of relevant samples
        :param min_sa_cov: minimum median coverage per sample
        :param min_sa_maf: minimum median mutant allele frequency per sample
        :param fpr: false positive rate of the used sequencing technology
        :param false_discovery_rate: control the false-positives with the benjamini-hochberg procedure
        :param min_absent_cov: minimum coverage at a non-significantly mutated position for absence classification
        :return: processed samples
        """

        # gene names and mutation pathways are not available in VCF files
        self.gene_names = None
        self.mut_functions = None

        # raw sequencing data
        self.mut_reads = defaultdict(dict)
        self.phred_coverage = defaultdict(dict)

        # calculate median coverage and median MAF of all confirmed present mutations in each sample
        self.sample_phred_coverages = defaultdict(list)
        self.sample_mafs = defaultdict(list)

        # Calculate p-values for confidence in presence and absence
        self.present_p_values = defaultdict(dict)
        self.absent_p_values = defaultdict(dict)

        heap = []
        tmp_vars = []

        self.sample_names = []
        self.mut_keys = []

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
                        self._add_variant(tmp_vars, fpr)

                else:       # variant at different position
                    # process old identical variants and remove them from the temporal list
                    self._add_variant(tmp_vars, fpr)

                    # add the new variant to the list to process
                    tmp_vars.append((variant, sample_name))

            else:
                # add new variant
                tmp_vars.append((variant, sample_name))

        logger.info("{} samples have been processed. ".format(len(self.sample_names)))
        logger.info("Completed reading of allele frequencies at {} mutated positions. ".format(len(self.mut_keys)))

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
        non_log_p0 = math.log(int_sets.NO_DATA_P0)
        non_log_p1 = math.log(1.0 - int_sets.NO_DATA_P0)

        # - - - - - - - - CLASSIFY MUTATIONS - - - - - - - - -
        for mut_id, mut_key in enumerate(self.mut_keys):

            for sample_name in self.sample_names:

                # logger.debug('Confidence p-value values {:.2e}, {:.2e} for {} mut-reads at {}x coverage'.format(
                #     self.present_p_values[sample_name][mut_key], self.absent_p_values[sample_name][mut_key],
                #     self.mut_reads[mut_key][sample_name], self.phred_coverage[mut_key][sample_name]))

                # was the null hypothesis rejected (declared as significant)?
                if (sample_name, mut_key) in sig_muts:

                    maf = float(self.mut_reads[mut_key][sample_name]) / self.phred_coverage[mut_key][sample_name]

                    self.data[mut_id].append(maf)
                    self.positives[sample_name] += 1

                # is there enough coverage supporting a conclusion
                elif min_absent_cov == 0 or sample_name not in self.phred_coverage[mut_key] \
                        or self.phred_coverage[mut_key][sample_name] >= min_absent_cov:

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

                # calculate posterior: log probability that VAF = 0
                if self.phred_coverage[mut_key][sample_name] < 0:   # no sequencing data in this sample
                    self.log_p01[mut_id].append([non_log_p0, non_log_p1])

                else:                                               # calculate posterior according to prior and data
                    p0, p1 = get_log_p0(self.phred_coverage[mut_key][sample_name], self.mut_reads[mut_key][sample_name],
                                        self.bi_error_rate, self.bi_c0)
                    self.log_p01[mut_id].append([p0, p1])

                # logger.debug('p0: {:.2e}, k: {}, n: {}.'.format(
                #     self.log_p0[mut_id][len(self.log_p0[mut_id])-1], self.mut_reads[mut_key][sample_name],
                #     self.phred_coverage[mut_key][sample_name]))

        for sample_name in self.sample_names:
            logger.info('Sample {} classifications: '.format(sample_name)
                        + '{} positives; {} negatives; {} unknowns;'.format(
                        self.positives[sample_name], self.negatives[sample_name],
                        self.unknowns[0][sample_name]+self.unknowns[1][sample_name]))

        for sample_name in self.sample_names:
            logger.debug('Sample {} classifications: '.format(sample_name)
                         + '{} positives; {} negatives; {} positive unknowns, {} negative unknowns;'.format(
                         self.positives[sample_name], self.negatives[sample_name],
                         self.unknowns[0][sample_name], self.unknowns[1][sample_name]))

        # calculate median coverage across all samples
        coverages = []
        for mut_key in self.mut_reads.keys():
            for sample_name in self.mut_reads[mut_key].keys():
                if self.phred_coverage[mut_key][sample_name] >= 0:
                    coverages.append(self.phred_coverage[mut_key][sample_name])

        logger.info('Median coverage in patient {}: {} (mean: {:.1f})'.format(self.name, np.median(coverages),
                                                                              np.mean(coverages)))

        return len(self.discarded_samples) + self.n

    def _add_variant(self, tmp_vars, fpr):
        """
        Process identical variants and remove them from the given list
        :param tmp_vars: list of identical variants in different samples
        :param fpr: false positive rate of the used sequencing technology
        """

        # process identical variants
        mut_key = '{}_{}_{}>{}'.format(tmp_vars[0][0].CHROM, tmp_vars[0][0].POS,
                                       tmp_vars[0][0].REF, tmp_vars[0][0].ALT[0])
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
                self.phred_coverage[mut_key][sample_name] = var.DP

                self.sample_phred_coverages[sample_name].append(var.DP)
                if var.AD[1] > 0:
                    self.sample_mafs[sample_name].append(var.BAF)

                self.present_p_values[sample_name][mut_key] = calculate_present_pvalue(var.AD[1], var.DP, fpr)

                # get variant in the next called sample
                if len(tmp_vars) > 0:
                    var, var_sample_name = tmp_vars.pop()
                else:
                    var_sample_name = -1

            else:
                # sequencing data information was not provided in this sample for this variant
                self.mut_reads[mut_key][sample_name] = -1       # -1 = unknown which is different from 0
                self.phred_coverage[mut_key][sample_name] = -1   # -1 = unknown which is different from 0

    def analyze_data(self):
        """
        Process data and apply various filters
        """

        # Possibility to implement filters!!!
        # Look at the average/median frequency of founder mutations
        # Filter mutations out with less than half of the average founder mutation frequency
        # data_utils.remove_contradicting_mutations(self.data)

        # determine in which samples each mutation is present (positives)
        for mut in range(len(self.data)):
            for sa_idx, maf in enumerate(self.data[mut]):
                if 0 < maf:
                    self.samples[sa_idx].add(mut)
                    self.mutations[mut].add(sa_idx)

        avg_mutations = 0
        for sa_idx, sa_name in enumerate(self.sample_names):
                logger.debug("Mutations present in sample {}: {}, {}".format(
                    sa_name, len(self.samples[sa_idx]),
                    str(self.samples[sa_idx]) if len(self.samples[sa_idx]) < 200 else ''))
                avg_mutations += len(self.samples[sa_idx])

        avg_mutations /= float(len(self.sample_names))
        logger.info('The average number of mutations per sample in patient {} is {:.1f}.'.format(
            self.name, avg_mutations))

        # total number of exonic mutations if information is available
        if self.mut_functions is not None and len(self.mut_functions):
            assert len(self.mut_functions) >= len(self.mutations), \
                'Each mutation has an associated function: {} >= {}'.format(len(self.mut_functions),
                                                                            len(self.mutations))
            no_exonic_muts = sum(1 for mut_func in self.mut_functions if mut_func in settings.EXONIC_FILTER)
            logger.info('Number of exonic mutations {}/{} (={:.3f}).'.format(no_exonic_muts, len(self.mutations),
                        (float(no_exonic_muts) / len(self.mutations))))

        for mut_idx in range(len(self.mut_keys)):
            self.mutations[mut_idx] = frozenset(self.mutations[mut_idx])
            # logger.debug("Mutation {} ({}) is present in: {}".format(self.mut_names[mut_idx],
            #            self.gene_names[mut_idx],
            #            (','.join(self.sample_names[sa_idx] for sa_idx in self.mutations[mut_idx]))))

        self.determine_sharing_status()

        # keep list of mutations present in some of the used samples
        # excluding unknowns
        self.present_mutations = get_present_mutations(self.data, include_unknowns=False)
        self.shared_mutations = get_shared_mutations(self.data, include_unknowns=False)

        # calculate homogeneity index (Jaccard similarity coefficient)
        self._calculate_genetic_similarity()

    def determine_sharing_status(self):
        """
        Determine which samples share the various mutations
        (allele frequencies are ignored - as long they are greater than 0)
        """

        for pos, fres in self.data.items():

            if all([(fre > 0.0) for fre in fres]):
                self.founders.add(pos)
            else:
                self.shared_muts[sum(1 for fre in fres if fre > 0.0)].add(pos)

        logger.info("{:.1%} ({}/{}) of all distinct mutations are founders.".format(
            float(len(self.founders))/len(self.mutations), len(self.founders), len(self.mutations)))
        logger.info('In average {:.1f} ({:.1%}) mutations are unique (private) per sample.'.format(
            float(len(self.shared_muts[1])) / len(self.sample_names),
            (float(len(self.shared_muts[1])) / len(self.sample_names)) /
            (sum(len(muts) for sa_idx, muts in self.samples.items()) / len(self.sample_names))))

        for shared in range(len(self.sample_names)-1, -1, -1):
            if len(self.shared_muts[shared]) > 0:
                logger.info('Mutations shared in {} samples ({} of {} = {:.3f}): {}'.format(shared,
                            len(self.shared_muts[shared]), len(self.mutations),
                            float(len(self.shared_muts[shared])) / len(self.mutations),
                            self.shared_muts[shared] if len(self.shared_muts[shared]) < 200 else ''))

        # compute the number of shared and additional mutations among the sample pairs
        # similar to a distance matrix
        self.common_muts = [[set() for _ in range(self.n)] for _ in range(self.n)]
        self.add_muts = [[0 for _ in range(self.n)] for _ in range(self.n)]

        # helper to save calculation time:
        cmuts = [[0 for _ in range(self.n)] for _ in range(self.n)]

        for s1 in range(self.n):
            for s2 in range(self.n):

                self.common_muts[s1][s2] = self.samples[s1].intersection(self.samples[s2])
                cmuts[s1][s2] = len(self.common_muts[s1][s2])
                self.add_muts[s1][s2] = self.samples[s1].difference(self.samples[s2])

                # logger.debug('Sample {} has {} mutations in common with sample {} and {} in addition. '.format(
                #    self.sample_names[s1], cmuts[s1][s2], self.sample_names[s2], len(self.add_muts[s1][s2])))

            logger.info(('Sample {} has most ({}) mutations in common with sample(s) {}'
                         ' and least ({}) with sample(s) {}. ').format(self.sample_names[s1],
                        max(cms for idx, cms in enumerate(cmuts[s1]) if idx != s1),
                        str([self.sample_names[idx] for idx, cms in enumerate(cmuts[s1])
                            if cms == max(cms for idx, cms in enumerate(cmuts[s1]) if idx != s1)]),
                        min(cmuts[s1]),
                        str([self.sample_names[idx] for idx, cms in enumerate(cmuts[s1]) if cms == min(cmuts[s1])])))

        # displays common mutations among the samples excluding founder mutations
        # similar to a distance matrix
        logger.debug('Parsimony-informative mutations among the samples: ')
        for s1 in range(self.n):
            logger.debug(self.sample_names[s1]+': '
                         + ' '.join((str((muts-len(self.founders)) if muts != -1 else -1)) for muts in cmuts[s1]))

        # mutation patterns are sharing its mutations in exactly the same samples (basis for the weighting scheme)
        self.mps = defaultdict(set)

        # Compute all clones sharing mutations in exactly the same samples
        # Clones with mutations in all or only one sample are uninformative for
        # the creation of the most parsimonious tree
        for mut_idx, samples in self.mutations.items():
            # if 1 < len(samples) < len(self.sample_names):
            if 0 < len(samples):
                self.mps[samples].add(mut_idx)

        # show the 10 clones supported by the most mutations
        for key, value in islice(sorted(self.mps.items(), key=lambda x: len(x[1]), reverse=True), 10):
            logger.debug('Mutation pattern {} shares mutations: {} '.format(str(key), value))

        logger.info('Total number of distinct mutation patterns: {}'.format(len(self.mps)))

    def _filter_samples(self, min_sa_cov, min_sa_maf):

        # filter low median coverage and low median MAF samples
        discarded_samples = []
        # set passed sample names
        self.sample_names = []

        for sample_name in sorted(self.sample_phred_coverages.keys(),
                                  key=lambda item: (item.split('_')[0], int(item.split('_')[1]))
                                  if len(item.split('_')) > 1 and item.split('_')[1].isdigit() else (item, item)):

            # discard a sample if its median coverage is lower than the threshold
            if np.median(self.sample_phred_coverages[sample_name]) < min_sa_cov:

                logger.warn('Sample {} is discarded! Median phred coverage of {} is below the threshold {}.'.format(
                    sample_name, np.median(self.sample_phred_coverages[sample_name]),
                    min_sa_cov))

                self.sample_phred_coverages.pop(sample_name)
                if self.sample_dis_phred_coverages is not None:
                    self.sample_dis_phred_coverages.pop(sample_name)
                if sample_name in self.sample_mafs.keys():
                    self.sample_mafs.pop(sample_name)

                self.present_p_values.pop(sample_name)
                self.absent_p_values.pop(sample_name)

                discarded_samples.append(sample_name)

            # discard a sample if its median MAF is lower than the threshold
            elif np.median(self.sample_mafs[sample_name]) < min_sa_maf:

                logger.warn('Sample {} is discarded!'.format(sample_name)
                            + ' Median mutant allele frequency (MAF) of {:.3f} is below the threshold {}.'.format(
                            np.median(self.sample_mafs[sample_name]), min_sa_maf))

                self.sample_phred_coverages.pop(sample_name)
                if self.sample_dis_phred_coverages is not None:
                    self.sample_dis_phred_coverages.pop(sample_name)
                self.sample_mafs.pop(sample_name)

                self.present_p_values.pop(sample_name)
                self.absent_p_values.pop(sample_name)

                discarded_samples.append(sample_name)

            # sample passed filtering
            else:
                self.sample_names.append(sample_name)
                logger.info('Sample {}: median phred coverage {:.1f}, median MAF {:.3f}'.format(
                    sample_name, np.median(self.sample_phred_coverages[sample_name]),
                    np.median(self.sample_mafs[sample_name])))
                # logger.info('Median distinct phred coverage in sample {}: {}'.format(
                #     sample_name, np.median(self.sample_dis_phred_coverages[sample_name])))

        return discarded_samples

    def _calculate_genetic_similarity(self):
        """
        Calculate the genetic similarity among all samples in this patient
        (1) calculate the genetic distance between all pairs
        (2) calculate the Jaccard similarity coefficient between all pairs
        """

        self.gen_dis = [[0 for _ in range(self.n)] for _ in range(self.n)]
        self.sim_coff = [[0 for _ in range(self.n)] for _ in range(self.n)]

        for s1_idx in range(len(self.data[0])):
            for s2_idx in range(len(self.data[0])):

                disagree = 0
                present_agree = 0
                no_known_variants = 0

                for mut_idx in range(len(self.data)):

                    if self.data[mut_idx][s1_idx] > 0:          # mutation present
                        if self.data[mut_idx][s2_idx] == 0:     # disagree
                            no_known_variants += 1
                            disagree += 1

                        # mutation present in both
                        elif self.data[mut_idx][s2_idx] > 0:    # agree
                            no_known_variants += 1
                            present_agree += 1

                        elif self.data[mut_idx][s2_idx] == NEG_UNKNOWN:
                            # no indication of a mutation at this position but low coverage
                            pass

                    elif self.data[mut_idx][s1_idx] == 0:     # mutation absent
                        if self.data[mut_idx][s2_idx] > 0:
                            no_known_variants += 1
                            disagree += 1

                        elif self.data[mut_idx][s2_idx] == 0:   # agree
                            pass

                        elif self.data[mut_idx][s2_idx] == POS_UNKNOWN:
                            # indication of a mutation at this position but insufficient mutant reads
                            # maybe penalize distance with two epsilon
                            pass

                    # mutation believed to be present but insufficient mutant reads
                    elif self.data[mut_idx][s1_idx] == POS_UNKNOWN:
                        if self.data[mut_idx][s2_idx] == 0:
                            pass

                    # mutation believed to be absent but insufficient coverage
                    elif self.data[mut_idx][s1_idx] == NEG_UNKNOWN:
                        if self.data[mut_idx][s2_idx] > 0:
                            pass

                self.gen_dis[s1_idx][s2_idx] = disagree
                self.sim_coff[s1_idx][s2_idx] = 1.0 if no_known_variants == 0 \
                    else (float(present_agree) / no_known_variants)

        # Produce latex table with the genetic distance between samples
        # print('Genetic distance across the samples:')
        # print('Sample & '+' & '.join(self.sample_names[sa_idx].replace('_', ' ') for sa_idx in range(self.n))+' \\\\')
        # print('\\hline ')
        # for s1_idx in range(self.n):
        #     print('{} & '.format(self.sample_names[s1_idx].replace('_', ' '))
        #           + ' & '.join('${}$'.format(self.gen_dis[s1_idx][s2_idx]) for s2_idx in range(self.n))+' \\\\')
        #
        # print('Similarity index based on the fraction of shared mutations:')
        # print('Sample & '+' & '.join(self.sample_names[sa_idx].replace('_', ' ') for sa_idx in range(self.n))+' \\\\')
        # print('\\hline ')
        # for s1_idx in range(self.n):
        #     print('{} & '.format(self.sample_names[s1_idx].replace('_', ' '))
        #           + ' & '.join('${:.2f}$'.format(self.sim_coff[s1_idx][s2_idx]) for s2_idx in range(self.n))+' \\\\')
