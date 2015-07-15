"""Parsing and processing of VCF files"""
__author__ = 'Johannes REITER'
__date__ = 'July 21, 2014'

import logging
import csv
import re
import os
from collections import namedtuple
from utils.sample import Sample
from utils.sample import Variant

# create logger for application
logger = logging.getLogger('vcf-parser')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('vcf-parser.log')
fh.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s: %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


class VCFParser(object):
    """
    Read and process VCF files
    """

    def __init__(self, filename):

        # dictionary holding all samples of the VCF files
        # key=sample name, value=sample object
        self.samples = None

        self._parse_vcf_file(filename)

    def _parse_vcf_file(self, filename, filter_zero_maf=False):
        """
        Read in and process the given VCF file
        :param filename: path to the VCF file
        """

        with open(filename) as vcf_file:
            logger.info('Reading VCF file {}'.format(filename))
            f_tsv = csv.reader(vcf_file, delimiter='\t')

            # process rows in VCF file
            named_row = None

            # regex patterns for making column names in VCF files valid identifiers
            p_replace = re.compile(r'( |/)')
            p_remove = re.compile(r'#')

            self.samples = dict()

            for row in f_tsv:
                if row[0].startswith('##'):
                    # skip the meta-information
                    continue
                elif row[0].startswith('#'):                # process VCF header
                    headers = [p_replace.sub('_', p_remove.sub('', e)) for e in row]

                    logger.debug('Header: {}'.format(headers))
                    named_row = namedtuple('variant', headers)

                    if len(row) > 9:        # samples are present and hence their format has to specified first
                        logger.info('Found data for {} samples: {}'.format(len(row)-9, headers[9:]))

                        # add identified samples and generate sample objects
                        for sample_name in headers[9:]:
                            self.samples[sample_name] = Sample(sample_name)

                    else:
                        raise ValueError('No data is found in the provided VCF file: {}'.format(filename))

                else:                                       # process variants
                    var = named_row(*row)

                    for filter_name in var.FILTER.split(';'):

                        if filter_name == 'PASS':
                            continue
                        elif filter_name == 'StrandBiasFilter' \
                                or filter_name == 'Mask' \
                                or filter_name == 'SnpCluster' \
                                or filter_name == 'HARD_TO_VALIDATE':     # Variant data did not pass general filtering
                            # logger.debug('Mutation at chr {} and pos {} did not pass the filtering.'
                            #               .format(r.CHROM, r.POS))
                            break
                        elif filter_name == 'mf1':               # variant did not pass MuTect filtering
                            break
                        elif filter_name == 'GATKStandardFilter':               # variant did not pass GATK filtering
                            break
                        else:
                            logger.warn('Unrecognized filter value: {}'.format(filter_name))
                            break

                    else:
                        # variant passed filter and is called for all provided samples

                        named_format = namedtuple('sample', var.FORMAT.split(':'))

                        # generate separate variant directories for each sample
                        for sa_idx, sa in enumerate(var[9:], 9):

                            # Standard cancer format: VCF files contains two samples named NORMAL and PRIMARY
                            # if 'NORMAL' in headers and 'PRIMARY' in headers, NORMAL could be skipped
                            # if headers[sa_idx] == 'NORMAL':
                            #     continue

                            try:
                                sample = named_format(*(sa.split(':')))

                                variant = generate_variant(var, sample)

                                # the reason for these are that multiple samples have been merged
                                # in a single VCF file, but almost all point mutations occur only in one patient
                                if filter_zero_maf and variant.BAF == 0:
                                    # logger.warn('Excluded variant {} since its BAF is 0.'.format(str(variant)))
                                    # self.variants[headers[sa_idx]][(variant.CHROM, variant.POS)] = variant
                                    pass
                                else:
                                    # add variant to dictionary
                                    # self.variants[headers[sa_idx]][(variant.CHROM, variant.POS)] = variant
                                    self.samples[headers[sa_idx]].add_variant(variant)

                                # logger.debug(row)

                            except TypeError:
                                if sa == './.':
                                    logger.debug('No data ({}) for sample {} at chr {} and pos {}'
                                                 .format(sa, headers[sa_idx], var.CHROM, var.POS))
                                else:
                                    logger.warn('Could not parse data {} for sample {} at chr {} and pos {}'
                                                .format(sa, headers[sa_idx], var.CHROM, var.POS))
                                    logger.info('Row {}'.format(row))

            for sample_name in headers[9:]:
                logger.debug('{} variants were detected in sample {}.'.format(
                    len(self.samples[sample_name].variants), sample_name))


def generate_variant(var, sample):
    """
    Create variant object with all the relevant information
    :param var: holds the VCF required information (chrom, pos, id, etc.)
    :param sample: holds all the information provided in each sample column in the VCF file
    :return variant:
    """
    # set VCF provided general information about the variant
    variant = Variant(var.CHROM,    # chromosome
                      var.POS,      # 1-based position of the start of the variant
                      var.ID,       # unique identifiers of the variant
                      var.REF,      # reference allele
                      var.ALT,      # comma separated list of alternate non-reference alleles
                      var.QUAL,     # phred-scaled quality score: -10log_10 p(no variant)
                      var.FILTER,   # site filtering information
                      var.INFO      # semicolon separated list of additional annotations
                      )

    # set read sequencing data
    variant.set_allelic_depth(sample.AD)
    variant.set_total_depth(sample.DP)

    # check if B allele frequency is provided otherwise calculate it from the allelic depth counts
    try:
        variant.set_baf(sample.FA)
    except AttributeError:
        if variant.AD[1] == 0:      # mutant allele has zero coverage
            baf = 0.0
        elif variant.AD[0] == 0:    # reference allele has zero coverage
            baf = 1.0
        else:
            baf = float(variant.AD[1]) / (variant.AD[0] + variant.AD[1])

        variant.set_baf(baf)

    return variant


def read_vcf_files(directory_name, excluded_samples=None):
    """
    Read all VCF files in the given directory and return a list of
    the samples including their variants
    :param directory_name: path to directory with VCF files
    :param excluded_samples: exclude variants in samples of this name (e.g. normal samples)
    :return: dictionary of relevant samples
    """

    samples = dict()    # relevant samples in directory

    for filename in os.listdir(directory_name):
        if filename.endswith('.vcf'):

            # parse VCF file
            s = read_vcf_file(os.path.join(directory_name, filename), excluded_samples)
            for sample_name, sample in s.items():
                samples[sample_name] = sample

    if len(samples) > 0:
        return samples
    else:
        logger.warn('No samples found in directory {}.'.format(directory_name))
        return samples


def read_vcf_file(vcf_file, excluded_samples=None):
    """
    Read the given VCF file return a list of the samples including their variants
    :param vcf_file: path to VCF file
    :param excluded_samples: exclude variants in samples of this name (e.g. normal samples)
    :return: dictionary of relevant samples
    """

    samples = dict()    # relevant samples in directory

    # parse VCF file
    vcf = VCFParser(vcf_file)
    for sample in vcf.samples.values():

        # exclude given samples (e.g. normal samples)
        if sample.name not in excluded_samples:
            samples[sample.name] = sample
            logger.info('Read sample {} with {} variants from file {}.'.format(
                sample.name, len(sample.variants), vcf_file))
        else:
            logger.info('Excluded sample {} with {} variants from file {}.'.format(
                sample.name, len(sample.variants), vcf_file))

    if len(samples) > 0:
        return samples
    else:
        logger.warn('No samples found in VCF file {}.'.format(vcf_file))
        return samples
