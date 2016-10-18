"""Read and write data from and to tab-separated-values files """
import logging
import csv
import re
import math
import numpy as np
from collections import namedtuple, defaultdict
from itertools import chain
import utils.int_settings as def_sets


__author__ = 'Johannes REITER'
__date__ = 'Sept 10, 2014'


# python 2, 3 compatibility
try:                    # python 2
    from itertools import izip_longest as zip_longest
except ImportError:     # python 3
    from itertools import zip_longest

# get logger for application
logger = logging.getLogger('treeomics')


def read_mutation_table(filename, normal_sample=None, excluded_columns=set(), exclude_chr_arm=False):
    """
    Reads TSV file with sequencing data of variants in multiple samples (columns)
    :param normal_sample: name of normal sample; if provided, subsequent filtering on these data is performed
    :param filename: path to TSV file
    :param excluded_columns: samples (column names) which should not be returned
    :return: dictionary of the data per variant and sample, dictionary of the gene names per variant
    """

    with open(filename, 'rU') as data_file:

        logger.debug('Reading data file {}'.format(filename))
        f_tsv = csv.reader(data_file, delimiter='\t')

        # process rows in TSF file
        named_row = None

        # regex patterns for making column names in TSV files valid identifiers
        p_replace = re.compile(r'(/)| ')
        # p_remove = re.compile(r'#| ')
        p_remove = re.compile(r'#')

        headers = None
        data = defaultdict(dict)
        if normal_sample is not None:
            normal_data = dict()
        else:
            normal_data = None
        gene_names = dict()
        sample_names = []

        for row in f_tsv:
            if row[0].startswith('#') or not len(row[0]):
                # skip comments
                continue
            elif row[0].startswith('Chr') or row[0].startswith('Gene'):                # process data table header
                headers = [p_replace.sub('_', p_remove.sub('', e)) for e in row]

                logger.debug('Header: {}'.format(headers))
                named_row = namedtuple('variant', headers)

                # determine where the sample columns start
                max_info_idx = max(headers.index('Change'), headers.index('Gene'), headers.index('Position'),
                                   headers.index('Driver') if 'Driver' in headers else -1,
                                   headers.index('EndPosition') if 'EndPosition' in headers else -1)
                first_sa_col = max_info_idx + 1

                if len(row) > first_sa_col:        # samples are present
                    logger.debug('Found data for {} samples: {}'.format(len(row)-4, headers[4:]))

                    # add identified samples
                    for sample_name in headers[first_sa_col:]:
                        if sample_name not in excluded_columns and sample_name != normal_sample:
                            sample_names.append(sample_name)
                        elif sample_name == sample_name:
                            logger.info('Normal sample {}'.format(sample_name))
                        else:
                            logger.info('Exclude sample {}'.format(sample_name))

                else:
                    raise ValueError('No data is found in the provided file: {}'.format(filename))

            else:                                       # process variants
                var = named_row(*row)

                for sa_idx, sample_name in enumerate(headers[first_sa_col:], first_sa_col):

                    # remove chromosome arm information from the key if parameter is set
                    if exclude_chr_arm and var.Chromosome.find('p') != -1:
                        key = var.Chromosome[:var.Chromosome.find('p')]+'__'+var.Position+'__'+var.Change
                    elif exclude_chr_arm and var.Chromosome.find('q') != -1:
                        key = var.Chromosome[:var.Chromosome.find('q')]+'__'+var.Position+'__'+var.Change
                    else:
                        key = var.Chromosome+'__'+var.Position+'__'+var.Change
                    gene_names[key] = var.Gene

                    if sample_name not in excluded_columns and sample_name != normal_sample:
                        if var[sa_idx].lower() != 'n/a':
                            data[key][sample_name] = int(var[sa_idx])
                        else:
                            data[key][sample_name] = -1
                    elif sample_name == normal_sample:
                        normal_data[key] = int(var[sa_idx])

        logger.info("Read {} entries in file {}. ".format(len(data), filename))

        return data, gene_names, normal_data


def read_csv_file(filename, normal_sample=None, excluded_columns=set()):
    """
    Reads CSV file with sequencing data of variants (rows) across multiple samples (columns)
    Reference allele count column names have to end with '_ref'
    Alternate allele count column names have to end with '_alt'
    :param filename: path to CSV file
    :param normal_sample: name of normal sample
    :param excluded_columns: name of samples to exclude
    :return: dictionary with coverage, dict with mut reads, dict of gene names,
    dict of coverage in normal, dict of mut reads in normal
    """

    with open(filename) as data_file:

        logger.debug('Reading data file {}'.format(filename))
        f_csv = csv.reader(data_file)

        # regex patterns for making column names in TSV files valid identifiers
        # p_replace = re.compile(r'(/)|=')
        # p_remove = re.compile(r'#| |\?|,|\(|\)')

        mut_reads = defaultdict(dict)
        coverage = defaultdict(dict)

        if normal_sample is not None:
            norm_coverage = dict()
            norm_mut_reads = dict()
        else:
            norm_coverage = None
            norm_mut_reads = None

        gene_names = dict()
        sample_names = []
        headers = reference_cols = alternate_cols = None
        chrom_col_idx = None
        pos_col_idx = None
        ref_allele_idx = alt_allele_idx = None
        gene_col_idx = None
        ref_norm_idx = alt_norm_idx = None

        for row in f_csv:
            if row[0].startswith('#') or not len(row[0]):
                # skip comments
                continue

                # process file header
            elif row[0].startswith('Chr') or row[0].startswith('Gene') or row[0].startswith('Hugo_Symbol'):
                headers = row

                logger.info('Read file {} with header: {}'.format(filename, headers))

                reference_cols = [col_idx for col_idx, h in enumerate(headers) if h.endswith('_ref')]
                alternate_cols = [col_idx for col_idx, h in enumerate(headers) if h.endswith('_alt')]
                assert len(reference_cols) == len(alternate_cols), \
                    'CSV data columns were incorrectly parsed! Make sure that data columns end with "_ref" and "_alt"'

                chrom_col_idx = headers.index('Chromosome')
                pos_col_idx = headers.index('Pos')
                ref_allele_idx = headers.index('Ref')
                alt_allele_idx = headers.index('Alt')
                if 'Hugo_Symbol' in headers:
                    gene_col_idx = headers.index('Hugo_Symbol')

                if len(reference_cols) > 0:        # samples are present

                    # add identified samples
                    excluded_sample_ids = set()
                    for idx, ref_col_idx in enumerate(reference_cols):
                        sa_name = headers[ref_col_idx][0:-4]
                        if sa_name not in excluded_columns and sa_name != normal_sample:
                            sample_names.append(sa_name)
                        elif sa_name == normal_sample:
                            ref_norm_idx = headers.index(normal_sample+'_ref')
                            alt_norm_idx = headers.index(normal_sample+'_alt')
                            logger.info('Normal sample {}'.format(sa_name))
                            excluded_sample_ids.add(idx)
                        else:
                            excluded_sample_ids.add(idx)
                            logger.info('Exclude sample {}'.format(sa_name))
                    for idx in excluded_sample_ids:
                        del reference_cols[idx]
                        del alternate_cols[idx]

                    logger.debug('Found data for {} samples: {}'.format(len(sample_names), sample_names))
                else:
                    raise ValueError('No data is found in the provided file: {}'.format(filename))

            else:                                       # process variants
                if headers is None:
                    raise RuntimeError('Header of CSV file needs to be provided before the body!')

                key = '{}__{}__{}>{}'.format(row[chrom_col_idx], row[pos_col_idx],
                                             row[ref_allele_idx], row[alt_allele_idx])
                if gene_col_idx is not None:
                    gene_names[key] = row[gene_col_idx]

                for ref_col_idx, alt_col_idx, sa_name in zip(reference_cols, alternate_cols, sample_names):
                    if row[ref_col_idx].lower() != 'n/a':
                        coverage[key][sa_name] = int(row[ref_col_idx]) + int(row[alt_col_idx])
                        mut_reads[key][sa_name] = int(row[alt_col_idx])
                    else:
                        coverage[key][sa_name] = -1
                        mut_reads[key][sa_name] = -1

                if ref_norm_idx is not None:
                    if row[ref_norm_idx].lower() != 'n/a':
                        norm_coverage[key] = int(row[ref_norm_idx]) + int(row[alt_norm_idx])
                        norm_mut_reads[key] = int(row[alt_norm_idx])
                    else:
                        norm_coverage[key] = -1
                        norm_mut_reads[key] = -1

        assert len(norm_coverage) == len(mut_reads), \
            'Length of provided coverage vector and variant reads vector do not match!'

        logger.info("Read {} entries in file {}. ".format(len(norm_coverage), filename))

        return coverage, mut_reads, gene_names, norm_coverage, norm_mut_reads


def read_table(filename, variant_key_column_names, variant_key_pattern, data_column_names):
    """
    Reads TSV file with possibly multiple id columns and possibly multiple data columns
    :param filename: path to TSV file
    :param variant_key_column_names: list of ordered columns which form the key in the return dictionary
    :param variant_key_pattern: python format string describing how the dictionary key is created
    :param data_column_names: samples (column names) which should be returned
    :return: dictionary of the data per variant and sample, dictionary of the gene names per variant
    """

    with open(filename) as data_file:

        logger.debug('Reading data file {}'.format(filename))
        f_tsv = csv.reader(data_file, delimiter='\t')

        # regex patterns for making column names in TSV files valid identifiers
        p_replace = re.compile(r'(/)')
        p_remove = re.compile(r'#| ')

        # process data table header
        row = next(f_tsv)
        headers = [p_replace.sub('_', p_remove.sub('', e)) for e in row]
        print('Read file {} with header: {}'.format(filename, headers))

        key_columns = []
        data_columns = []
        for name in chain(variant_key_column_names, data_column_names):
            for column_idx, column_name in enumerate(headers):
                if name == column_name:
                    if name in variant_key_column_names:
                        key_columns.append(column_idx)
                    else:
                        data_columns.append(column_idx)
                    break
            else:
                raise ValueError('Column {} is not found in the provided file: {}'.format(
                    name, filename))
        for column_idx, column_name in enumerate(headers):
            if column_name == 'Sample':
                sample_column = column_idx
                break
        else:
            raise ValueError('Sample column is not found in the provided file: {}'.format(filename))

        data = defaultdict(dict)

        for row in f_tsv:
            if row[0].startswith('#') or not len(row[0]):
                # skip comments and blank lines
                continue
            else:                                       # process data
                # generate custom key string
                key = ''
                for key_idx, pat in zip_longest(key_columns, variant_key_pattern, fillvalue=''):
                    key += row[key_idx] + pat

                data[key][row[sample_column]] = [row[data_idx] for data_idx in data_columns]

        print("Completed reading {} entries from file {}. ".format(len(data), filename))

        return data


def write_posterior_table(filepath, sample_names, estimated_purities, sample_vafs, mut_positions,
                          gene_names, log_p01, betas):
    """
    Write file with posterior probabilities and parameter details to given file path
    :param filepath: path to output file with posterior probabilities
    :param sample_names: ordered list of sample names
    :param estimated_purities: estimated sample purities
    :param sample_vafs: dictionary with variant allele frequencies per sample
    :param mut_positions: data tuples about the mutation: chromosome, start position, and end position
    :param gene_names: arrary with gene names
    :param betas: beta values for the beta distribution used in the prior calculation
    :param log_p01: 2-dimensional array with log probabilities
    """

    with open(filepath, 'w') as post_file:

        post_writer = csv.writer(post_file, delimiter='\t')

        # write header
        header = ['Chromosome', 'StartPosition', 'EndPosition', 'Gene']
        header += [sample_name for sample_name in sample_names]
        post_writer.writerow(header)

        # write settings and sample properties
        post_writer.writerow(['#EstPurities', '', '', ''] +
                             [('{:.3%}'.format(estimated_purities[sa_name]) if sa_name in estimated_purities else '')
                              for sa_name in sample_names])
        post_writer.writerow(['#MedianVAFs', '', '', ''] +
                             ['{:.3%}'.format(np.median(sample_vafs[sa_name])) for sa_name in sample_names])
        post_writer.writerow(['#PriorAlpha', '', '', ''] +
                             ['{:.3f}'.format(def_sets.PSEUDO_ALPHA) for _ in range(len(sample_names))])
        post_writer.writerow(['#PriorBeta', '', '', ''] +
                             ['{:.3f}'.format(betas[sa_name]) for sa_name in sample_names])

        # run through all variants and write their posterior probabilities to the file
        for mut_idx, mut_pos in enumerate(mut_positions):

            row = list()
            row.append(mut_pos[0])
            row.append(mut_pos[1])
            row.append(mut_pos[2])

            if gene_names is not None:
                row.append(gene_names[mut_idx])
            else:
                row.append('')

            row += ['{:.9g}'.format(math.exp(log_p01[mut_idx][sa_idx][1])) for sa_idx in range(len(sample_names))]
            post_writer.writerow(row)


def write_mutation_patterns(phylogeny, filepath):
    """
    Generate file with the identified most reliable mutation and evolutionarily compatible mutation patterns
    :param phylogeny: data structure around phylogeny
    :param filepath: path to output file
    """

    with open(filepath, 'w') as mp_file:

        logger.debug('Write mutation pattern file {}'.format(filepath))
        mp_writer = csv.writer(mp_file, delimiter='\t')

        # write header
        mp_writer.writerow(['# Identified most reliable and evolutionarily compatible mutation patterns'])

        for mp in sorted(phylogeny.compatible_nodes):
            # if len(samples) == len(self.patient.sample_names) + len(self.sc_sample_ids):          # founder mut.
            #     self.mlh_founders.add(mut_idx)
            # elif 1 < len(samples) < len(self.patient.sample_names) + len(self.sc_sample_ids):     # shared mut.
            #     self.shared_mlh_mps[frozenset(samples)].add(mut_idx)
            # elif len(samples) == 1:                                                          # unique mut.
            #     for sa_idx in samples:
            #         self.mlh_unique_mutations[sa_idx].add(mut_idx)
            # else:
            #     self.mlh_absent_mutations.add(mut_idx)

            if len(mp) == 0:
                continue

            mp_writer.writerow(sorted([phylogeny.patient.sample_names[sa_idx] if phylogeny.patient.sc_names is None
                                       else phylogeny.patient.sc_names[sa_idx] for sa_idx in mp]))
