"""Read data from tab-separated-values files"""
__author__ = 'jreiter'
__date__ = 'Sept 10, 2014'


import logging
import csv
import re
from collections import namedtuple, defaultdict
from itertools import chain

# python 2, 3 compatibility
try:                    # python 2
    from itertools import izip_longest as zip_longest
except ImportError:     # python 3
    from itertools import zip_longest

# get logger for application
logger = logging.getLogger('treeomics')


def read_mutation_table(filename, excluded_columns=set(), exclude_chr_arm=False):
    """
    Reads TSV file with sequencing data of variants in multiple samples (columns)
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
        gene_names = dict()
        sample_names = []

        for row in f_tsv:
            if row[0].startswith('#') or not len(row[0]):
                # skip comments
                continue
            elif row[0].startswith('Chr'):                # process data table header
                headers = [p_replace.sub('_', p_remove.sub('', e)) for e in row]

                logger.debug('Header: {}'.format(headers))
                named_row = namedtuple('variant', headers)

                if len(row) > 4:        # samples are present
                    logger.debug('Found data for {} samples: {}'.format(len(row)-4, headers[4:]))

                    # add identified samples
                    for sample_name in headers[4:]:
                        if sample_name not in excluded_columns:
                            sample_names.append(sample_name)
                        else:
                            logger.debug('Exclude sample {}'.format(sample_name))

                else:
                    raise ValueError('No data is found in the provided file: {}'.format(filename))

            else:                                       # process variants
                var = named_row(*row)

                for sa_idx, sample_name in enumerate(headers[4:], 4):

                    # remove chromosome arm information from the key if parameter is set
                    if exclude_chr_arm and var.Chromosome.find('p') != -1:
                        key = var.Chromosome[:var.Chromosome.find('p')]+'__'+var.Position+'__'+var.Change
                    elif exclude_chr_arm and var.Chromosome.find('q') != -1:
                        key = var.Chromosome[:var.Chromosome.find('q')]+'__'+var.Position+'__'+var.Change
                    else:
                        key = var.Chromosome+'__'+var.Position+'__'+var.Change
                    gene_names[key] = var.Gene

                    if sample_name not in excluded_columns:
                        if var[sa_idx].lower() != 'n/a':
                            data[key][sample_name] = int(var[sa_idx])
                        else:
                            data[key][sample_name] = -1

        logger.info("Read {} entries in file {}. ".format(len(data), filename))

        return data, gene_names


def read_table(filename, variant_key_column_names, variant_key_pattern, data_column_names):
    """
    Reads TSV file with possibly multiple id columns and possibly multiple data columns
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
