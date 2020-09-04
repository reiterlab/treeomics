"""Classify observed variants as cancer """

import logging
import csv
import sys
import re
import os
from collections import namedtuple
import seaborn as sns

try:    # check if varcode and pyensembl is available (necessary for Windows)
    from treeomics.utils.mutation_effects import get_top_effect_name
    VARCODE = True

except ImportError:
    # mutation effect prediction will not be performed since VarCode is not avilable
    VARCODE = False
    get_top_effect_name = None

__author__ = 'Johannes REITER'
__date__ = 'Feb 22, 2017'


# get logger
logger = logging.getLogger(__name__)


class Driver:

    MaxSourceSupport = 0    # maximum number of sources supporting a driver
    Colors = None           # driver color depending on support

    def __init__(self, gene_name, mutation_effect=None, cgc_driver=None, sources=None, variant=None, cgc_drivers=None):
        """
        Initialize driver
        :param gene_name: name of gene
        :param mutation_effect:
        :param cgc_driver: Is the variant in the CGC?
        :param sources: how was the gene identified as putative driver?
        :param variant: instance of VarCode variant is only needed if variant should be checked for CGC
        :param cgc_drivers: List of CGC
        """

        self.gene_name = gene_name
        self._sources = None
        self.sources = set() if sources is None else sources
        self.mutation_effect = mutation_effect
        self.genomic_location = None
        self.mutation_location = None
        self.mutation_change = None

        if cgc_driver is not None:
            self.cgc_driver = cgc_driver

        elif cgc_driver is None and cgc_drivers is not None:
            if gene_name in cgc_drivers:
                dri_pos = cgc_drivers[gene_name].genomic_location
                self.genomic_location = dri_pos

                if dri_pos is None or variant is None:
                    # no positions provided => assume it's a driver
                    self.cgc_driver = True

                # check if variant is at the same chromosome and is within given region
                elif (dri_pos[0] == variant.contig and
                        (dri_pos[1] <= variant.start <= dri_pos[2] or dri_pos[1] <= int(variant.end) <= dri_pos[2])):

                    # variant is among CGC region
                    self.cgc_driver = True

                else:
                    self.cgc_driver = False
            else:
                self.cgc_driver = False

        else:           # don't know if in CGC
            self.cgc_driver = None

    @property
    def sources(self):
        return self._sources

    @sources.setter
    def sources(self, sources):
        if len(sources) > Driver.MaxSourceSupport:
            Driver.MaxSourceSupport = len(sources)
        self._sources = sources

    @staticmethod
    def colors():

        if Driver.Colors is None or len(Driver.Colors) != Driver.MaxSourceSupport:

            cp = sns.light_palette("darkred", n_colors=Driver.MaxSourceSupport+1)
            # for r, g, b, _ in cp:
            #     cs.append("#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b)))
            Driver.Colors = cp

        return Driver.Colors


def _check_gene_match(patient_mut_location, user_drivers):
    pat_chrom, pat_pos, _ = patient_mut_location
    for driver in user_drivers.values():
        driver_chrom, driver_pos = driver.mutation_location
        if pat_chrom == driver_chrom and pat_pos == driver_pos:
            return True
    return False



def potential_driver(gene_name, patient_mut_location, user_drivers, variant=None, cgc_drivers=None):
    """
    Checks if gene name is in set among potential driver genes
    Moreover, if VarCode is installed, it checks for non-synonymous mutations and variants that affect splicing
    Last, if a dictionary with CGC is provided, it checks, if the gene name is in the list and if the position
    is in the right place
    :param gene_name: gene name where the variant occurred
    :param user_drivers: set with potential driver genes defined by user
    :param variant: instance VarCode variant
    :param cgc_drivers: dictionary of Cancer Gene Census
    :return: variant in a driver gene, potential mutation effect, in CGC, mutation effect
    """

    # driver_gene = gene_name in user_drivers
    driver_gene = _check_gene_match(patient_mut_location, user_drivers)

    if driver_gene and variant is not None and VARCODE:
        mut_effect = get_top_effect_name(variant)
        put_driver = is_functional(mut_effect)

        if not put_driver:
            return driver_gene, put_driver, None, mut_effect

    else:   # we can't predict the mutation effect and hence has to assume there is one
        put_driver = driver_gene
        mut_effect = None

    if cgc_drivers is None:
        return driver_gene, put_driver, None, mut_effect

    # are there known positions for this driver gene?
    if cgc_drivers is not None:
        if gene_name in cgc_drivers:
            dri_pos = cgc_drivers[gene_name].genomic_location

            if dri_pos is None or variant is None:
                # no positions provided => assume it's a driver
                cgc_driver = True

            # check if variant is at the same chromosome and is within given region
            elif (dri_pos[0] == variant.contig and
                    (dri_pos[1] <= variant.start <= dri_pos[2] or dri_pos[1] <= int(variant.end) <= dri_pos[2])):
                # variant is among CGC region
                cgc_driver = True

            else:
                cgc_driver = False

            return driver_gene, put_driver, cgc_driver, mut_effect

        else:   # gene name is not CGC
            return driver_gene, put_driver, False, mut_effect

    else:       # no CGC provided, hence, we don't know if the variant is in the CGC
        return driver_gene, put_driver, None, mut_effect


def is_functional(effect_name):
    """
    Can this variant effect the expression of a gene?
    :param effect_name: varcode variant top priority effect, see https://github.com/hammerlab/varcode
    :return: True for non-synonymous mutations and variants that affect splicing
    """

    # AlternateStartCodon: Replace annotated start codon with alternative start codon (e.g. "ATG>CAG")
    if effect_name == 'AlternateStartCodon':
        return True

    # ComplexSubstitution: Insertion and deletion of multiple amino acids
    elif effect_name == 'ComplexSubstitution':
        return True

    # Deletion:	Coding mutation which causes deletion of amino acid(s)
    elif effect_name == 'Deletion':
        return True

    # ExonLoss:	Deletion of entire exon, significantly disrupts protein
    elif effect_name == 'ExonLoss':
        return True

    # ExonicSpliceSite:	Mutation at the beginning or end of an exon, may affect splicing
    elif effect_name == 'ExonicSpliceSite':
        return True

    # FivePrimeUTR:	Variant affects 5' untranslated region before start codon
    elif effect_name == 'FivePrimeUTR':
        return False

    # FrameShiftTruncation:	A frameshift which leads immediately to a stop codon (no novel amino acids created)
    elif effect_name == 'FrameShiftTruncation':
        return True

    # FrameShift: Out-of-frame insertion or deletion of nucleotides, causes novel protein sequence
    # and often premature stop codon
    elif effect_name == 'FrameShift':
        return True

    # IncompleteTranscript:	Can't determine effect since transcript annotation is incomplete
    # (often missing either the start or stop codon)
    elif effect_name == 'IncompleteTranscript':
        return False

    # Insertion: Coding mutation which causes insertion of amino acid(s)
    elif effect_name == 'Insertion':
        return True

    # Intergenic: Occurs outside of any annotated gene
    elif effect_name == 'Intergenic':
        return False

    # Intragenic: Within the annotated boundaries of a gene but not in a region that's transcribed into pre-mRNA
    elif effect_name == 'Intragenic':
        return False

    # IntronicSpliceSite: Mutation near the beginning or end of an intron but less likely
    # to affect splicing than donor/acceptor mutations
    elif effect_name == 'IntronicSpliceSite':
        return True

    # Intronic:	Variant occurs between exons and is unlikely to affect splicing
    elif effect_name == 'Intronic':
        return False

    # NoncodingTranscript: Transcript doesn't code for a protein
    elif effect_name == 'NoncodingTranscript':
        return False

    # PrematureStop: Insertion of stop codon, truncates protein
    elif effect_name == 'PrematureStop':
        return True

    # Silent: Mutation in coding sequence which does not change the amino acid sequence of the translated protein
    elif effect_name == 'Silent':
        return False

    # SpliceAcceptor: Mutation in the last two nucleotides of an intron, likely to affect splicing
    elif effect_name == 'SpliceAcceptor':
        return True

    # SpliceDonor: Mutation in the first two nucleotides of an intron, likely to affect splicing
    elif effect_name == 'SpliceDonor':
        return True

    # StartLoss: Mutation causes loss of start codon, likely result is that
    # an alternate start codon will be used down-stream (possibly in a different frame)
    elif effect_name == 'StartLoss':
        return True

    # StopLoss:	Loss of stop codon, causes extension of protein by translation of nucleotides from 3' UTR
    elif effect_name == 'StopLoss':
        return True

    # Substitution:	Coding mutation which causes simple substitution of one amino acid for another
    elif effect_name == 'Substitution':
        return True

    # ThreePrimeUTR: Variant affects 3' untranslated region after stop codon of mRNA
    elif effect_name == 'ThreePrimeUTR':
        return False

    elif effect_name == 'unknown':
        return True

    else:
        logger.warning('Variant effect {} is not yet supported.'.format(effect_name))
        return True


def get_drivers(cgc_path, user_driver_path, reference_genome):
    """
    Read provided CSV files and extract likely driver genes and hotspot mutations
    :param cgc_path:
    :param user_driver_path:
    :param reference_genome:
    :return: dictionary of driver genes with instances of class Driver
    """

    # is path to file with cancer census gene set provided?
    if cgc_path is not None and os.path.isfile(cgc_path):
        if reference_genome not in cgc_path:
            logger.error('Given the name of the Cancer Gene Consensus file , the genomic coordinates '
                         'may not match the given reference genome {}: {}'.format(reference_genome, cgc_path))
        cgc_drivers = read_driver_file(cgc_path)
    else:
        cgc_drivers = dict()

    # is a CSV file with user defined driver gene names provided?
    if user_driver_path is not None and os.path.isfile(user_driver_path):
        # positions can only be provided through the cancer gene census csv file
        # see settings.py
        user_drivers = read_driver_file(user_driver_path)
    else:
        user_drivers = dict()

    # merge the lists of drivers
    # return merge_driver_lists(cgc_drivers, user_drivers)
    return cgc_drivers, user_drivers


def read_driver_file(driver_list_path, cancer_type=None):
    """
    Path to CSV file with cancer drivers
    Column "Gene_Symbol" with the gene name is required,
    columns "Genome_Location" and "CancerType" are optional
    :param driver_list_path: path to CSV file with driver gene
    :param cancer_type: only read driver genes for the given cancer type
    :return: dictionary with driver genes as keys and instances of class Driver as values
    """

    with open(driver_list_path, 'rU') as driver_file:

        logger.debug('Reading cancer driver list file {}'.format(driver_list_path))
        f_csv = csv.reader(driver_file)

        # regex patterns for making column names in CSV files valid identifiers
        p_replace = re.compile(r'(/)| |[(]|[)]')
        # p_remove = re.compile(r'#| ')
        p_remove = re.compile(r'#')

        headers = None
        driver_dict = dict()

        for row in f_csv:
            if row[0].startswith('#') or not len(row[0]):
                # skip comments
                continue
            elif headers is None:  # process data table header
                headers = [p_replace.sub('_', p_remove.sub('', e)) for e in row]

                logger.debug('Header: {}'.format(headers))
                # process rows in CSV file
                DriverEntry = namedtuple('DriverEntry', headers)

            else:  # process entries in given list of drivers
                driver = DriverEntry(*row)

                if cancer_type is not None and hasattr(DriverEntry, 'CancerType'):
                    if driver.CancerType != cancer_type:
                        continue

                if driver.Gene_Symbol in driver_dict:
                    d = driver_dict[driver.Gene_Symbol]
                else:
                    d = Driver(driver.Gene_Symbol)
                    driver_dict[driver.Gene_Symbol] = d

                if hasattr(DriverEntry, 'Genome_Location'):
                    # extract genome location
                    chrom, position = driver.Genome_Location.split(':', 1)
                    if position == '' or position == '-':
                        start_pos = 0
                        end_pos = sys.maxsize
                    else:
                        start_pos, end_pos = position.split('-', 1)

                    if (driver.Gene_Symbol in driver_dict and
                            driver_dict[driver.Gene_Symbol].genomic_location is not None and
                            driver_dict[driver.Gene_Symbol].genomic_location != (chrom, int(start_pos), int(end_pos))):
                        raise RuntimeError('Same driver gene name with different genomic location has been provided!')

                    d.genomic_location = (chrom, int(start_pos), int(end_pos))

                if hasattr(DriverEntry, 'Sources'):
                    sources = driver.Sources.split(';')
                    if len(sources) > 0 and (len(sources) > 1 or sources[0] != ''):
                        d.sources = set(sources)

                if hasattr(DriverEntry, 'Mutation_Location'):
                    chrom, position = driver.Mutation_Location.split(':', 1)
                    d.mutation_location = (chrom, int(position))

                if hasattr(DriverEntry, 'Change'):
                    d.mutation_change = driver.Change

        logger.info("Read {} entries in driver list file {}{}. ".format(
            len(driver_dict), driver_list_path,
            ' of cancer type '+cancer_type if (cancer_type is not None
                                               and hasattr(DriverEntry, 'CancerType')) else ''))

        return driver_dict


def merge_driver_lists(*driver_dicts):
    """
    Merge multiple dictionaries of driver genes, values genomic locations encoded as tuple (chromosome, start, end)
    :param driver_dicts:
    :return: merged dictionary
    """
    merged_driver_dict = dict()
    for driver_dict in driver_dicts:
        for name, driver in driver_dict.items():
            if name in merged_driver_dict:
                if merged_driver_dict[name] is None:
                    merged_driver_dict[name] = driver
                elif (merged_driver_dict[name].genomic_location is not None and
                        merged_driver_dict[name].genomic_location != driver.genomic_location):
                    raise RuntimeError('Same driver gene {} with different genomic location has been provided!'.format(
                        name))
                else:
                    merged_driver_dict[name].genomic_location = driver.genomic_location
                    for s in driver.sources:
                        merged_driver_dict[name].sources.add(s)
            else:
                merged_driver_dict[name] = driver

    logger.debug('Generated list with {} unique driver genes.'.format(len(merged_driver_dict)))

    return merged_driver_dict
