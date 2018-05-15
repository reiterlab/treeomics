#!/usr/bin/python
"""Data structure around sequencing sample"""
import logging
import heapq
import math

__author__ = 'Johannes REITER'

# create logger for application
logger = logging.getLogger('vcf-parser')
logger.setLevel(logging.DEBUG)


class Sample(object):
    """
    Data structure class for a DNA sequencing sample
    """

    def __init__(self, name):

        # holds all variants of the sample in the VCF file in a minheap
        self.variants = []
        # name of the sample
        self.name = name

    def add_variant(self, variant):
        """
        Push variant to heap of variants
        :param variant:
        """
        heapq.heappush(self.variants, variant)

        # mut_key = '{}_{}_{}>{}'.format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])


class Variant(object):
    """
    Data structure class for a DNA variant (SNV or short indel)
    """

    def __init__(self, chrom, pos, identifier, ref, alt, qual=None, filter_info=None, info=None,
                 gene_name=None, ccf=None, var_type=None):
        if chrom.startswith('chr'):
            self.CHROM = chrom[3:]          # chromosome
        else:
            self.CHROM = chrom              # chromosome
        self.POS = pos              # 1-based position of the start of the variant
        self.ID = identifier        # unique identifiers of the variant
        self.REF = ref              # reference allele
        self.ALT = alt.split(',')   # comma separated list of alternate non-reference alleles
        self.QUAL = qual            # phred-scaled quality score: -10log_10 p(no variant)
        self.FILTER = filter_info   # site filtering information
        self.INFO = info            # semicolon separated list of additional, user extensible annotations

        self.AD = None              # allelic depths for the ref and alt alleles (in ordered list)
        self.DP = None              # total read depth
        self.BAF = None             # B-allele frequency
        self.GENE_NAME = gene_name  # name of gene where variant occurred
        # functional type of mutation, eg. missense; instance of class VarType (utils.var_type.py)
        self.VAR_TYPE = var_type
        self.CCF = ccf              # Cancer cell fraction

        # self.mut_key = '{}_{}_{}>{}'.format(self.CHROM, self.POS, self.REF, self.ALT[0])

        if len(self.ALT) > 1:
            logger.warning('Multiple alternate alleles are given.')

    def set_allelic_depth(self, ad):
        """
        Set allelic depths for the ref and alt alleles given in order lists
        :param ad: list of allele counts for reference and alternates
        """

        self.AD = [int(ad) if int(ad) >= 0 else float('nan') for ad in ad.split(',')]

    def set_total_depth(self, dp):
        """
        Set the total read depth per sample
        :param dp: total reads in the tumor sample at this position
        """

        self.DP = float('nan') if math.isnan(float(dp)) or int(dp) < 0 else int(dp)

    def set_baf(self, fa):
        """
        Set the B Allele Frequency
        Fractions of reads (excluding MQ0 from both ref and alt) supporting each
        reported alternative allele, per sample
        :param fa: fraction of reads supporting the B allele (could be a list)
        """

        if len(self.ALT) > 1:
            logger.warning("List of BAFs should be calculated for multiple alternate alleles.")

        self.BAF = float(fa)

    def set_ccf(self, ccf):
        """
        Set cancer cell fraction for the ref and alt alleles in the order listed
        :param ccf: list of cancer cell fractions for reference and alternates
        """

        # self.CCF = [float(ccf) for ccf in ccf.split(',')]
        self.CCF = float(ccf)

    def __lt__(self, other):    # called if x < y

        if self.CHROM < other.CHROM:
            return True
        elif self.CHROM > other.CHROM:
            return False
        elif self.POS < other.POS:
            return True
        elif self.POS > other.POS:
            return False
        elif self.REF < other.REF:
            return True
        elif self.REF > other.REF:
            return False
        elif self.ALT < other.ALT:
            return True
        elif self.ALT > other.ALT:
            return False
        elif self.DP and other.DP and self.DP < other.DP:
            return True
        elif self.DP and other.DP and self.DP > other.DP:
            return False
        elif self.BAF and other.BAF and self.BAF < other.BAF:
            return True
        elif self.BAF and other.BAF and self.BAF > other.BAF:
            return False
        else:
            return False

    def __le__(self, other):    # called if x <= y

        if self.__lt__(other) or self.__eq__(other):
            return True
        else:
            return False

    def __gt__(self, other):    # called if x > y

        if self.CHROM > other.CHROM:
            return True
        elif self.CHROM < other.CHROM:
            return False
        elif self.POS > other.POS:
            return True
        elif self.POS < other.POS:
            return False
        elif self.REF > other.REF:
            return True
        elif self.REF < other.REF:
            return False
        elif self.ALT > other.ALT:
            return True
        elif self.ALT < other.ALT:
            return False
        elif self.DP and other.DP and self.DP > other.DP:
            return True
        elif self.DP and other.DP and self.DP < other.DP:
            return False
        elif self.BAF and other.BAF and self.BAF > other.BAF:
            return True
        elif self.BAF and other.BAF and self.BAF < other.BAF:
            return False
        else:
            return False

    def __ge__(self, other):    # called if x >= y

        if self.__gt__(other) or self.__eq__(other):
            return True
        else:
            return False

    def __eq__(self, other):    # called if x == y

        if self.CHROM == other.CHROM \
                and self.POS == other.POS \
                and self.REF == other.REF \
                and self.ALT == other.ALT \
                and self.DP and other.DP and self.DP == other.DP \
                and self.BAF and other.BAF and self.BAF == other.BAF:
            return True
        else:
            return False

    def __ne__(self, other):    # called if x != y or x<>y

        if self.__eq__(other):
            return False
        else:
            return True

    def __str__(self):
        if self.BAF:
            return 'Chr {}, pos {}, ref {}, alt {}, baf {:.3f} (info: {})'.format(
                self.CHROM, self.POS, self.REF, str(self.ALT), self.BAF, self.INFO)
        else:
            return 'Chr {}, pos {}, ref {}, alt {} (info: {})'.format(
                self.CHROM, self.POS, self.REF, str(self.ALT), self.INFO)
