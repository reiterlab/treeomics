"""Statistical calculations"""
import logging
import math
from scipy.stats import binom
from scipy.special import logsumexp
from scipy.special import betainc
from scipy.special import gammaln
import numbers
import utils.int_settings as def_sets


__author__ = 'jreiter, jgerold'
__date__ = 'July 21, 2015'

# get logger for application
logger = logging.getLogger('treeomics')


def calculate_present_pvalue(mut_reads, coverage, fpr):
    """
    Calculate the p-value for the given FPR and number of mutant reads and coverage at this position
    :param mut_reads: number of mutant reads
    :param coverage: coverage at this position
    :param fpr: false-positive rate
    :return: p-value
    """

    return 1.0 - binom.cdf(mut_reads-1, coverage, fpr)


def calculate_absent_pvalue(mut_reads, coverage, min_maf):
    """
    Calculate the p-value that no mutation exists with a larger MAF
    for the number of mutant reads and coverage at this position
    P(X<=k|H_0); H_0: MAF > min_freq
    min frequency already accounts for false-positives
    :param mut_reads: number of mutant reads
    :param coverage: coverage at this position
    :param min_maf: minimal mutant allele frequency
    :return: p-value
    """
    return binom.cdf(mut_reads, coverage, min_maf)


def find_significant_mutations(p_values, fdr):
    """
    Find the significantly mutated positions using the Benjamini-Hochberg procedure (BH step-up)
    to control the false discovery rate at the given level fdr
    :param p_values: unsorted dictionary of p_values
    :param fdr: false discovery rate
    :return: set of keys which were determined to be significant
    """

    sig_muts = set()

    for i, (key, p_value) in enumerate(sorted(p_values.items(), key=lambda x: x[1]), 1):

        # reject null hypothesis (declare as significantly mutated)
        if p_value <= fdr * i / len(p_values):
            sig_muts.add(key)
        else:
            logger.debug('Conventional classification: ' +
                         'All variants with a p-value greater than {:.3e} (threshold: >{:.5e}) '.format(
                           p_value, fdr * i / len(p_values)) + 'are not significantly mutated.')
            break

    # logger.debug('{} variants out of {} are significantly mutated.'.format(len(sig_muts), len(p_values)))

    return sig_muts


def loglp(n, k, p, e):
    """
    The log likelihood of p, the true variant fraction
    :param n: total coverage
    :param k: variant count
    :param p: variant allele frequency
    :param e: sequencing error rate
    :return: log likelihood of p
    """

    return k*math.log((p*(1-e) + (1-p)*e)) + (n-k)*math.log((p*e+(1-p)*(1-e)))


def get_log_p0(n, k, e, c0, pseudo_alpha=None, pseudo_beta=None, cutoff_f=None):
    """
    Returns the (log) probability that the variant is absent and present
    :param n: coverage (number of reads)
    :param k: observed number of reads reporting the variant
    :param e: sequencing error rate
    :param c0: prior mixture parameter of delta function and uniform distribution
    :param pseudo_alpha: alpha parameter for the beta distributed part of the prior
    :param pseudo_beta: beta parameter for the beta distributed part of the prior
    :param cutoff_f: cutoff frequency for variants being absent
    :return: tuple (log probability that variant is absent, log probability that variant is present)
    """

    # cutoff frequency
    if cutoff_f is None:
        cutoff_f = 0.05
        logger.warning('Cutoff absent frequency was not set in the Bayesian inference model! Assumed {:1e}.'.format(
            cutoff_f))

    if not isinstance(n, numbers.Real) or not isinstance(k, numbers.Real):
        raise(RuntimeError('Sequencing read counts need to be numbers: n={}, k={}!'.format(n, k)))

    if math.isnan(n) or math.isnan(k):
        logger.warning('Sequencing read counts should be numbers: n={}, k={}!'.format(n, k))

    if e > cutoff_f:
        raise RuntimeError('Error rate e={} can not be higher than the calculated cutoff absent frequency {}'.format(
            e, cutoff_f))

    if k > n:
        raise RuntimeError('Number of variant reads cannot be higher than the sequencing depth: {} <= {}'.format(k, n))

    if pseudo_alpha is None:
        pseudo_alpha = def_sets.PSEUDO_ALPHA
    if pseudo_beta is None:
        pseudo_beta = def_sets.PSEUDO_BETA

    # pseudocounts are added to n and k, but removed for the computation of p0 at the end
    n_new = n + pseudo_alpha - 1 + pseudo_beta - 1
    k_new = k + pseudo_alpha - 1

    # overall weight assigned to the delta spike
    delta_value = math.log(c0) + loglp(n, k, 0, e)
    # whole integral normalizing constant so that c0 is recovered without data when cutoff = 0
    beta_norm_const = math.log(1-2*e) + (math.log(1.0 - c0) + gammaln(pseudo_alpha + pseudo_beta) -
                                         gammaln(pseudo_alpha) - gammaln(pseudo_beta))
    # correct the cutoff for change of variables
    new_cutoff = cutoff_f*(1-2*e) + e
    # compute the integral of allele frequencies below the cutoff, given we are in the beta function
    tmp_beta_inc = betainc(k+pseudo_alpha, n-k+pseudo_beta, new_cutoff)
    if tmp_beta_inc == 0.0:
        fraction_below_cutoff = -1e10   # very small number in log space
    else:
        fraction_below_cutoff = math.log(tmp_beta_inc)

    # compute the total weight of the beta distribution
    total_weight_beta = (-1*math.log(1-2*e) + gammaln(k_new+1) + gammaln(n_new-k_new+1) - gammaln(n_new+2) +
                         beta_norm_const)

    posterior_term_1 = -logsumexp([-fraction_below_cutoff, delta_value - fraction_below_cutoff - total_weight_beta])
    posterior_term_2 = -logsumexp([0, total_weight_beta - delta_value])
    # print("Posterior term 1 ", posterior_term_1)
    # print("Posterior term 2 ", posterior_term_2)
    p0 = logsumexp([posterior_term_1, posterior_term_2])
    try:
        if p0 >= 0.0:
            p1 = -1e10
        elif p0 > -1e-10:
            p1 = math.log(-p0)
        else:
            p1 = math.log(-math.expm1(p0))
    except ValueError:
        logger.error('ERROR: {}'.format(p0))
        raise RuntimeError('Posterior probability could not be calculated!')

    return p0, p1
