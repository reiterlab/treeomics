"""Statistical calculations"""
__author__ = 'jreiter, jgerold'
__date__ = 'July 21, 2014'

import logging
import math
from scipy.stats import binom
from scipy.misc import logsumexp
from scipy.integrate import quad
from scipy.special import gammaln
import utils.int_settings as def_sets


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
            logger.info('All variants with a p-value greater than {:.3e} (threshold: >{:.5e}) '.format(
                p_value, fdr * i / len(p_values)) + 'are not significantly mutated.')
            break

    logger.info('{} variants out of {} are significantly mutated.'.format(len(sig_muts), len(p_values)))

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


def get_log_p0(n, k, e, c0, pseudo_alpha=def_sets.PSEUDO_ALPHA, pseudo_beta=def_sets.PSEUDO_BETA):
    """
    Returns the (log) probability that the variant is absent
    :param n: coverage (number of reads)
    :param k: observed number of reads reporting the variant
    :param e: sequencing error rate
    :param c0: prior mixture parameter of delta function and uniform distribution
    :param pseudo_alpha: alpha parameter for the beta prior
    :param pseudo_beta: beta parameter for the beta prior
    :return: tuple (log probability that variant is absent, log probability that variant is present)
    """

    # pseudocounts are added to n and k, but removed for the computation of p0 at the end
    n = n + pseudo_alpha - 1 + pseudo_beta - 1
    k = k + pseudo_alpha - 1

    def exploglp(p):
        return math.exp(loglp(n, k, p, e))

    nonzero_integral = 0	    # The part of the posterior not at zero

    # Integration tolerances
    abs_tol = 0
    rel_tol = math.exp(-15)

    if n <= 300:    # coverage is below 0, integration is easy
        integral = quad(exploglp, a=0, b=1, epsrel=rel_tol, epsabs=0)
        nonzero_integral = math.log(integral[0])
        error = integral[1]

    elif float(k)/n >= e:     # approximation is easy: gamma approximation
        # logger.debug("Using gamma approx")
        # approximate the integral to be performed as a beta function after a small change of variables
        nonzero_integral = gammaln(k+1) + gammaln(n-k+1) - gammaln(n+2) - math.log(1-2*e)

    else:       # coverage high and VAF is below error rate, approximation is hard
        # Possible source of error for coverage larger than current dataset
        integral = quad(exploglp, a=0, b=1, epsrel=rel_tol, epsabs=0)
        if integral[0] == 0:
            # We will hit this limit when n -> large with k fixed, and then we want this section
            # of the likelihood to be very small
            # It may also occur at other times and give difficult (ie extremely small) behavior
            logger.warn("Warning! Exceeded performance of numerical integration, using approximation")
            # nonzeroInt = gammaln(k+1) + gammaln(n-k+1) - gammaln(n+2) - math.log(1-2*e)
            # LINEAR APPROX nonzeroInt = loglp(n, k, 0, e) - math.log((2*e-1)*k/e+(n-k)/(1-e)*(1-2*e))
            # SIMPSONS APPROX
            # This method approximates the integral using the first and second derivatives at 0
            # by estimating when a second order approximation crosses zero and integrating from 0
            # up to there (called xzero)
            # dldf = k/e*(1-2*e) + (n-k)/(1-e)*(2*e-1)
            # d2ldf2 = -1*k/(math.pow(e, 2))*math.pow(1-2*e,2) - (n-k)/(path.pow(1-e,2))*math.pow(2*e-1,2)
            # f0 = exploglp(0)
            # xzero = (-1*f0*dldf - math.sqrt(math.pow(f0*dldf,2)-2*f0*(f0*d2ldf2 + dldf)))/(d0*d2ldf2 + dldf)
            # area = f0*xzero + f0*dldf*math.pow(xzero,2)/2 + (f0*d2ldf2+dldf)/6*math.pow(xzero,3)
            # dldf = k/e*(1-2*e) + (n-k)/(1-e)*(2*e-1)
            #
            # d2ldf2 = -1*k/(math.pow(e, 2))*math.pow(1-2*e, 2) - (n-k)/(math.pow(1-e, 2))*math.pow(2*e-1, 2)
            # f0 = exploglp(0)
            # lf0 = loglp(n, k, 0, e)
            # # xzero = (-1*f0*dldf - math.sqrt(math.pow(f0*dldf,2)-2*f0*(f0*d2ldf2 + dldf)))/(f0*d2ldf2 + dldf)
            # x1 = lf0+math.log(-1*dldf)
            # x2 = 1.0/2*math.log(f0*math.pow(dldf, 2)-2*(f0*d2ldf2 + dldf)) + lf0/2
            # lx0 = logsumexp([x1, x2])
            # area = lf0 + lx0    # + f0*dldf*math.pow(xzero,2)/2 + (f0*d2ldf2+dldf)/6*math.pow(xzero,3))
            # nonzero_integral = area

            #Normal approx see supplement

            normalApprox = -1*(e-k/n)*(e-k/n)*n/(e*(1-e)) - math.log(e-k/n) - 1/2*math.log(n) + 1/2*math.log(2*math.pi) + math.log(e*(1-e))
            correction = loglp(n, k, k/n, 0) + 1/2*math.log(2*math.pi)
            nonzero_integral = normalApprox + correction

        else:
            nonzero_integral = math.log(integral[0])



    # carry along the beta function normalizing constant
    nonzero_integral += (math.log(1-c0) + gammaln(pseudo_alpha + pseudo_beta)
                         - gammaln(pseudo_alpha) - gammaln(pseudo_beta))
    # remove the pseudocounts for the delta function
    zero_integral = (math.log(c0) + loglp(n-pseudo_alpha + 1 - pseudo_beta + 1,
                                          k - pseudo_alpha + 1, 0, e))

    denom = logsumexp([zero_integral, nonzero_integral])

    # return p0, 1-p0 in log space
    return zero_integral - denom, nonzero_integral - denom
