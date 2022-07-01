import numpy as np
import math


def gumbel_r_mom(x):
    scale = np.sqrt(6) / np.pi * np.std(x)
    loc = np.mean(x) - np.euler_gamma * scale
    return loc,scale

def gamma_mom(x):
    var = np.var(x)
    scale = var / np.mean(x)
    loc = (np.mean(x)*np.mean(x))/(var)
    shape = np.mean(x)/scale
    return loc, scale, shape

def stat(x):
    df1 = x.describe().transpose()
    df1["rho"] = x.corr(method="spearman")
    df1["skew"] = x.skew()
    return df1


def gev_lmom(x):
    """
    Determine parameters GEV distribution using L-moments
    https://link.springer.com/content/pdf/10.1007%2F978-3-319-44234-1.pdf
    """
    b0 = np.mean(x)
    lengtereeks = len(x)
    values_sorted = sorted(x)
    f1 = 1 / (lengtereeks*(lengtereeks-1))
    a = 0
    for j in range(1, lengtereeks):
        ab = (j) * values_sorted[j]
        a = ab + a
    b1 = f1 * a
    f2 = 1 / (len(x) * (len(x) - 1) * (len(x) - 2))
    a2 = 0
    for k in range(2, len(x)):
        ab2 = (k) * (k - 1) * values_sorted[k]
        a2 = ab2 + a2
    b2 = f2 * a2
    labda1 = b0
    labda2 = (2 * b1) - b0
    labda3 = (6 * b2) - (6 * b1) + b0

    t2 = labda2 / labda1
    t3 = labda3 / labda2
    c = (2 / (t3 + 3)) - (np.log(2) / np.log(3))
    from scipy.special import gamma, factorial
    kappa = (7.8590 * c) + (2.9554 * (math.pow(c,2)))
    alpha = (kappa * labda2) / (gamma(kappa + 1) * (1 - 2**-kappa))
    epsilon = labda1 + ((alpha / kappa) * (gamma(kappa + 1) - 1))
    loc = epsilon
    scale = alpha
    c = kappa
    return loc, scale, c

def plotpositions(i, N, plot_position_name = 'Benard'):
    """
    Determines how to plot observations on for the return period plots, based on a given plotposition type.

    Parameters:
    -----------
    i: float
        Original observation for which to find the plot position.
    N: int
        Total length of dataset.
    plot_position_name: str
        Chosen plot position. Choose from: Beard, Gringorten, Hazen, Weibull, Benard

    Returns:
    --------
    pp: float
        The plot position for the observation
    """
    if plot_position_name == 'Beard':
        pp = 1 / ((i - 0.31) / (N + 0.38))
        return pp
    elif plot_position_name == 'Gringorten':
        pp = 1 / ((i - 0.44) / (N + 0.12))
        return pp
    elif plot_position_name == 'Hazen':
        pp = 1 / ((i - 0.5) / N)
        return pp
    elif plot_position_name == 'Weibull':
        pp = 1 / ((i / (N + 1)))
        return pp
    elif plot_position_name == 'Benard':
        pp = 1 / ((i - 0.3) / (N + 0.4))
        return pp