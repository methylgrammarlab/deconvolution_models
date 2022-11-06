'''
portions of this code are copied from CelFiE:
Christa Caggiano, Barbara Celona, Fleur Garton, Joel Mefford, Brian Black, Catherine Lomen-Hoerth,
Andrew Dahl, Noah Zaitlen, "Comprehensive cell type decomposition of circulating cell-free DNA with CelFiE",
Nature Communications, May 2021

https://github.com/christacaggiano/celfie

Portions of the code are licensed under:

                    GNU AFFERO GENERAL PUBLIC LICENSE
                       Version 3, 19 November 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
 
'''


import bottleneck as bn  # substantially speeds up calculations with nan's
import numpy as np
import pandas as pd
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *

#%%
#CelFiE code

def add_pseudocounts(value, array, meth, meth_depths):

    """finds values of gamma where logll cannot be computed, adds pseudo-counts to make
    computation possible
    value: checks for a value that will prevent computation; either 0 or 1
    array: gamma array to check for inproper value
    meth: np array of methylation counts
    meth_depths: np array of total number of reads (meth counts + unmethylated counts)
    """

    axis0, axis1 = np.where(
        array == value  # find indices where value isn't able to be computed
    )

    meth[axis0, axis1] += 1  # add one read to methylated counts
    meth_depths[axis0, axis1] += 2  # adds two reads to total counts


def check_gamma(array):
    """checks for values of gamma where log likelihood cannot be computed, returns
    true if can be computed
    array: np array to check
    """

    return (0 in array) or (1 in array)


def expectation(gamma, alpha):
    """calculates the components needed for loglikelihood for each iteration of gamma and alpha
    gamma: np matrix of the estimated 'true' methylation proportions
    alpha: np matrix of estimated mixing proportions
    """

    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]

    p0 = (1.0 - gamma) * alpha
    p1 = gamma * alpha

    p0 /= np.nansum(p0, axis=0)[np.newaxis, ...]
    p1 /= np.nansum(p1, axis=0)[np.newaxis, ...]

    return p0, p1


def log_likelihood(p0, p1, x_depths, x, y_depths, y, gamma, alpha):
    """calculates the log likelihood P(X, Z, Y | alpha, gamma)
    p0: probability that read is methylated
    p1: probability read is unmethylated
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts
    gamma: estimated true methylation proportions
    alpha: estimated mixing proportions
    """

    tissues, sites, individuals = p0.shape[0], p0.shape[1], p0.shape[2]

    # Reshape arrays for faster computation
    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]

    y = y[..., np.newaxis]
    y_depths = y_depths[..., np.newaxis]

    x = x.T[np.newaxis, ...]
    x_depths = x_depths.T[np.newaxis, ...]

    ll = 0
    ll += np.sum((y + p1 * x) * np.log(gamma))
    ll += np.sum((y_depths - y + p0 * (x_depths - x)) * np.log(1.0 - gamma))
    ll += np.sum((p1 * x + (x_depths - x) * p0) * np.log(alpha))

    return ll


def maximization(p0, p1, x, x_depths, y, y_depths):

    """maximizes log-likelihood, calculated in the expectation step
    calculates new alpha and gamma given these new parameters
    p0: probability that read is methylated
    p1: probability read is unmethylated
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts
    """

    individuals = p0.shape[2]

    # initialize vector
    ones_vector = np.ones(shape=(y.shape[0]))
    new_alpha = np.zeros((x.shape[0], y.shape[0]))

    # in case of overflow or error, transform nans to 0 and inf to large float
    p0 = np.nan_to_num(p0)
    p1 = np.nan_to_num(p1)
    x = np.nan_to_num(x)
    x_depths = np.nan_to_num(x_depths)

    # break up calculation into two terms
    term0 = 0
    term1 = 0

    for n in range(individuals):

        new_alpha[n, :] = np.dot(p1[:, :, n], x[n, :]) + np.matmul(
            p0[:, :, n], (x_depths[n, :] - x[n, :])
        )

        term1 += p1[:, :, n] * (np.outer(ones_vector, x[n, :]))
        term0 += p0[:, :, n] * (np.outer(ones_vector, x_depths[n, :] - x[n, :]))

    gamma = (term1 + y) / (term0 + term1 + y_depths)  # calculate new gamma

    # check if gamma goes out of bounds, if so add psuedocounts to misbehaving y values
    if check_gamma(gamma):
        add_pseudocounts(1, gamma, y, y_depths)
        add_pseudocounts(0, gamma, y, y_depths)
        gamma = (term1 + y) / (term0 + term1 + y_depths)  # recalculate gamma

    # return alpha to be normalized to sum to 1
    normalized_new_alpha = new_alpha / np.sum(new_alpha, axis=1)[:, np.newaxis]
    return normalized_new_alpha, gamma



def em(x, x_depths, y, y_depths, num_iterations, convergence_criteria):
    """take in the input cfdna matrices and the reference data and
    runs the EM for the specified number of iterations, or stops once the
    convergence_criteria is reached
    x: methylated cfDNA read counts
    x_depths: depth of cfDNA
    y: methylated reference counts
    y_depths: depth of cfDNA
    convergence_criteria: difference between alpha + gamma before stopping
    """

    # randomly intialize alpha for each iteration
    alpha = np.random.uniform(size=(x.shape[0], y.shape[0]))
    alpha /= np.sum(alpha, axis=1)[:, np.newaxis]  # make alpha sum to 1
    # begin by checking for instances where there are no counts for y or y_depths
    add_pseudocounts(1, np.nan_to_num(y / y_depths), y, y_depths)
    add_pseudocounts(0, np.nan_to_num(y / y_depths), y, y_depths)

    # intialize gamma to reference values
    gamma = y / y_depths
    i = 0
    # perform EM for a given number of iterations
    for i in range(num_iterations):
        p0, p1 = expectation(gamma, alpha)
        a, g = maximization(p0, p1, x, x_depths, y, y_depths)

        # check convergence of alpha and gamma
        alpha_diff = np.mean(abs(a - alpha)) / np.mean(abs(alpha))
        gamma_diff = np.nanmean(abs(g - gamma)) / np.nanmean(abs(gamma)) #IU: changed to nanmean

        if i and ( #I added this
            alpha_diff + gamma_diff < convergence_criteria
        ):  # if convergence criteria, break
            break

        else:  # set current evaluation of alpha and gamma
            alpha = a
            # print(i, alpha)
            gamma = g

    ll = log_likelihood(
        p0, p1, x_depths, x, y_depths, y, gamma, alpha
    )  # print ll for random restarts
    return alpha, gamma, ll, i

# import matplotlib.pyplot as plt
# import seaborn as sns
# fig, ax = plt.subplots()
# # plt.scatter(np.arange(1, i + 2), all_l, label="Q function")
# plt.scatter(np.arange(1, i + 2), [a[0][1] for a in alphas], label="t1")
# # plt.scatter(np.arange(1, i + 2), [a[1] for a in alphas], label="t2")
# plt.xlabel("iteration")
# plt.show()

