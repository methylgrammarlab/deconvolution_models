
import numpy as np
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *
from itertools import compress
from scipy.special import logsumexp

def length_one(x):
    nvals = np.count_nonzero(x)
    res = np.zeros((nvals, x.shape[1]))
    row_ind, col_ind = np.where(x>0)
    for i in range(nvals):
        res[i,col_ind[i]] = x[row_ind[i], col_ind[i]]
    return res

def shuffle(x):
    res = np.zeros(x.shape)
    row_ind, col_ind = np.where(x>0)
    new_ind = row_ind.copy()
    np.random.seed(123)
    np.random.shuffle(new_ind)
    res[new_ind, col_ind] = x[row_ind, col_ind]
    return res

class CelfiePlus:
    '''
    Read-based EM Algorithm for Deconvolution of Methylation sequencing
    '''

    def __init__(self, mixtures, beta, num_iterations=50, convergence_criteria=0.001, alpha=None):
        '''
        :param mixtures: data for deconvolution. c reads by m cpg sites
        :param beta: methylation probability for reference atlas
        :param num_iterations: maximum iterations for em
        :param convergence_criteria: stopping criteria for em
        :param alpha: override random initialization
        '''
        self.x = self.filter_empty_rows(mixtures)
        self.beta = [self.add_pseudocounts(x) for x in beta]
        self.filter_no_coverage()
        # self.x = [length_one(self.x[0])]
        # self.x = [shuffle(self.x[0])]
        self.num_iterations = num_iterations
        self.convergence_criteria = convergence_criteria
        self.x_c_m = [(x == METHYLATED) for x in self.x]
        self.x_c_u = [(x == UNMETHYLATED) for x in self.x]
        self.x_c_v = [~(x == NOVAL) for x in self.x]
        self.t = self.beta[0].shape[0]
        self.alpha = alpha

        self.log_beta = [np.log(x) for x in self.beta]
        self.log_one_minus_beta = [np.log(1-x) for x in self.beta]
        self.log_term1 = self.calc_term1()

    def calc_term1(self):
        '''
        since beta is constant, so is the
        first term of the likelihood
        :return: log term1
        '''
        t1 = []
        for window in range(len(self.x)):
            x_c_m =  self.x_c_m[window].astype(int)
            x_c_u = self.x_c_u[window].astype(int)
            log_beta = np.nan_to_num(self.log_beta[window].T)
            log_one_minus_beta = np.nan_to_num(self.log_one_minus_beta[window].T)
            t1.append((np.matmul(x_c_m, log_beta) + np.matmul(x_c_u, log_one_minus_beta)).T)
        return t1

    def add_pseudocounts(self, arr):
        '''
        move slightly away from 1 and 0
        :param arr: numpy array
        :return: changes array inplace
        '''
        pc = 1e-10
        arr[arr==0] += pc
        arr[arr==1] -= pc
        return arr

    def filter_empty_rows(self, reads):
        '''
        remove mixture rows with no data
        :param reads: mixture
        :return: filtered mixture
        '''
        filtered = []
        for region in reads:
            if region.size > 0:
                filtered_region = region[~(region==NOVAL).all(axis=1),:]
                filtered.append(filtered_region)
            else: #empty region, we'll keep it for alignment later
                filtered.append(region)
        return filtered

    def filter_no_coverage(self):
        '''
        filter cpg sites with no information in atlas
        :return:
        '''
        ref_cov = [~np.isnan(a).all(axis=0) for a in self.beta] #no data in ref, remove
        self.beta =[self.beta[i][:,ref_cov[i]] for i in range(len(self.x))]
        self.x = [self.x[i][:,ref_cov[i]] for i in range(len(self.x))]
        has_cov = np.array([(~(x == NOVAL)).any() for x in self.x]) #empty regions, remove
        self.beta = list(compress(self.beta, has_cov))
        self.x = list(np.array(self.x)[has_cov])

    def log_likelihood(self, alpha):
        '''
        logP(x|alpha, beta)
        :param alpha: cell type proportions
        :return: log likelihood
        '''
        ll = 0
        for window in range(len(self.x)):
            log_t1 = self.log_term1[window]
            T, C = log_t1.shape
            ll += np.sum(logsumexp(np.tile(np.log(alpha), (C, 1)).T + log_t1, axis=0))
        return ll

    def log_expectation(self, alpha):
        '''
        P(z=1|x, alpha_i)
        :param alpha: cell type proportions
        :return: probability of z
        '''
        z = []
        for window in range(len(self.x)):
            T, C = self.log_term1[window].shape
            a = np.tile(np.log(alpha), (C, 1)).T + self.log_term1[window]
            b = logsumexp(a, axis=0)
            log_z = a - np.tile(b, (T,1))
            z.append(np.exp(log_z))
        return z

    def maximization(self, z):
        '''
        argmax value of cell type proportions
        :param z: probability of cell type indicator
        :return: alpha
        '''
        all_z = np.hstack(z)
        new_alpha = np.sum(all_z, axis=1)
        new_alpha /= np.sum(new_alpha)
        return new_alpha

    def test_convergence(self, new_alpha):
        alpha_diff = np.mean(abs(new_alpha - self.alpha)) / np.mean(abs(self.alpha))
        return alpha_diff < self.convergence_criteria

    def init_alpha(self):
        alpha = np.random.uniform(size=(self.t))
        alpha /= np.sum(alpha)
        self.alpha = np.array([0.6784, 0.3216]) ###

    def two_step(self):
        '''
        perform EM for a given number of iterations
        :return: cell type proportions, log-likelihood
        '''
        if not self.alpha:
            self.init_alpha()
        ll = []
        for i in range(self.num_iterations):
            ll.append(self.get_ll())
            z = self.log_expectation(self.alpha)
            new_alpha = self.maximization(z)
            if i and self.test_convergence(new_alpha):
                break

            else:  # set current evaluation of alpha and gamma
                self.alpha = new_alpha


        return self.alpha, i

    def get_ll(self):
        return self.log_likelihood(self.alpha)