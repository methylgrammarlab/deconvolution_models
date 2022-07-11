
import numpy as np
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *

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
        '''
        self.x = self.filter_empty_rows(mixtures)
        self.beta = [self.add_pseudocounts(x) for x in beta]
        self.filter_no_coverage()

        self.num_iterations = num_iterations
        self.convergence_criteria = convergence_criteria
        self.x_c_m = [(x == METHYLATED) for x in self.x]
        self.x_c_u = [(x == UNMETHYLATED) for x in self.x]
        self.x_c_v = [~(x == NOVAL) for x in self.x]
        self.t = self.beta[0].shape[0]
        c, m = [arr.shape[0] for arr in self.x], [arr.shape[1] for arr in self.x]
        self.c, self.m = np.sum(c), np.sum(m)
        self.alpha = alpha

    def add_pseudocounts(self, arr):
        arr[arr==0] += pseudocount
        arr[arr==1] -= pseudocount
        return arr

    def filter_no_coverage(self):
        has_cov = np.array([(~(x == NOVAL)).any() for x in self.x])
        self.beta = list(np.array(self.beta)[has_cov])
        self.x = list(np.array(self.x)[has_cov])

    def filter_empty_rows(self, reads):
        filtered = []
        for region in reads:
            if region.size > 0:
                filtered_region = region[~(region==NOVAL).all(axis=1),:]
                filtered.append(filtered_region)
            else: #empty region, we'll keep it for alignment later
                filtered.append(region)
        return filtered

    def simplified_expectation(self, alpha):
        z = []
        for window in range(len(self.x)):
            beta_tm = self.beta[window]
            a = self.x_c_m[window][np.newaxis, :, :] * beta_tm[:,np.newaxis, :]
            b = self.x_c_u[window][np.newaxis, :, :] * (1 - beta_tm)[:,np.newaxis, :]
            vals = np.nan_to_num(a+b, nan=1)
            z_window = np.prod(vals, axis=2, where=self.x_c_v[window]) #only where data
            z_window *= alpha[:,np.newaxis]
            z_window /= np.sum(z_window, axis=0)
            z.append(z_window)
        return z

    def maximization(self, z):
        '''
        argmax value of cel type proportions
        :param z: cell type indicator
        :return: alpha
        '''
        all_z = np.hstack(z)
        new_alpha = np.sum(all_z, axis=1) / self.c
        new_alpha /= np.sum(new_alpha)
        return new_alpha

    def test_convergence(self, new_alpha):
        alpha_diff = np.mean(abs(new_alpha - self.alpha)) / np.mean(abs(self.alpha))
        return alpha_diff < self.convergence_criteria

    def init_z(self):
        self.z = []
        for window in range(len(self.x)):
            C = self.x[window].shape[0]
            z = np.random.uniform(size=(self.t, C))
            z /= np.sum(z, axis =0) #each read has to come from somewhere
            self.z.append(z)

    def init_alpha(self):
        alpha = np.random.uniform(size=(self.t))
        alpha /= np.sum(alpha)
        self.alpha = alpha

    def two_step(self):
        '''
        perform EM for a given number of iterations
        :return: cell type proportions, log-likelihood
        '''
        self.init_alpha()

        for i in range(self.num_iterations):
            z = self.simplified_expectation(self.alpha)
            new_alpha = self.maximization(z)
            if i and self.test_convergence(new_alpha):
                break

            else:  # set current evaluation of alpha and gamma
                self.alpha = new_alpha
        return self.alpha, i
