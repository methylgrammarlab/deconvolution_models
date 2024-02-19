
import numpy as np
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *

class READMeth:
    '''
    Read-based EM Algorithm for Deconvolution of Methylation sequencing
    '''
    pseudocount = 1e-10

    def __init__(self, mixtures, lambda_t, theta_high, theta_low, num_iterations=50, convergence_criteria=0.001, alpha=None):
        '''
        :param mixtures: data for deconvolution. c reads by m cpg sites
        :param lambda_t: prob of epistateH per cell type
        :param theta_high: methylation prob under epistateH per cpg
        :param theta_low: methylation prob under epistateL per cpg
        :param num_iterations: maximum iterations for em
        :param convergence_criteria: stopping criteria for em
        '''
        self.x = mixtures
        self.x_c_v = [(~(x == NOVAL)).any() for x in self.x]
        self.Lt = self.add_pseudocounts(lambda_t)
        self.thetaH = theta_high
        self.thetaL = theta_low
        self.filter_no_coverage()
        self.log_Lt =  [np.log(t) for t in self.Lt]
        self.log_one_minus_Lt = [np.log(1-t) for t in self.Lt]
        self.log_thetaH = [np.log(t) for t in self.thetaH]
        self.log_thetaL = [np.log(t) for t in self.thetaL]
        self.log_one_minus_thetaH = [np.log(1-t) for t in self.thetaH]
        self.log_one_minus_thetaL = [np.log(1-t) for t in self.thetaL]
        self.num_iterations = num_iterations
        self.convergence_criteria = convergence_criteria
        self.x_c_m = [(x == METHYLATED) for x in self.x]
        self.x_c_u = [(x == UNMETHYLATED) for x in self.x]
        self.t = self.Lt[0].shape[0]
        c, m = [arr.shape[0] if arr.any() else NOVAL for arr in self.x], [arr.shape[1] if arr.any() else NOVAL for arr in self.x]
        self.c, self.m = np.sum(c), np.sum(m)
        self.alpha = alpha

    def filter_no_coverage(self):
        self.x = self.x[self.x_c_v]
        self.Lt = self.Lt[self.x_c_v]
        self.thetaH = self.thetaH[self.x_c_v]
        self.thetaL = self.thetaL[self.x_c_v]


    def add_pseudocounts(self, list_of_arrays):
        '''
        avoid prob 0 and 1 for logarithm
        :param arr: array of probability
        :return: array without 0 and 1
        '''
        for arr in list_of_arrays:
            arr[arr==1] -= self.pseudocount
            arr[arr==0] += self.pseudocount
        return list_of_arrays

    def simplified_expectation(self, alpha):
        z = []
        for window in range(len(self.x)):
            beta_tm = self.thetaH[window]*self.Lt[window][:, np.newaxis] + self.thetaL[window]*(1-self.Lt[window][:, np.newaxis])
            a = self.x_c_m[window][np.newaxis, :, :] * beta_tm[:,np.newaxis, :]
            b = self.x_c_u[window][np.newaxis, :, :] * (1 - beta_tm)[:,np.newaxis, :]
            z_window = np.prod(a+b, axis=2)
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
        i=0
        for i in range(self.num_iterations):
            z = self.simplified_expectation(self.alpha)
            new_alpha = self.maximization(z)
            if i and self.test_convergence(new_alpha):
                break

            else:  # set current evaluation of alpha and gamma
                self.alpha = new_alpha
        return self.alpha, i



