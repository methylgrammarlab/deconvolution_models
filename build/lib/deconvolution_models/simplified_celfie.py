
import numpy as np
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *
from itertools import compress
from scipy.special import logsumexp


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
        # self.x = [self.x[0]]
        # self.beta = [self.beta[0]]
        self.num_iterations = num_iterations
        self.convergence_criteria = convergence_criteria
        self.x_c_m = [(x == METHYLATED) for x in self.x]
        self.x_c_u = [(x == UNMETHYLATED) for x in self.x]
        self.x_c_v = [~(x == NOVAL) for x in self.x]
        self.t = self.beta[0].shape[0]
        c, m = [arr.shape[0] for arr in self.x], [arr.shape[1] for arr in self.x]
        self.c, self.m = np.sum(c), np.sum(m)
        self.alpha = alpha

        self.log_beta = [np.log(x) for x in self.beta]
        self.log_one_minus_beta = [np.log(1-x) for x in self.beta]


    def add_pseudocounts(self, arr):
        arr[arr==0] += 100*pseudocount
        arr[arr==1] -= 100*pseudocount
        return arr

    def filter_no_coverage(self):
        ref_cov = [~np.isnan(a).all(axis=0) for a in self.beta] #no data in ref, remove
        self.beta =[self.beta[i][:,ref_cov[i]] for i in range(len(self.x))]
        self.x = [self.x[i][:,ref_cov[i]] for i in range(len(self.x))]
        has_cov = np.array([(~(x == NOVAL)).any() for x in self.x]) #empty regions, remove
        self.beta = list(compress(self.beta, has_cov))
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

    def q_function(self, z, alpha):
        ll = 0
        for window in range(len(self.x)):
            p_tilde = z[window]
            beta_tm = self.beta[window]
            x_c_m =  self.x_c_m[window].astype(int)
            x_c_u = self.x_c_u[window].astype(int)
            T, C, M = p_tilde.shape
            for t in range(T):
                for c in range(C):
                    for m in range(M):
                        ll += p_tilde[t,c,m]*(x_c_m[c,m]*np.log(beta_tm[t,m]) + x_c_u[c,m]*np.log(1-beta_tm[t,m]))
                        ll += p_tilde[t,c,m]*np.log(alpha[t])
        return ll

    def simplified_expectation(self, alpha):
        z = []
        for window in range(len(self.x)):
            beta_tm = self.beta[window]
            x_c_m =  self.x_c_m[window].astype(int)
            x_c_u = self.x_c_u[window].astype(int)
            T, M = beta_tm.shape
            C = x_c_u.shape[0]
            z_window = np.zeros((T, C, M))
            for t in range(T):
                for c in range(C):
                    for m in range(M):
                        z_window[t,c,m] = alpha[t] * (beta_tm[t,m]**x_c_m[c,m])*((1-beta_tm[t,m])**x_c_u[c,m])
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
        new_alpha = np.sum(all_z, axis=1).sum(axis=1)
        new_alpha /= np.sum(new_alpha)
        return new_alpha

    def test_convergence(self, new_alpha):
        alpha_diff = np.mean(abs(new_alpha - self.alpha)) / np.mean(abs(self.alpha))
        return alpha_diff < self.convergence_criteria

    def init_alpha(self):
        # alpha = np.random.uniform(size=(self.t))
        # alpha /= np.sum(alpha)
        # self.alpha = alpha
        # self.alpha = np.array([0.8, 0.2])
        self.alpha = np.array([0.6784, 0.3216])

    def two_step(self):
        '''
        perform EM for a given number of iterations
        :return: cell type proportions, log-likelihood
        '''
        self.init_alpha()
        q = []
        for i in range(self.num_iterations):
            z = self.simplified_expectation(self.alpha)
            new_alpha = self.maximization(z)
            q.append(self.q_function(z, new_alpha))

            if i and self.test_convergence(new_alpha):
                break

            else:  # set current evaluation of alpha and gamma
                self.alpha = new_alpha


        return self.alpha, i, q
#%%
# beta_tm = [np.zeros((2,2)), np.ones((2,3))]
# beta_tm[0][0,:], beta_tm[1][0,:] = 1,0
# mixture = [np.zeros((100,2)),np.zeros((100,3))]
# [x.fill(UNMETHYLATED) for x in mixture]
# mixture[0][:20,:], mixture[1][:80,:] = METHYLATED, METHYLATED
# r = CelfiePlus(mixture, beta_tm, num_iterations=1000,
#                 convergence_criteria=0.001)
# alpha, i ,q = r.two_step()