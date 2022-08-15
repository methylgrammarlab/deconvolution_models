
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

    def __init__(self, mixtures, y, y_depths, num_iterations=50, convergence_criteria=0.001, alpha=None):
        '''
        :param mixtures: data for deconvolution. c reads by m cpg sites
        :param beta: methylation probability for reference atlas
        :param num_iterations: maximum iterations for em
        :param convergence_criteria: stopping criteria for em
        '''
        self.x = self.filter_empty_rows(mixtures)
        self.y, self.y_depths = y, y_depths
        # self.x,  self.y, self.y_depths  = [self.x[150]],[self.y[150]],[self.y_depths[150]]
        # self.filter_no_coverage()

        self.beta = self.add_pseudocounts(self.y, self.y_depths)

        self.num_iterations = num_iterations
        self.convergence_criteria = convergence_criteria
        self.x_c_m = [(x == METHYLATED) for x in self.x]
        self.x_c_u = [(x == UNMETHYLATED) for x in self.x]
        self.x_c_v = [~(x == NOVAL) for x in self.x]
        self.t = self.beta[0].shape[0]
        c, m = [arr.shape[0] for arr in self.x], [arr.shape[1] for arr in self.x]
        self.c, self.m = np.sum(c), np.sum(m)
        self.alpha = alpha

    def filter_no_coverage(self):
        ref_cov = [self.y_depths[a].sum(axis=0) != NOVAL for a in range(len(self.y_depths))] #no data in ref, remove
        self.y = [self.y[i][:,ref_cov[i]] for i in range(len(ref_cov))]
        self.y_depths = [self.y_depths[i][:,ref_cov[i]] for i in range(len(ref_cov))]
        self.x = [self.x[i][:,ref_cov[i]] for i in range(len(self.x))]

        has_cov = np.array([(~(x == NOVAL)).any() for x in self.x]) #empty regions, remove
        self.y = list(compress(self.y, has_cov))
        self.y_depths = list(compress(self.y_depths, has_cov))
        self.x = list(np.array(self.x)[has_cov])


    def add_pseudocounts(self, a, a_depths):
        res = []
        for i in range(len(a)):
            arr = a[i]/a_depths[i]
            arr[a[i]==0] = ((a[i]+1)/(a_depths[i]+2))[a[i]==0]
            arr[a[i]==a_depths[i]] = ((a[i]+1)/(a_depths[i]+2))[a[i]==a_depths[i]]
            res.append(arr)
        return res

    def calc_term1(self):
        log_beta = [np.log(x) for x in self.beta]
        log_one_minus_beta = [np.log(1 - x) for x in self.beta]

        T, C, M = 0, 1, 2  # cell_types, reads, CpG sites
        t1 = []
        for window in range(len(self.x)):
            t1.append(np.nansum(self.x_c_m[window][np.newaxis, :, :] * log_beta[window][:, np.newaxis, :] +
                                self.x_c_u[window][np.newaxis, :, :] * log_one_minus_beta[window][:, np.newaxis, :],
                                axis=M))
        return t1

    def q_function(self, z, alpha):
        term1 = self.calc_term1()
        ll = 0
        for window in range(len(self.x)):
            ll += np.sum(z[window]*term1[window])
            ll += np.sum(z[window]*np.log(alpha)[:, np.newaxis])
            ll += np.sum(self.y[window]*np.log(self.beta[window]) +
                         (self.y_depths[window]-self.y[window])*np.log(1-self.beta[window]))
        return ll


    def filter_empty_rows(self, reads):
        filtered = []
        for region in reads:
            if region.size > 0:
                filtered_region = region[~(region==NOVAL).all(axis=1),:]
                filtered.append(filtered_region)
            else: #empty region, we'll keep it for alignment later
                filtered.append(region)
        return filtered

    def log_expectation(self, alpha):
        z = []
        term1 = self.calc_term1()
        for window in range(len(self.x)):
            T, C = term1[window].shape
            a = np.tile(np.log(alpha), (C, 1)).T + term1[window]
            b = logsumexp(a, axis=0)
            log_z = a - np.tile(b, (T,1))
            z.append(np.exp(log_z))
        return z


    def maximization(self, z):
        '''
        argmax value of cel type proportions
        :param z: cell type indicator
        :return: alpha
        '''
        all_z = np.hstack(z)
        new_alpha = np.sum(all_z, axis=1)
        new_alpha /= np.sum(new_alpha)

        new_y = []
        new_y_depths = []
        for i in range(len(self.beta)):
            new_y.append(self.y[i] + (self.x_c_m[i][np.newaxis,:,:]*z[i][:,:,np.newaxis]).sum(axis=1))
            new_y_depths.append(self.y_depths[i] + (self.x_c_v[i][np.newaxis,:,:]*z[i][:,:,np.newaxis]).sum(axis=1))
        self.y, self.y_depths = new_y, new_y_depths
        self.beta = self.add_pseudocounts(new_y, new_y_depths)
        return new_alpha

    def test_convergence(self, new_alpha):
        alpha_diff = np.mean(abs(new_alpha - self.alpha)) / np.mean(abs(self.alpha))
        return alpha_diff < self.convergence_criteria

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
        # self.alpha = np.array([0.01,0.99])
        q = []
        for i in range(self.num_iterations):
            z = self.log_expectation(self.alpha)
            q.append(self.q_function(z, self.alpha))
            new_alpha = self.maximization(z)
            if i and self.test_convergence(new_alpha):
                break

            else:  # set current evaluation of alpha and gamma
                self.alpha = new_alpha
        return self.alpha, i



#%%

# mixture = [np.zeros((100,2)),np.zeros((100,3))]
# [x.fill(UNMETHYLATED) for x in mixture]
# mixture[0][:20,:], mixture[1][:80,:] = METHYLATED, METHYLATED
# beta_tm = [np.zeros((2, 2)), np.ones((2, 3))]
# beta_tm[0][0, :], beta_tm[1][0, :] = 1, 0
# y = [x*100 for x in beta_tm]
# y_depths = [np.zeros((2,2)),np.zeros((2,3))]
# [t.fill(100) for t in y_depths]
# r = CelfiePlus(mixture, y, y_depths, num_iterations=1000,
#                 convergence_criteria=0.001)
# alpha, _ = r.two_step()