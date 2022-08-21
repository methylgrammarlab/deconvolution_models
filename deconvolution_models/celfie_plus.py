
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
        self.y, self.y_depths = self.add_pseudocounts(y, y_depths)
        self.x,  self.y, self.y_depths  = [self.x[150]],[self.y[150]],[self.y_depths[150]]
        self.filter_no_coverage()
        self.beta = [self.y[i]/self.y_depths[i] for i in range(len(self.y))]

        self.num_iterations = num_iterations
        self.convergence_criteria = convergence_criteria
        self.x_c_m = [(x == METHYLATED) for x in self.x]
        self.x_c_u = [(x == UNMETHYLATED) for x in self.x]
        self.x_c_v = [~(x == NOVAL) for x in self.x]
        self.t = self.beta[0].shape[0]
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

    def filter_empty_rows(self, reads):
        filtered = []
        for region in reads:
            if region.size > 0:
                filtered_region = region[~(region==NOVAL).all(axis=1),:]
                filtered.append(filtered_region)
            else: #empty region, we'll keep it for alignment later
                filtered.append(region)
        return filtered

    def add_pseudocounts(self, a, a_depths):
        res = []
        res_depths = []
        for i in range(len(a)):
            new_a = a[i].copy()
            new_a_depths = a_depths[i].copy()
            new_a[new_a==0] += 1
            new_a_depths[new_a==0] += 2
            new_a[new_a==new_a_depths] += 1
            new_a_depths[new_a == new_a_depths] += 2
            res.append(new_a)
            res_depths.append(new_a_depths)
        return res, res_depths

    def calc_term1(self):
        '''
        since beta is constant, so is the
        first term of the likelihood
        :return: log term1
        '''
        t1 = []
        log_beta = [np.log(x) for x in self.beta]
        log_one_minus_beta = [np.log(1 - x) for x in self.beta]
        for window in range(len(self.x)):
            x_c_m =  self.x_c_m[window].astype(int)
            x_c_u = self.x_c_u[window].astype(int)
            log_beta_win = np.nan_to_num(log_beta[window]).T
            log_one_minus_beta_win = np.nan_to_num(log_one_minus_beta[window]).T
            t1.append((np.matmul(x_c_m, log_beta_win) + np.matmul(x_c_u, log_one_minus_beta_win)).T)
        return t1

    def log_expectation(self, alpha):
        '''
        P(z=1|x, alpha_i)
        :param alpha: cell type proportions
        :return: probability of z
        '''
        z = []
        term1 = self.calc_term1()

        for window in range(len(self.x)):
            T, C = term1[window].shape
            a = np.tile(np.log(alpha), (C, 1)).T + term1[window]
            b = logsumexp(a, axis=0)
            log_z = a - np.tile(b, (T,1))
            z.append(np.exp(log_z))
        return z

    def log_likelihood(self, alpha):
        '''
        logP(x|alpha, beta)
        :param alpha: cell type proportions
        :return: log likelihood
        '''
        ll = 0
        term1 = self.calc_term1()
        for window in range(len(self.x)):
            log_t1 = term1[window]
            T, C = log_t1.shape
            ll += np.sum(logsumexp(np.tile(np.log(alpha), (C, 1)).T + log_t1, axis=0))
            ll += np.nansum(self.y[window]*np.log(self.beta[window]) + (self.y_depths[window] - self.y[window])*np.log(1-self.beta[window]))
        return ll

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
            z_win = np.nan_to_num(z[i])
            x_c_m = np.nan_to_num(self.x_c_m[i])
            x_c_v = np.nan_to_num(self.x_c_v[i])
            new_y.append(np.matmul(z_win, x_c_m) + self.y[i])
            new_y_depths.append(np.matmul(z_win, x_c_v) + self.y_depths[i])
        y, y_depths = self.add_pseudocounts(new_y, new_y_depths)
        new_beta = [y[a]/ y_depths[a] for a in range(len(y))]
        return new_alpha, new_beta

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
        if not self.alpha:
            self.init_alpha()
        self.alpha = np.array([0.6784, 0.3216])
        ll = []
        for i in range(self.num_iterations):
            ll.append(self.get_ll())
            z = self.log_expectation(self.alpha)
            new_alpha, new_beta = self.maximization(z)
            if i and self.test_convergence(new_alpha):
                break

            else:  # set current evaluation of alpha and gamma
                self.alpha = new_alpha
                self.beta = new_beta
        return self.alpha, i

    def get_ll(self):
        return self.log_likelihood(self.alpha)



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
# import matplotlib.pyplot as plt
#
# plt.scatter(np.arange(1,len(ll)+1),ll)
# plt.show()

