import numpy as np
from scipy.special import logsumexp
METHYLATED = 1
UNMETHYLATED = 0


class READMe:
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
        self.Lt = self.add_pseudocounts(lambda_t)
        self.log_Lt =  [np.log(t) for t in self.Lt]
        self.log_one_minus_Lt = [np.log(1-t) for t in self.Lt]
        self.thetaH = theta_high
        self.thetaL = theta_low
        self.log_thetaH = [np.log(t) for t in self.thetaH]
        self.log_thetaL = [np.log(t) for t in self.thetaL]
        self.log_one_minus_thetaH = [np.log(1-t) for t in self.thetaH]
        self.log_one_minus_thetaL = [np.log(1-t) for t in self.thetaL]
        self.num_iterations = num_iterations
        self.convergence_criteria = convergence_criteria
        self.x_c_m = [(x == METHYLATED) for x in self.x]
        self.x_c_u = [(x == UNMETHYLATED) for x in self.x]
        self.t = self.Lt[0].shape[0]
        c, m = [arr.shape[0] for arr in self.x], [arr.shape[1] for arr in self.x]
        self.c, self.m = np.sum(c), np.sum(m)
        self.alpha = alpha


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

    def one_window_mu(self, x_c_m, x_c_u, Lt, log_thetaH, log_thetaL,
                      log_one_minus_thetaH, log_one_minus_thetaL, z):
        '''

        :param x_c_m: true where x is methylated
        :param x_c_u: true where x is unmethylated
        :param Lt: lambda values
        :param log_thetaH: array of high state for cpgs in window
        :param log_thetaL: array of low state for cpgs in window
        :param z: probability of cell origin per read
        :return: mu per read
        '''
        log_prior_high = np.log(np.sum(np.multiply(z, Lt[:,np.newaxis]), axis=0))
        log_prior_low = np.log(np.sum(np.multiply(z, (1-Lt)[:,np.newaxis]), axis=0))
        log_posterior_high = np.sum(x_c_m*(log_thetaH[np.newaxis,:]) + x_c_u*log_one_minus_thetaH[np.newaxis, :], axis=1)
        log_posterior_low = np.sum(x_c_m*(log_thetaL)[np.newaxis, :] + x_c_u*log_one_minus_thetaL[np.newaxis, :], axis=1)
        total_high = log_prior_high + log_posterior_high
        total_low = log_prior_low + log_posterior_low
        log_mu = total_high - logsumexp([total_high, total_low], axis=0)
        return np.exp(log_mu)

    def init_one_window_mu(self, x_c_m, x_c_u, Lt, thetaH, thetaL, alpha):
        prior_high = np.repeat(np.sum(alpha*Lt), x_c_m.shape[0])
        prior_low = 1-prior_high
        posterior_high = np.prod((thetaH[np.newaxis, :]*x_c_m) + ((1-thetaH)[np.newaxis, :]*x_c_u), axis=1)
        posterior_low = np.prod((thetaL[np.newaxis, :]*x_c_m) + ((1-thetaL)[np.newaxis, :]*x_c_u), axis=1)
        total_high = prior_high*posterior_high
        total_low = prior_low*posterior_low
        mu = total_high/(total_low+total_high)
        return mu

    def one_window_mu_no_log(self, x_c_m, x_c_u, Lt, thetaH, thetaL, z):
        prior_high = np.sum(z*Lt[:,np.newaxis], axis=0)
        prior_low = 1-prior_high
        posterior_high = np.prod((thetaH[np.newaxis, :]*x_c_m) + ((1-thetaH)[np.newaxis, :]*x_c_u), axis=1)
        posterior_low = np.prod((thetaL[np.newaxis, :]*x_c_m) + ((1-thetaL)[np.newaxis, :]*x_c_u), axis=1)
        total_high = prior_high*posterior_high
        total_low = prior_low*posterior_low
        mu = total_high/(total_low+total_high)
        return mu

    def calc_mu(self, z):
        mu = []
        for window in range(len(z)):
            # mu_win = self.one_window_mu(self.x_c_m[window], self.x_c_u[window], self.Lt[window], self.log_thetaH[window],
            #                             self.log_thetaL[window], self.log_one_minus_thetaH[window],
            #                             self.log_one_minus_thetaL[window], z[window])
            mu_win = self.one_window_mu_no_log(self.x_c_m[window], self.x_c_u[window], self.Lt[window],
                                               self.thetaH[window], self.thetaL[window], z[window])
            mu.append(mu_win)
        return mu

    def init_alpha(self):
        alpha = np.random.uniform(size=(self.t))
        alpha /= np.sum(alpha)
        self.alpha = alpha

    def init_mu(self):
        mu = []
        for window in range(len(self.thetaH)):
            mu_win = self.init_one_window_mu(self.x_c_m[window], self.x_c_u[window], self.Lt[window],
                                               self.thetaH[window], self.thetaL[window], self.alpha)
            mu.append(mu_win)
        return mu

    def one_window_z(self, log_mu, log_one_minus_mu, log_Lt, log_one_minus_Lt, log_alpha):
        log_high = log_alpha[:, np.newaxis] + log_Lt[:, np.newaxis] + log_mu[np.newaxis, :]
        log_low = log_alpha[:, np.newaxis] + log_one_minus_Lt[:, np.newaxis] + log_one_minus_mu[np.newaxis, :]
        total = logsumexp([log_high, log_low], axis = 0)
        log_z = total - logsumexp(total, axis = 0)
        return np.exp(log_z)

    def one_window_z_no_log(self, mu, Lt, alpha):
        high = mu*Lt[:,np.newaxis]
        low = (1-mu)*(1-Lt)[:,np.newaxis]
        total = alpha[:,np.newaxis]*(high+low)
        z = total/np.sum(total, axis=0)
        return z

    def calc_z(self, new_mu, alpha):
        z = []
        # log_alpha = np.log(alpha)
        for window in range(len(new_mu)):
            # z_t_c = self.one_window_z(np.log(new_mu[window]), np.log(1-new_mu[window]),
            #                         self.log_Lt[window], self.log_one_minus_Lt[window],
            #                         log_alpha)
            z_t_c = self.one_window_z_no_log(new_mu[window], self.Lt[window], alpha)
            z.append(z_t_c)
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

    def em(self):
        '''
        perform EM for a given number of iterations
        :return: cell type proportions, log-likelihood
        '''
        self.init_alpha()
        mu = self.init_mu()

        for i in range(self.num_iterations):
            # print(i)
            self.z = self.calc_z(mu, self.alpha)
            new_alpha = self.maximization(self.z)
            mu = self.calc_mu(self.z)

            if i and self.test_convergence(new_alpha):
                print(i, "readme break")
                break

            else:  # set current evaluation of alpha and gamma
                self.alpha = new_alpha
        print("readme")
        return self.alpha, i
#%%
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("data_file", type=str)
parser.add_argument("metadata_file", type=str)
parser.add_argument("t", type=int)
parser.add_argument("depth", type=int)
parser.add_argument("m_per_window", type=int)
parser.add_argument("windows_per_t", type=int)
parser.add_argument("num_iterations", type=int)
parser.add_argument("outpath", type=str)

args = parser.parse_args()
reads = np.load(args.data_file, allow_pickle=True)
thetaH, thetaL, lambdas = np.load(args.metadata_file, allow_pickle=True)
r = READMe(reads, lambdas, thetaH,
           thetaL, num_iterations=args.num_iterations,
           convergence_criteria=0.001)
estimated_alpha, i = r.em()
np.save(args.outpath, [estimated_alpha, np.array([i])], allow_pickle=True)
