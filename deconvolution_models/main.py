'''
MIT License

Copyright (c) 2022 irene unterman and ben berman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import click
import json
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
# sys.path.append("/Users/ireneu/PycharmProjects/deconvolution_models/deconvolution_models")
sys.path.append("/Users/ireneu/PycharmProjects/deconvolution_models/tests")
from deconvolution_models.celfie import em as celfie
from deconvolution_models.celfie_ish import CelfieISH as celfie_ish
from deconvolution_models.celfie_ish_reatlas import CelfieISHReatlas as reatlas
from deconvolution_models.epistate import Epistate as epistate
from deconvolution_models.UXM import uxm
import numpy as np
from epiread_tools.epiparser import EpireadReader, CoordsEpiread,AtlasReader, EpiAtlasReader, UXMAtlasReader
from epiread_tools.epiformat import epiformat_to_reader
from epiread_tools.naming_conventions import *
from epiread_tools.em_utils import calc_coverage, calc_methylated, calc_percent_U


class NotImplementedError(Exception):
    """Exception raised for calling nonexistent methods.

    Attributes:
        name -- class name
        method --method name
    """
    def __init__(self, method, name):
        self.message = "%s is not implemented for model %s"%(method, name)
        super().__init__(self.message)


class DECmodel:
    '''
    general class for deconvolution
    '''
    def __init__(self, config):
        self.config = config
        if "epiformat" in self.config:
            self.reader = epiformat_to_reader[self.config["epiformat"]]
        self.name, self.alpha, self.i = None, None, None

    def read_mixture(self):
        raise NotImplementedError("read_mixture", self.name)

    def read_atlas(self):
        raise NotImplementedError("read_atlas", self.name)

    def load_npy(self):
        raise NotImplementedError("load_npy", self.name)

    def write_npy(self):
        '''
        convenient for pipelines (easier to load)
        :return:
        '''
        np.save(self.config['outfile'], [self.alpha, np.array([self.i])], allow_pickle=True)

    def write_output(self):
        alpha = (self.alpha*100).reshape(1, -1)
        np.savetxt(self.config['outfile'], alpha, fmt='%.3f', delimiter='\t')

    def deconvolute(self):
        raise NotImplementedError("deconvolute", self.name)

    def run_from_npy(self):
        self.load_npy()
        self.deconvolute()
        self.write_npy()

    def run_model(self):
        self.read_mixture()
        self.read_atlas()
        self.deconvolute()
        if self.config["npy"]:
            self.write_npy()
        else:
            self.write_output()

class Celfie(DECmodel):
    '''
    wrapper for https://github.com/christacaggiano/celfie
    '''

    def __init__(self, config):
        super().__init__(config)

    def read_mixture(self):
        reader = self.reader(self.config)
        self.interval_order, self.matrices, self.cpgs, self.origins = reader.get_matrices_for_intervals()
        self.methylation, self.coverage = calc_methylated(self.matrices), calc_coverage(self.matrices)
        self.x, self.x_depths = np.hstack(self.methylation), np.hstack(self.coverage) #no summing
        if self.config["summing"]:
            self.x, self.x_depths = np.array([np.sum(x) for x in self.methylation]), \
                                    np.array([np.sum(x) for x in self.coverage])
    def read_atlas(self):
        reader = AtlasReader(self.config)
        atlas_intervals, y, y_depths = reader.meth_cov_to_meth_cov() #no summing
        interval_to_atlas = dict(zip([str(x) for x in atlas_intervals], np.arange(len(atlas_intervals))))
        self.y = [y[interval_to_atlas[str(x)]] for x in self.interval_order]
        self.y_depths = [y_depths[interval_to_atlas[str(x)]] for x in self.interval_order]

        if self.config["summing"]:
            self.y, self.y_depths = np.vstack([np.sum(x, axis=1) for x in self.y]).T, \
                                    np.vstack([np.sum(x, axis=1) for x in self.y_depths]).T
        else:
            self.y, self.y_depths = np.hstack(self.y), np.hstack(self.y_depths)


    def load_npy(self):
        self.matrices = np.load(self.config["data_file"], allow_pickle=True)
        self.y, self.y_depths = np.load(self.config["metadata_file"], allow_pickle=True)
        self.methylation, self.coverage = calc_methylated(self.matrices), calc_coverage(self.matrices)
        self.x, self.x_depths = np.hstack(self.methylation), np.hstack(self.coverage) #might be problematic with missing data
        if self.config["summing"]:
            self.x, self.x_depths = np.array([np.sum(x) for x in self.methylation]), \
                                    np.array([np.sum(x) for x in self.coverage])
            shapes = np.array([x.reshape(-1,1).shape[0] for x in self.methylation])
            shape_indices = np.cumsum(shapes)
            shape_indices = np.array([0]+list(shape_indices))
            cumsum_y = np.hstack([np.zeros((self.y.shape[0],1)), np.cumsum(self.y, axis = 1)])
            cumsum_y_depths = np.hstack([np.zeros((self.y_depths.shape[0],1)), np.cumsum(self.y_depths, axis = 1)])
            self.y = cumsum_y[:,shape_indices[1:]] - cumsum_y[:,shape_indices[:-1]]
            self.y_depths = cumsum_y_depths[:,shape_indices[1:]] - cumsum_y_depths[:,shape_indices[:-1]]

    def deconvolute(self):
        restarts = []
        for r in range(self.config['random_restarts']):
            estimated_alpha, estimated_gamma, ll, i = celfie(
                self.x.reshape((1,-1)), self.x_depths.reshape((1,-1)), self.y.astype(int), self.y_depths,
                self.config['num_iterations'], self.config['stop_criterion']
            )
            restarts.append((ll, estimated_alpha, i))
        ll_max, alpha_max, i_max = max(restarts)
        self.alpha, self.i = alpha_max.flatten(), i_max

class UXM(DECmodel):
    '''
    wrapper for https://github.com/nloyfer/UXM_deconv
    '''
    def __init__(self, config):
        super().__init__(config)
        self.U = [] #percent u
        self.N = [] #number of fragments = weights
        self.i = 0 #iteration, for compatibility with pipeline, means nothing
        self.min_length = self.config["min_length"]

    def read_mixture(self):
        reader = self.reader(self.config)
        self.interval_order, self.matrices, self.cpgs, self.origins = reader.get_matrices_for_intervals()
        self.matrices = [x.todense() for x in self.matrices]
        self.calc_u()

    def calc_u(self):
        for mat in self.matrices:
            x_c_v = np.array(mat != NOVAL)
            # filter short reads
            self.len_filt = (np.sum(x_c_v, axis=1) >= self.min_length).ravel()
            if not np.sum(self.len_filt):  # empty region
                self.U.append(0)
                self.N.append(0)
            else:
                small = mat[self.len_filt, :]
                self.U.append(calc_percent_U(small, self.config["u_threshold"]))
                self.N.append(small.shape[0])

    def filter_empty(self):
        empty = np.array(self.N) == 0
        self.U = list(np.array(self.U)[~empty])
        self.atlas=self.atlas[~empty,:]
        self.N = list(np.array(self.N)[~empty])
        #TODO: filter cpgs,intervals

    def read_atlas(self):
        reader = UXMAtlasReader(self.config)
        self.atlas_intervals, self.atlas = reader.read()
        self.sort_intervals()

    def sort_intervals(self):
        atlas_to_index = {str(v): k for k, v in dict(enumerate(self.atlas_intervals)).items()}
        atlas_order = np.array([atlas_to_index[str(x)] for x in self.interval_order])
        self.atlas = self.atlas[atlas_order, :]
        reorder = [self.atlas_intervals[i] for i in atlas_order]
        self.atlas_intervals = reorder

    def load_npy(self):
        self.matrices = list(np.load(self.config["data_file"], allow_pickle=True))
        self.calc_u()
        self.atlas = np.vstack(np.load(self.config["metadata_file"], allow_pickle=True))

    def deconvolute(self):
        self.filter_empty() #Note: this messes up the interval/cpgs
        self.i = np.sum(self.N)
        if self.config["weights"]:
            self.alpha = uxm(self.atlas, self.U, self.N)
        else:
            self.alpha = uxm(self.atlas, self.U)


class CelfieISH(DECmodel):

    def __init__(self, config):
        super().__init__(config)
        self.name = "celfie-ISH"
        if self.config['random_restarts'] > 1:
            raise NotImplementedError("random_restarts", self.name)

    def read_mixture(self):
        reader = self.reader(self.config)
        self.interval_order, self.matrices, self.cpgs, self.origins = reader.get_matrices_for_intervals()
        self.matrices = [x.toarray() for x in self.matrices]


    def read_atlas(self):
        reader = AtlasReader(self.config)
        atlas_intervals, atlas_matrices = reader.meth_cov_to_beta_matrices()
        interval_to_mat = dict(zip([str(x) for x in atlas_intervals], atlas_matrices))
        self.atlas_matrices = [interval_to_mat[str(x)] for x in self.interval_order]

    def load_npy(self):
        self.matrices = list(np.load(self.config["data_file"], allow_pickle=True))
        self.atlas_matrices = np.load(self.config["metadata_file"], allow_pickle=True)

    def deconvolute(self):
        r = celfie_ish(self.matrices, self.atlas_matrices, origins=None, num_iterations=self.config['num_iterations'],
                       convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()

class ReAtlas(CelfieISH):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "CelFiE-ISH ReAtlas"

    def read_atlas(self):
        reader = AtlasReader(self.config)
        atlas_intervals, y, y_depths = reader.meth_cov_to_meth_cov()
        interval_to_atlas = dict(zip([str(x) for x in atlas_intervals], np.arange(len(atlas_intervals))))
        self.y = [y[interval_to_atlas[str(x)]] for x in self.interval_order]
        self.y_depths = [y_depths[interval_to_atlas[str(x)]] for x in self.interval_order]


    def load_npy(self):
        self.matrices = list(np.load(self.config["data_file"], allow_pickle=True))
        self.y, self.y_depths = np.load(self.config["metadata_file"], allow_pickle=True)

    def deconvolute(self):
        r = reatlas(self.matrices, self.y, self.y_depths, num_iterations=self.config['num_iterations'],
                        convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()


class Epistate(CelfieISH):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "Epistate"

    def load_npy(self):
        self.matrices = np.load(self.config["data_file"], allow_pickle=True)
        self.thetaH, self.thetaL, self.lambdas = np.load(self.config["metadata_file"], allow_pickle=True)

    def read_atlas(self):
        #load lambdas and thetas
        reader = EpiAtlasReader(self.config)
        self.lambda_intervals, self.lambdas = reader.read_lambdas()
        self.theta_intervals, self.thetaH, self.thetaL = reader.read_thetas()
        self.sort_intervals()

    def sort_intervals(self):
        '''
        make sure everything is in the same order
        mostly useful if data missing in atlas
        :return:
        '''
        lambda_order = dict(zip([str(x) for x in self.lambda_intervals], self.lambdas))
        thetaH_order = dict(zip([str(x) for x in self.theta_intervals], self.thetaH))
        thetaL_order = dict(zip([str(x) for x in self.theta_intervals], self.thetaL))
        interval_order, lambdas, thetaH, thetaL, cpgs, matrices, origins = [], [], [], [],[],[], []
        for i, interval in enumerate(self.interval_order):
            if str(interval) in lambda_order: #exists in atlas
                interval_order.append(interval)
                lambdas.append(lambda_order[str(interval)])
                thetaH.append(thetaH_order[str(interval)])
                thetaL.append(thetaL_order[str(interval)])
                cpgs.append(self.cpgs[i])
                origins.append(self.origins[i])
                matrices.append(self.matrices[i])
        self.interval_order, self.lambdas, self.thetaH, self.thetaL, self.cpgs, self.matrices, self.origins =\
        interval_order, lambdas, thetaH, thetaL, cpgs, matrices, origins


    def deconvolute(self):
        r = epistate(self.matrices, self.lambdas, self.thetaH, self.thetaL, origins=None,
                     num_iterations=self.config["num_iterations"],
                     convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.em()
        self.runner = r

#%%

@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option('--model',
              type=click.Choice(['uxm','celfie','sum-celfie', 'celfie-ish','reatlas', 'epistate'], case_sensitive=False))
@click.option('-m', '--mixture', help='mixture file to deconvolute')
@click.option('-a', '--atlas_file', help='atlas with beta values')
@click.option('--minimal_cpg_per_read', type=int, help='only reads with at least n cpgs will be considered', default=1)
@click.option('-j', '--json', help='run from json config file')
@click.option('--cpg_coordinates', help='sorted cpg bed file')
@click.option('--outfile', help='output file path')
@click.option('-j', '--json', help='run from json config file')
@click.option('-i', '--genomic_intervals', help='interval(s) to process. formatted chrN:start-end, separated by commas')
@click.option('-b', '--bedfile', help='the intervals are in a bedfile',
              is_flag=True, default=False)
@click.option('--header', is_flag=True, default=False, help="bedgraph with regions to process has header")
@click.option('-A', '--coords', is_flag=True, help='epiread files contain coords', default=False)
@click.option('--epiformat',
              type=click.Choice(['old_epiread','old_epiread_A','pat'], case_sensitive=False))
@click.option('--num_iterations', type=int,  help="maximal iterations")
@click.option('--stop_criterion', type=float,  help="minimal improvement required to continue")
@click.option('--random_restarts', type=int,  help="number of initializations (only one returned)")
@click.option('--data_file', help='mixture file (for simulated data only)')
@click.option('--metadata_file', help='atlas file (for simulated data only)')
@click.option('--lambdas', help='lambda estimates per region (specific to epistate)')
@click.option('--thetas', help='theta estimates per region (specific to epistate)')
@click.option('--percent_u', help='atlas file with %U values (specific to UXM)')
@click.option('--weights', help='weights per marker region (specific to UXM)')
@click.option('--u_threshold',type=float, help='maximal methylation to be considered U (specific to UXM)')
@click.option('--min_length',type=int, help='only reads with at least n cpgs will be considered (specific to UXM). same as minimal_cpg_per_read but applied at the deconvolution level')


@click.option('-s', '--summing', help='sum each marker region (CelFiE sum)',
              is_flag=True, default=False)
@click.option('--npy', help='output .npy file instead of text (for pipelines)',
              is_flag=True, default=False)
@click.version_option()
@click.pass_context
def main(ctx, **kwargs):
    """deconvolute epiread file using atlas"""
    config = {}
    config.update(kwargs)
    config.update(dict([item.strip('--').split('=') for item in ctx.args]))
    if kwargs["json"] is not None:
        with open(kwargs["json"], "r") as jconfig:
            config.update(json.load(jconfig))

    if "epiread_files" not in config:
        config["epiread_files"] = [config["mixture"]]

    if config["model"]=='celfie':
        model=Celfie
        config["summing"]=False
    elif config["model"]=="sum-celfie":
        model=Celfie
        config["summing"]=True
    elif config["model"]=='celfie-ish':
        model = CelfieISH
    elif config["model"]=='reatlas':
        model = ReAtlas
    elif config["model"]=='uxm':
        model = UXM
    else: # config["model"]=='epistate':
        model=Epistate
    em_model = model(config)
    if config["data_file"] is not None and config["metadata_file"] is not None:
        em_model.run_from_npy()
    else:
        em_model.run_model()

if __name__ == '__main__':
    main()

#%%
# import os
# os.chdir("/Users/ireneu/PycharmProjects/deconvolution_models")
# config = {"cpg_coordinates": "demo/hg19.CpG.bed.sorted.gz", "bedfile":True,
#           "genomic_intervals":"demo/U250.tsv",
#           "outfile":"/Users/ireneu/berman_lab/ALS/test.bedgraph",
#           "epiformat":"old_epiread_A", "header":False, "epiread_files":["demo/mixture.epiread.gz"],
#           "atlas_file": "demo/beta_atlas.txt",
#           "data_file":"/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/2_rep25_data.npy",
#           "metadata_file":"/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/2_rep25_metadata_reatlas.npy",
#   "num_iterations": 10, "stop_criterion": 1e-07, "random_restarts": 1, "summing":False,
#           "min_length":4, "u_threshold":0.25,
#
#           }

# config = {"cpg_coordinates": "demo/hg19.CpG.bed.sorted.gz", "bedfile":True,
#           "genomic_intervals":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/sensitivity_200723_U250_merged_regions_file.bed",
#           "outfile":"/Users/ireneu/berman_lab/ALS/test.bedgraph",
#           "epiformat":"old_epiread_A", "header":False, "epiread_files":["/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/sensitivity_200723_U250_4_rep15_mixture.epiread.gz"],
#           "atlas_file": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/sensitivity_200723_U250_atlas_over_regions.txt",
#             "percent_u": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/sensitivity_200723_U250_percent_U.bedgraph",
#   "num_iterations": 10, "stop_criterion": 1e-05, "random_restarts": 1, "summing":False,
#           "min_length":1, "u_threshold":0.25, "npy":False, "weights":False, "minimal_cpg_per_read":1
#           }
#
# config = {"bedfile": True, "header": False, "cpg_coordinates": "tests/data/hg38_pat_cpg_from_netanel.bed.gz",
#         "npy": True, "depth": 4.5, "num_iterations": 30000, "random_restarts": 1,
#           "true_alpha": "[0.00201613,0.00403226,0.00604839,0.00806452,0.01008065,0.01209677,0.0141129 ,0.01612903,0.01814516,0.02016129,0.02217742,0.02419355,0.02620968,0.02822581,0.03024194,0.03225806,0.03427419,0.03629032,0.03830645,0.04032258,0.04233871,0.04435484,0.04637097,0.0483871 ,0.05040323,0.05241935,0.05443548,0.05645161,0.05846774,0.06048387,0.0625]",
#           "stop_criterion": 1e-07, "min_length": 4, "u_threshold": 0.25,
#           "epiread_files": ["tests/data/HU.10.filtered.pat.gz"],
#           "epiformat": "pat",
#           "atlas_file": "tests/data/Filippo_atlas_over_regions.txt",
#           "genomic_intervals": "tests/data//Filippo_merged_regions_file.bed",
#           "cell_types": ["Adipocytes", "Endothel", "Bladder-Ep", "Blood-B", "Blood-Granul", "Blood-Mono+Macro",
#                          "Blood-NK", "Blood-T", "Eryth-prog", "Breast-Basal-Ep", "Breast-Luminal-Ep", "Neuron",
#                          "Oligodend", "Head-Neck-Ep", "Gastric-Ep", "Small-Int-Ep", "Colon-Ep", "Heart-Cardio",
#                          "Fallopian-Ep", "Kidney-Ep", "Liver-Hep", "Lung-Ep-Alveo", "Lung-Ep-Bron", "Ovary-Ep",
#                          "Pancreas-Acinar", "Pancreas-Duct", "Pancreas-Alpha", "Pancreas-Beta", "Pancreas-Delta",
#                          "Prostate-Ep", "Thyroid-Ep", "Megakaryocyte"],
#           "lambdas": "", "percent_u": "tests/data/atlas_U1000_32cellTypes_hg38_for_irene.tsv", "weights": False,
#           "summing": False,"thetas": ""}
#
# em_model = Celfie(config)
# em_model.run_model()
