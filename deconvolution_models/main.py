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
from deconvolution_models.celfie_plus import CelfiePlus as celfie_plus
from deconvolution_models.celfie_plus_reatlas import CelfiePlus as reatlas
from deconvolution_models.epistate import CelfiePlus as epistate
from deconvolution_models.epistate_plus import READMeth as epistate_plus
from deconvolution_models.UXM import uxm
# from epistate_plus_simplified import READMeth as epistate_plus
import numpy as np
from epiread_tools.epiparser import EpireadReader, CoordsEpiread, epiformat_to_reader,AtlasReader, EpiAtlasReader, UXMAtlasReader
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


class EMmodel:
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

    def write_output(self):
        np.save(self.config['outfile'], [self.alpha, np.array([self.i])], allow_pickle=True)

    def deconvolute(self):
        raise NotImplementedError("deconvolute", self.name)

    def run_from_npy(self):
        self.load_npy()
        self.deconvolute()
        self.write_output()

    def run_model(self):
        self.read_mixture()
        self.read_atlas()
        self.deconvolute()
        self.write_output()

class Celfie(EMmodel):

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

class UXM(EMmodel):
    def __init__(self, config):
        super().__init__(config)
        self.U = [] #percent u
        self.N = [] #number of fragments = weights
        self.i = 0
        self.min_length = self.config["min_length"]

    def read_mixture(self):
        reader = self.reader(self.config)
        self.interval_order, self.matrices, self.cpgs, self.origins = reader.get_matrices_for_intervals()
        self.calc_u()

    def calc_u(self):
        for mat in self.matrices:
            x_c_v = np.array(mat.todense() != NOVAL)
            # filter short reads
            len_filt = (np.sum(x_c_v, axis=1) >= self.min_length).ravel()
            if not np.sum(len_filt):  # empty region
                self.U.append(0)
                self.N.append(0)
            else:
                small = mat[len_filt, :]
                self.U.append(calc_percent_U(small, self.config["u_threshold"]))
                self.N.append(small.shape[0])

    def read_atlas(self):
        reader = UXMAtlasReader(self.config)
        self.atlas_intervals, self.atlas_mat = reader.read()
        self.sort_intervals()

    def sort_intervals(self):
        atlas_to_index = {str(v): k for k, v in dict(enumerate(self.atlas_intervals)).items()}
        atlas_order = np.array([atlas_to_index[str(x)] for x in self.interval_order])
        self.atlas_mat = self.atlas_mat[atlas_order,:]
        self.atlas_intervals = self.atlas_intervals[atlas_order]

    def load_npy(self):
        self.matrices = list(np.load(self.config["data_file"], allow_pickle=True))
        self.calc_u()
        self.atlas = np.load(self.config["metadata_file"], allow_pickle=True)

    def deconvolute(self):
        if self.config["weights"]:
            self.alpha = uxm(np.vstack(self.atlas), self.U, self.N)
        else:
            self.alpha = uxm(np.vstack(self.atlas), self.U)


class CelfiePlus(EMmodel):

    def __init__(self, config):
        super().__init__(config)
        self.name = "celfie-plus"
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
        r = celfie_plus(self.matrices, self.atlas_matrices, origins=None, num_iterations=self.config['num_iterations'],
                        convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()

class ReAtlas(CelfiePlus):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "plus-reatlas"

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


class Epistate(CelfiePlus): #TODO: load lambdas and thetas from fil
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "epistate"

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
        r = epistate(self.matrices, self.lambdas, self.thetaH, self.thetaL,
                     num_iterations=self.config["num_iterations"],
                   convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()

class EpistatePlus(Epistate):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "epistate-plus"

    def deconvolute(self):
        r = epistate_plus(self.matrices, self.lambdas, self.thetaH, self.thetaL,origins=None,
                     num_iterations=self.config["num_iterations"],
                   convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.em()
        self.runner = r
#%%

@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option('--model',
              type=click.Choice(['celfie','sum-celfie', 'celfie-plus','reatlas', 'epistate', 'epistate-plus'], case_sensitive=False))
@click.option('-j', '--json', help='run from json config file')
@click.version_option()
@click.pass_context
def main(ctx, **kwargs):
    """deconvolute epiread file using atlas"""
    with open(kwargs["json"], "r") as jconfig:
        config = json.load(jconfig)
    config.update(kwargs)
    config.update(dict([item.strip('--').split('=') for item in ctx.args]))

    if config["model"]=='celfie':
        model=Celfie
        config["summing"]=False
    elif config["model"]=="sum-celfie":
        model=Celfie
        config["summing"]=True
    elif config["model"]=='celfie-plus':
        model = CelfiePlus
    elif config["model"]=='reatlas':
        model = ReAtlas
    elif config["model"]=='epistate':
        model=Epistate
    else:
        model=EpistatePlus

    em_model = model(config)
    em_model.run_model()

if __name__ == '__main__':
    main()

#%%
#
# config = {"bedfile": True, "header": False, "cpg_coordinates": "/Users/ireneu/PycharmProjects/old_in-silico_deconvolution/debugging/hg19.CpG.bed.sorted.gz",
#           "num_iterations": 1000, "random_restarts": 1,
#           "stop_criterion":0.001, "summing":False,
#           "epiread_files":["/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/271122_U250_pancreatic_under_2_1_rep4_mixture.epiread.gz"],
#           "atlas_file":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/271122_U250_pancreatic_under_2_atlas_over_regions.txt",
#           "lambdas":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/271122_U250_pancreatic_under_2_lambdas.bedgraph",
#           "thetas":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/271122_U250_pancreatic_under_2_thetas.bedgraph",
#           "genomic_intervals":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/271122_U250_pancreatic_under_2_merged_regions_file.bed",
#           "epiformat":"old_epiread_A", "slop":0}
# #
# config = {"bedfile": True, "header": False,"cpg_coordinates": "/Users/ireneu/PycharmProjects/old_in-silico_deconvolution/debugging/hg19.CpG.bed.sorted.gz",
#           "depth": 0.2, "num_iterations": 1000, "random_restarts": 1, "true_alpha": "[0.04761905,0.0952381 ,0.14285714,0.19047619,0.23809524,0.28571429]",
#           "stop_criterion": 0.001, "epiread_files": ["/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25_1_rep46_mixture.epiread.gz"],
#           "epiformat": "old_epiread_A", "atlas_file": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25_atlas_over_regions.txt",
#           "genomic_intervals": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25_merged_regions_file.bed",
#           "lambdas": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25__lambdas.bedgraph",
#           "thetas": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25__thetas.bedgraph",
#           "summing":True}
# #%%
# model = CelfiePlus(config)
# model.run_model()
# import numpy as  np
# a = np.arange(1,26)
# b=(a/np.sum(a))
# config = {"bedfile": True, "header": False,"cpg_coordinates": "/Users/ireneu/PycharmProjects/old_in-silico_deconvolution/debugging/hg19.CpG.bed.sorted.gz",
#           "depth": 10, "num_iterations": 1000, "random_restarts": 1, "true_alpha": np.array2string(b, max_line_width=np.inf, separator=","),
#           "stop_criterion": 0.001, "epiread_files": ["/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25_1_rep46_mixture.epiread.gz"],
#           "epiformat": "old_epiread_A", "atlas_file": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25_atlas_over_regions.txt",
#           "percent_u":"/Users/ireneu/PycharmProjects/bimodal_detector/results/test_uxm_percent_U.bedgraph",
#           "genomic_intervals": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25_merged_regions_file.bed",
#           "lambdas": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25__lambdas.bedgraph",
#           "thetas": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/060223_pancreatic_U25__thetas.bedgraph",
#           "data_file": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/10_rep0_data.npy",
#           "u_threshold":0.25,"min_length":4,"weights":True, "cell_types": list(range(6)),
#           "metadata_file":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/10_rep0_metadata_uxm.npy",
#           "summing":False}
# #%%
# model = UXM(config)
# model.run_model()