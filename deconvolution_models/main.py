import click
import json
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
# sys.path.append("/Users/ireneu/PycharmProjects/deconvolution_models/deconvolution_models")
sys.path.append("/Users/ireneu/PycharmProjects/deconvolution_models/tests")
from deconvolution_models.celfie import em as celfie
from deconvolution_models.celfie_plus import CelfiePlus as celfie_plus
from deconvolution_models.celfie_plus_reatlas import CelfiePlus as reatlas
from deconvolution_models.epistate import READMeth as epistate
from deconvolution_models.epistate_plus import READMeth as epistate_plus
# from epistate_plus_simplified import READMeth as epistate_plus
import numpy as np
from epiread_tools.epiparser import EpireadReader, CoordsEpiread, epiformat_to_reader,AtlasReader, EpiAtlasReader
from epiread_tools.naming_conventions import *
from epiread_tools.em_utils import calc_coverage, calc_methylated


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
        self.interval_order, self.matrices, self.cpgs = reader.get_matrices_for_intervals()
        self.methylation, self.coverage = calc_methylated(self.matrices), calc_coverage(self.matrices)
        self.x, self.x_depths = np.hstack(self.methylation), np.hstack(self.coverage) #no summing
        if self.config["summing"]:
            self.x, self.x_depths = np.array([np.sum(x) for x in self.methylation]), \
                                    np.array([np.sum(x) for x in self.coverage])
    def read_atlas(self):
        reader = AtlasReader(self.config)
        self.y, self.y_depths = reader.meth_cov_to_meth_cov() #no summing
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


class CelfiePlus(EMmodel):

    def __init__(self, config):
        super().__init__(config)
        self.name = "celfie-plus"
        if self.config['random_restarts'] > 1:
            raise NotImplementedError("random_restarts", self.name)

    def read_mixture(self):
        reader = self.reader(self.config)
        self.interval_order, self.matrices, self.cpgs = reader.get_matrices_for_intervals()
        self.matrices = [x.toarray() for x in self.matrices]

    def read_atlas(self):
        reader = AtlasReader(self.config)
        self.atlas_matrices = reader.meth_cov_to_beta_matrices()

    def load_npy(self):
        self.matrices = list(np.load(self.config["data_file"], allow_pickle=True))
        self.atlas_matrices = np.load(self.config["metadata_file"], allow_pickle=True)

    def deconvolute(self):
        r = celfie_plus(self.matrices, self.atlas_matrices, num_iterations=self.config['num_iterations'],
                        convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()

class ReAtlas(CelfiePlus):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "plus-reatlas"

    def read_atlas(self):
        reader = AtlasReader(self.config)
        self.y, self.y_depths = reader.meth_cov_to_meth_cov()


    def load_npy(self):
        self.matrices = list(np.load(self.config["data_file"], allow_pickle=True))
        self.y, self.y_depths = np.load(self.config["metadata_file"], allow_pickle=True)

    def deconvolute(self):
        r = reatlas(self.matrices, self.y, self.y_depths, num_iterations=self.config['num_iterations'],
                        convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()


class Epistate(CelfiePlus): #TODO: load lambdas and thetas from file
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "epistate"

    def load_npy(self):
        self.matrices = np.load(self.config["data_file"], allow_pickle=True)
        self.thetaH, self.thetaL, self.lambdas = np.load(self.config["metadata_file"], allow_pickle=True)

    def read_atlas(self):
        #load lambdas and thetas
        reader = EpiAtlasReader(self.config)
        self.lambdas = reader.read_lambdas()
        self.thetaH, self.thetaL = reader.read_thetas()

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
        r = epistate_plus(self.matrices, self.lambdas, self.thetaH, self.thetaL,
                     num_iterations=self.config["num_iterations"],
                   convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.em()
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
#           "stop_criterion":0.001, "summing":True,
#           "epiread_files":["/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/all_models_linear_6_rep9_mixture.epiread.gz"],
#           "atlas_file":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/all_models_linear_atlas_over_regions.txt",
#           "lambdas":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/all_models_linear_lambdas.bedgraph",
#           "thetas":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/all_models_linear_thetas.bedgraph",
#           "genomic_intervals":"/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/all_models_linear_merged_regions_file.bed",
#           "epiformat":"old_epiread_A", "slop":0}
#
# config = {"data_file":"/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/18_rep5_data.npy",
#           "metadata_file":"/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/18_rep5_metadata_celfie.npy",
#                     "num_iterations": 1000, "random_restarts": 1,
#                     "stop_criterion":0.001, "summing":False
#           }
# r = EpistatePlus(config)
# r.run_model()

# true_alpha = np.array([0.00307692, 0.00615385, 0.00923077, 0.01230769, 0.01538462,
#        0.01846154, 0.02153846, 0.02461538, 0.02769231, 0.03076923,
#        0.03384615, 0.03692308, 0.04      , 0.04307692, 0.04615385,
#        0.04923077, 0.05230769, 0.05538462, 0.05846154, 0.06153846,
#        0.06461538, 0.06769231, 0.07076923, 0.07384615, 0.07692308])
# [self.q[i] >= self.q[i-1] for i in range(1, len(self.q))]