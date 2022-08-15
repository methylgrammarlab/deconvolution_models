import click
import json
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
# sys.path.append("/Users/ireneu/PycharmProjects/deconvolution_models/deconvolution_models")
from deconvolution_models.celfie import em as celfie
from deconvolution_models.celfie_plus import CelfiePlus as celfie_plus
from deconvolution_models.epistate import READMeth as epistate
from deconvolution_models.epistate_plus import READMeth as epistate_plus
import numpy as np
from epiread_tools.epiparser import EpireadReader, CoordsEpiread, epiformat_to_reader,AtlasReader
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
        self.x, self.x_depths = np.array([np.sum(x) for x in self.methylation]),\
                                np.array([np.sum(x) for x in self.coverage])

    def read_atlas(self):
        reader = AtlasReader(self.config)
        self.y, self.y_depths = reader.meth_cov_to_sum()
        self.y, self.y_depths = self.y.T, self.y_depths.T

    def load_npy(self):
        self.matrices = np.load(self.config["data_file"], allow_pickle=True)
        self.y, self.y_depths = np.load(self.config["metadata_file"], allow_pickle=True)
        self.methylation, self.coverage = calc_methylated(self.matrices), calc_coverage(self.matrices)
        # self.x, self.x_depths = np.array([np.sum(x) for x in self.methylation]),\
        #                         np.array([np.sum(x) for x in self.coverage])
        self.x, self.x_depths = np.hstack(self.methylation),np.hstack(self.coverage)

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
        self.y, self.y_depths = reader.meth_cov_to_meth_cov()


    def load_npy(self): #broken
        self.matrices = list(np.load(self.config["data_file"], allow_pickle=True))
        self.atlas_matrices = np.load(self.config["metadata_file"], allow_pickle=True)
        self.y = [x*100 for x in self.atlas_matrices]
        self.y_depths = [np.ones(x.shape) for x in self.y]
        [x.fill(100) for x in self.y_depths]

    def deconvolute(self):
        r = celfie_plus(self.matrices, self.y, self.y_depths, num_iterations=self.config['num_iterations'],
                        convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()

class Epistate(CelfiePlus): #TODO: load lambdas and thetas from file
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "epistate"

    def load_npy(self):
        self.matrices = np.load(self.config["data_file"], allow_pickle=True)
        self.thetaH, self.thetaL, self.lambdas = np.load(self.config["metadata_file"], allow_pickle=True)

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
              type=click.Choice(['celfie', 'celfie-plus', 'epistate', 'epistate-plus'], case_sensitive=False))
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
    elif config["model"]=='celfie-plus':
        model = CelfiePlus
    elif config["model"]=='epistate':
        model=Epistate
    else:
        model=EpistatePlus

    em_model = model(config)
    em_model.run_model()

# if __name__ == '__main__':
#     main()

#%%
#
config = {"bedfile": True, "header": False, "cpg_coordinates": "/Users/ireneu/PycharmProjects/old_in-silico_deconvolution/debugging/hg19.CpG.bed.sorted.gz",
          "depth": 30, "num_iterations": 1000, "random_restarts": 1,
          # "epiread_files": ["/Users/ireneu/PycharmProjects/deconvolution_in_silico_pipeline/data/mixtures/acin_endo_5_rep9_mixture.epiread.gz"],
          # "atlas_file": "/Users/ireneu/PycharmProjects/deconvolution_in_silico_pipeline/data/mixtures/acin_endo_atlas_over_regions.txt",
          "stop_criterion":0.0001,
          "genomic_intervals": "/Users/ireneu/PycharmProjects/deconvolution_in_silico_pipeline/data/netanel_regions.txt", "epiformat":"old_epiread",
          "data_file": "/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/debugging/test_rep0_data.npy",
"metadata_file":"/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/debugging/test_rep0_metadata_celfie.npy"}

r = Celfie(config)
r.run_from_npy()


