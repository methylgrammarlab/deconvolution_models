import click
from celfie import em as celfie
from celfie_plus import CelfiePlus as celfie_plus
import numpy as np
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
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
        self.reader = epiformat_to_reader[self.config["epiformat"]]
        self.name, self.alpha, self.i = None, None, None

    def read_mixture(self):
        raise NotImplementedError("read_mixture", self.name)

    def read_atlas(self):
        raise NotImplementedError("read_atlas", self.name)

    def write_output(self):
        np.save(self.outfile, [self.alpha, np.array([self.i])], allow_pickle=True)

    def deconvolute(self):
        raise NotImplementedError("deconvolute", self.name)

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

    def deconvolute(self):
        restarts = []
        for r in range(self.config['random_restarts']):
            estimated_alpha, estimated_gamma, ll, i = celfie(
                self.x.reshape((1,-1)), self.x_depths.reshape((1,-1)), self.y.T, self.y_depths.T,
                self.config['max_iterations'], self.config['stop_criterion']
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

    def deconvolute(self):
        assert [x.shape[1] for x in self.matrices] == [x.shape[1] for x in self.atlas_matrices]
        r = celfie_plus(self.matrices, self.atlas_matrices, num_iterations=self.config['max_iterations'],
                        convergence_criteria=self.config['stop_criterion'])
        self.alpha, self.i = r.two_step()

class Epistate(CelfiePlus): #TODO: load lambdas and thetas
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "epistate"

    def deconvolute(self):
        # r = epistate(reads, lambdas, thetaH, thetaL, num_iterations=args.num_iterations,
        #            convergence_criteria=self.stop_criterion)
        # self.alpha, self.i = r.two_step()
        pass

class EpistatePlus(CelfiePlus):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "epistate-plus"

    def deconvolute(self):
        # r = epistate_plus(reads, lambdas, thetaH, thetaL, num_iterations=args.num_iterations,
        #            convergence_criteria=self.stop_criterion)
        # estimated_alpha, i = r.em()
        pass
#%%

# @click.command()
# @click.argument('mixture')
# @click.argument('atlas')
# @click.argument('regions')
# @click.argument('cpg_coordinates')
# @click.argument('outfile')
# @click.option('--header', is_flag=True, help="bedgraph with regions to process has header")
# @click.option('-A', '--coords', is_flag=True, help='epiread files contain coords', default=False)
# @click.option('--model',
#               type=click.Choice(['celfie', 'celfie-plus', 'epistate', 'epistate-plus'], case_sensitive=False))
# @click.option('--max_iterations', default=1000)
# @click.option('--random_restarts', default=1)
# @click.option('--stop_criterion', default=0.001)
# @click.version_option()
# def main(mixture,atlas,regions, cpg_coorfinates, outfile, **kwargs):
#     """deconvolute epiread file using atlas"""
#     if kwargs["model"]=='celfie':
#         model=Celfie
#     elif kwargs["model"]=='celfie-plus':
#         model = CelfiePlus
#     elif kwargs["model"]=='epistate':
#         model=Epistate
#     else:
#         model=EpistatePlus
#     epiformat = "old_epiread"
#     if kwargs["coords"]:
#         epiformat = "old_epiread_A"
#
#     em_model = model(mixture,atlas,regions, cpg_coorfinates, outfile,epiformat,**kwargs)
#     em_model.run_model()
#
# if __name__ == '__main__':
#     main()

#%%

config = {"epiread_files": ["/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_6_rep9_mixture.epiread.gz"],
          "atlas_file": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_atlas_over_tims.txt",
          "genomic_intervals": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_processed_tims.txt",
          "cpg_coordinates":  "/Users/ireneu/PycharmProjects/in-silico_deconvolution/debugging/hg19.CpG.bed.sorted.gz",
           "outfile":"sample_output.something",
            "epiformat" : "old_epiread",
          "bedfile":True,
          "header":True,
          "max_iterations":1000,
          "random_restarts":1,
          "stop_criterion":0.001
          }
em_model = CelfiePlus(config)
em_model.run_model()