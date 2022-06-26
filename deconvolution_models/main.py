import click
from parsing import MixtureReader, AtlasReader
from celfie import em as celfie
from celfie_plus import READMeth as celfie_plus
from epistate import READMeth as epistate
from epistate_plus import READMeth as epistate_plus
import numpy as np

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
    def __init__(self, mixture_file, atlas_file, regions, CpG_coordinates, outfile, epiformat,
                 max_iterations, random_restarts, stop_criterion, header):
        self.mixture = mixture_file
        self.atlas = atlas_file
        self.regions = regions
        self.cpg_file = CpG_coordinates
        self.outfile = outfile
        self.epiformat = epiformat
        self.max_iterations = max_iterations
        self.random_restarts = random_restarts
        self.stop_criterion = stop_criterion
        self.header=header
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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def read_mixture(self):
        reader = MixtureReader(self.regions, self.mixture, self.cpg_file)
        self.x, self.x_depths = reader.sum_regions()

    def read_atlas(self):
        reader = AtlasReader(self.regions, self.mixture, self.cpg_file, header=self.header)
        self.y, self.y_depths = reader.sum_regions()

    def deconvolute(self):
        restarts = []
        for r in range(self.random_restarts):
            estimated_alpha, estimated_gamma, ll, i = celfie(
                self.x, self.x_depths, self.y, self.y_depths, self.max_iterations, self.stop_criterion
            )
            restarts.append((ll, estimated_alpha, i))
        ll_max, alpha_max, i_max = max(restarts)
        self.alpha, self.i = alpha_max.flatten(), i_max


class CelfiePlus(EMmodel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "celfie-plus"
        if self.random_restarts > 1:
            raise NotImplementedError("random_restarts", self.name)

    def read_mixture(self):
        reader = MixtureReader(self.regions, self.mixture, self.cpg_file)
        self.mixture_matrices = reader.get_matrices()

    def read_atlas(self):
        reader = AtlasReader(self.regions, self.mixture, self.cpg_file, header=self.header)
        self.atlas_matrics = reader.get_matrices()

    def deconvolute(self):
        assert [x.shape[1] for x in self.mixture_matrices] == [x.shape[1] for x in self.atlas_matrics]
        r = celfie_plus(self.mixture_matrices, self.atlas_matrics, num_iterations=self.max_iterations,
                        convergence_criteria=self.stop_criterion)
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
mixture = "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_6_rep9_mixture.epiread.gz"
atlas = "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_atlas_over_tims.txt"
regions = "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_processed_tims.txt"
cpg_coordinates = "/Users/ireneu/PycharmProjects/in-silico_deconvolution/debugging/hg19.CpG.bed.sorted.gz"
outfile="sample_output.something"
epiformat = "old_epiread"
em_model = Celfie(mixture,atlas,regions, cpg_coordinates, outfile,epiformat,
                  max_iterations=1000,random_restarts=1,stop_criterion=0.001,
                  header=True)
em_model.run_model()