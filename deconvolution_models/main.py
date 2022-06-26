import click
from parsing import MixtureReader, AtlasReader

class EMmodel:
    def __init__(self, mixture_file, atlas_file,outfile, max_iterations, random_restarts, stop_criterion):
        pass

    def read_mixture(self):
        pass

    def read_atlas(self):
        pass

    def write_output(self):
        pass

    def deconvolute(self):
        pass

    def run_model(self):
        self.read_mixture()
        self.read_atlas()
        self.deconvolute()
        self.write_output()

class Celfie(EMmodel):

    def __init__(self, mixture_file, atlas_file,outfile, max_iterations, random_restarts):
        super().__init__(mixture_file, atlas_file,outfile, max_iterations, random_restarts)

    def read_mixture(self):
        reader = MixtureReader(regions, mixture, CpG_file, outfile)
        self.x, self.x_depths = reader.sum_regions()

    def read_atlas(self):
        reader = AtlasReader(tims, atlas_file, cpg_file, outfile)
        self.y, self.y_depths = reader.sum_regions()

    def deconvolute(self):
        restarts = []
        for r in range(self.random_restarts):
            estimated_alpha, estimated_gamma, ll, i = em(
                self.x, self.x_depths, self.y, self.y_depths, self.max_iterations, self.stop_criterion
            )
            restarts.append((ll, estimated_alpha, i))
        ll_max, alpha_max, i_max = max(restarts)
        alpha_max = alpha_max.flatten()

    def write_output(self):
        np.save(args.outpath, [alpha_max, np.array([i_max])], allow_pickle=True)

class CelfiePlus(EMmodel):

    def __init__(self, mixture_file, atlas_file,outfile, max_iterations, random_restarts):
        super().__init__(mixture_file, atlas_file,outfile, max_iterations, random_restarts)
        assert [x.shape[1] for x in reads] == [x.shape[1] for x in atlas]

    def read_mixture(self):
        reader = MixtureReader(regions, mixture, CpG_file, outfile)
        self.mixture_matrics = reader.get_matrices()

    def read_atlas(self):
        reader = AtlasReader(tims, atlas_file, cpg_file, outfile)
        self.atlas_matrics = reader.get_matrices()

    def deconvolute(self):
        r = READMe(self.mixture_matrics, self.atlas_matrics, num_iterations=self.max_iterations,
                   convergence_criteria=0.001)
        estimated_alpha, i = r.two_step()

class Epistate(EMmodel):
    def __init__(self, mixture_file, atlas_file,outfile, max_iterations, random_restarts):
        super().__init__(mixture_file, atlas_file,outfile, max_iterations, random_restarts)
        assert [x.shape[1] for x in reads] == [x.shape[1] for x in atlas]

    def deconvolute(self):
        r = READMe(reads, lambdas, thetaH, thetaL, num_iterations=args.num_iterations,
                   convergence_criteria=0.001)
        estimated_alpha, i = r.two_step()

class EpistatePlus(Epistate):
    def __init__(self, mixture_file, atlas_file, outfile, max_iterations, random_restarts):
        super().__init__(mixture_file, atlas_file, outfile, max_iterations, random_restarts)
        assert [x.shape[1] for x in reads] == [x.shape[1] for x in atlas]



@click.command()
@click.argument('cpg_coordinates')
@click.argument('epireads')
@click.argument('outfile')
@click.option('-i', '--intervals', help='interval(s) to process. formatted chrN:start-end, separated by commas')
@click.option('-b', '--bedfile', help='bed file chrom start end with interval(s) to process. tab delimited',
              is_flag=True, default=False)
@click.option('--header', is_flag=True, help="bedgraph with regions to process has header")
@click.option('-A', '--coords', is_flag=True, help='epiread files contain coords', default=False)
@click.option('--model',
              type=click.Choice(['celfie', 'celfie-plus', 'epistate', 'epistate-plus'], case_sensitive=False))
def main(cpg_coordinates, epireads, outfile ,*args, **kwargs):
    """deconvolute epiread file using atlas"""
    if kwargs["model"]=='celfie':
        model=Celfie
    elif kwargs["model"]=='celfie-plus':
        model = CelfiePlus
    elif kwargs["model"]=='epistate':
        model=Epistate
    else:
        model=EpistatePlus
    model(mixture_file, atlas_file,outfile, max_iterations, random_restarts, stop_criterion)
    model.run_model()

if __name__ == '__main__':
    main()