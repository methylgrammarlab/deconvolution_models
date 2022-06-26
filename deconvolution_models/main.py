import click

#TODO:
#read user args
#parse mixture and atlas file for model
#run model


@click.command()
@click.argument('cpg_coordinates')
@click.argument('epireads')
@click.argument('outfile')
@click.option('-i', '--intervals', help='interval(s) to process. formatted chrN:start-end, separated by commas')
@click.option('-b', '--bedfile', help='bed file chrom start end with interval(s) to process. tab delimited',
              is_flag=True, default=False)
@click.option('--header', is_flag=True, help="bedgraph with regions to process has header")
@click.option('-A', '--coords', is_flag=True, help='epiread files contain coords', default=False)
def main(cpg_coordinates, epireads, outfile ,*args, **kwargs):
    """ biscuit epiread to bedgraph converter"""
    epiread_files = epireads.split(",")

    if kwargs["intervals"]:
        genomic_intervals = kwargs["intervals"].split(",")
    elif kwargs["bedfile"]:
        genomic_intervals = kwargs["bedfile"]
    else:
        raise ValueError("either specify intervals or add bed file. For whole genome use -b with chrom sizes")
    epiformat = "old_epiread"
    if kwargs["coords"]:
        epiformat = "old_epiread_A"
    runner = EpiRunner(genomic_intervals, cpg_coordinates, epiread_files, outfile, epiformat,
                       kwargs["header"], kwargs["bedfile"])
    runner.tobedgraph()


if __name__ == '__main__':
    main()