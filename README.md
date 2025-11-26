
# Package: deconvolution_models

This package contains models for methylation deconvolution.

## Installation

To install the package, you can use `pip` with the following commands:

Basic install:
```shell
pip install git+https://github.com/methylgrammarlab/deconvolution_models
```

If you encounter any issues with the basic install, you may need to use a [personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) as follows:
```shell
pip install git+https://<PERSONAL ACCESS TOKEN>@github.com/methylgrammarlab/deconvolution_models.git
```

Please note that this package requires [epiread-tools](https://github.com/methylgrammarlab/epiread-tools) to run, which should be downloaded automatically during the installation process.

## Reference files
A record of this repository along with reference files in epiread and pat formats can be found at: https://zenodo.org/records/10799835

## Getting Started

To view the full list of input arguments and usage instructions, run the following command:
```shell
deconvolution --help
```

- `--model`: Specify the deconvolution model to use. Available options are 'uxm', 'celfie', 'sum-celfie', 'celfie-ish', 'reatlas', and 'epistate'.
- `-m`, `--mixture`: Path to the mixture file to deconvolute.
- `-a`, `--atlas_file`: Path to the atlas file containing methylation and coverage for each cell type.
- `--minimal_cpg_per_read`: Set the minimum number of CpGs required for a read to be considered. Default is 1.
- `-j`, `--json`: Run the deconvolution using a JSON config file.
- `--cpg_coordinates`: Path to the sorted CpG bed file. Should contain coordinates for all CpG sites in bed format. 
- `--outfile`: Path to the output file (to be generated).
- `-i`, `--genomic_intervals`: Specify the genomic intervals to process. Use the format chrN:start-end, separated by commas.
- `-b`, `--bedfile`: The intervals (-i) are in a bedfile, rather than a comma delimited list.
- `--header`: The bedgraph file with regions to process has a header.
- `--epiformat`: Specify the format of the epiread files. Available options are 'old_epiread', 'old_epiread_A', and 'pat'.
- `--num_iterations`: Set the maximum number of iterations.
- `--stop_criterion`: Set the minimal improvement required to continue deconvolution. 
- `--random_restarts`: Set the number of initializations (only one will be returned).
- `--data_file`: Specify the mixture file (for simulated data only).
- `--metadata_file`: Specify the atlas file (for simulated data only).
- `--lambdas`: Specify the lambda estimates per region (specific to epistate).
- `--thetas`: Specify the theta estimates per region (specific to epistate).
- `--percent_u`: Specify the atlas file with %U values (specific to UXM).
- `--weights`: Specify the weights per marker region (specific to UXM).
- `--u_threshold`: Set the maximal methylation value to be considered as U (specific to UXM).
- `--min_length`: Set the minimum number of CpGs required for a read to be considered at the deconvolution level (specific to UXM). Same as `--minimal_cpg_per_read` but applied at the deconvolution level.
- `-s`, `--summing`: Perform summing for each marker region (CelFiE sum).


Alternatively, if you want to deconvolute using a JSON configuration file, you can use the following command:
```shell
deconvolution -j <config.json>
```
JSON will override command line arguments provided. 

### Example Code

Below are example commands for running deconvolution using different models. All the necessary files for the example code can be found in the "demo" directory.

#### CelFiE-ISH:
```shell
deconvolution --model celfie-ish -a demo/beta_atlas.txt -m demo/mixture.epiread.gz \
--cpg_coordinates demo/sample_cpg_file.bed.gz -i demo/U250.tsv -b \
--epiformat old_epiread_A --num_iterations 10000 --stop_criterion 0.0000001 \
--random_restarts 1 --outfile demo/sample_output.txt
```

#### CelFiE-ISH ReAtlas:
```shell
deconvolution --model reatlas -a demo/beta_atlas.txt -m demo/mixture.epiread.gz \
--cpg_coordinates demo/sample_cpg_file.bed.gz -i demo/U250.tsv -b \
--epiformat old_epiread_A --num_iterations 10000 --stop_criterion 0.0000001 \
--random_restarts 1 --outfile demo/sample_output.txt
```

#### Epistate:
```shell
deconvolution --model epistate --lambdas demo/lambdas.txt --thetas demo/thetas.txt \
-m demo/mixture.epiread.gz --cpg_coordinates demo/sample_cpg_file.bed.gz \
-i demo/U250.tsv -b --epiformat old_epiread_A --random_restarts 1 \
--num_iterations 10000 --stop_criterion 0.0000001 --outfile demo/sample_output.txt
```
The Epistate atlas files can be generated with the [bimodal-detector](https://github.com/methylgrammarlab/bimodal_detector) package.

#### Comparing to other deconvolution models:

To compare with the [CelFiE](https://github.com/christacaggiano/celfie) model, you can use the following command:
```shell
deconvolution --model celfie -a demo/beta_atlas.txt -m demo/mixture.epiread.gz \
--cpg_coordinates demo/sample_cpg_file.bed.gz -i demo/U250.tsv -b \
--epiformat old_epiread_A --num_iterations 10000 --stop_criterion 0.0000001 \
--random_restarts 1 --outfile demo/sample_output.txt
```

To compare with the [UXM](https://github.com/nloyfer/UXM_deconv) model, you can use the following command:
```shell
deconvolution --model uxm --percent_u demo/U_atlas.txt -m demo/mixture.epiread.gz \
--cpg_coordinates demo/sample_cpg_file.bed.gz -i demo/U250.tsv -b \
--epiformat old_epiread_A --min_length 5 --u_threshold 0.1 \
--random_restarts 1 --outfile demo/sample_output.txt
```
## Cite the paper
Unterman, I., Avrahami, D., Katsman, E. Triche, T.J.Jr., Glaser, B. & Berman, B.P. CelFiE-ISH: a probabilistic model for multi-cell type deconvolution from single-molecule DNA methylation haplotypes. Genome Biol 25, 151 (2024). https://doi.org/10.1186/s13059-024-03275-x


