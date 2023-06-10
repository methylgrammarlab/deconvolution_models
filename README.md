
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

## Getting Started

To view the full list of input arguments and usage instructions, run the following command:
```shell
deconvolution --help
```

Alternatively, if you want to deconvolute using a JSON configuration file, you can use the following command:
```shell
deconvolution -j <config.json>
```

### Example Code

Below are example commands for running deconvolution using different models. All the necessary files for the example code can be found in the "demo" directory.

#### CelFiE-ISH:
```shell
deconvolution --model celfie-ish -a demo/beta_atlas.txt -m demo/mixture.epiread.gz \
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



