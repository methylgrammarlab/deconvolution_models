# Package: deconvolution_models

Contains models for methylation deconvolution


## Usage

basic install
```
pip install git+https://github.com/methylgrammarlab/deconvolution_models
```
alternatively, you may need to use an [access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
```
pip install git+https://<PERSONAL ACCESS TOKEN>@github.com/methylgrammarlab/deconvolution_models.git
```

to deconvolute:
```
deconvolution -j <config.json>

```
#### Configuration Parameters

The following parameters can be specified in the `config` dictionary or overwritten with specific arguments:
1. Any and all parameters for epiread parsing (see https://github.com/methylgrammarlab/epiread-tools)
2. 
- `model` (str): one of 'uxm','celfie','sum-celfie', 'celfie-plus','reatlas', 'epistate-plus'
- `num_iterations` (int): maximal iterations.
- `stop_criterion` (float): minimal improvement required
- `random_restarts` (int): number of initializations (only one returned)
- `atlas_file` (str): atlas with beta values
- `lambdas` (str): lambda estimates per region (specific to epistate)
- `thetas` (str): theta estimates per region (specific to epistate)
- `percent_u` (str): atlas file with %U values (specific to uxm)
- `summing` (bool): sum each marker region (CelFiE sum)
- `data_file` - mixture file (for simulated data only)
- `metadata_file` - atlas file (for simulated data only)
