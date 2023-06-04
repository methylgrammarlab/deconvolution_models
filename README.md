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
config.json should contain the following parameters:

1. Any and all parameters for epiread parsing (see https://github.com/methylgrammarlab/epiread-tools)
2. "model": one of 'uxm','celfie','sum-celfie', 'celfie-plus','reatlas', 'epistate-plus'
3. Parameters common to multiple models:
* "num_iterations" - maximal iterations
* "stop_criterion" - minimal improvement required
* "random_restarts" - number of initializations (only one returned)
* "atlas_file" - atlas with beta values

3. Model specific:
Epistate: 
* "lambdas" - lambda estimates per region
* "thetas" - theta estimates per region
UXM:
* "percent_u" - atlas file with %U values
CelFiE:
* "summing" - True if sum each marker region (CelFiE sum)
4. For simulated data only: (shortcut instead of simulating proper epiread files):
* "data_file" - mixture
* "metadata_file" - atlas