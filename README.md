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
you may override any param via the command line 
```
deconvolution --model <model> -j <config.json> outfile=<output_filepath>
```