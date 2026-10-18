# Architecture

## Purpose

`deconvolution_models` is a Python 3 package for estimating cell-type proportions from read-level DNA methylation mixtures and reference atlases. It exposes several deconvolution methods through one command-line interface.

## Runtime Flow

1. The installed `deconvolution` console command invokes `deconvolution_models.main:main`.
2. The Click command builds a configuration from command-line options and optional extra key/value arguments. A supplied JSON configuration is applied last and therefore overrides command-line values.
3. `main.py` selects a model wrapper (`Celfie`, `CelfieISH`, `ReAtlas`, `UXM`, or `Epistate`). `sum-celfie` uses the `Celfie` wrapper with regional summing enabled.
4. The wrapper uses `epiread-tools` readers to load and align mixture reads, CpG coordinates, genomic intervals, and the model-specific atlas.
5. The wrapper calls the numerical implementation in the corresponding package module and writes cell-type proportions, per-read probabilities, or a NumPy pipeline artifact.

## Module Map

- `deconvolution_models/main.py`: Click CLI, configuration assembly, shared `DECmodel` lifecycle, model-specific I/O wrappers, and output handling.
- `deconvolution_models/celfie.py`: CelFiE expectation-maximization implementation.
- `deconvolution_models/celfie_ish.py`: CelFiE-ISH read-level mixture model and probability/clipping support.
- `deconvolution_models/celfie_ish_reatlas.py`: CelFiE-ISH variant that incorporates reference-atlas methylation and coverage.
- `deconvolution_models/epistate.py`: Epistate mixture model using per-region lambda and high/low methylation parameters.
- `deconvolution_models/UXM.py`: UXM optimization using unmethylated-read proportions and optional region weights.
- `tests/test_models.py`: Pytest fixtures and deterministic unit checks for the CelFiE and CelFiE-ISH numerical cores; three atlas test placeholders currently contain no assertions.
- `demo/`: Small example inputs and expected/generated outputs used by commands in `README.md`.
- `supplementary/`: Larger reference interval resources.

## Interfaces and File Formats

- Supported model names: `uxm`, `celfie`, `sum-celfie`, `celfie-ish`, `reatlas`, and `epistate`.
- Mixtures are read-level methylation data in `old_epiread`, `old_epiread_A`, or `pat` format and are normally bgzip-compressed and tabix-indexed.
- Genomic intervals are supplied as a comma-delimited `chr:start-end` list or a BED-like file. Coordinate conventions must remain consistent among the mixture, CpG-coordinate file, intervals, and atlas.
- Atlas shape and required fields vary by model: methylation/coverage pairs for CelFiE-family models, percent-unmethylated values for UXM, and lambda/theta files for Epistate.
- Default text output is a tab-delimited row of percentage estimates. `--probs` writes per-read origins and probabilities; `--npy` writes a NumPy object containing estimates and iteration metadata.
- `--minimal_cpg_per_read` defaults to 1. CelFiE-ISH currently permits only one random restart. Probability clipping uses `--read_max_clipping_cell_type_probabilty`, whose current default is `0.05`.

## Packaging and Dependencies

- Packaging uses `setuptools` through `setup.py`; the current package version there is `0.0.2`.
- The console entry point is `deconvolution = deconvolution_models.main:main`.
- Core numerical dependencies include NumPy, pandas, SciPy, and Bottleneck; Click provides the CLI.
- `epiread-tools` is installed directly from its GitHub repository and owns the parsers, naming conventions, coverage calculations, and atlas readers used by the wrappers.
- `requirements.txt` is a pinned legacy environment specification, while `setup.py` contains the smaller runtime dependency set.

## Current Constraints

- Reference atlases are loaded into memory and must cover and align with the requested marker regions/CpGs as documented in `TUTORIAL.md`.
- `main.py` and `tests/test_models.py` contain historical absolute `sys.path` additions for a former developer environment. Normal installed-package operation is intended to resolve `epiread-tools` through the declared dependency.
- The repository currently uses `setup.py` rather than a `pyproject.toml` build configuration.
