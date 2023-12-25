# deconvolution_models Tutorial

## Environment & Installation

**deconvolution_models** is a Python 3 package. To install, use the following command:

```shell
pip install git+https://github.com/methylgrammarlab/deconvolution_models
```

Alternatively, set up your Python 3 environment with the requirements.txt file and clone this repository.

In this tutorial, we will focus on running CelFiE-ISH to deconvolute WGBS data.

## Input Files

### Marker Regions

Marker regions represent intervals where methylation differences between cell types are expected. For the Loyfer atlas, we recommend using the U250 regions available in their supplementary table 4b. If you have your own reference dataset, utilize [UXM](https://github.com/nloyfer/UXM_deconv) or [TIM](https://github.com/christacaggiano/celfie) to select a marker region set. Marker regions can be provided either in a tab-delimited file format:

```shell
chr1    1500    1770
chr3    1350    1400
```

Or as a comma-delimited list to the shell:

```shell
chr1:1500-1770,chr3:1350-1400
```

Options:

- `-i`, `--genomic_intervals`: Specify the genomic intervals to process in the format chrN:start-end, separated by commas.
- `-b`, `--bedfile`: The intervals (-i) are in a bedfile instead of a comma-delimited list.
- `--header`: Specifies that the bedgraph file with regions to process has a header (e.g., CHROM START END).

### Reference Atlas

The reference atlas for CelFiE-ISH should contain one row *per CpG* within each marker region, with methylation and coverage information for each reference cell type:

```shell
chr start   end alpha_meth  alpha_cov   beta_meth   beta_cov
chr1    1500    1501    3   10  9   10
chr1    1505    1506    29  30  4   5
```

Options:

- `-a`, `--atlas_file`: Path to the atlas file containing methylation and coverage for each cell type.

Notes:

* Each reference cell type is represented by two columns (methylation and coverage). If you have multiple samples of each cell type, merge them.
* CelFiE-ISH requires one row for each CpG, so a marker region with 100 CpGs would be represented by 100 rows.
* The actual coverage is only used in the ReAtlas model. If you have array data with beta values, you can arrange them in a similar format; a beta of 0.7 can be represented as methylation 700 coverage 1000.
* Other models (UXM, Epistate) require atlas files with different formats (see below).
* The reference atlas may contain the entire genome or only the marker regions. Only the marker regions will be taken into account.
* This file should be tabix-indexed for access:

```shell
bgzip atlas.file
tabix -p bed atlas.file.gz
```

### CpG File

- `--cpg_coordinates`: Path to the sorted CpG bed file. Should contain coordinates for all CpG sites in bed format.

This file has the coordinates of all CpGs in the genome:

```shell
chr1    10468   10469   CpG1
chr1    10470   10471   CpG2
chr1    10483   10484   CpG3
chr1    10488   10489   CpG4
chr1    10492   10493   CpG5
```

All files to deconvolute should match this file. If a CpG is assumed but not present in the CpG file:

```shell
chr1    10468 10493 A00222:175:HGM72DMXX:2:1107:11089:34882 1       -       10468,10470,10472 -CC     .       .
```

Note: 10472 is not in the CpG file! This will cause an error.

The CpG file should not be limited just to the marker regions to ensure proper alignment of partially overlapping reads. It should also be zipped and indexed.

### Mixture File

This is the file you want to deconvolute, and it can have one of several formats:

- `-m`, `--mixture`: Path to the mixture file to deconvolute.
- `--epiformat`: Specify the format of the epiread files. Available options are 'old_epiread', 'old_epiread_A', and 'pat'.

#### [old_epiread](https://huishenlab.github.io/biscuit/old_epiread_format/):

```shell
chr19    read_456    1    +    3040315    CCCCTCCC    .    .
chr19    read_789    1    +    3078472    CC    3078510    T
```

#### old_epiread_A:

Generated with the -A flag in [biscuit](https://huishenlab.github.io/biscuit), it contains all the CpG coordinates:

```shell
chr1    1045697 1045787 A00222:175:HGM72DMXX:2:1107:11089:34882 1       -       1045697,1045700,1045787 -CC  
```

Ensure these are generated correctly. The coordinate length (1045697,1045700,1045787) should always equal the pattern length (-CC), and all coordinates should be listed in the CpG file.

#### [pat](https://github.com/nloyfer/wgbs_tools/blob/master/docs/pat_format.md)

```shell
chr1    46       CT    1
chr1    47       CC..TC  1
chr1    47       T       13
```

If you use pat, ALL coordinates should be in "pat," meaning the CpG number is used instead of genomic. The CpG file will read something like:

```shell
chr1    1   2   CpG1
chr1    2   3   CpG2
chr1    3   4   CpG3
```

The marker regions should also be translated to "pat." To do that, keep a "rosetta" file with both genomic and "pat" coordinates:

```shell
chr1    10468   10469   1
chr1    10470   10471   2
chr1    10483   10484   3
chr1    10488   10489   4
```

Use [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/map.html) to find the interval borders:

```shell
bedtools map -a genomic.sorted.tsv -b hg19.pat_rosetta.sorted.bed -c 4,4 -o min,max
```

The atlas should also be in "pat."

### Additional Information

- **Can I run this from a bed file? Bedgraph?

**
  - No, you need read-level information for CelFIE-ISH. The methylation and coverage are fine for the atlas, but the mixture needs to be read-based.

- **Can I run directly from BAM?**
  - BAM is currently not directly supported. You can generate epiread/pat files from BAM.
 
## Run Parameters
CelFiE-ISH runs Expectation-Maximization to estimate the cell type proportions. 
- `--num_iterations`: Set the maximum number of iterations.
- `--stop_criterion`: Set the minimal improvement required to continue deconvolution. 
- `--random_restarts`: Set the number of initializations (only one will be returned).


- `--model`: Specify the deconvolution model to use. Available options are 'uxm', 'celfie', 'sum-celfie', 'celfie-ish', 'reatlas', and 'epistate'.
- `--minimal_cpg_per_read`: Set the minimum number of CpGs required for a read to be considered. Default is 1.
- `-j`, `--json`: Run the deconvolution using a JSON config file.
- `--outfile`: Path to the output file (to be generated).



## Other models
- `--lambdas`: Specify the lambda estimates per region (specific to epistate).
- `--thetas`: Specify the theta estimates per region (specific to epistate).
- 
- `--percent_u`: Specify the atlas file with %U values (specific to UXM).
- `--weights`: Specify the weights per marker region (specific to UXM).
- `--u_threshold`: Set the maximal methylation value to be considered as U (specific to UXM).
- `--min_length`: Set the minimum number of CpGs required for a read to be considered at the deconvolution level (specific to UXM). Same as `--minimal_cpg_per_read` but applied at the deconvolution level.
- 
- `-s`, `--summing`: Perform summing for each marker region (CelFiE sum).