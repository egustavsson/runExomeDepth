# runExomeDepth

Here are some scripts to help call CNVs from WES data using [`ExomeDepth`](https://cran.r-project.org/web/packages/ExomeDepth/index.html). 

# Getting Started

## Input

To run this analysis you just need following input:

  - a set of BAM files - one sample per BAM file 
  - indexed BAM files (.bai) for the above BAM files
  - a BED file of the target region of your exome or targeted sequencing data. If this is not supplied hg19 will be used.

It is advisable ro do the indexing of the BAM files prior to running the pipeline as if the index dates are not newer than the BAM file ExomeDepth will throw an error. This can be done by:
```bash
samtools index input.bam # for  single sample or;
samtools index *.bam # multiple samples
```

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies are installed via conda through the `environment.yml` file

## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/egustavsson/runExomeDepth.git
```

## Analysis steps

### 1. Activate the conda environment
Create the conda environment for the pipeline which will install all the dependencies:

```bash
cd runExomeDepth
conda env create -f environment.yml
```
This step only needs to be done once. The conda environment needs to be activated prior to running any analysis. Therefore, make sure to first activate the conda environment using the command `conda activate runExomeDepth`.

This can be done like this:
```bash
cd runExomeDepth
conda activate runExomeDepth
```
To deactivate the conda environment:
```bash
conda deactivate
```

### 2. Install ExomeDepth and required R packages
To make sure that required R packages, including `ExomeDepth` are installed simply run the script `install-packages.R`.

From the `./runExomDepth` directory, the script can be ran like this:
```bash
Rscript install-packages.R
```
R packages also only needs to be installed once and are saved within the conda environment.

### 2. Call CNVs with ExomeDepth

#### Input data
To run the ExomeDepth.R you need the following input data:

| Parameter | Description |
| --- | --- |
| `--targets` | bed file with exon targets. This is optiona and if none is given hg19 will be used |
| `--test-sample` | `TSV` file with  paths to the BAM files to call CNVs for, one per line |
| `--baseline-samples` | `TSV` file with  paths to the BAM files used for the baseline, one per line |
| `--output-directory` | path to output directory |

Example of `--test-sample` and `--baseline-samples` required `TSV` files:

```bash
/path/to/test_sample1.bam
/path/to/test_sample2.bam
/path/to/test_sample3.bam
```

#### Run the script
Once you have the required input data, follow these steps to run the `ExomeDepth.R` script:

```
Rscript ExomeDepth.R \
        --targets /path/to/targets.bed \
        --test-sample test_samples.tsv \
        --baseline-samples baseline_samples.tsv \
        --output-directory /path/to/output_folder/
```