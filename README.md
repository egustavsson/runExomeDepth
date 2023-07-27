# runExomeDepth

This repo contains code to help call CNVs from WES data using [`ExomeDepth`](https://cran.r-project.org/web/packages/ExomeDepth/index.html). The main script `ExomeDepth.R` allows multiple samples to be analysed at the same time.

# Getting Started

## Input

To run this analysis you need the following input:

  - a set of BAM files for which to call CNVs - one sample per BAM file 
  - a set of BAM files to use as the baseline - one sample per BAM file
  - indexed BAM files (.bai) for the above BAM files
  - a BED file of the target region of your exome or targeted sequencing data. If this is not supplied hg19 will be used.
  - annotation file (GTF/GFF). If this is not supplied ensembl version 71 (hg19) will be used.

It is advisable to do the indexing of the BAM files prior to running the pipeline. If the dates of the index files are older than the BAM files, ExomeDepth will throw an error. Indexing of BAM files can be done by:
```bash
samtools index input.bam # for  single sample or;
samtools index -M *.bam # multiple samples
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
First you neeed to create the conda environment which will install all the dependencies. This step only needs to be done once:

```bash
cd runExomeDepth
conda env create -f environment.yml
```
After the conda environment has been created it needs to be activated prior to running any analysis. Therefore, make sure to first activate the conda environment using the command `conda activate runExomeDepth`.

It can be done like this:
```bash
cd runExomeDepth
conda activate runExomeDepth
```
To deactivate the conda environment simply run:
```bash
conda deactivate
```

### 2. Install ExomeDepth and required R packages
While base R and R essentials are installed throught the conda environment, required R packages, including `ExomeDepth` are installed by running the script `install-packages.R`.

From the `./runExomDepth` directory, the script can be ran like this:
```bash
Rscript install-packages.R
```
R packages also only needs to be installed once and are saved within the conda environment.

### 2. Call CNVs with ExomeDepth

#### Input data
To run the `ExomeDepth.R` you need the following input data:

| Parameter | Description |
| --- | --- |
| `--targets` | bed file with exon targets. This is optional and if none is given hg19 will be used |
| `--annotation` | GTF/GFF file with gene coordinates. This is optional and if none is given ensembl version 71 (hg19) will be used |
| `--test-sample` | TSV file with  paths to the BAM files to call CNVs for, one per line |
| `--baseline-samples` | TSV file with  paths to the BAM files used for the baseline, one per line |
| `--output-directory` | path to output directory |

Example of `--test-sample` and `--baseline-samples` required TSV files:

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
        --annotation /path/to/annotation.gff \
        --test-sample test_samples.tsv \
        --baseline-samples baseline_samples.tsv \
        --output-directory /path/to/output_folder/
```

#### Output

The output is a CSV file per sample and a log file.