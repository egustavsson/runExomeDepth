# runExomeDepth

<!-- badges: start -->
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

This repo contains code to help call CNVs from WES, or targeted sequencing data, using [`ExomeDepth`](https://github.com/vplagnol/ExomeDepth). It allows multiple samples to be analysed at the same time with a set of baseline samples.
# Getting Started

## Prerequisites

To run this you will need to have `Docker` installed on your system. `Docker` is a platform that allows you to create, deploy, and run applications in containers. More information about `Docker` can be found at https://www.docker.com/. Follow the instructions below to install `Docker` if you have not already:

### Installing Docker

#### For Linux:

1. **Open a Terminal:** Use your Linux distribution's package manager to install `Docker`. The command varies depending on the distribution:
   
   **Ubuntu/Debian**: 
   ```
   sudo apt-get install docker-ce docker-ce-cli containerd.io
   ```
   **Fedora**: 
   ```
   sudo dnf -y install docker-ce
   ```
   **CentOS**: 
   ```
   sudo yum install docker-ce docker-ce-cli containerd.io
   ```
2. **Start the Docker Service:** Use `sudo systemctl start docker` to start the Docker service.

3. **Verify Installation:** Check if Docker is installed correctly by typing `docker --version` in the terminal.



## Input

To run this analysis you will need the following input:

  - a set of BAM files for which to call CNVs - one sample per BAM file 
  - a set of BAM files to use as the baseline - one sample per BAM file
  - indexed BAM files (.bai) for all the above BAM files
  - a BED file of the target region of your exome capture or targeted sequencing data. If this is not supplied hg19 will be used.
  - an annotation file (GTF/GFF). This needs to match the build of your targets. If this is not supplied ensembl version 71 (hg19) will be used.

It is advisable to do the indexing of the BAM files prior to running the pipeline. If the dates of the index files are older than the BAM files, ExomeDepth will throw an error. Indexing of BAM files can be done by:

```bash
samtools index input.bam # for a single sample or;
samtools index -M *.bam # multiple samples
```

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies are installed via conda through the `environment.yml` file

## Installation

Clone the directory. From command line simply run the following command from the directory where you wish to install the repo:

```bash
git clone --recursive https://github.com/egustavsson/runExomeDepth.git
```

## Analysis steps

### 1. Create the conda environment
First you neeed to create the conda environment which will install all the dependencies. This step only needs to be done once:

```bash
cd runExomeDepth
conda env create -f environment.yml
```

### 2. Activate the conda environment
After the conda environment has been created it needs to be activated prior to running the analysis. Therefore, make sure to first activate the conda environment using the command `conda activate runExomeDepth`.

It can be done like this:
```bash
cd runExomeDepth
conda activate runExomeDepth
```
After you are done with the analysis you can deactivate the conda environment by:
```bash
conda deactivate
```

### 3. Install ExomeDepth and required R packages
While `R` and `R-essentials` are installed throught the conda environment, other required R packages, including `ExomeDepth` are installed by running the script `install-packages.R`.

From the `./runExomDepth` directory, the script can be ran like this:

```bash
Rscript install-packages.R
```

Installing R packages using this script only needs to be done once and are saved within the conda environment.

### 2. Call CNVs with ExomeDepth

#### Input data
The main script to call CNVs with is called `ExomeDepth.R`. Make sure you have the following input data prior to running it:

| Parameter | Description |
| --- | --- |
| `--targets` | bed file with exon targets. This is optional and if none is given hg19 will be used |
| `--annotation` | GTF/GFF file with gene coordinates. This is optional and if none is given ensembl version 71 (hg19) will be used |
| `--test-samples` | TSV file with  paths to the BAM files to call CNVs for, one per line |
| `--baseline-samples` | TSV file with  paths to the BAM files used for the baseline, one per line |
| `--output-directory` | path to output directory |

Example of `--test-samples` and `--baseline-samples` required TSV files looks like this:

```bash
/path/to/test_sample1.bam
/path/to/test_sample2.bam
/path/to/test_sample3.bam
```

Example files can also be found here [test_samples.tsv](./test_samples.tsv) and here [baseline_samples.tsv](./baseline-samples.tsv)

#### Run the script
Once you have the required input data, follow these steps to run the `ExomeDepth.R` script:

```
Rscript ExomeDepth.R \
        --targets /path/to/targets.bed \
        --annotation /path/to/annotation.gff \
        --test-samples test_samples.tsv \
        --baseline-samples baseline_samples.tsv \
        --output-directory /path/to/output_folder/
```

#### Output
```
working directory  
|--- {sample}_CNV.csv      # A CSV file with called CNVs per sample  
|--- {sample}_stderr.log   # stdout to a sample-specific output log file
|--- {sample}_stdout.log   # stderr to a sample-specific output log file
