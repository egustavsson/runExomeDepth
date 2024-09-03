# runExomeDepth withut docker

While we strongly advise against it, due to issues with version and comptability of dependencies, this analysis can be ran withpout the use of docker.

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