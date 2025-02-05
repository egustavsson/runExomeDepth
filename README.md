# runExomeDepth

<!-- badges: start -->
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/669066146.svg)](https://doi.org/10.5281/zenodo.14008789)
<!-- badges: end -->

This repository provides a tool to call copy number variations (CNVs) from Whole Exome Sequencing (WES) or targeted sequencing data using [`ExomeDepth`](https://github.com/vplagnol/ExomeDepth). The tool is containerised using `Docker` to eliminate compatibility issues and allows for analysis of multiple samples simultaneously with a set of baseline samples.

# Getting Started

Follow these steps to set up and run the pipeline using `Docker`. While using `Docker` is recommended, you can run it without `Docker`. 

## Prerequisites

You will need to have `Docker` installed on your system. For installation instructions, visit Docker's official website at https://www.docker.com/.

> [!WARNING]
> While it is possible to run this without Docker, we strongly recommend against it due to potential compatibility issues with dependencies. If you still wish to proceed without Docker, you can find the necessary information here: [without docker](without_docker/README.md)

To run this you will need to have `Docker` installed on your system. `Docker` is a platform that allows you to create, deploy, and run applications in containers. More information about `Docker` can be found at https://www.docker.com/. Follow the instructions below to install `Docker` if you have not already:

### Docker Installation (Linux)

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

### Running the Docker Container

#### Step 1: Pull the Docker Image

The first step is to pull the Docker image from Docker Hub. This ensures you have the latest version of the `runexomedepth` image.

```bash
docker pull murphydaviducl/runexomedepth:latest
```

This command downloads the latest version of the `runexomedepth` Docker image to your local machine.

#### Step 2: Prepare Directories

You need to have specific directories on your host machine that will be mounted into the Docker container. These directories will be used for input data and output results.

For example:
- `/mnt/example`: Directory for general data.
- `/mnt/example-2`: Another directory for additional data.
- `/mnt/test_data`: A third directory for more data.
- `/data/working/`: Directory where output files will be stored.

Make sure these directories exist or replace them with paths that suit your environment.

> [!NOTE]
>***Mounting Folders in Docker***
>
>Mounting folders in Docker allows you to share a directory from your computer (the host) with a specific location inside a Docker container. This makes the files accessible to the software running inside the container.
>
>***Why mount a folder?***
> 
>To allow the container to access, read, or modify files on your computer. By default, the container does not have access to your files unless you explicitly share them.
>
>***How to mount a folder***
> 
>Use the `-v` option in the Docker command to specify:
>- The folder on your computer (host) that you want to share.
>- The location inside the container where you want the folder to appear.

#### Step 3: Run the Docker Container

Use the following command to run the Docker container with the required directories mounted. This command will start an interactive shell inside the container.

```bash
docker run -v /host/data1:/container/data1 \
           -v /host/data2:/container/data2 \
           -v /host/test_data:/container/test_data \
           -v /host/working_dir:/container/working_dir \
           -it murphydaviducl/runexomedepth:latest /bin/bash
```

#### Explanation:
- `v /host/data1:/container/data1`: Maps `/host/data1` from your machine to `/container/data1` in the container.
- `v /host/working_dir:/container/working_dir`: Maps `/host/working_dir` on the host to `/container/working_dir` in the container.
- `it`: Runs the container interactively, allowing you to interact with the shell.
- `murphydaviducl/runexomedepth:latest`: Specifies the Docker image.
- `/bin/bash`: Starts a Bash shell inside the container.

#### Step 4: Navigate to the Data Directory

Change the directory to the location where your data and scripts are stored inside the container.

```bash
cd /data/
```

This command navigates to the `/data/` directory, which was mounted from the host machine.

# Input Requirements

To run the analysis, you need:

  - a set of `BAM` files for which to call CNVs - one sample per `BAM` file 
  - a set of `BAM` files to use as the baseline - one sample per `BAM` file
  - indexed `BAM` files (`BAI`) for all the above `BAM` files
  - a `BED` file of the target region of your exome capture or targeted sequencing data. If this is not supplied hg19 will be used.
  - an annotation file (`GTF`/`GFF`). This needs to match the build of your targets. If this is not supplied ensembl version 71 (hg19) will be used. These can be downloaded from [Ensembl](https://www.ensembl.org/index.html) or [GENCODE](https://www.ensembl.org/index.html). 

The latest GENCODE annotation can be found here: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/

The latest ensembl annotation can be found here: https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/

It is advisable to do the indexing of the BAM files prior to running the pipeline. If the dates of the index files are older than the BAM files, ExomeDepth will throw an error. Indexing of BAM files can be done by:

```bash
samtools index input.bam # for a single sample or;
samtools index -M *.bam # multiple samples
```

#  Call CNVs with ExomeDepth

## Define input data
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

Example files can also be found here [test_samples.tsv](./example_files/test_samples.tsv) and here [baseline_samples.tsv](./example_files/baseline-samples.tsv)

## Run the script
Once you have the required input data, follow these steps to run the `ExomeDepth.R` script:

```
Rscript ExomeDepth.R \
        --targets /path/to/targets.bed \
        --annotation /path/to/annotation.gff \
        --test-samples test_samples.tsv \
        --baseline-samples baseline_samples.tsv \
        --output-directory /path/to/output_folder/
```

# Output

The output will be a CSV file with called CNVs, a stderr and a stdout log file per sample stored in the working directory.

```
working directory  
|--- {sample}_CNV.csv      # A CSV file with called CNVs per sample  
|--- {sample}_stderr.log   # stdout to a sample-specific output log file
|--- {sample}_stdout.log   # stderr to a sample-specific output log file

```

# Citation
If using this tool please cite the original ExomeDepth manuscript

Plagnol, Vincent, et al. "A robust model for read count data in exome sequencing experiments and implications for copy number variant calling." *Bioinformatics* 28.21 (2012). https://doi.org/10.1093/bioinformatics/bts526
