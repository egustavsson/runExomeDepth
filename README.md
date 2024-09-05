# runExomeDepth

<!-- badges: start -->
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

This repository provides a pipeline to call copy number variations (CNVs) from Whole Exome Sequencing (WES) or targeted sequencing data using [`ExomeDepth`](https://github.com/vplagnol/ExomeDepth). The tool allows for analysis of multiple samples simultaneously with a set of baseline samples.

# Getting Started

Follow these steps to set up and run the pipeline using Docker. While [`Docker`](https://www.docker.com/) is recommended, you can run it without Docker. 

## Prerequisites

You will need to have Docker installed on your system. For installation instructions, visit Docker's official website at https://www.docker.com/.

> Note: While it is possible to run this without Docker, we strongly recommend against it due to potential compatibility issues with dependencies. If you still wish to proceed without Docker, you can find the necessary information here: [without docker](without_docker/README.md)

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

##### Understanding Mounting Folders in Docker

Mounting folders in Docker is similar to connecting a USB drive to your computer or mapping a network drive. When you "mount" a folder, you're making a directory from your computer (the host) accessible to a specific location inside the Docker container (a virtual environment where the software runs). For example, if you have important data saved in a folder on your computer, you can mount that folder into the Docker container so that the software inside the container can access, read, or modify the files directly. This is crucial because the container itself doesn't have access to your computer's files unless you explicitly share (or "mount") them. The -v option in the Docker command is used to specify which folders on your computer should be made available inside the container and where they should appear. It’s like telling Docker, "Here’s where the files are on my computer, and here’s where I want you to put them inside your environment."

#### Step 3: Run the Docker Container

Use the following command to run the Docker container with the required directories mounted. This command will start an interactive shell inside the container.

```bash
docker run -v /mnt/example:/mnt/example \
           -v /mnt/example-2:/mnt/example-2 \
           -v /mnt/test_data:/mnt/test_data_1 \
           -v /data/working/:/data \
           -it murphydaviducl/runexomedepth:latest /bin/bash
```

#### Explanation:
- `-v /mnt/example:/mnt/example`: This maps the directory `/mnt/example` on the host to `/mnt/example` inside the container. 
- `-it`: Runs the container in interactive mode with a pseudo-TTY, allowing you to interact with the shell.
- `murphydaviducl/runexomedepth:latest`: Specifies the Docker image to use.
- `/bin/bash`: The command to run inside the container, which in this case is starting a Bash shell.

#### Step 4: Activate the Conda Environment

Once inside the container, activate the `runExomeDepth` Conda environment. This environment contains all the necessary dependencies.

```bash
conda activate runExomeDepth
```

This command activates the Conda environment named `runExomeDepth` inside the Docker container.

#### Step 5: Navigate to the Data Directory

Change the directory to the location where your data and scripts are stored inside the container.

```bash
cd /data/
```

This command navigates to the `/data/` directory, which was mounted from the host machine.

# Input Requirements

To run the analysis, you need:

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

Example files can also be found here [test_samples.tsv](./test_samples.tsv) and here [baseline_samples.tsv](./baseline-samples.tsv)

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
