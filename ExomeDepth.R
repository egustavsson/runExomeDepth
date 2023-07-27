# Script    : ExomeDepth.R
# Objective : To call CNVs using ExomeDepth
# Written by: egustavsson

library(optparse)
library(ExomeDepth)
library(GenomicRanges)
library(tidyverse)
library(rtracklayer)

# Function Definitions -----------------------------------------------------

callCNVs <- function(targets, annotation, test_sample, baseline_samples, output_directory) {

  # Check if targets are provided; if not, generate exons.hg19 object
  if (missing(targets) || is.null(targets)) {
    data("exons.hg19")
    targets <- exons.hg19
  } else {
    targets <- read.table(targets, header = FALSE, col.names = c("chrom", "start", "end", "info"))
  }
  
  # Check if annotations are provided; if not, generate genes.hg19 object
  if (missing(annotation) || is.null(annotation)) {
    data("genes.hg19")
    annotation <- %>% dplyr::rename(gene_name = name) %>% GRanges()

    # add chr
    annotation$chromosome <- ifelse(startsWith(annotation$chromosome, "chr"),
                                                 annotation$chromosome, 
                                                 paste0("chr", annotation$chromosome))                                                 
  } else {
    annotation <- rtracklayer::import(opt$annotation) %>% .[.$type == "gene"] %>% unique()
  }

  Counts <- getBamCounts(bed.frame = targets,
                         bam.files = c(test_sample, baseline_samples),
                         include.chr = TRUE) %>%
    setNames(gsub("^X(\\d+)", "\\1", names(.))) # Remove 'X' from column names starting with a number. 
  
  Counts.df <- as.data.frame(Counts)
  
  my.reference.set <- as.matrix(Counts.df[, basename(baseline_samples)])
  
  my.test <- Counts[, basename(test_sample)]
    
  my.choice <- select.reference.set(test.counts = my.test,
                                    reference.counts = my.reference.set,
                                    bin.length = (Counts.df$end - Counts.df$start) / 1000,
                                    n.bins.reduced = 10000)
    
  my.matrix <- as.matrix(Counts.df[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix, MAR = 1, FUN = sum)
    
  all.exons <- new('ExomeDepth', 
                   test = my.test,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = Counts$chromosome,
                        start = Counts$start,
                        end = Counts$end,
                        name = Counts$exon)

  # sort by BF value and annotate
  CNV_calls <- all.exons@CNV.calls %>% arrange(desc(BF)) %>% GRanges()

  # Find overlaps using "any" method to handle partial overlaps
overlap_hits <- findOverlaps(CNV_calls, annotation, type = "any")

# Combine gene names for overlapping ranges
combine_gene_names <- function(gene_names) {
  return(paste(unique(gene_names), collapse = ","))
}

# Initialize a list to store gene names for each CNV_calls entry
gene_names_list <- vector("list", length(CNV_calls))

# Loop through the overlaps and update the gene_names_list
for (i in seq_along(gene_names_list)) {
  # Find all overlaps for the current CNV_calls entry
  current_overlap_indices <- subjectHits(overlap_hits)[queryHits(overlap_hits) == i]
  
  if (length(current_overlap_indices) > 0) {
    # Extract the gene names from annotation for the current overlaps
    overlapping_genes <- mcols(annotation[current_overlap_indices])$gene_name
    
    # Combine gene names with commas and store them in the gene_names_list
    gene_names_list[[i]] <- combine_gene_names(overlapping_genes)
  }
}

# Update the "CNV_calls" object with the combined gene names, preserving NA values
CNV_calls$gene_name <- unlist(lapply(gene_names_list, function(x) if (is.null(x)) NA else x))
  
  # Generate the output filename based on the test sample name
  sample_name <- gsub("\\.bam$", "", basename(test_sample))
  output_file <- file.path(output_directory, paste0(sample_name, "_CNV.csv"))
    
  write.csv(file = output_file,
            x = CNV_calls,
            row.names = FALSE)
  
  # Print completion message for the test sample
  cat("Analysis completed for", sample_name, "\n")
}

# Parse command-line options
option_list <- list(
  make_option("--test-sample", dest="test_sample", type="character"),
  make_option("--baseline-samples", dest="baseline_samples", type="character"),
  make_option("--targets", dest="targets", type="character"),
  make_option("--output-directory", dest="output_directory", type="character")  # Updated option name
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Redirect output to a log file to suppress warnings and messages
log_file <- file.path(opt$output_directory, "exomedepth_log.txt")
sink(log_file, append = FALSE)

# Read test samples from TSV
test_samples <- read.table(opt$test_sample, header = FALSE, col.names = "test_sample_path")

# Read baseline samples from TSV
baseline_samples <- read.table(opt$baseline_samples, header = FALSE, col.names = "baseline_sample_path")

# Run the analysis for each test sample
for (test_sample_path in test_samples$test_sample_path) {
  callCNVs(
    targets = opt$targets,
    test_sample = test_sample_path,
    baseline_samples = baseline_samples$baseline_sample_path,
    output_directory = opt$output_directory  # Updated argument name
  )
}

# Close the sink to restore the standard output
sink()
