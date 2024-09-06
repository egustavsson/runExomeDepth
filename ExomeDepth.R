# Script    : ExomeDepth.R
# Objective : To call CNVs using ExomeDepth
# Written by: egustavsson

library(optparse)
library(ExomeDepth)
library(GenomicRanges)
library(tidyverse)
library(rtracklayer)

# Function Definitions -----------------------------------------------------

callCNVs <- function(targets, annotation, test_sample, output_directory) {

  test_Counts <- getBamCounts(bed.frame = targets,
                            bam.files = test_sample,
                            include.chr = TRUE) %>%
  setNames(gsub("^X(\\d+)", "\\1", names(.))) # Remove 'X' from column names starting with a number.

  Counts.df <- as.data.frame(
      dplyr::left_join(test_Counts, base_Counts)
      )
  
  my.reference.set <- as.matrix(Counts.df[, basename(baseline_samples$baseline_sample_path)])
  
  my.test <- test_Counts[, basename(test_sample)]
    
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
  
  # Check if any CNVs were called
  if (nrow(all.exons@CNV.calls) == 0) {
    # Generate the output filename based on the test sample name
    sample_name <- gsub("\\.bam$", "", basename(test_sample))
    output_file <- file.path(output_directory, paste0(sample_name, "_CNV.csv"))
    
    # Create a dataframe with the specified headers and no data
    empty_df <- data.frame(
      seqnames = character(),
      start = integer(),
      end = integer(),
      width = integer(),
      strand = character(),
      start.p = integer(),
      end.p = integer(),
      type = character(),
      nexons = integer(),
      id = character(),
      BF = numeric(),
      reads.expected = numeric(),
      reads.observed = numeric(),
      reads.ratio = numeric(),
      gene_name = character()
    )
    
    # Write the empty dataframe to CSV
    write.csv(empty_df, file = output_file, row.names = FALSE)
    cat("No CNVs called for", sample_name, "using the baseline samples provided.\n")
  } else {
    # Continue processing if there are CNVs
    CNV_calls <- all.exons@CNV.calls %>% GRanges()

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
}

# Parse command-line options
option_list <- list(
  make_option("--test-samples", dest="test_samples", type="character"),
  make_option("--baseline-samples", dest="baseline_samples", type="character"),
  make_option("--targets", dest="targets", type="character"),
  make_option("--annotation", dest="annotation", type="character"),
  make_option("--output-directory", dest="output_directory", type="character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input options
if (!file.exists(opt$test_samples)) {
  stop("Error: The file specified for --test-samples does not exist: ", opt$test_samples)
}

if (!file.exists(opt$baseline_samples)) {
  stop("Error: The file specified for --baseline-samples does not exist: ", opt$baseline_samples)
}

if (!dir.exists(opt$output_directory)) {
  stop("Error: The directory specified for --output-directory does not exist: ", opt$output_directory)
}

# Read test samples from TSV
test_samples <- read_tsv(opt$test_samples, col_names = "test_sample_path", show_col_types = F)

# Read baseline samples from TSV
baseline_samples <- read_tsv(opt$baseline_samples, col_names = "baseline_sample_path", show_col_types = F)

# Initialize targets and annotation once outside the loop
if (is.null(opt$targets)) {
  data("exons.hg19")
  targets <- exons.hg19
} else {
  targets <- read.table(opt$targets, header = FALSE, col.names = c("chrom", "start", "end", "info"))
}

if (is.null(opt$annotation)) {
  data("genes.hg19")
  annotation <- genes.hg19 %>%
    dplyr::rename(gene_name = name) %>%
    mutate(chromosome = paste0("chr", chromosome)) %>%
    GRanges()
} else {
  annotation <- rtracklayer::import(opt$annotation) %>%
    .[.$type == "gene"] %>% unique() # Ensure "chr" within seqnames if needed
}

# run getBamCounts() for baseline samples first so that they do not need to be generated for each iteration
base_Counts <- getBamCounts(bed.frame = targets,
                            bam.files = baseline_samples$baseline_sample_path,
                            include.chr = TRUE) %>%
  setNames(gsub("^X(\\d+)", "\\1", names(.)))

# Run the analysis for each test sample
for (test_sample_path in test_samples$test_sample_path) {
  # Generate the log filenames based on the test sample name
  sample_name <- gsub("\\.bam$", "", basename(test_sample_path))
  output_log_file <- file.path(opt$output_directory, paste0(sample_name, "_stdout.log"))
  message_log_file <- file.path(opt$output_directory, paste0(sample_name, "_stderr.log"))

  # Try-Catch block to ensure sinks are closed properly
  tryCatch({
    # Check if the output directory exists
    if (!dir.exists(opt$output_directory)) {
      stop("Output directory does not exist or is not accessible: ", opt$output_directory)
    }
    
    # Redirect stdout to the sample-specific output log file
    sink(output_log_file, append = FALSE)
    
    # Open connection for stderr redirection
    message_log_conn <- file(message_log_file, open = "wt")
    sink(message_log_conn, append = FALSE, type = "message")
    
    # Call the function for each test sample
    callCNVs(
      targets = targets,
      annotation = annotation,
      test_sample = test_sample_path,
      output_directory = opt$output_directory
    ) 
    gc() # garbage collection
  }, finally = {
    # Ensure sinks are closed in the proper order
    sink(type = "message")  # Stop redirecting stderr
    close(message_log_conn) # Close the stderr connection
    sink()  # Stop redirecting stdout
  })
}
