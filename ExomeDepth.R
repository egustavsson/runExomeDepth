# Updated ExomeDepth.R

library(optparse)
library(ExomeDepth)
library(tidyverse)

# Function Definitions -----------------------------------------------------

callCNVs <- function(targets, test_sample, baseline_samples, output_directory) {

  # Check if targets are provided; if not, generate exons.hg19 object
  if (missing(targets) || is.null(targets)) {
    data("exons.hg19")
    targets <- exons.hg19
  } else {
    # Read the bed file into a data frame
    # Replace "your_bed_file.bed" with the actual path to your bed file
    targets <- read.table(targets, header = FALSE, col.names = c("chrom", "start", "end", "info"))
  }
  
  Counts <- getBamCounts(bed.frame = targets,
                         bam.files = c(test_sample, baseline_samples),
                         include.chr = TRUE) %>%
    setNames(gsub("^X(\\d+)", "\\1", names(.))) # Remove 'X' from column names starting with a number. Make sure magrittr is installed
  
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
  
  # Generate the output filename based on the test sample name
  sample_name <- gsub("\\.bam$", "", basename(test_sample))
  output_file <- file.path(output_directory, paste0(sample_name, "_CNV.csv"))
    
  write.csv(file = output_file,
            x = all.exons@CNV.calls,
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
