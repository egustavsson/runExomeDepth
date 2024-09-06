# Required for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.ma.imperial.ac.uk/")
}

# Bioconductor packages with specific versions
bioc_packages <- list(
  "rtracklayer" = "1.60.1",
  "GenomicRanges" = "1.52.0",
  "GenomeInfoDb" = "1.36.3",
  "IRanges" = "2.34.1",
  "S4Vectors" = "0.38.2",
  "BiocGenerics" = "0.46.0",
  "Rsamtools" = "2.16.0",
  "Biostrings" = "2.68.1",
  "GenomicAlignments" = "1.36.0"
)

# Function to install a Bioconductor package with a specific version
install_bioc_with_version <- function(package, version) {
  if (!requireNamespace(package, quietly = TRUE) || packageVersion(package) != version) {
    BiocManager::install(package, version = version)
  }
}

# Install Bioconductor packages
for (package in names(bioc_packages)) {
  install_bioc_with_version(package, bioc_packages[[package]])
}

# CRAN packages with specific versions
cran_packages <- list(
  "lubridate" = "1.9.2",
  "forcats" = "1.0.0",
  "stringr" = "1.4.0",
  "dplyr" = "1.1.3",
  "purrr" = "1.0.2",
  "readr" = "2.1.4",
  "tidyr" = "1.3.0",
  "tibble" = "3.2.1",
  "ggplot2" = "3.4.4",
  "tidyverse" = "2.0.0",
  "optparse" = "1.7.3",
  "ExomeDepth" = "1.1.16"
)

# Function to install a CRAN package with a specific version
install_cran_with_version <- function(package, version) {
  if (!requireNamespace(package, quietly = TRUE) || packageVersion(package) != version) {
    install.packages(package, repos = "https://cran.ma.imperial.ac.uk/", version = version)
  }
}

# Install CRAN packages
for (package in names(cran_packages)) {
  install_cran_with_version(package, cran_packages[[package]])
}

# Check and print installed packages
installed_packages <- c(names(bioc_packages), names(cran_packages))
for (package in installed_packages) {
  if (requireNamespace(package, quietly = TRUE)) {
    cat("Package", package, "version", packageVersion(package), "is installed.\n")
  } else {
    cat("Package", package, "is not installed.\n")
  }
}
