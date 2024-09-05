# Required for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.ma.imperial.ac.uk/")
}

# Bioconductor packages
bioc_packages <- c("Biostrings", "IRanges", "Rsamtools", "GenomicRanges", "GenomicAlignments", "rtracklayer")

# Function to install a Bioconductor package if it is not already installed
install_bioc_if_not_installed <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

# Install Bioconductor packages
for (package in bioc_packages) {
  install_bioc_if_not_installed(package)
}

# Function to install a CRAN package if it is not already installed
install_cran_if_not_installed <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = "https://cran.ma.imperial.ac.uk/")
  }
}

# Install CRAN packages
cran_packages <- c("optparse", "ExomeDepth")
for (package in cran_packages) {
  install_cran_if_not_installed(package)
}

# Check and print installed packages
installed_packages <- c(bioc_packages, cran_packages)
for (package in installed_packages) {
  if (requireNamespace(package, quietly = TRUE)) {
    cat("Package", package, "is installed.\n")
  } else {
    cat("Package", package, "is not installed.\n")
  }
}
