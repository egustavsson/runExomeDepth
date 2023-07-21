# Rrequired for Bioconductor packages
install.packages("BiocManager", repos = "https://cran.ma.imperial.ac.uk/")

# Bioconductor packages
bioc_packages <- c("Biostrings", "IRanges", "Rsamtools", "GenomicRanges", "GenomicAlignments")

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

# Cunction to install a CRAN package if it is not already installed
install_cran_if_not_installed <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = "https://cran.ma.imperial.ac.uk/")
  }
}

# Install CRAN packages
cran_packages <- c("optparse", "ExomeDepth", "magrittr")
for (package in cran_packages) {
  install_cran_if_not_installed(package)
}



