#!/usr/bin/env Rscript
#
# dependencies.R
#
# One-shot installer for everything the LineagePhenotyping R toolchain
# needs. Re-running is safe: already-installed packages are skipped.
#
# Usage:
#   Rscript dependencies.R
#
# If you want a hermetic environment instead, switch to renv:
#   renv::init()
#   renv::install(<list below>)
#   renv::snapshot()

required_pkgs <- c(
    # Core analysis
    "ggplot2",
    "ggrepel",
    "gridExtra",
    "tidyr",
    "dplyr",
    "stringr",
    "reshape2",
    "scales",

    # Pipeline runner
    "yaml",

    # Visualization helpers
    "plotly",
    "patchwork",
    "scatterplot3d",
    "plotrix",
    "rgl",

    # Lineage / phenotyping specific
    "umap",
    "gdata"
)

# Bioconductor-only packages (different installer)
bioc_pkgs <- c(
    "ggtree"
)

# CRAN install ------------------------------------------------------------

missing_cran <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) == 0) {
    message("All CRAN dependencies already installed.")
} else {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# Bioconductor install ----------------------------------------------------

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        message("Installing BiocManager (required for: ", paste(missing_bioc, collapse = ", "), ")")
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# Final summary -----------------------------------------------------------

still_missing <- c(
    required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)],
    bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
)

if (length(still_missing) > 0) {
    stop("Could not install: ", paste(still_missing, collapse = ", "),
         "\nInstall manually and re-run dependencies.R.")
}

message("dependencies.R: all required packages are available.")
