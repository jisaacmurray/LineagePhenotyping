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
    "gdata",

    # Phase 5 — defect-tree plotting
    "ape",
    "remotes"
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

# GitHub install ---------------------------------------------------------
# LIVEtools is the johnmurraylab tree-/3D-/trajectory-plotting package.
# Phase 5 (defect-colored lineage trees) calls LIVEtools::plot_lineage_tree()
# and LIVEtools::paginate_tree_plots(). We install from the canonical
# johnmurraylab/LIVE_tools repo (NOT LIVEtools-paper, which is a static
# archive of the original publication).

if (!requireNamespace("LIVEtools", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    message("Installing LIVEtools from johnmurraylab/LIVE_tools (GitHub)...")
    remotes::install_github("johnmurraylab/LIVE_tools",
                            upgrade = "never", quiet = TRUE)
}

if (!requireNamespace("LIVEtools", quietly = TRUE)) {
    stop("Could not install LIVEtools.\n",
         "  Try manually:\n",
         "    remotes::install_github(\"johnmurraylab/LIVE_tools\")\n",
         "  If you see auth errors, ensure GITHUB_PAT is set or the repo is public.")
}

message("dependencies.R: all required packages are available.")
