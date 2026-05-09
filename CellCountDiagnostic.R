#!/usr/bin/env Rscript
#
# CellCountDiagnostic.R
#
# Phase 5.3H1: characterize the "cliff" pattern where mutant embryo
# tracking degrades sharply at late timepoints — fewer cells survive
# in the recorded trace, the rotation alignment is biased by the
# remaining (often mispositioned) cells, and per-timepoint position
# deviation trees show artefactual spikes for the surviving cells.
#
# Provides two helpers:
#   analyze_cell_counts_per_t(positions_path)  — long frame of n_cells
#                                                 alive per (time, embryo)
#   plot_cell_counts_per_t(...)                — overlay plot for spotting
#                                                 cliffs visually + a per-
#                                                 embryo small-multiples
#                                                 page that's easier to
#                                                 read on dense datasets
#
# Run as a standalone script: produces
#   <output_dir>/plots/_diagnostics/<name>_cell_counts_vs_t.pdf
#   <output_dir>/<name>_cell_counts_vs_t.csv
#
# Usage (from R):
#   source("CellCountDiagnostic.R")
#   df <- analyze_cell_counts_per_t(
#             "/path/to/JIM721/JIM721positions.txt",
#             dataset_name = "JIM721")
#   plot_cell_counts_per_t(df, output_dir = "/path/to/JIM721_phase2_check")

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})


#' Per-(time, embryo) count of cells with non-NA position data.
#'
#' Reads the wide positions.txt format (columns are X_<emb>, Y_<emb>,
#' Z_<emb> per embryo) and returns a long frame:
#'   data.frame(time, embryo, n_cells, dataset)
#'
#' "Alive" here = an embryo's X column for that (cell, time) is not NA.
#' Matches what the rest of the pipeline considers a tracked cell.
#'
#' @param positions_path path to <name>positions.txt
#' @param dataset_name used to populate the `dataset` column (helpful
#'   when overlaying multiple datasets on one plot).
analyze_cell_counts_per_t <- function(positions_path, dataset_name = NULL) {
    if (!file.exists(positions_path)) {
        stop("positions file not found: ", positions_path)
    }
    raw <- read.delim(positions_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("cell", "time") %in% colnames(raw))) {
        stop("positions file missing 'cell' and/or 'time' columns: ", positions_path)
    }
    cn <- colnames(raw)
    x_cols <- grep("^X_", cn, value = TRUE)
    if (length(x_cols) == 0) {
        stop("positions file has no X_<embryo> columns: ", positions_path)
    }
    embryos <- sub("^X_", "", x_cols)

    out_list <- vector("list", length(embryos))
    for (i in seq_along(embryos)) {
        e <- embryos[i]
        xc <- paste0("X_", e)
        if (!(xc %in% cn)) next
        present <- !is.na(raw[[xc]])
        # Count cells present per timepoint for this embryo.
        agg <- aggregate(present, by = list(time = raw$time), FUN = sum)
        colnames(agg)[2] <- "n_cells"
        agg$embryo <- e
        out_list[[i]] <- agg
    }
    df <- do.call(rbind, out_list)
    df$dataset <- if (is.null(dataset_name)) basename(positions_path) else dataset_name
    df[order(df$dataset, df$embryo, df$time), c("dataset", "embryo", "time", "n_cells"), drop = FALSE]
}


#' Render a cell-counts-vs-time diagnostic PDF.
#'
#' Produces two pages:
#'   1. Overlay — every embryo as a thin line, coloured by dataset.
#'   2. Per-embryo small-multiples (one panel per embryo) for inspecting
#'      individual embryos for sudden cliffs.
#'
#' Also writes the long frame as a CSV next to the PDF for downstream
#' filtering / threshold tuning.
#'
#' @param df long frame from analyze_cell_counts_per_t().
#' @param output_dir base output directory; PDF lands in
#'   `<output_dir>/plots/_diagnostics/`.
#' @param name filename prefix (default = first dataset name).
plot_cell_counts_per_t <- function(df, output_dir, name = NULL) {
    if (nrow(df) == 0) {
        warning("Empty cell-count frame; nothing to plot")
        return(invisible(NULL))
    }
    if (is.null(name)) name <- df$dataset[1]

    diag_dir <- if (exists(".plots_dir", mode = "function")) {
        .plots_dir(output_dir, "_diagnostics")
    } else {
        d <- file.path(output_dir, "plots", "_diagnostics")
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        d
    }
    pdf_path <- file.path(diag_dir, paste0(name, "_cell_counts_vs_t.pdf"))
    csv_path <- file.path(output_dir, paste0(name, "_cell_counts_vs_t.csv"))

    write.csv(df, csv_path, row.names = FALSE)
    message("[CellCountDiagnostic] wrote ", csv_path)

    # Overlay page
    p_overlay <- ggplot(df, aes(x = time, y = n_cells,
                                  group = interaction(dataset, embryo),
                                  color = dataset)) +
        geom_line(alpha = 0.6, linewidth = 0.4) +
        theme_bw() +
        labs(title = "Cells alive per timepoint (overlay)",
             subtitle = sprintf("%d embryos across %d dataset(s)",
                                 length(unique(df$embryo)),
                                 length(unique(df$dataset))),
             x = "Time (min)", y = "n cells with tracked position")

    # Small-multiples page
    p_facets <- ggplot(df, aes(x = time, y = n_cells)) +
        geom_line(linewidth = 0.4) +
        facet_wrap(~ embryo, scales = "free_x") +
        theme_bw() +
        theme(strip.text = element_text(size = 7)) +
        labs(title = "Cells alive per timepoint (per-embryo)",
             x = "Time (min)", y = "n cells")

    pdf(pdf_path, width = 11, height = 8.5)
    print(p_overlay)
    print(p_facets)
    dev.off()
    message("[CellCountDiagnostic] wrote ", pdf_path)

    invisible(c(csv_path, pdf_path))
}


# Allow running as a script: Rscript CellCountDiagnostic.R <positions.txt> <output_dir> [name]
if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        stop("Usage: Rscript CellCountDiagnostic.R <positions.txt> <output_dir> [name]")
    }
    positions_path <- args[1]
    output_dir <- args[2]
    name <- if (length(args) >= 3) args[3] else NULL

    df <- analyze_cell_counts_per_t(positions_path, dataset_name = name)
    plot_cell_counts_per_t(df, output_dir = output_dir, name = name)
}
