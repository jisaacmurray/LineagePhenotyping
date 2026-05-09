#!/usr/bin/env Rscript
#
# PlotDefectTrees.R
#
# Orchestrates rendering of defect-colored lineage trees for one mutant
# dataset using LIVEtools::plot_lineage_tree() + the per-defect adapters
# in DefectScoreFrames.R. Produces ~80 PDFs across:
#
#   * cell-cycle deviation (mean/max across embryos + per-embryo)
#   * position deviation per-cell mean/max + per-timepoint magnitude +
#     per-timepoint AP-signed + per-timepoint radial-signed
#   * division-orientation deviation
#   * (optional) any user-supplied arbitrary phenotype CSV
#
# When `subdir_layout = "by_kind"` (default), outputs land under
# `<output_dir>/plots/<cc|position/cell|position/per_t|position/arrows|
# angle|expression|boxplots|summary>/[per_embryo/]`. See plot_paths.R.

suppressPackageStartupMessages({
    library(LIVEtools)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
})


# -----------------------------------------------------------------------------
# Defect-specific color schemes (Phase 5.1 + 5.2)
# -----------------------------------------------------------------------------

# `defect_scheme(kind, observed)` returns a list of args spliced into
# plot_lineage_tree() to control color_low/high, colors, color_values,
# value_min, value_max, truncate_to_last_data, drop_empty_cells per defect.
#
# Hard caps on outer endpoints:
# - cc_dev: ±20 minutes (large negatives are typically technical errors)
# - cc_dev_z: ±10
# - position_dev_AP: ±10 µm
# - position_dev_radial_signed: ±8 µm
#
# Outlier-quantile (99%) clip for unsigned kinds (position, dot) so a
# single outlier doesn't blow out the upper end of the scale.
.outlier_clip_hi <- function(observed, q = 0.99) {
    if (length(observed) == 0) return(1)
    v <- observed[!is.na(observed) & is.finite(observed)]
    if (length(v) == 0) return(1)
    as.numeric(quantile(v, q))
}

defect_scheme <- function(kind = c("position_dev", "position_dev_z",
                                    "position_dev_AP",
                                    "position_dev_radial_signed",
                                    "cc_dev", "cc_dev_z",
                                    "dot_dev", "generic"),
                          observed = NULL,
                          neutral_band = NULL) {
    kind <- match.arg(kind)

    switch(kind,
        position_dev = list(
            value_min = 3,
            value_max = max(10, .outlier_clip_hi(observed, 0.99)),
            colors = c("grey90", "blue", "yellow", "red"),
            truncate_to_last_data = TRUE
        ),
        position_dev_z = list(
            value_min = -5, value_max = 5,
            colors = c("green", "grey90", "grey90", "red"),
            color_values = c(-5, -2, 2, 5),
            truncate_to_last_data = TRUE
        ),
        position_dev_AP = list(
            value_min = -10, value_max = 10,
            colors = c("green", "grey90", "grey90", "red"),
            color_values = c(-10,
                             if (!is.null(neutral_band)) -neutral_band else -3,
                             if (!is.null(neutral_band))  neutral_band else  3,
                             10),
            truncate_to_last_data = TRUE
        ),
        position_dev_radial_signed = list(
            value_min = -8, value_max = 8,
            colors = c("green", "grey90", "grey90", "red"),
            color_values = c(-8,
                             if (!is.null(neutral_band)) -neutral_band else -3,
                             if (!is.null(neutral_band))  neutral_band else  3,
                             8),
            truncate_to_last_data = TRUE
        ),
        cc_dev = list(
            value_min = -20, value_max = 20,
            colors = c("green", "grey90", "grey90", "red"),
            color_values = c(-20,
                             if (!is.null(neutral_band)) -neutral_band else -5,
                             if (!is.null(neutral_band))  neutral_band else  5,
                             20),
            truncate_to_last_data = TRUE,
            drop_empty_cells = TRUE
        ),
        cc_dev_z = list(
            value_min = -10, value_max = 10,
            colors = c("green", "grey90", "grey90", "red"),
            color_values = c(-10,
                             if (!is.null(neutral_band)) -neutral_band else -3,
                             if (!is.null(neutral_band))  neutral_band else  3,
                             10),
            truncate_to_last_data = TRUE,
            drop_empty_cells = TRUE
        ),
        dot_dev = list(
            value_min = 0,
            value_max = max(0.5, .outlier_clip_hi(observed, 0.99)),
            colors = c("grey90", "blue", "yellow", "red"),
            truncate_to_last_data = TRUE
        ),
        generic = list(
            colors = c("grey90", "blue", "yellow", "red"),
            percentile_bounds = c(0.05, 0.95),
            truncate_to_last_data = TRUE
        )
    )
}


# Wrap a named list of schemes so .tree_to_pdf knows to render a row per
# scheme. Phase 5.1 used this for the B=5/B=3 CC double-tree; Phase 5.2
# settled on B=5 only and dropped the double-tree, but the helper stays
# as a building block for future use.
multi_scheme <- function(...) {
    x <- list(...)
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
        stop("multi_scheme entries must be named")
    }
    structure(x, class = c("defect_multi_scheme", "list"))
}


# Build a single tree (ggplot) for one defect score.
.build_tree <- function(adapter_result, title,
                        root, end_time_milestone,
                        scheme,
                        na_color = "grey80",
                        smoothing_w = 0,
                        branch_width = 1.4,
                        plot_height_in = 12,
                        lineage_subset = NULL) {
    if (nrow(adapter_result$df) == 0) {
        warning("Empty defect-score frame: ", title)
        return(NULL)
    }
    args <- c(list(
        CD = adapter_result$df,
        root = root,
        value_col = "value",
        resolution = adapter_result$resolution,
        smoothing_w = smoothing_w,
        lineage_subset = lineage_subset,
        end_time_milestone = end_time_milestone,
        na_color = na_color,
        branch_width = branch_width,
        plot_height_in = plot_height_in,
        title = title
    ), scheme)
    do.call(LIVEtools::plot_lineage_tree, args)
}


# Render a defect tree as one PDF — single plot, multi-lineage row, or
# multi-row grid (when `scheme` is a `defect_multi_scheme`, one row per
# sub-scheme; combined with `split_roots` columns this gives an
# (n_schemes × n_lineages) grid).
.tree_to_pdf <- function(adapter_result, file, title,
                         root, end_time_milestone,
                         scheme,
                         na_color = "grey80",
                         smoothing_w = 0,
                         branch_width = 1.4,
                         split_roots = NULL,
                         split_ncol = 3,
                         width = 11, height = 8.5,
                         split_width = 24, split_height = 14) {
    if (nrow(adapter_result$df) == 0) {
        warning("Empty defect-score frame; skipping ", file)
        return(invisible(NULL))
    }

    is_multi <- inherits(scheme, "defect_multi_scheme")
    schemes <- if (is_multi) scheme else list(. = scheme)

    if ((is.null(split_roots) || length(split_roots) == 0) && !is_multi) {
        p <- .build_tree(adapter_result, title, root, end_time_milestone,
                         scheme = scheme,
                         na_color = na_color, smoothing_w = smoothing_w,
                         branch_width = branch_width,
                         plot_height_in = height)
        ggplot2::ggsave(file, plot = p, width = width, height = height,
                        units = "in", device = "pdf")
        return(invisible(file))
    }

    if (is.null(split_roots) || length(split_roots) == 0) {
        lineage_list <- list(list(label = "", subset = NULL, root_for_panel = root))
    } else {
        lineage_list <- lapply(split_roots, function(sr) {
            list(label = sr, subset = sr, root_for_panel = sr)
        })
    }
    n_cols <- length(lineage_list)
    n_rows <- length(schemes)
    panel_h <- split_height / n_rows

    panels <- vector("list", n_rows * n_cols)
    panel_idx <- 0
    for (s_name in names(schemes)) {
        sub_scheme <- schemes[[s_name]]
        s_label <- if (is_multi) paste0(" [", s_name, "]") else ""
        for (lp in lineage_list) {
            panel_idx <- panel_idx + 1
            label <- if (nzchar(lp$label)) {
                paste0(lp$label, " lineage", s_label)
            } else {
                trimws(s_label)
            }
            panels[[panel_idx]] <- .build_tree(
                adapter_result, label,
                root = lp$root_for_panel, end_time_milestone = end_time_milestone,
                scheme = sub_scheme,
                na_color = na_color, smoothing_w = smoothing_w,
                branch_width = branch_width,
                plot_height_in = panel_h,
                lineage_subset = lp$subset)
        }
    }
    panels <- panels[!vapply(panels, is.null, logical(1))]
    if (length(panels) == 0) {
        warning("No panels produced: ", title)
        return(invisible(NULL))
    }
    grid <- patchwork::wrap_plots(panels, ncol = n_cols, nrow = n_rows) +
        patchwork::plot_annotation(
            title = title,
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(size = 18, face = "bold")))
    ggplot2::ggsave(file, plot = grid,
                    width = split_width, height = split_height,
                    units = "in", device = "pdf", limitsize = FALSE)
    invisible(file)
}


# -----------------------------------------------------------------------------
# Top-level orchestrator
# -----------------------------------------------------------------------------

#' Render a full set of defect-colored lineage trees for a mutant dataset.
#'
#' @param name dataset name (used as filename prefix).
#' @param devs list returned by `AnalyzePositions()`. Pass `NULL` to skip
#'   the per-timepoint position-deviation tree.
#' @param ccdevs_path path to `<Name>_ccDevs.csv`. `NULL` = derived from
#'   `output_dir` + `name`.
#' @param dots_path path to `<Name>_dots.csv`. `NULL` = derived.
#' @param positions_path path to `<Name>positions.txt`. Required.
#' @param mean_pos_dev_csv path to `<Name>_CellMeanPositionDevs.csv`. `NULL` = derived.
#' @param max_pos_dev_csv path to `<Name>_CellMaxPositionDevs.csv`. `NULL` = derived.
#' @param data_dir,output_dir path roots.
#' @param wt_ref_dir directory containing the WT rotated CSVs (needed
#'   for AP-signed and radial-signed per-timepoint trees).
#' @param root tree root for the aggregated plots.
#' @param end_time_milestone optional `list(lineage = <name>, n_cells = <int>)`
#'   to crop plots at a developmental milestone.
#' @param split_roots character vector of per-panel roots (default ABa/ABp/P1).
#' @param arbitrary_scores optional list of user-supplied phenotype tables.
#' @param per_embryo emit one PDF per mutant embryo for each defect kind.
#' @return character vector of paths to the new PDFs.
PlotDefectTrees <- function(name,
                            devs = NULL,
                            ccdevs_path = NULL,
                            dots_path = NULL,
                            positions_path = NULL,
                            mean_pos_dev_csv = NULL,
                            max_pos_dev_csv = NULL,
                            data_dir = ".",
                            output_dir = ".",
                            wt_ref_dir = NULL,
                            root = "P0",
                            end_time_milestone = NULL,
                            split_roots = c("ABa", "ABp", "P1"),
                            arbitrary_scores = NULL,
                            per_embryo = TRUE,
                            min_cells_per_timepoint = NULL,
                            na_color = "grey80",
                            split_ncol = 3,
                            split_width = 24, split_height = 14,
                            single_width = 14, single_height = 8.5) {
    # Phase 5.3H2(a): cliff filter for the per-timepoint trees only.
    # Read from option if not passed explicitly so YAML configs can set it.
    if (is.null(min_cells_per_timepoint)) {
        min_cells_per_timepoint <- getOption("LineagePhenotyping.min_cells_per_timepoint", NULL)
    }
    .default <- function(x, def) if (is.null(x)) def else x

    # --- resolve default paths ---
    ccdevs_path      <- .default(ccdevs_path,      file.path(output_dir, paste0(name, "_ccDevs.csv")))
    dots_path        <- .default(dots_path,        file.path(output_dir, paste0(name, "_dots.csv")))
    mean_pos_dev_csv <- .default(mean_pos_dev_csv, file.path(output_dir, paste0(name, "_CellMeanPositionDevs.csv")))
    max_pos_dev_csv  <- .default(max_pos_dev_csv,  file.path(output_dir, paste0(name, "_CellMaxPositionDevs.csv")))
    positions_path   <- .default(positions_path,   file.path(data_dir, name, paste0(name, "positions.txt")))

    if (!file.exists(positions_path)) {
        stop("PlotDefectTrees: positions file not found at ", positions_path,
             " — pass `positions_path` explicitly")
    }

    out_files <- c()

    # Local dispatcher.
    .pdf_one <- function(adapter, file, title, scheme, smoothing_w = 0) {
        .tree_to_pdf(
            adapter, file, title,
            root = root,
            end_time_milestone = end_time_milestone,
            scheme = scheme,
            na_color = na_color,
            smoothing_w = smoothing_w,
            split_roots = split_roots,
            split_ncol = split_ncol,
            width = single_width, height = single_height,
            split_width = split_width, split_height = split_height
        )
    }

    # Render aggregated + per-embryo PDFs for one defect kind.
    .render_kind <- function(kind, agg_adapter_fn,
                              base_label, file_prefix,
                              modes = c("mean", "max"),
                              scheme_fn = NULL,
                              subdir = NULL) {
        if (is.null(scheme_fn)) {
            scheme_fn <- function(observed) defect_scheme(kind, observed = observed)
        }
        agg_dir <- .plots_dir(output_dir, subdir)
        emb_dir <- .plots_dir(output_dir, subdir, per_embryo = TRUE)
        for (mode in modes) {
            adapter <- tryCatch(agg_adapter_fn(mode),
                error = function(e) { message("  error: ", conditionMessage(e)); NULL })
            if (is.null(adapter) || nrow(adapter$df) == 0) next
            obs <- adapter$df$value[!is.na(adapter$df$value)]
            scheme <- scheme_fn(obs)
            f <- file.path(agg_dir, paste0(name, "_", file_prefix, "_", mode, ".pdf"))
            title <- paste0(name, ": ", base_label, " (", mode, " across embryos)")
            .pdf_one(adapter, f, title, scheme = scheme)
            out_files <<- c(out_files, f)
        }

        if (per_embryo) {
            adapter_emb <- tryCatch(agg_adapter_fn("none"),
                error = function(e) { message("  error: ", conditionMessage(e)); NULL })
            if (!is.null(adapter_emb) && nrow(adapter_emb$df) > 0
                && "embryo" %in% colnames(adapter_emb$df)) {
                obs <- adapter_emb$df$value[!is.na(adapter_emb$df$value)]
                scheme <- scheme_fn(obs)
                emb_ids <- unique(adapter_emb$df$embryo)
                emb_ids <- emb_ids[!is.na(emb_ids)]
                for (emb in emb_ids) {
                    sub <- list(df = adapter_emb$df[adapter_emb$df$embryo == emb, ],
                                 resolution = adapter_emb$resolution)
                    if (sum(!is.na(sub$df$value)) == 0) next
                    safe_emb <- gsub("[^A-Za-z0-9._-]", "_", emb)
                    f <- file.path(emb_dir,
                                   paste0(name, "_", file_prefix, "_per_embryo_", safe_emb, ".pdf"))
                    title <- paste0(name, ": ", base_label, " (", emb, ")")
                    .pdf_one(sub, f, title, scheme = scheme)
                    out_files <<- c(out_files, f)
                }
            }
        }
    }

    # --- CC deviation (single ±5 grey band, hard ±20 cap) ---
    if (file.exists(ccdevs_path)) {
        message("[PlotDefectTrees] CC deviation (±5 grey band) ...")
        .render_kind(
            kind = "cc_dev",
            agg_adapter_fn = function(mode) {
                cc_dev_long(ccdevs_path, positions_path, mode = "deviation",
                            aggregate_embryos = mode)
            },
            base_label = "CC deviation",
            file_prefix = "DefectTrees_CCDev",
            modes = c("mean", "max"),
            subdir = "cc"
        )
    } else {
        message("[PlotDefectTrees] skipping CC trees (file missing): ", ccdevs_path)
    }

    # --- Position deviation (per-cell, from CSVs) ---
    if (file.exists(mean_pos_dev_csv)) {
        message("[PlotDefectTrees] Position deviation (CellMean*.csv) ...")
        .render_kind(
            kind = "position_dev",
            agg_adapter_fn = function(mode) {
                position_dev_long_from_csv(mean_pos_dev_csv, positions_path,
                                           aggregate_embryos = mode)
            },
            base_label = "mean position deviation (um)",
            file_prefix = "DefectTrees_PositionDevMean",
            modes = c("mean"),
            subdir = "position/cell"
        )
    } else {
        message("[PlotDefectTrees] skipping mean PositionDev tree (file missing): ",
                mean_pos_dev_csv)
    }
    if (file.exists(max_pos_dev_csv)) {
        message("[PlotDefectTrees] Position deviation (CellMax*.csv) ...")
        .render_kind(
            kind = "position_dev",
            agg_adapter_fn = function(mode) {
                position_dev_long_from_csv(max_pos_dev_csv, positions_path,
                                           aggregate_embryos = mode)
            },
            base_label = "max position deviation (um)",
            file_prefix = "DefectTrees_PositionDevMax",
            modes = c("max"),
            subdir = "position/cell"
        )
    } else {
        message("[PlotDefectTrees] skipping max PositionDev tree (file missing): ",
                max_pos_dev_csv)
    }

    # --- Position deviation per-timepoint magnitude ---
    posdevs_path <- file.path(output_dir, paste0(name, "_PositionDevs.csv"))
    if (file.exists(posdevs_path)) {
        if (!is.null(min_cells_per_timepoint)) {
            message("[PlotDefectTrees] Position deviation per-timepoint ",
                    "(from PositionDevs.csv, min_cells=",
                    min_cells_per_timepoint, ") ...")
        } else {
            message("[PlotDefectTrees] Position deviation per-timepoint (from PositionDevs.csv) ...")
        }
        .render_kind(
            kind = "position_dev",
            agg_adapter_fn = function(mode) {
                position_dev_per_t_long(posdevs_path, aggregate_embryos = mode,
                                         min_cells_per_timepoint = min_cells_per_timepoint)
            },
            base_label = "instantaneous position deviation (um)",
            file_prefix = "DefectTrees_PositionDev_per_t",
            modes = c("mean"),
            subdir = "position/per_t"
        )
    } else {
        message("[PlotDefectTrees] skipping per-timepoint PositionDev (",
                posdevs_path, " missing)")
    }

    # --- Direction-encoded position deviation (signed AP + signed radial) ---
    rotated_x <- file.path(output_dir, paste0(name, "_rotatedX.csv"))
    if (file.exists(rotated_x) && !is.null(wt_ref_dir) && dir.exists(wt_ref_dir)) {
        wt_basename <- basename(wt_ref_dir)
        message("[PlotDefectTrees] Position deviation per-timepoint AP_signed ...")
        .render_kind(
            kind = "position_dev_AP",
            agg_adapter_fn = function(mode) {
                position_dev_per_t_components_long(
                    mut_rotated_dir = output_dir, mut_name = name,
                    wt_rotated_dir = wt_ref_dir, wt_name = wt_basename,
                    kind = "AP_signed", aggregate_embryos = mode,
                    min_cells_per_timepoint = min_cells_per_timepoint)
            },
            base_label = "AP-signed position deviation (um, +post / -ant)",
            file_prefix = "DefectTrees_PositionDev_AP_per_t",
            modes = c("mean"),
            subdir = "position/per_t"
        )

        message("[PlotDefectTrees] Position deviation per-timepoint radial_signed ...")
        .render_kind(
            kind = "position_dev_radial_signed",
            agg_adapter_fn = function(mode) {
                position_dev_per_t_components_long(
                    mut_rotated_dir = output_dir, mut_name = name,
                    wt_rotated_dir = wt_ref_dir, wt_name = wt_basename,
                    kind = "radial_signed", aggregate_embryos = mode,
                    min_cells_per_timepoint = min_cells_per_timepoint)
            },
            base_label = "radial-signed position deviation (um, +out / -in)",
            file_prefix = "DefectTrees_PositionDev_radial_per_t",
            modes = c("mean"),
            subdir = "position/per_t"
        )
    } else {
        message("[PlotDefectTrees] skipping AP/radial signed trees (rotated CSVs or wt_ref_dir missing)")
    }

    # --- Division-orientation deviation ---
    if (file.exists(dots_path)) {
        message("[PlotDefectTrees] Division-orientation deviation ...")
        .render_kind(
            kind = "dot_dev",
            agg_adapter_fn = function(mode) {
                dot_dev_long(dots_path, positions_path,
                             mode = "abs_diff",
                             aggregate_embryos = mode)
            },
            base_label = "|Mutant_Mean - WT_Mean| division dot",
            file_prefix = "DefectTrees_DotDev",
            modes = c("mean"),
            subdir = "angle"
        )
    } else {
        message("[PlotDefectTrees] skipping DotDev tree (file missing): ", dots_path)
    }

    # --- arbitrary user-supplied scores ---
    if (!is.null(arbitrary_scores) && length(arbitrary_scores) > 0) {
        for (sc in arbitrary_scores) {
            label <- sc$label %||% basename(tools::file_path_sans_ext(sc$file))
            message("[PlotDefectTrees] Arbitrary score: ", label)
            adapter <- tryCatch(
                arbitrary_long(
                    sc$file,
                    cell_col  = sc$cell_col  %||% "cell",
                    time_col  = sc$time_col  %||% "time",
                    value_col = sc$value_col %||% "value",
                    positions_path = positions_path),
                error = function(e) { message("  error: ", conditionMessage(e)); NULL })
            if (!is.null(adapter)) {
                f <- file.path(.plots_dir(output_dir, "arbitrary"),
                               paste0(name, "_DefectTrees_Arbitrary_", label, ".pdf"))
                .pdf_one(adapter, f, paste0(name, ": ", label),
                         scheme = defect_scheme("generic"))
                out_files <- c(out_files, f)
            }
        }
    }

    message("[PlotDefectTrees] Done. ", length(out_files), " files written.")
    invisible(out_files)
}


# `%||%` for callers; defined here in case the runner hasn't already.
if (!exists("%||%", mode = "function")) {
    `%||%` <- function(a, b) if (is.null(a)) b else a
}


#' Re-render defect trees for one or more sublineages without re-running the pipeline.
#'
#' Reads the defect-score CSVs already on disk and calls `PlotDefectTrees()`
#' with `split_roots = sublineage`, writing into
#' `<output_dir>/plots/_subsets/<sublineage_safe>/`.
#'
#' @param name dataset name.
#' @param output_dir directory containing the existing defect CSVs.
#' @param sublineage character vector — each entry produces one subset.
#' @param kinds optional defect-kind filter (NULL = all kinds with data).
#' @param modes "mean", "max", or both.
#' @param per_embryo emit per-embryo PDFs in addition to aggregated.
#' @param wt_ref_dir,positions_path optional overrides.
#' @param ... forwarded to PlotDefectTrees().
replot_defect_subset <- function(name, output_dir, sublineage,
                                  kinds = NULL,
                                  modes = c("mean"),
                                  per_embryo = FALSE,
                                  wt_ref_dir = NULL,
                                  positions_path = NULL,
                                  ...) {
    if (!dir.exists(output_dir)) {
        stop("output_dir not found: ", output_dir)
    }
    if (length(sublineage) == 0) stop("sublineage is empty")

    out_files <- character()
    for (sl in sublineage) {
        safe_sl <- gsub("[^A-Za-z0-9._-]", "_", sl)
        subset_root <- file.path(output_dir, "plots", "_subsets", safe_sl)
        dir.create(subset_root, recursive = TRUE, showWarnings = FALSE)

        # PlotDefectTrees() reads <output_dir>/<name>_PositionDevs.csv,
        # <name>_rotatedX/Y/Z.csv etc. via paths anchored on `output_dir`.
        # Because we redirect `output_dir` to the subset folder so writes
        # land there, link the source data files in so the existence
        # checks for per-t / AP / radial blocks succeed without copying
        # the (potentially large) CSVs.
        for (data_file in c(paste0(name, "_PositionDevs.csv"),
                             paste0(name, "_rotatedX.csv"),
                             paste0(name, "_rotatedY.csv"),
                             paste0(name, "_rotatedZ.csv"))) {
            src <- file.path(output_dir, data_file)
            dst <- file.path(subset_root, data_file)
            if (file.exists(src) && !file.exists(dst)) {
                tryCatch(file.symlink(src, dst),
                         error = function(e) {
                             file.copy(src, dst, overwrite = FALSE)
                         })
            }
        }

        message("[replot_defect_subset] ", name, " / ", sl,
                " -> ", subset_root)

        ccdevs_path      <- file.path(output_dir, paste0(name, "_ccDevs.csv"))
        dots_path        <- file.path(output_dir, paste0(name, "_dots.csv"))
        mean_pos_dev_csv <- file.path(output_dir, paste0(name, "_CellMeanPositionDevs.csv"))
        max_pos_dev_csv  <- file.path(output_dir, paste0(name, "_CellMaxPositionDevs.csv"))

        if (is.null(positions_path)) {
            candidates <- c(
                file.path(dirname(output_dir), name, paste0(name, "positions.txt")),
                file.path(output_dir, paste0(name, "positions.txt"))
            )
            hit <- candidates[file.exists(candidates)]
            if (length(hit) == 0) {
                stop("Couldn't find <name>positions.txt; pass positions_path explicitly")
            }
            positions_path_use <- hit[1]
        } else {
            positions_path_use <- positions_path
        }

        if (!is.null(kinds)) {
            if (!("cc_dev" %in% kinds))             ccdevs_path      <- ""
            if (!("dot_dev" %in% kinds))            dots_path        <- ""
            if (!any(c("position_dev") %in% kinds)) mean_pos_dev_csv <- max_pos_dev_csv <- ""
        }

        prev_layout <- getOption("LineagePhenotyping.subdir_layout", "by_kind")
        options(LineagePhenotyping.subdir_layout = "by_kind")
        on.exit(options(LineagePhenotyping.subdir_layout = prev_layout), add = TRUE)

        new_files <- PlotDefectTrees(
            name = name,
            ccdevs_path = if (nzchar(ccdevs_path)) ccdevs_path else NULL,
            dots_path = if (nzchar(dots_path)) dots_path else NULL,
            positions_path = positions_path_use,
            mean_pos_dev_csv = if (nzchar(mean_pos_dev_csv)) mean_pos_dev_csv else NULL,
            max_pos_dev_csv = if (nzchar(max_pos_dev_csv)) max_pos_dev_csv else NULL,
            data_dir = output_dir,
            output_dir = subset_root,
            wt_ref_dir = wt_ref_dir,
            split_roots = sl,
            per_embryo = per_embryo,
            ...
        )
        out_files <- c(out_files, new_files)
    }

    invisible(out_files)
}
