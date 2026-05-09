#!/usr/bin/env Rscript
#
# DefectScoreFrames.R
#
# Adapters that turn LineagePhenotyping's defect-score outputs into the
# canonical input format consumed by `LIVEtools::plot_lineage_tree()`:
#
#     data.frame(cell, time, value, [embryo])
#
# Every adapter returns `list(df = <frame>, resolution = <"cell"|"timepoint">)`.
# `resolution` is the recommended `plot_lineage_tree(..., resolution = ...)`
# value for the score's natural granularity (CC + dot are per-cell aggregates;
# position deviations have per-timepoint resolution when sourced from the
# `devs` list returned by `AnalyzePositions()`).
#
# Most adapters need a `<Name>positions.txt` file to know each cell's
# observed birth/end times. Per-cell aggregate scores (CC, dot, position-dev
# CSVs) are *broadcast* across each cell's lifespan so downstream
# `plot_lineage_tree` can render them with either resolution.
#
# All adapters return rows in the order produced by `read_positions_to_cd()`:
# one row per (cell, time, embryo) tuple. Use `aggregate_embryos = "mean"`
# to collapse to one value per (cell, time).

suppressPackageStartupMessages({
    library(dplyr)
})


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

#' Read a `<Name>positions.txt` file into a long-format CD-like frame.
#'
#' The on-disk format is wide: rows are (cell, time), columns are
#' `X_<embryo>`, `Y_<embryo>`, `Z_<embryo>` triples. This function
#' un-pivots the per-embryo position triples into long format, dropping
#' rows where the embryo's X position is `NA` (cell not present in that
#' embryo at that time).
#'
#' @return data.frame with columns `cell, time, x, y, z, embryo`.
read_positions_to_cd <- function(positions_path) {
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

    out <- vector("list", length(embryos))
    for (i in seq_along(embryos)) {
        e <- embryos[i]
        xc <- paste0("X_", e); yc <- paste0("Y_", e); zc <- paste0("Z_", e)
        if (!all(c(xc, yc, zc) %in% cn)) next
        sub <- raw[, c("cell", "time", xc, yc, zc), drop = FALSE]
        colnames(sub) <- c("cell", "time", "x", "y", "z")
        sub$embryo <- e
        sub <- sub[!is.na(sub$x), , drop = FALSE]
        out[[i]] <- sub
    }
    do.call(rbind, out)
}


# Internal helper: broadcast a (cell -> value) lookup across each cell's
# observed timepoints, optionally keeping `embryo` as a column.
.broadcast_per_cell <- function(per_cell, cd, keep_embryo = FALSE) {
    if (keep_embryo && "embryo" %in% colnames(cd)) {
        out <- merge(cd[, c("cell", "time", "embryo")], per_cell,
                     by = "cell", all.x = TRUE)
    } else {
        out <- merge(unique(cd[, c("cell", "time")]), per_cell,
                     by = "cell", all.x = TRUE)
    }
    out[order(out$cell, out$time), , drop = FALSE]
}


# Internal helper: aggregate a numeric vector with NA-safe handling.
.agg_fn <- function(mode = c("mean", "max", "median", "sum", "first", "last")) {
    mode <- match.arg(mode)
    switch(mode,
        mean   = function(x) mean(x, na.rm = TRUE),
        max    = function(x) {
                     y <- x[!is.na(x)]
                     if (length(y) == 0) NA_real_ else max(y)
                 },
        median = function(x) median(x, na.rm = TRUE),
        sum    = function(x) sum(x, na.rm = TRUE),
        first  = function(x) { y <- x[!is.na(x)]; if (length(y) == 0) NA_real_ else y[1] },
        last   = function(x) { y <- x[!is.na(x)]; if (length(y) == 0) NA_real_ else y[length(y)] }
    )
}


# -----------------------------------------------------------------------------
# Cell-cycle deviation (per-cell)
# -----------------------------------------------------------------------------

#' Build a long defect frame from `<Name>_ccDevs.csv`.
#'
#' `<Name>_ccDevs.csv` (written by `AnalyzeDivTimes.R:208`) has one row per
#' `(cell, embryo)` with columns including `Cell, Embryo, CC Deviation,
#' CC Z score`. This adapter aggregates across embryos and broadcasts each
#' cell's value across its observed timepoints from the positions file.
#'
#' Handles both the space-separated column-name variant (`"CC Deviation"`)
#' and the underscore-separated variant (`"CC_Deviation"`) that
#' AnalyzeDivTimes has produced over time.
#'
#' @param ccdevs_path path to `<Name>_ccDevs.csv`.
#' @param positions_path path to `<Name>positions.txt`.
#' @param mode which column to use: `"deviation"` (CC Deviation, default),
#'   `"z"` (CC Z score), or `"cc_diff"` (CC - WT CC).
#' @param aggregate_embryos `"mean"` (default), `"max"`, or `"none"`.
#'   `"none"` keeps the embryo dimension (one row per cell-time-embryo).
#' @return list(df = data.frame(cell, time, value, [embryo]),
#'              resolution = "cell")
cc_dev_long <- function(ccdevs_path, positions_path,
                        mode = c("deviation", "z", "cc_diff"),
                        aggregate_embryos = c("mean", "max", "none")) {
    mode <- match.arg(mode)
    aggregate_embryos <- match.arg(aggregate_embryos)
    if (!file.exists(ccdevs_path)) {
        stop("ccdevs file not found: ", ccdevs_path)
    }
    ccd <- read.csv(ccdevs_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!("Cell" %in% colnames(ccd))) {
        stop("ccdevs file lacks `Cell` column: ", ccdevs_path)
    }

    # Column-name variants: AnalyzeDivTimes.R has produced both
    # space-separated ("CC Deviation") and underscore ("CC_Deviation")
    # variants over time. Pick whichever is present.
    .pick_col <- function(cands) {
        hit <- cands[cands %in% colnames(ccd)]
        if (length(hit) == 0) NULL else hit[1]
    }
    cc_col   <- .pick_col(c("CC",            "CC"))
    wtcc_col <- .pick_col(c("WT CC",         "WT_CC"))
    dev_col  <- .pick_col(c("CC Deviation",  "CC_Deviation"))
    z_col    <- .pick_col(c("CC Z score",    "CC_Z_score"))

    if (mode == "deviation") {
        if (is.null(dev_col)) stop("ccdevs file lacks `CC Deviation`/`CC_Deviation`")
        ccd$.value <- ccd[[dev_col]]
    } else if (mode == "z") {
        if (is.null(z_col)) stop("ccdevs file lacks `CC Z score`/`CC_Z_score`")
        ccd$.value <- ccd[[z_col]]
    } else if (mode == "cc_diff") {
        if (is.null(cc_col) || is.null(wtcc_col)) {
            stop("ccdevs file lacks CC and/or WT CC columns")
        }
        ccd$.value <- ccd[[cc_col]] - ccd[[wtcc_col]]
    }

    cd <- read_positions_to_cd(positions_path)

    if (aggregate_embryos == "none") {
        # Keep per-embryo resolution; broadcast each (cell, embryo)'s value
        # across that embryo's observed timepoints.
        per_cell_emb <- ccd[, c("Cell", "Embryo", ".value"), drop = FALSE]
        colnames(per_cell_emb) <- c("cell", "embryo", "value")
        # ccdevs uses periods; positions.txt may use periods OR underscores
        # (e.g. "X_20230307_JIM721_tab.1_L4" vs "X20230307_JIM721_tab.1_L2").
        # Normalize: strip leading "X" + remove dots, underscores, hyphens
        # before joining. Try strict join first; if no rows match, try the
        # normalized fallback.
        strict <- merge(cd[, c("cell", "time", "embryo")], per_cell_emb,
                        by = c("cell", "embryo"), all.x = TRUE)
        if (sum(!is.na(strict$value)) == 0) {
            cd$.embryo_norm <- gsub("[._-]", "", sub("^X", "", cd$embryo))
            per_cell_emb$.embryo_norm <- gsub("[._-]", "", sub("^X", "", per_cell_emb$embryo))
            merged <- merge(cd[, c("cell", "time", "embryo", ".embryo_norm")],
                            per_cell_emb[, c("cell", ".embryo_norm", "value")],
                            by = c("cell", ".embryo_norm"), all.x = TRUE)
            merged$.embryo_norm <- NULL
            return(list(df = merged[order(merged$cell, merged$time), ],
                        resolution = "cell"))
        }
        return(list(df = strict[order(strict$cell, strict$time), ],
                    resolution = "cell"))
    }

    fn <- .agg_fn(aggregate_embryos)
    per_cell <- aggregate(ccd$.value, by = list(cell = ccd$Cell), FUN = fn)
    colnames(per_cell)[2] <- "value"
    list(df = .broadcast_per_cell(per_cell, cd, keep_embryo = FALSE),
         resolution = "cell")
}


# -----------------------------------------------------------------------------
# Position deviation (per-cell, from CellMean/Max CSVs)
# -----------------------------------------------------------------------------

#' Build a long defect frame from a CellMean/Max/MaxPositionDevs CSV.
#'
#' `<Name>_CellMeanPositionDevs.csv` etc. have rows = cells, columns =
#' embryos. This adapter aggregates across embryos and broadcasts each
#' cell's value across its observed timepoints from the positions file.
#'
#' @param csv_path path to `<Name>_CellMeanPositionDevs.csv` (or `_CellMax...`).
#' @param positions_path path to `<Name>positions.txt`.
#' @param aggregate_embryos `"mean"` (default), `"max"`, `"median"`, or `"none"`.
#' @return list(df = data.frame(cell, time, value, [embryo]),
#'              resolution = "cell")
position_dev_long_from_csv <- function(csv_path, positions_path,
                                        aggregate_embryos = c("mean", "max",
                                                              "median", "none")) {
    aggregate_embryos <- match.arg(aggregate_embryos)
    if (!file.exists(csv_path)) {
        stop("position-dev CSV not found: ", csv_path)
    }
    pd <- read.csv(csv_path, stringsAsFactors = FALSE,
                   check.names = FALSE, row.names = 1)
    # Some CSVs prepend a duplicate "Cells" column; drop if present.
    if ("Cells" %in% colnames(pd)) pd$Cells <- NULL

    cd <- read_positions_to_cd(positions_path)

    if (aggregate_embryos == "none") {
        # Keep per-embryo: long-pivot the wide matrix.
        long <- tidyr::pivot_longer(
            data.frame(cell = rownames(pd), pd, check.names = FALSE),
            cols = -.data$cell, names_to = "embryo", values_to = "value")
        # Embryo-name normalization (strip "X" prefix, remove ./_/-)
        cd$.embryo_norm  <- gsub("[._-]", "", sub("^X", "", cd$embryo))
        long$.embryo_norm <- gsub("[._-]", "", sub("^X", "", long$embryo))
        merged <- merge(cd[, c("cell", "time", "embryo", ".embryo_norm")],
                        long[, c("cell", ".embryo_norm", "value")],
                        by = c("cell", ".embryo_norm"), all.x = TRUE)
        merged$.embryo_norm <- NULL
        return(list(df = merged[order(merged$cell, merged$time), ],
                    resolution = "cell"))
    }

    fn <- .agg_fn(aggregate_embryos)
    per_cell_vals <- apply(as.matrix(pd), 1, fn)
    per_cell <- data.frame(cell = rownames(pd), value = per_cell_vals,
                            stringsAsFactors = FALSE)
    list(df = .broadcast_per_cell(per_cell, cd, keep_embryo = FALSE),
         resolution = "cell")
}


# -----------------------------------------------------------------------------
# Position deviation (per-timepoint, from PositionDevs.csv)
# -----------------------------------------------------------------------------

#' Build a per-timepoint position-deviation frame from `<Name>_PositionDevs.csv`.
#'
#' AnalyzePositions writes a `<Name>_PositionDevs.csv` with one row per
#' (cell, time, embryo) tuple containing the d-magnitude position
#' deviation. This adapter optionally aggregates across embryos.
#'
#' @param posdevs_path path to `<Name>_PositionDevs.csv`.
#' @param aggregate_embryos `"mean"` (default), `"max"`, `"median"`, or `"none"`.
#' @return list(df = data.frame(cell, time, value, [embryo]),
#'              resolution = "timepoint")
position_dev_per_t_long <- function(posdevs_path,
                                     aggregate_embryos = c("mean", "max",
                                                           "median", "none")) {
    aggregate_embryos <- match.arg(aggregate_embryos)
    if (!file.exists(posdevs_path)) {
        stop("PositionDevs.csv not found: ", posdevs_path)
    }
    pd <- read.csv(posdevs_path, stringsAsFactors = FALSE,
                   check.names = FALSE, row.names = 1)
    # Decode rownames into (cell, time): rows are like "ABala:120"
    rn <- rownames(pd)
    parts <- strsplit(rn, ":", fixed = TRUE)
    cells <- vapply(parts, `[`, character(1), 1)
    times <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2)))

    # Drop any leading "Cell"/"Time"/"WT" non-embryo columns
    summary_cols <- c("Cell", "Time", "WT D mean", "WT D SD", "WT D count",
                      "WT_D_mean", "WT_D_SD", "WT_D_count")
    emb_cols <- setdiff(colnames(pd), summary_cols)
    pd_emb <- pd[, emb_cols, drop = FALSE]

    if (aggregate_embryos == "none") {
        emb_names <- colnames(pd_emb)
        out_list <- vector("list", length(emb_names))
        for (j in seq_along(emb_names)) {
            out_list[[j]] <- data.frame(
                cell = cells, time = times, embryo = emb_names[j],
                value = as.numeric(pd_emb[, j]),
                stringsAsFactors = FALSE)
        }
        df <- do.call(rbind, out_list)
        df <- df[!is.na(df$time), , drop = FALSE]
        return(list(df = df[order(df$cell, df$time), ],
                    resolution = "timepoint"))
    }

    fn <- .agg_fn(aggregate_embryos)
    agg_vals <- apply(as.matrix(pd_emb), 1, fn)
    df <- data.frame(cell = cells, time = times,
                     value = as.numeric(agg_vals), stringsAsFactors = FALSE)
    df <- df[!is.na(df$time), , drop = FALSE]
    list(df = df[order(df$cell, df$time), ], resolution = "timepoint")
}


#' Build a per-timepoint signed-component position-deviation frame.
#'
#' Decomposes the d-magnitude position deviation into directional
#' components by reading the on-disk rotated position CSVs:
#' `<mut_name>_rotatedX/Y/Z.csv` (mutant in WT-aligned frame, columns =
#' mutant embryos) and `<wt_name>_rotatedX/Y/Z.csv` (WT reference).
#'
#' Axis convention (`PlotPositionDevs.R:249-251`):
#'   x = AP (anterior/posterior), y = LR (left/right), z = DV (dorsal/ventral).
#'
#' Two `kind` options:
#' - `"AP_signed"` = mutant_X − mean(WT_X). Negative = anterior shift,
#'   positive = posterior shift.
#' - `"radial_signed"` = mutant_R − WT_R, where `R = sqrt(Y² + Z²)` is
#'   the cell's distance from the embryo's AP axis (the centerline).
#'   Negative = mutant cell is closer to the centerline than WT,
#'   positive = farther. Collapses LR (mirror-symmetric) with DV but
#'   preserves the in/out direction.
#'
#' @param mut_rotated_dir directory containing `<mut_name>_rotatedX/Y/Z.csv`.
#' @param mut_name basename used as the mutant rotated-file prefix.
#' @param wt_rotated_dir directory containing the WT rotated CSVs.
#' @param wt_name basename used as the WT-file prefix.
#' @param kind `"AP_signed"` or `"radial_signed"`.
#' @param aggregate_embryos `"mean"` (default), `"max"`, `"median"`, or `"none"`.
#' @return list(df = data.frame(cell, time, value, [embryo]),
#'              resolution = "timepoint")
position_dev_per_t_components_long <- function(mut_rotated_dir, mut_name,
                                                wt_rotated_dir, wt_name,
                                                kind = c("AP_signed",
                                                          "radial_signed"),
                                                aggregate_embryos = c("mean",
                                                                       "max",
                                                                       "median",
                                                                       "none")) {
    kind <- match.arg(kind)
    aggregate_embryos <- match.arg(aggregate_embryos)

    .read_rot <- function(dir, name, axis) {
        f <- file.path(dir, paste0(name, "_rotated", axis, ".csv"))
        if (!file.exists(f)) stop("rotated file not found: ", f)
        as.matrix(read.csv(f, stringsAsFactors = FALSE,
                            check.names = FALSE, row.names = 1))
    }

    mx <- .read_rot(mut_rotated_dir, mut_name, "X")
    if (kind == "radial_signed") {
        my <- .read_rot(mut_rotated_dir, mut_name, "Y")
        mz <- .read_rot(mut_rotated_dir, mut_name, "Z")
    }
    wx <- .read_rot(wt_rotated_dir,  wt_name,  "X")
    if (kind == "radial_signed") {
        wy <- .read_rot(wt_rotated_dir,  wt_name,  "Y")
        wz <- .read_rot(wt_rotated_dir,  wt_name,  "Z")
    }

    # WT means per (cell, time) across WT embryos, aligned to mutant rownames.
    rn <- rownames(mx)
    wt_x_mean <- rowMeans(wx, na.rm = TRUE)[rn]
    if (kind == "radial_signed") {
        wt_y_mean <- rowMeans(wy, na.rm = TRUE)[rn]
        wt_z_mean <- rowMeans(wz, na.rm = TRUE)[rn]
        wt_R <- sqrt(wt_y_mean^2 + wt_z_mean^2)
    }

    if (kind == "AP_signed") {
        dev_mat <- mx - wt_x_mean
    } else {
        mut_R <- sqrt(my^2 + mz^2)
        dev_mat <- mut_R - wt_R
    }

    parts <- strsplit(rn, ":", fixed = TRUE)
    cells <- vapply(parts, `[`, character(1), 1)
    times <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2)))

    if (aggregate_embryos == "none") {
        emb_names <- colnames(dev_mat)
        out_list <- vector("list", length(emb_names))
        for (j in seq_along(emb_names)) {
            out_list[[j]] <- data.frame(
                cell = cells, time = times, embryo = emb_names[j],
                value = as.numeric(dev_mat[, j]),
                stringsAsFactors = FALSE)
        }
        df <- do.call(rbind, out_list)
        df <- df[!is.na(df$time), , drop = FALSE]
        return(list(df = df[order(df$cell, df$time), ],
                    resolution = "timepoint"))
    }

    fn <- .agg_fn(aggregate_embryos)
    agg_vals <- apply(dev_mat, 1, fn)
    df <- data.frame(cell = cells, time = times,
                     value = as.numeric(agg_vals), stringsAsFactors = FALSE)
    df <- df[!is.na(df$time), , drop = FALSE]
    list(df = df[order(df$cell, df$time), ], resolution = "timepoint")
}


#' Build a long defect frame from `AnalyzePositions()`'s return list.
#'
#' AnalyzePositions returns `list(dDevs = ..., xDevs = ..., yDevs = ..., zDevs = ...)`
#' where each `Devs` is a vector indexed by `cell:time` keys (one per timepoint, see
#' `AnalyzePositions.R:524`). This adapter decodes the keys into `(cell, time)`
#' and returns a per-timepoint long frame — the natural high-resolution
#' source for position-deviation tree coloring.
#'
#' Note: `xDevs` etc. as returned by `AnalyzePositions` may be from the
#' last embryo processed (the function loops over embryos and the inner
#' return at line 458 is per-embryo). For multi-embryo aggregation use
#' `position_dev_long_from_csv()` against the `_CellMean/Max...` CSVs.
#'
#' @param devs the return list from `AnalyzePositions()`.
#' @param kind `"d"` (default, magnitude), `"x"`, `"y"`, or `"z"`.
#' @return list(df = data.frame(cell, time, value), resolution = "timepoint")
position_dev_long_from_devs <- function(devs, kind = c("d", "x", "y", "z")) {
    kind <- match.arg(kind)
    fld <- paste0(kind, "Devs")
    if (!(fld %in% names(devs))) {
        stop("devs list lacks field: ", fld)
    }
    v <- devs[[fld]]
    if (is.null(v) || length(v) == 0) {
        return(list(df = data.frame(cell = character(),
                                     time = numeric(),
                                     value = numeric()),
                    resolution = "timepoint"))
    }
    keys <- names(v)
    parts <- strsplit(keys, ":", fixed = TRUE)
    cells <- vapply(parts, `[`, character(1), 1)
    times <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2)))
    df <- data.frame(cell = cells, time = times,
                     value = as.numeric(v), stringsAsFactors = FALSE)
    df <- df[!is.na(df$time), , drop = FALSE]
    list(df = df[order(df$cell, df$time), ], resolution = "timepoint")
}


# -----------------------------------------------------------------------------
# Division-orientation deviation (per-cell)
# -----------------------------------------------------------------------------

#' Build a long defect frame from `<Name>_dots.csv`.
#'
#' `<Name>_dots.csv` is wide: one row per cell, columns include `WT_Mean,
#' WT_SD, Mutant_Mean, Mutant_SD, pVals, peak.names.WT_Mean..` plus one
#' column per mutant embryo with that embryo's dot value.
#'
#' Modes:
#'   "mutant_mean"  = Mutant_Mean (or per-embryo dot if aggregate_embryos="none")
#'   "abs_diff"     = |Mutant_Mean - WT_Mean|
#'   "signed_diff"  = Mutant_Mean - WT_Mean
#'   "angle_dev"    = 1 - Mutant_Mean (0 = same as WT, larger = more deviated)
#'
#' @return list(df = data.frame(cell, time, value, [embryo]),
#'              resolution = "cell")
dot_dev_long <- function(dots_path, positions_path,
                         mode = c("mutant_mean", "abs_diff", "signed_diff", "angle_dev"),
                         aggregate_embryos = c("mean", "max", "median", "none")) {
    mode <- match.arg(mode)
    aggregate_embryos <- match.arg(aggregate_embryos)
    if (!file.exists(dots_path)) {
        stop("dots file not found: ", dots_path)
    }
    dots <- read.csv(dots_path, stringsAsFactors = FALSE,
                     check.names = FALSE, row.names = 1)
    needed <- c("Mutant_Mean", "WT_Mean")
    missing_cols <- setdiff(needed, colnames(dots))
    if (length(missing_cols) > 0 && mode != "mutant_mean") {
        stop("dots file lacks columns: ", paste(missing_cols, collapse = ", "))
    }

    cd <- read_positions_to_cd(positions_path)

    if (aggregate_embryos == "none") {
        summary_cols <- c("WT_Mean", "WT_SD", "Mutant_Mean", "Mutant_SD",
                          "pVals", "peak.names.WT_Mean..", "Peak", "peak")
        emb_cols <- setdiff(colnames(dots), summary_cols)
        if (length(emb_cols) == 0) {
            stop("dots file has no per-embryo columns to aggregate from")
        }
        long <- tidyr::pivot_longer(
            data.frame(cell = rownames(dots), dots[, emb_cols, drop = FALSE],
                       check.names = FALSE),
            cols = -.data$cell, names_to = "embryo", values_to = "mut_dot")
        wt_mean <- setNames(dots$WT_Mean, rownames(dots))
        long$value <- switch(mode,
            mutant_mean  = long$mut_dot,
            abs_diff     = abs(long$mut_dot - wt_mean[long$cell]),
            signed_diff  = long$mut_dot - wt_mean[long$cell],
            angle_dev    = 1 - long$mut_dot
        )
        long <- long[, c("cell", "embryo", "value")]
        cd$.embryo_norm  <- gsub("[._-]", "", sub("^X", "", cd$embryo))
        long$.embryo_norm <- gsub("[._-]", "", sub("^X", "", long$embryo))
        merged <- merge(cd[, c("cell", "time", "embryo", ".embryo_norm")],
                        long[, c("cell", ".embryo_norm", "value")],
                        by = c("cell", ".embryo_norm"), all.x = TRUE)
        merged$.embryo_norm <- NULL
        return(list(df = merged[order(merged$cell, merged$time), ],
                    resolution = "cell"))
    }

    if (aggregate_embryos == "max") {
        summary_cols <- c("WT_Mean", "WT_SD", "Mutant_Mean", "Mutant_SD",
                          "pVals", "peak.names.WT_Mean..", "Peak", "peak")
        emb_cols <- setdiff(colnames(dots), summary_cols)
        emb_mat <- as.matrix(dots[, emb_cols, drop = FALSE])
        per_cell_mut <- apply(emb_mat, 1, function(x) {
            y <- x[!is.na(x)]; if (length(y) == 0) NA_real_ else max(y)
        })
    } else {
        if (aggregate_embryos == "median") {
            summary_cols <- c("WT_Mean", "WT_SD", "Mutant_Mean", "Mutant_SD",
                              "pVals", "peak.names.WT_Mean..", "Peak", "peak")
            emb_cols <- setdiff(colnames(dots), summary_cols)
            emb_mat <- as.matrix(dots[, emb_cols, drop = FALSE])
            per_cell_mut <- apply(emb_mat, 1, median, na.rm = TRUE)
        } else {
            per_cell_mut <- dots$Mutant_Mean
        }
    }

    vals <- switch(mode,
        mutant_mean  = per_cell_mut,
        abs_diff     = abs(per_cell_mut - dots$WT_Mean),
        signed_diff  = per_cell_mut - dots$WT_Mean,
        angle_dev    = 1 - per_cell_mut
    )
    per_cell <- data.frame(cell = rownames(dots), value = vals,
                           stringsAsFactors = FALSE)
    list(df = .broadcast_per_cell(per_cell, cd, keep_embryo = FALSE),
         resolution = "cell")
}


# -----------------------------------------------------------------------------
# Arbitrary user-supplied long frame (CSV/TSV passthrough)
# -----------------------------------------------------------------------------

#' Convert a user-supplied long-format CSV/TSV into the canonical adapter
#' shape so it can be fed to `plot_lineage_tree()` directly.
#'
#' If `time_col` is present and non-NULL, returns `resolution = "timepoint"`;
#' otherwise broadcasts via `positions_path` and returns `"cell"`.
arbitrary_long <- function(path_or_df, cell_col = "cell", time_col = "time",
                            value_col = "value", positions_path = NULL,
                            sep = "auto") {
    if (is.data.frame(path_or_df)) {
        df <- path_or_df
    } else {
        path <- path_or_df
        if (sep == "auto") {
            sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
        }
        df <- read.delim(path, sep = sep, stringsAsFactors = FALSE,
                         check.names = FALSE)
    }
    if (!(cell_col %in% colnames(df))) {
        stop("cell column not in file: ", cell_col)
    }
    if (!(value_col %in% colnames(df))) {
        stop("value column not in file: ", value_col)
    }
    out <- data.frame(cell = df[[cell_col]],
                       value = df[[value_col]],
                       stringsAsFactors = FALSE)
    if (!is.null(time_col) && time_col %in% colnames(df)) {
        out$time <- df[[time_col]]
        return(list(df = out[order(out$cell, out$time), c("cell", "time", "value")],
                    resolution = "timepoint"))
    }
    if (is.null(positions_path)) {
        stop("arbitrary_long needs positions_path to broadcast cell-level values")
    }
    cd <- read_positions_to_cd(positions_path)
    list(df = .broadcast_per_cell(out, cd, keep_embryo = FALSE),
         resolution = "cell")
}
