#!/usr/bin/env Rscript
#
# tests/test_defect_trees.R
#
# Smoke test for the Phase-5 defect-tree plotting:
#   1. LIVEtools::plot_lineage_tree() builds a valid ggplot in both
#      "cell" and "timepoint" resolution modes.
#   2. NA values render in `na_color` (cells with no defect data).
#   3. The Sulston E/MS daughter-order convention is enforced (E above MS
#      in the rendered y-axis, even though the underlying tree is
#      alphabetical at every other node).
#   4. The DefectScoreFrames.R adapters round-trip a synthetic input.
#
# Run via:  Rscript tests/test_defect_trees.R
# Exits 0 on success, 1 on any failure.

suppressPackageStartupMessages({
    library(LIVEtools)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

# Locate runner-relative paths so this test can be run from any cwd.
.find_runner_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", args, value = TRUE)
    if (length(f) > 0) return(normalizePath(dirname(sub("^--file=", "", f[1]))))
    getwd()
}
TEST_DIR <- .find_runner_dir()
PKG_DIR  <- normalizePath(file.path(TEST_DIR, ".."))

source(file.path(PKG_DIR, "DefectScoreFrames.R"))
source(file.path(PKG_DIR, "PlotDefectTrees.R"))

# -----------------------------------------------------------------------------
# Synthetic 15-cell tree covering both AB and P1 sublineages — every cell has
# its sister so MakeNewick is happy.
# -----------------------------------------------------------------------------
mk <- function(cell, t_start, t_end, value) {
    data.frame(cell = cell, time = seq(t_start, t_end), value = value,
               stringsAsFactors = FALSE)
}
CD <- do.call(rbind, list(
    mk("P0",  1,  5, 10),
    mk("AB",  5, 12, 20),
    mk("P1",  5, 12,  5),
    mk("ABa", 12, 25, 30),
    mk("ABp", 12, 25, 25),
    mk("EMS", 12, 20, 15),
    mk("P2",  12, 25,  8),
    mk("MS",  20, 30, 12),
    mk("E",   20, 30, 18),
    mk("MSa", 30, 40, 14),
    mk("MSp", 30, 40, 14),
    mk("Ea",  30, 40, 11),
    mk("Ep",  30, 40, 11),
    mk("MSaa", 40, 50, 16),
    mk("MSap", 40, 50, NA),  # missing data -> na_color
    mk("MSpa", 40, 50, 16),
    mk("MSpp", 40, 50, 16),
    mk("Eal", 40, 50, 13),
    mk("Ear", 40, 50, 13),
    mk("Epl", 40, 50, 13),
    mk("Epr", 40, 50, 13)
))

failures <- character()
.check <- function(cond, msg) {
    if (!isTRUE(cond)) {
        cat("  FAIL:", msg, "\n")
        failures <<- c(failures, msg)
    } else {
        cat("  ok:  ", msg, "\n")
    }
}

# -----------------------------------------------------------------------------
# Test 1: timepoint-resolution returns a ggplot with the expected layers
# -----------------------------------------------------------------------------
cat("Test 1: resolution = 'timepoint'\n")
p_t <- plot_lineage_tree(CD, root = "P0", value_col = "value",
                         resolution = "timepoint",
                         percentile_bounds = c(0.1, 0.9))
.check(inherits(p_t, "ggplot"), "timepoint mode returns a ggplot")
geoms <- vapply(p_t$layers, function(L) class(L$geom)[1], character(1))
.check(sum(geoms == "GeomSegment") >= 2,
       "timepoint mode draws >= 2 GeomSegment layers (branches + divisions)")

# -----------------------------------------------------------------------------
# Test 2: cell-resolution
# -----------------------------------------------------------------------------
cat("\nTest 2: resolution = 'cell'\n")
p_c <- plot_lineage_tree(CD, root = "P0", value_col = "value",
                         resolution = "cell", cell_aggregate = "mean")
.check(inherits(p_c, "ggplot"), "cell mode returns a ggplot")
geoms_c <- vapply(p_c$layers, function(L) class(L$geom)[1], character(1))
.check(sum(geoms_c == "GeomSegment") >= 2,
       "cell mode draws >= 2 GeomSegment layers")

# -----------------------------------------------------------------------------
# Test 3: NA -> na_color renders without error
# -----------------------------------------------------------------------------
cat("\nTest 3: na_color path\n")
tmp_pdf <- tempfile(fileext = ".pdf")
ggplot2::ggsave(tmp_pdf, plot = p_t, width = 8, height = 6, units = "in")
.check(file.exists(tmp_pdf) && file.info(tmp_pdf)$size > 1000,
       sprintf("PDF written with NA cells (%.1fk bytes)",
               file.info(tmp_pdf)$size / 1024))
unlink(tmp_pdf)

# -----------------------------------------------------------------------------
# Test 4: E/MS daughter-order convention. With overrides applied, the y of
# any E descendant should exceed the y of any MS descendant.
# -----------------------------------------------------------------------------
cat("\nTest 4: Sulston E/MS daughter ordering\n")
build_data <- p_t$layers[[2]]$data  # CD_plot rows
e_cells   <- c("E", "Ea", "Ep", "Eal", "Ear", "Epl", "Epr")
ms_cells  <- c("MS", "MSa", "MSp", "MSaa", "MSap", "MSpa", "MSpp")
e_y  <- build_data$tree_y[build_data$cell %in% e_cells]
ms_y <- build_data$tree_y[build_data$cell %in% ms_cells]
.check(length(e_y) > 0 && length(ms_y) > 0,
       "E and MS cells both present in plot data")
.check(min(e_y) > max(ms_y),
       sprintf("All E cells y (%.1f-%.1f) above all MS cells y (%.1f-%.1f)",
               min(e_y), max(e_y), min(ms_y), max(ms_y)))

# -----------------------------------------------------------------------------
# Test 5: Within-lineage alphabetical order. MSaa < MSap < MSpa < MSpp by y.
# -----------------------------------------------------------------------------
cat("\nTest 5: within-MS alphabetical order\n")
ms_grand <- c("MSaa", "MSap", "MSpa", "MSpp")
y_ms <- vapply(ms_grand,
               function(c) build_data$tree_y[build_data$cell == c][1],
               numeric(1))
.check(all(diff(y_ms) > 0),
       sprintf("MS granddaughters in alphabetical order (bottom->top): %s",
               paste(sprintf("%s=%.1f", names(y_ms), y_ms), collapse = ", ")))

# -----------------------------------------------------------------------------
# Test 6: arbitrary_long passthrough
# -----------------------------------------------------------------------------
cat("\nTest 6: arbitrary_long passthrough\n")
adapter <- arbitrary_long(CD, cell_col = "cell", time_col = "time",
                           value_col = "value")
.check(adapter$resolution == "timepoint",
       "arbitrary_long detects timepoint resolution when time_col present")
.check(all(c("cell", "time", "value") %in% colnames(adapter$df)),
       "arbitrary_long produces (cell, time, value) frame")

# -----------------------------------------------------------------------------
# Test 7: paginate_tree_plots writes a multi-page PDF
# -----------------------------------------------------------------------------
cat("\nTest 7: paginate_tree_plots\n")
tmp_pdf <- tempfile(fileext = ".pdf")
paginate_tree_plots(list(p_t, p_c, p_t, p_c, p_t),
                    tmp_pdf, ncol = 2, plots_per_page = 4,
                    title = "smoke test", width = 12, height = 9)
.check(file.exists(tmp_pdf) && file.info(tmp_pdf)$size > 5000,
       sprintf("paginated PDF written (%.1fk bytes)",
               file.info(tmp_pdf)$size / 1024))
unlink(tmp_pdf)

# -----------------------------------------------------------------------------
# Test 8: truncate_to_last_data trims trailing-NA branches
# -----------------------------------------------------------------------------
cat("\nTest 8: truncate_to_last_data\n")
CD_trim <- rbind(
    CD,
    # Add 15 trailing-NA rows on ABa (lifespan was 12-25; now extended to 40 with NAs)
    data.frame(cell = "ABa", time = 26:40, value = NA)
)
p_trunc <- plot_lineage_tree(CD_trim, root = "P0", value_col = "value",
                              resolution = "timepoint",
                              truncate_to_last_data = TRUE)
seg_data <- p_trunc$layers[[2]]$data
aba_max_t <- max(seg_data$time[seg_data$cell == "ABa"], na.rm = TRUE)
.check(aba_max_t == 25,
       sprintf("truncate_to_last_data clipped ABa at last non-NA time (got %d, expected 25)",
               aba_max_t))

p_no_trunc <- plot_lineage_tree(CD_trim, root = "P0", value_col = "value",
                                 resolution = "timepoint",
                                 truncate_to_last_data = FALSE)
aba_max_t_full <- max(p_no_trunc$layers[[2]]$data$time[p_no_trunc$layers[[2]]$data$cell == "ABa"], na.rm = TRUE)
.check(aba_max_t_full == 40,
       sprintf("truncate_to_last_data = FALSE keeps trailing rows (got %d, expected 40)",
               aba_max_t_full))

# -----------------------------------------------------------------------------
# Test 9: drop_empty_cells prunes all-NA leaves
# -----------------------------------------------------------------------------
cat("\nTest 9: drop_empty_cells\n")
CD_empty <- CD
CD_empty$value[CD_empty$cell == "ABp"] <- NA   # ABp now fully NA
p_drop <- plot_lineage_tree(CD_empty, root = "P0", value_col = "value",
                             resolution = "timepoint",
                             drop_empty_cells = TRUE)
cells_kept <- unique(p_drop$layers[[2]]$data$cell)
.check(!("ABp" %in% cells_kept),
       sprintf("drop_empty_cells removed ABp (cells kept: %s)",
               paste(cells_kept, collapse = ",")))
.check("ABa" %in% cells_kept,
       "drop_empty_cells kept ABa (which has data)")

p_keep <- plot_lineage_tree(CD_empty, root = "P0", value_col = "value",
                             resolution = "timepoint",
                             drop_empty_cells = FALSE)
.check("ABp" %in% unique(p_keep$layers[[2]]$data$cell),
       "drop_empty_cells = FALSE retains all-NA ABp")

# -----------------------------------------------------------------------------
# Test 10: plot_lineage_tree accepts color_values + value_min/max for piecewise
# diverging gradient (CC-style).
# -----------------------------------------------------------------------------
cat("\nTest 10b: .plots_dir() routes by_kind layout\n")
source(file.path(PKG_DIR, "plot_paths.R"))
tmp_root <- tempfile("plots_dir_smoke_")
dir.create(tmp_root, recursive = TRUE)
options(LineagePhenotyping.subdir_layout = "by_kind")
.check(.plots_dir(tmp_root, "cc") == file.path(tmp_root, "plots", "cc"),
       "by_kind: cc -> <root>/plots/cc")
.check(.plots_dir(tmp_root, "position/per_t", per_embryo = TRUE) ==
        file.path(tmp_root, "plots", "position/per_t", "per_embryo"),
       "by_kind + per_embryo nests under per_embryo/")
options(LineagePhenotyping.subdir_layout = "flat")
.check(.plots_dir(tmp_root, "cc") == tmp_root,
       "flat layout collapses to <root>")
options(LineagePhenotyping.subdir_layout = "by_kind")  # restore
unlink(tmp_root, recursive = TRUE)

cat("\nTest 10: piecewise diverging color scale\n")
CD_signed <- CD
set.seed(42)
CD_signed$value <- runif(nrow(CD_signed), -25, 25)
p_div <- plot_lineage_tree(CD_signed, root = "P0", value_col = "value",
                            resolution = "timepoint",
                            value_min = -20, value_max = 20,
                            colors = c("green", "grey90", "grey90", "red"),
                            color_values = c(-20, -5, 5, 20))
.check(inherits(p_div, "ggplot"),
       "piecewise diverging gradient renders to a ggplot")
gradn_present <- any(vapply(p_div$scales$scales,
                             function(s) inherits(s, "ScaleContinuous"),
                             logical(1)))
.check(gradn_present, "ScaleContinuous (gradientn) attached")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n----\n")
if (length(failures) == 0) {
    cat("All tests passed.\n")
    quit(status = 0)
} else {
    cat(length(failures), "test(s) failed:\n")
    for (f in failures) cat("  - ", f, "\n", sep = "")
    quit(status = 1)
}
