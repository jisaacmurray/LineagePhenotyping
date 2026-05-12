#!/usr/bin/env Rscript
#
# plot_paths.R — shared helper for routing display outputs into per-kind
# subdirectories under <output_dir>/plots/. Phase 5.2B.
#
# Display files (.pdf, .png, .jpg) are not consumed by any downstream R
# code in this pipeline (verified by grep), so they can be relocated
# aggressively without breaking analysis steps. Data files (.csv, .tsv,
# .txt, .RData) stay at top level.
#
# Layout modes:
#   "by_kind" (default): output_dir/plots/<kind>/[per_embryo/]
#   "flat"             : output_dir (legacy / back-compat for byte-diff
#                        verification against the Phase-2 baseline).
#
# Layout is read from `getOption("LineagePhenotyping.subdir_layout")`
# and defaults to "by_kind" when unset. `run_pipeline.R` sets this
# option from the YAML config at startup; ad-hoc callers can override
# with `options(LineagePhenotyping.subdir_layout = "flat")`.
#
# Known `kind` values (canonical):
#   "cc"          — cell-cycle defect plots
#   "position"    — fallback for misc position outputs
#   "position/cell"     — per-cell position aggregates
#   "position/per_t"    — per-timepoint position trees
#   "position/arrows"   — 3D arrow plots
#   "angle"       — division-orientation defect plots
#   "expression"  — ExpVsDev scatter plots
#   "boxplots"    — comparative boxplots
#   "summary"     — multi-panel summary grids (Spatial/Lineage/WT_stats)
#
# Use `per_embryo = TRUE` to nest under `per_embryo/`.

.plots_dir <- function(output_dir, kind = NULL, per_embryo = FALSE) {
    layout <- getOption("LineagePhenotyping.subdir_layout", "by_kind")
    if (!(layout %in% c("by_kind", "flat"))) {
        warning("Unknown LineagePhenotyping.subdir_layout: ", layout,
                " — falling back to 'by_kind'")
        layout <- "by_kind"
    }
    if (layout == "flat" || is.null(kind) || !nzchar(kind)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        return(output_dir)
    }
    # by_kind
    base <- file.path(output_dir, "plots", kind)
    if (per_embryo) base <- file.path(base, "per_embryo")
    dir.create(base, recursive = TRUE, showWarnings = FALSE)
    base
}
