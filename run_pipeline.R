#!/usr/bin/env Rscript

# run_pipeline.R
#
# Single config-driven entry point for the LineagePhenotyping analysis pipeline.
# Replaces the per-dataset wrappers (LineagePhenotyping_<dataset>.R) with one
# generic runner that reads its parameters from a YAML config file.
#
# Usage:
#   Rscript run_pipeline.R path/to/config.yaml
#
# The runner resolves paths in the following way:
#   - Absolute paths in the config are used as-is.
#   - Relative paths are resolved against `data_dir` (which is itself resolved
#     against the config file's directory if given as relative).
#   - This means configs are portable: copying a config + its data_dir to a
#     new machine and re-running just works.
#
# See configs/template.yaml for the full schema.

suppressPackageStartupMessages({
    library(yaml)
    library(LIVEtools)
})

# ------------------------------------------------------------------
# Locate the runner's own directory so we can source siblings
# ------------------------------------------------------------------
.find_runner_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(normalizePath(dirname(sub("^--file=", "", file_arg[1]))))
    }
    # Sourced interactively - fall back to CWD
    return(getwd())
}
RUNNER_DIR <- .find_runner_dir()

# ------------------------------------------------------------------
# Parse config
# ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript run_pipeline.R <config.yaml>")
}
config_path <- normalizePath(args[1])
if (!file.exists(config_path)) stop("Config file not found: ", config_path)

config_dir <- dirname(config_path)
cfg <- yaml.load_file(config_path)

# Required keys
required <- c("name", "data_dir")
missing <- setdiff(required, names(cfg))
if (length(missing) > 0) stop("Config missing required keys: ", paste(missing, collapse = ", "))

# Resolve a path: absolute -> as-is; relative -> against `base`.
.resolve <- function(path, base) {
    if (is.null(path) || is.na(path)) return(NULL)
    if (path == "") return(NULL)
    if (substr(path, 1, 1) == "/" || grepl("^[A-Za-z]:[\\\\/]", path)) {
        return(normalizePath(path, mustWork = FALSE))
    }
    return(normalizePath(file.path(base, path), mustWork = FALSE))
}

name        <- cfg$name
data_dir    <- .resolve(cfg$data_dir, config_dir)
output_dir  <- if (!is.null(cfg$output_dir)) .resolve(cfg$output_dir, config_dir) else file.path(data_dir, name)
wt_ref_dir  <- if (!is.null(cfg$wt_ref_dir)) .resolve(cfg$wt_ref_dir, data_dir) else file.path(data_dir, "Richard_et_al_plus_comma_WT")
exp_file    <- if (!is.null(cfg$expression_file)) .resolve(cfg$expression_file, data_dir) else NULL
emb_meta    <- if (!is.null(cfg$embryo_metadata_file)) .resolve(cfg$embryo_metadata_file, data_dir) else file.path(data_dir, "embryo_metadata.csv")
cell_names  <- if (!is.null(cfg$cell_names_file)) .resolve(cfg$cell_names_file, RUNNER_DIR) else file.path(RUNNER_DIR, "CellNames.csv")
cell_order  <- if (!is.null(cfg$cell_lineage_order_file)) .resolve(cfg$cell_lineage_order_file, RUNNER_DIR) else file.path(RUNNER_DIR, "Cells_350min_lineageOrder.csv")

# Numeric params (with defaults that match LineagePhenotyping_ceh76.R)
params      <- if (!is.null(cfg$params)) cfg$params else list()
expCutoff   <- if (!is.null(params$expCutoff))   as.numeric(params$expCutoff)   else 200
sig         <- if (!is.null(params$sig))         as.numeric(params$sig)         else 3
microns     <- if (!is.null(params$microns))     as.numeric(params$microns)     else 5
posDevTime  <- if (!is.null(params$posDevTime))  as.numeric(params$posDevTime)  else 250
minDivTime  <- if (!is.null(params$minDivTime))  as.numeric(params$minDivTime)  else 70
peak_recalc <- if (!is.null(params$peak_recalc))   isTRUE(params$peak_recalc)   else TRUE
dim_reduction <- if (!is.null(params$dim_reduction)) isTRUE(params$dim_reduction) else TRUE
defect_trees  <- if (!is.null(params$defect_trees))  isTRUE(params$defect_trees)  else TRUE
defect_tree_root      <- if (!is.null(params$defect_tree_root))      params$defect_tree_root      else "P0"
defect_tree_milestone <- params$defect_tree_milestone
defect_tree_split     <- if (!is.null(params$defect_tree_split_roots)) params$defect_tree_split_roots else c("ABa", "ABp", "P1")
defect_tree_arbitrary <- params$defect_tree_arbitrary

# Display-output layout (Phase 5.2). "by_kind" routes all PDFs/PNGs/JPGs
# under output_dir/plots/<kind>/. "flat" reproduces pre-5.2 byte-identical
# behavior. The .plots_dir() helper (sourced from plot_paths.R) reads
# this option.
subdir_layout <- if (!is.null(params$subdir_layout)) params$subdir_layout else "by_kind"
options(LineagePhenotyping.subdir_layout = subdir_layout)

# Phase 5.2D: arrow-frame jpg directories are temp scratch used to
# assemble the corresponding .pdf files. Default delete after the PDF
# is sealed.
keep_arrow_jpgs <- if (!is.null(params$keep_arrow_jpgs)) isTRUE(params$keep_arrow_jpgs) else FALSE
options(LineagePhenotyping.keep_arrow_jpgs = keep_arrow_jpgs)

# ------------------------------------------------------------------
# Set CWD to data_dir so any residual relative paths inside helper
# scripts resolve sensibly. The functions themselves now use
# absolute paths via the data_dir/output_dir/wt_ref_dir args, but
# this also lets the legacy DimensionalityReductionHelpers code find
# embryo_metadata.csv if a config didn't override it.
# ------------------------------------------------------------------
prev_wd <- getwd()
on.exit(setwd(prev_wd), add = TRUE)
setwd(data_dir)

# ------------------------------------------------------------------
# Source the analysis library
# ------------------------------------------------------------------
.source_runner <- function(name) source(file.path(RUNNER_DIR, name))
.source_runner("plot_paths.R")    # .plots_dir() helper, must come first
.source_runner("functions.R")
.source_runner("DimensionalityReductionHelpers.R")
.source_runner("AnalyzeDivTimes.R")
.source_runner("AnalyzePositions.R")
.source_runner("AnalyzeRotation.R")
.source_runner("PlotDefectSummaries.R")
.source_runner("PlotPositionDevs.R")
.source_runner("PlotComparisonBoxplots.R")
.source_runner("DefectScoreFrames.R")
.source_runner("PlotDefectTrees.R")

# ------------------------------------------------------------------
# Set up logging
# ------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(output_dir, "run_log.txt")

cat(sprintf("[%s] run_pipeline.R\n", format(Sys.time())), file = log_path)
cat(sprintf("  config:        %s\n", config_path), file = log_path, append = TRUE)
cat(sprintf("  name:          %s\n", name),       file = log_path, append = TRUE)
cat(sprintf("  data_dir:      %s\n", data_dir),   file = log_path, append = TRUE)
cat(sprintf("  output_dir:    %s\n", output_dir), file = log_path, append = TRUE)
cat(sprintf("  wt_ref_dir:    %s\n", wt_ref_dir), file = log_path, append = TRUE)
cat(sprintf("  expression:    %s\n", exp_file %||% "(none)"), file = log_path, append = TRUE)
cat(sprintf("  params: expCutoff=%g sig=%g microns=%g posDevTime=%g minDivTime=%g peak_recalc=%s dim_reduction=%s subdir_layout=%s keep_arrow_jpgs=%s\n",
            expCutoff, sig, microns, posDevTime, minDivTime, peak_recalc, dim_reduction,
            subdir_layout, keep_arrow_jpgs),
    file = log_path, append = TRUE)

message(sprintf("[run_pipeline] Analyzing %s", name))
message(sprintf("[run_pipeline] data_dir   = %s", data_dir))
message(sprintf("[run_pipeline] output_dir = %s", output_dir))

# ------------------------------------------------------------------
# Load peak expression
# ------------------------------------------------------------------
peak <- NULL
if (!is.null(exp_file) && file.exists(exp_file)) {
    peak <- ReadPeakExpression(exp_file)
} else if (file.exists(file.path(RUNNER_DIR, "ceh36_peak.csv"))) {
    message("[run_pipeline] No expression_file in config; falling back to ceh36_peak.csv")
    peak <- ReadPeakExpression(file.path(RUNNER_DIR, "ceh36_peak.csv"))
}

# ------------------------------------------------------------------
# Pipeline
# ------------------------------------------------------------------

# 1. Division Times
DivTimeResults <- AnalyzeDivTimes(
    name,
    Expression = peak,
    data_dir = data_dir,
    output_dir = output_dir,
    wt_ref_dir = wt_ref_dir,
    embryo_metadata_file = emb_meta,
    dim_reduction = dim_reduction
)
Cells <- DivTimeResults$Cells

# 1b. Optional peak recalculation using the discovered cell list
if (peak_recalc && !is.null(exp_file) && file.exists(exp_file) && !is.null(Cells)) {
    message(sprintf("[run_pipeline] Recalculating peak expression for all %d cells found in dataset", length(Cells)))
    peak <- ReadPeakExpression(exp_file, cells = Cells)
}

# 2. Positions
devs <- AnalyzePositions(
    name,
    Expression = peak,
    CalculateNeighbors = TRUE,
    expCutoff = expCutoff,
    data_dir = data_dir,
    output_dir = output_dir,
    wt_ref_dir = wt_ref_dir,
    embryo_metadata_file = emb_meta,
    dim_reduction = dim_reduction
)

# 3. Rotation
AnalyzeRotation(
    name,
    peak = peak,
    data_dir = data_dir,
    output_dir = output_dir,
    wt_ref_dir = wt_ref_dir
)

# 4. Defect summaries
PlotDefectSummaries(
    name, peak,
    sig = sig, microns = microns, expCutoff = expCutoff, minDivTime = minDivTime,
    data_dir = data_dir,
    output_dir = output_dir,
    wt_ref_dir = wt_ref_dir,
    cell_names_file = cell_names,
    cell_lineage_order_file = cell_order
)

# 5. Position deviation arrows
MeanPosDevs <- PlotDeviationsList(
    Name = name,
    exp = peak,
    t = posDevTime,
    data_dir = data_dir,
    output_dir = output_dir,
    wt_ref_dir = wt_ref_dir
)

# 6. Expression vs deviation
PlotExpVsDev(
    name,
    # Underscore separator (was "." pre-5.2) so the file looks like the
    # rest of the pipeline outputs and is easy to grep for.
    outfile = paste0(name, "_ExpVsDev.pdf"),
    exp = peak,
    data_dir = data_dir,
    output_dir = output_dir
)

# 7. Comparative boxplots
PlotComparisonBoxplots(
    name = name,
    data_dir = data_dir,
    output_dir = output_dir,
    wt_ref_dir = wt_ref_dir,
    embryo_metadata_file = emb_meta
)

# 8. Defect-colored lineage trees (Phase 5)
if (defect_trees) {
    message("[run_pipeline] Rendering defect-colored lineage trees...")
    tryCatch(
        PlotDefectTrees(
            name = name,
            devs = devs,
            ccdevs_path      = file.path(output_dir, paste0(name, "_ccDevs.csv")),
            dots_path        = file.path(output_dir, paste0(name, "_dots.csv")),
            mean_pos_dev_csv = file.path(output_dir, paste0(name, "_CellMeanPositionDevs.csv")),
            max_pos_dev_csv  = file.path(output_dir, paste0(name, "_CellMaxPositionDevs.csv")),
            positions_path   = file.path(data_dir, name, paste0(name, "positions.txt")),
            data_dir = data_dir, output_dir = output_dir, wt_ref_dir = wt_ref_dir,
            root = defect_tree_root,
            end_time_milestone = defect_tree_milestone,
            split_roots = defect_tree_split,
            arbitrary_scores = defect_tree_arbitrary
        ),
        error = function(e) {
            message("[run_pipeline] PlotDefectTrees failed: ", conditionMessage(e))
        }
    )
}

cat(sprintf("[%s] Done.\n", format(Sys.time())), file = log_path, append = TRUE)
message("[run_pipeline] Done.")
