source("functions.R")
source("AnalyzeDivTimes.R")
source("AnalyzePositions.R")
source("AnalyzeRotation.R")
source("PlotDefectSummaries.R")
source("PlotPositionDevs.R")

args <- commandArgs(trailingOnly = TRUE)
print(args)


moviesList <- args[1]
print(moviesList)
print(length(args))
exp_file <- args[2]
expCutoff <- as.integer(args[3])
sig <- as.integer(args[4])
microns <- as.integer(args[5])
posDevTime <- as.integer(args[6])

# Handle peak expression loading
peak <- NULL
if (!is.na(exp_file) && file.exists(exp_file)) {
    # If we don't have Cells yet, ReadPeakExpression will fallback to info in exp_file
    peak <- ReadPeakExpression(exp_file)
}

# 1. Division Times and Kinetics
# Check if ccDevs already exists to allow resuming
# 1. Division Times and Kinetics
# Check if ccDevs already exists to allow resuming
# RETURNS: list(cc_dev=cc_dev, Cells=Cells) or just cc_dev depending on implementation.
# We need to capture the Cells list to correctly compute peak expression for all cells.
DivTimeResults <- AnalyzeDivTimes(moviesList, Expression = peak)
cc_dev <- DivTimeResults$cc_dev
Cells <- DivTimeResults$Cells

# Re-calculate peak expression with the full cell list if expression file is available
# This ensures cells in the current dataset that are missing from the expression file
# get imputed values from their ancestors.
if (!is.na(exp_file) && file.exists(exp_file) && !is.null(Cells)) {
    message("Recalculating peak expression for all ", length(Cells), " cells found in dataset")
    peak <- ReadPeakExpression(exp_file, cells = Cells)
}

# 2. Positions and Spatial Deviations
# Check for one of the final position output files
devs <- AnalyzePositions(moviesList, Expression = peak, CalculateNeighbors = TRUE, expCutoff = expCutoff)

# 3. Rotation Analysis and Defect Summaries
# These are relatively fast or need to re-run if plotting thresholds changed
AnalyzeRotation(moviesList, peak = peak)
PlotDefectSummaries(moviesList, peak, sig = sig, microns = microns, expCutoff = expCutoff, minDivTime = 70)

# 4. Final Arrow Plots and Expression correlation
MeanPosDevs <- PlotDeviationsList(Name = moviesList, exp = peak, t = posDevTime)
PlotExpVsDev(moviesList, outfile = paste(moviesList, "ExpVsDev.pdf", sep = "."), exp = peak)
