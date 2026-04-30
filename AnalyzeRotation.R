source("functions.R")
library(rgl)
library(gdata)
library(plotrix)

# AnalyzeRotation
#
# Compares division orientations (sister-cell unit vectors) between a mutant
# dataset and a wild-type reference, computing per-cell t-test p-values.
#
# Parameters:
#   Name       (chr)  dataset identifier; mutant inputs and outputs are
#                     keyed by this name.
#   peak       (vec)  named vector of expression peaks; defaults to global
#                     ceh36_peak for back-compat (avoid relying on this).
#   data_dir   (chr)  root containing per-dataset input subdirectories.
#                     Default "." reproduces legacy CWD-relative behavior.
#   output_dir (chr)  where to write outputs; defaults to file.path(data_dir, Name).
#   wt_ref_dir (chr)  WT reference directory. The basename of this path is
#                     used as the file-name prefix for WT files
#                     (e.g. "Richard_et_al_plus_comma_WTDivTimeNorm.tsv").
#                     Defaults to file.path(data_dir, "Richard_et_al_plus_comma_WT").

AnalyzeRotation <- function(Name, peak = ceh36_peak,
                            data_dir = ".",
                            output_dir = NULL,
                            wt_ref_dir = NULL) {
    if (is.null(output_dir)) output_dir <- file.path(data_dir, Name)
    if (is.null(wt_ref_dir)) wt_ref_dir <- file.path(data_dir, "Richard_et_al_plus_comma_WT")
    wt_prefix <- basename(wt_ref_dir)

    message("Analyzing rotated division orientations", Name)

    # Load reference data
    WTDivTimes <- read.table(file.path(wt_ref_dir, paste0(wt_prefix, "DivTimeNorm.tsv")), header = T, row.names = 1, stringsAsFactors = F, check.names = FALSE)
    WTDivTimes[, 1] <- NULL

    # calculate mutant CC and div times to get cells union
    mutant_input_dir <- file.path(data_dir, Name)
    mutantRotatedX <- read.csv(file.path(output_dir, paste0(Name, "_rotatedX.csv")), row.names = 1, check.names = FALSE)

    Cells <- union(rownames(mutantRotatedX), rownames(WTDivTimes))

    WTRotatedX <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_rotatedX.csv")), row.names = 1, check.names = FALSE)
    WTRotatedY <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_rotatedY.csv")), row.names = 1, check.names = FALSE)
    WTRotatedZ <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_rotatedZ.csv")), row.names = 1, check.names = FALSE)

    # Ensure output directory exists
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    Times <- data.frame(do.call("rbind", strsplit(as.character(rownames(WTRotatedX)), ":", fixed = TRUE)))
    rownames(Times) <- rownames(WTRotatedX)
    minTimes <- tapply(as.numeric(as.character(Times[, 2])), Times[, 1], min)
    maxTimes <- tapply(as.numeric(as.character(Times[, 2])), Times[, 1], max)


    WTXMeans <- rowMeans(WTRotatedX, na.rm = T)
    WTYMeans <- rowMeans(WTRotatedY, na.rm = T)
    WTZMeans <- rowMeans(WTRotatedZ, na.rm = T)


    mutantRotatedX <- read.csv(file.path(output_dir, paste0(Name, "_rotatedX.csv")), row.names = 1, check.names = FALSE)
    mutantRotatedY <- read.csv(file.path(output_dir, paste0(Name, "_rotatedY.csv")), row.names = 1, check.names = FALSE)
    mutantRotatedZ <- read.csv(file.path(output_dir, paste0(Name, "_rotatedZ.csv")), row.names = 1, check.names = FALSE)

    mutantXMeans <- rowMeans(mutantRotatedX, na.rm = T)
    mutantYMeans <- rowMeans(mutantRotatedY, na.rm = T)
    mutantZMeans <- rowMeans(mutantRotatedZ, na.rm = T)

    # replace eventually with an sapply)
    #    Mutant_Dot=data.frame(col.names=colnames(mutantRotatedX))
    Mutant_Dot <- NULL
    #    WT_Dot=data.frame(col.names=colnames(WTRotatedX))
    WT_Dot <- NULL
    RowNames <- NULL
    for (thisCell in Cells) {
        thisSister <- GetSister(thisCell)
        thisParent <- GetParent(thisCell)

        Landmarks <- c("AB", "E", "MS", "C", "D")
        if (is.element(thisCell, Landmarks)) {
            print(paste(thisCell, thisSister, thisParent))
        }
        birthTime <- minTimes[thisCell]
        divTime <- maxTimes[thisParent]
        if (!is.na(thisSister) && !is.element(thisParent, RowNames[, 1]) && thisCell < thisSister && !is.na(WTXMeans[paste(thisCell, birthTime, sep = ":")]) && !is.na(mutantXMeans[paste(thisCell, birthTime, sep = ":")])) {
            # only use a, l and d daughters to avoid redundancy

            WT_Vector <- c(
                WTXMeans[paste(thisSister, birthTime, sep = ":")] - WTXMeans[paste(thisCell, birthTime, sep = ":")],
                WTYMeans[paste(thisSister, birthTime, sep = ":")] - WTYMeans[paste(thisCell, birthTime, sep = ":")],
                WTZMeans[paste(thisSister, birthTime, sep = ":")] - WTZMeans[paste(thisCell, birthTime, sep = ":")]
            )
            WTUnitVector <- WT_Vector / sqrt(sum(WT_Vector^2))

            #            rbind(Mutant_Dot,sapply(colnames(mutantRotatedX),function(X){
            #                    thisVector=c(mutantRotatedX[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedX[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedY[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedY[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedZ[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedZ[paste(thisCell,birthTime,sep=":"),X])
            #                    thisUnitVector=thisVector/sqrt(sum(thisVector^2))
            #                    as.double(WTUnitVector %*% thisUnitVector)
            #                    }))


            mags <- sqrt((mutantRotatedX[paste(thisSister, birthTime, sep = ":"), ] - mutantRotatedX[paste(thisCell, birthTime, sep = ":"), ])^2 +
                (mutantRotatedY[paste(thisSister, birthTime, sep = ":"), ] - mutantRotatedY[paste(thisCell, birthTime, sep = ":"), ])^2 +
                (mutantRotatedZ[paste(thisSister, birthTime, sep = ":"), ] - mutantRotatedZ[paste(thisCell, birthTime, sep = ":"), ])^2)


            Mutant_Dot <- rbind(
                Mutant_Dot,
                (mutantRotatedX[paste(thisSister, birthTime, sep = ":"), ] - mutantRotatedX[paste(thisCell, birthTime, sep = ":"), ]) / mags * WTUnitVector[1] +
                    (mutantRotatedY[paste(thisSister, birthTime, sep = ":"), ] - mutantRotatedY[paste(thisCell, birthTime, sep = ":"), ]) / mags * WTUnitVector[2] +
                    (mutantRotatedZ[paste(thisSister, birthTime, sep = ":"), ] - mutantRotatedZ[paste(thisCell, birthTime, sep = ":"), ]) / mags * WTUnitVector[3]
            )

            WTmags <- sqrt((WTRotatedX[paste(thisSister, birthTime, sep = ":"), ] - WTRotatedX[paste(thisCell, birthTime, sep = ":"), ])^2 +
                (WTRotatedY[paste(thisSister, birthTime, sep = ":"), ] - WTRotatedY[paste(thisCell, birthTime, sep = ":"), ])^2 +
                (WTRotatedZ[paste(thisSister, birthTime, sep = ":"), ] - WTRotatedZ[paste(thisCell, birthTime, sep = ":"), ])^2)

            WT_Dot <- rbind(
                WT_Dot,
                (WTRotatedX[paste(thisSister, birthTime, sep = ":"), ] - WTRotatedX[paste(thisCell, birthTime, sep = ":"), ]) / WTmags * WTUnitVector[1] +
                    (WTRotatedY[paste(thisSister, birthTime, sep = ":"), ] - WTRotatedY[paste(thisCell, birthTime, sep = ":"), ]) / WTmags * WTUnitVector[2] +
                    (WTRotatedZ[paste(thisSister, birthTime, sep = ":"), ] - WTRotatedZ[paste(thisCell, birthTime, sep = ":"), ]) / WTmags * WTUnitVector[3]
            )

            #            rbind(WT_Dot,sapply(colnames(WTRotatedX),function(X){
            #                    thisVector=c(mutantRotatedX[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedX[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedY[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedY[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedZ[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedZ[paste(thisCell,birthTime,sep=":"),X])
            #                    thisUnitVector=thisVector/sqrt(sum(thisVector^2))
            #                    as.double(WTUnitVector %*% thisUnitVector)
            #                    }))

            RowNames <- rbind(RowNames, thisParent)


            # Get XYZ positions of this cell and sister in each movie in WT
            # Get XYZ positions of this cell and sister in each mutant movie
            # Calculate mean vector in WT and angles (mean and SD)
            # Calculate mean vector in mutant and angles relative to WT (dot scores)
        }
    }

    rownames(WT_Dot) <- RowNames[, 1]
    rownames(Mutant_Dot) <- RowNames[, 1]

    WT_Mean <- rowMeans(WT_Dot, na.rm = T)
    WT_SD <- apply(WT_Dot, 1, sd, na.rm = T)

    Mutant_Mean <- rowMeans(Mutant_Dot, na.rm = T)
    Mutant_SD <- apply(Mutant_Dot, 1, sd, na.rm = T)
    pVals <- sapply(names(WT_Mean), function(X) {
        obj <- try(t.test(WT_Dot[X, ], Mutant_Dot[X, ]), silent = T)
        if (is(obj, "try-error")) {
            return(NA)
        } else {
            return(obj$p.value)
        }
    })


    write.csv(data.frame(WT_Mean, WT_SD, Mutant_Mean, Mutant_SD, pVals, Peak = peak[names(WT_Mean)], Mutant_Dot, check.names = FALSE), file = file.path(output_dir, paste0(Name, "_dots.csv")))
}
