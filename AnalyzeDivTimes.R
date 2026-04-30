source("functions.R")
source("DimensionalityReductionHelpers.R")
library(rgl)
library(gdata)
library(plotrix)

AnalyzeDivTimes <- function(Name, Expression = NULL,
                            data_dir = ".",
                            output_dir = NULL,
                            wt_ref_dir = NULL,
                            embryo_metadata_file = NULL,
                            dim_reduction = TRUE) {
    if (is.null(output_dir)) output_dir <- file.path(data_dir, Name)
    if (is.null(wt_ref_dir)) wt_ref_dir <- file.path(data_dir, "Richard_et_al_plus_comma_WT")
    wt_prefix <- basename(wt_ref_dir)
    mutant_input_dir <- file.path(data_dir, Name)
    if (is.null(embryo_metadata_file)) embryo_metadata_file <- "embryo_metadata.csv"

    message(paste("Analyzing division times", Name))

    # Load WT Data
    WTDivTimes <- read.table(file.path(wt_ref_dir, paste0(wt_prefix, "DivTimeNorm.tsv")), header = T, row.names = 1, stringsAsFactors = F, check.names = FALSE)
    WTDivTimes[, 1] <- NULL
    WTCellCycles <- read.table(file.path(wt_ref_dir, paste0(wt_prefix, "CCLengthNorm.tsv")), header = T, row.names = 1, stringsAsFactors = F, check.names = FALSE)
    WTCellCycles[, 1] <- NULL

    # calculate mutant CC and div times
    MutantDivTimes <- read.table(file.path(mutant_input_dir, paste0(Name, "DivTimeNorm.tsv")), header = T, row.names = 1, stringsAsFactors = F, check.names = FALSE)
    MutantDivTimes[, 1] <- NULL
    MutantCellCycles <- read.table(file.path(mutant_input_dir, paste0(Name, "CCLengthNorm.tsv")), header = T, row.names = 1, stringsAsFactors = F, check.names = FALSE)
    MutantCellCycles[, 1] <- NULL
    MutantTerminalLengths <- read.table(file.path(mutant_input_dir, paste0(Name, "CCLengthMinTerminal.tsv")), header = T, row.names = 1, stringsAsFactors = F, check.names = FALSE)

    # Define unified Cells list
    Cells <- union(rownames(MutantDivTimes), rownames(WTDivTimes))

    # calculate WT statistics for these cells
    WTDivMeans <- rowMeans(WTDivTimes[Cells, ], na.rm = T)
    WTDivSDs <- apply(WTDivTimes[Cells, ], 1, sd, na.rm = T)
    WTCCMeans <- rowMeans(WTCellCycles[Cells, ], na.rm = T)
    WTCCSDs <- apply(WTCellCycles[Cells, ], 1, sd, na.rm = T)

    WTCCCounts <- rowSums(!is.na(WTCellCycles[Cells, ]))
    WTDivCounts <- rowSums(!is.na(WTDivTimes[Cells, ]))

    # Require at least 2 observations in WT
    WTDivMeans[WTDivCounts < 2] <- NA
    WTDivSDs[WTDivCounts < 2] <- NA
    WTCCMeans[WTCCCounts < 2] <- NA
    WTCCSDs[WTCCCounts < 2] <- NA

    # Ensure directory exists
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # calculate mutant statistics
    MutantDivMeans <- rowMeans(MutantDivTimes[Cells, ], na.rm = T)
    MutantDivSDs <- apply(MutantDivTimes[Cells, ], 1, sd, na.rm = T)
    MutantCCMeans <- rowMeans(MutantCellCycles[Cells, ], na.rm = T)
    MutantCCSDs <- apply(MutantCellCycles[Cells, ], 1, sd, na.rm = T)

    # Fill Expression if NULL
    if (is.null(Expression)) {
        Expression <- rep(NA, length(Cells))
        names(Expression) <- Cells
    }

    # Ensure Expression vector matches Cells
    ExpVals <- Expression[Cells]

    library(ggplot2)

    # Define Lineage Groups and Colors
    CellData <- data.frame(Cell = Cells, stringsAsFactors = F)
    CellData$Lineage <- "Other"
    CellData$Lineage[grep("^ABa", Cells)] <- "ABa"
    CellData$Lineage[grep("^ABp", Cells)] <- "ABp"
    CellData$Lineage[grep("^MS", Cells)] <- "MS"
    CellData$Lineage[grep("^E", Cells)] <- "E"
    CellData$Lineage[grep("^C", Cells)] <- "C"
    CellData$Lineage[grep("^D", Cells)] <- "D"
    CellData$Lineage[grep("^P|^Z", Cells)] <- "P/Z"

    LineageColors <- c("ABa" = "red", "ABp" = "blue", "MS" = "green3", "E" = "magenta", "C" = "cyan3", "D" = "gold2", "P/Z" = "grey50", "Other" = "black")

    MutantEmbryos <- colnames(MutantCellCycles)
    pdf(file.path(output_dir, paste0(Name, "CC_plots.pdf")))

    # Plot helper
    PlotComparison <- function(mutantData, wtData, title, x_lab, y_lab) {
        df <- data.frame(Mutant = mutantData, WT = wtData, Lineage = CellData$Lineage)
        df <- df[!is.na(df$Mutant) & !is.na(df$WT), ]
        if (nrow(df) == 0) {
            return(NULL)
        }

        p <- ggplot(df, aes(x = Mutant, y = WT, color = Lineage)) +
            geom_point(alpha = 0.6, size = 1) +
            scale_color_manual(values = LineageColors) +
            theme_bw() +
            labs(title = title, x = x_lab, y = y_lab) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey")
        print(p)
    }

    # THIS PLOTS DIV TIME FOR EACH MUTANT SERIES VS WT AVERAGE
    for (i in MutantEmbryos) {
        PlotComparison(MutantDivTimes[Cells, i], WTDivMeans[Cells],
            title = paste("Division Time -", i),
            x_lab = paste("Mutant Division Time (min) -", i),
            y_lab = "Wild Type Mean Division Time (min)"
        )
    }

    # THIS PLOTS CC LENGTH FOR EACH MUTANT SERIES VS WT AVERAGE
    for (i in MutantEmbryos) {
        PlotComparison(MutantCellCycles[Cells, i], WTCCMeans[Cells],
            title = paste("Cell Cycle Length -", i),
            x_lab = paste("Mutant Cell Cycle Length (min) -", i),
            y_lab = "Wild Type Mean Cell Cycle Length (min)"
        )
    }

    # P-Value Plot
    CC_Dif_pValues <- sapply(Cells, function(X) {
        obj <- try(t.test(WTCellCycles[X, ], MutantCellCycles[X, ]), silent = T)
        if (is(obj, "try-error")) {
            return(NA)
        } else {
            return(obj$p.value)
        }
    })

    df_p <- data.frame(WTDiv = WTDivMeans[Cells], pVal = log10(CC_Dif_pValues), Lineage = CellData$Lineage)
    df_p <- df_p[!is.na(df_p$WTDiv) & !is.na(df_p$pVal), ]
    if (nrow(df_p) > 0) {
        p <- ggplot(df_p, aes(x = WTDiv, y = pVal, color = Lineage)) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = LineageColors) +
            theme_bw() +
            labs(title = "Cell Cycle Length P-Values", x = "Wild Type Mean Division Time (min)", y = "log10(p-value)") +
            geom_hline(yintercept = log10(0.05), linetype = "dotted", color = "red")
        print(p)
    }

    AlphOrder <- 1:length(Cells)
    names(AlphOrder) <- Cells

    df_alp <- data.frame(Order = AlphOrder[Cells], pVal = log10(CC_Dif_pValues), Lineage = CellData$Lineage)
    df_alp <- df_alp[!is.na(df_alp$pVal), ]
    if (nrow(df_alp) > 0) {
        p <- ggplot(df_alp, aes(x = Order, y = pVal, color = Lineage)) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = LineageColors) +
            theme_bw() +
            labs(title = "P-Values (Linear Order)", x = "Cells (Alphabetical Order Index)", y = "log10(p-value)")
        print(p)
    }


    # Calculate WT Birth Times
    WTBirthMeans <- WTDivMeans[Cells] - WTCCMeans[Cells]

    # Generate and output table of embryo/cell deviations
    CCdevs <- data.frame()
    for (i in MutantEmbryos) {
        CCdevs <- rbind(CCdevs, data.frame(Cells, i, ExpVals, WTBirthMeans, WTDivMeans[Cells],
            MutantDivTimes[Cells, i], WTDivMeans[Cells], WTDivSDs[Cells], MutantCellCycles[Cells, i], WTCCMeans[Cells], WTCCSDs[Cells],
            MutantCellCycles[Cells, i] - WTCCMeans[Cells], (MutantCellCycles[Cells, i] - WTCCMeans[Cells]) / WTCCSDs[Cells],
            MutantTerminalLengths[Cells, i], MutantTerminalLengths[Cells, i] - WTCCMeans[Cells], (MutantTerminalLengths[Cells, i] - WTCCMeans[Cells]) / WTCCSDs[Cells],
            row.names = paste(Cells, i, sep = "_")
        ))
    }
    # Note: WT Div Time included twice (Col 5 and Col 7) - cleaning up
    # Original Col 7 (WT Div time) kept for compat. New Col 5 (WTDivMeans) for request.
    # Actually, let's keep it clean.

    colnames(CCdevs) <- c(
        "Cell", "Embryo", "Expression", "WT_Birth_Time", "WT_Div_time",
        "Div_time", "WT_Div_time_dup", "WT_Div_time_SD", "CC", "WT_CC", "WT_CC_SD", "CC_Deviation", "CC_Z_score",
        "Terminal_Cell_Length", "Terminal_cell_delta", "Terminal_Cell_Z"
    )
    # Removing Dup col is cleaner but risky for indices.
    # Current indices used in PlotDefectSummaries:
    # 9: CC Deviaton -> Now 12
    # 10: CC Z Score -> Now 13
    # 12: Term Cell Delta -> Now 15
    # 13: Term Cell Z -> Now 16
    # 3: Div Time -> Now 6
    # 4: WT Div Time -> Now 7 (or 5)

    # Wait, simple insertion SHIFTS indices.
    # CCdevs users: PlotDefectSummaries use NAMES mostly:
    # theseCCDevs[,"CC.Deviation"] work by name.
    # theseCCDevs[,1] works by index (assumes Cell).
    # theseCCDevs[,9] used in AnalyzeDivTimes itself later! LINE 160.
    # I MUST update line 160 etc indices.

    # Correcting indices for internal usage below:
    # Cell=1, Embryo=2, Exp=3, Birth=4, Div=5, MDiv=6, WTDiv=7, ...
    # CC Deviation is now 12.
    # CC Z Score is now 13.
    # Terminal Delta is now 15.
    # Terminal Z is now 16.
    # Div Time is 6.
    # WT Div Time is 7.

    # Output full list, and filtered list containing just potentially interesting deviant cells
    write.csv(CCdevs, file = file.path(output_dir, paste0(Name, "_ccDevs.csv")))
    write.csv(
        subset(
            CCdevs,
            # Filter for CC deviation > 5 and z score > 3
            ((abs(CCdevs[, 12]) > 5 & abs(CCdevs[, 13]) > 3) |
                # Filter for cells that should but don't divide and have observed lengths > 5min more than WT CC and z score > 3
                (CCdevs[, 15] > 5 & CCdevs[, 16] > 3)) |
                # Filter for cells that divide inappropriately
                (CCdevs[, 6] > 0 & is.nan(CCdevs[, 7]))
        ),
        file = file.path(output_dir, paste0(Name, "5min_3sd_Devs.csv"))
    )

    # Plot CC deviation versus Z score
    df_dev <- CCdevs
    df_dev$Lineage <- CellData$Lineage[match(df_dev$Cell, CellData$Cell)]

    p <- ggplot(df_dev, aes(x = CC_Deviation, y = CC_Z_score, color = Lineage)) +
        geom_point(alpha = 0.4, size = 0.5) +
        scale_color_manual(values = LineageColors) +
        theme_bw() +
        labs(title = "Cell Cycle Deviation vs Z-Score", x = "CC Deviation (min)", y = "CC Z-Score")
    print(p)

    # Plot CC Z score vs Div Time
    # Ensure columns exist before calculation
    cols_to_avg <- c("CC_Z_score", "Terminal_Cell_Z")
    df_dev$CombinedZ <- rowMeans(df_dev[, cols_to_avg, drop = FALSE], na.rm = TRUE)

    p <- ggplot(df_dev, aes(x = Div_time, y = CombinedZ, color = Lineage)) +
        geom_point(alpha = 0.4, size = 0.5) +
        scale_color_manual(values = LineageColors) +
        theme_bw() +
        labs(title = "Combined CC Z-Score vs Division Time", x = "Division Time (min)", y = "Mean CC/Terminal Z-Score")
    print(p)

    # Plot CC Z score vs Alpha Order
    df_dev$Order <- AlphOrder[df_dev$Cell]

    p <- ggplot(df_dev, aes(x = Order, y = CombinedZ, color = Lineage)) +
        geom_point(alpha = 0.4, size = 0.5) +
        scale_color_manual(values = LineageColors) +
        theme_bw() +
        labs(title = "Combined CC Z-Score vs Alpha Order", x = "Cells (Alphabetical Order Index)", y = "Mean CC/Terminal Z-Score")
    print(p)

    dev.off()

    # OUTPUTS CLUSTERABLE DEVIATION AND CC LENGTH SPREADSHEETS
    # Modified to include Expression and separate Cell/P-val
    write.table(data.frame(Cells, round(CC_Dif_pValues, 4), ExpVals, round(WTBirthMeans, 1), round(WTDivMeans[Cells], 1), round(data.frame(WTCellCycles[Cells, ], MutantCellCycles[Cells, ]))), file = file.path(output_dir, paste0(Name, "CC.txt")), sep = "\t", na = "", quote = F, row.names = F)
    write.table(data.frame(Cells, round(CC_Dif_pValues, 4), ExpVals, round(WTBirthMeans, 1), round(WTDivMeans[Cells], 1), round(data.frame(WTCellCycles[Cells, ] - WTCCMeans[Cells], MutantCellCycles[Cells, ] - WTCCMeans[Cells]))), file = file.path(output_dir, paste0(Name, "CCdev.txt")), sep = "\t", na = "", quote = F, row.names = F)

    # Count and output table of outliers per cell
    # Count and output table of outliers per cell
    LateCalls <- sapply(Cells, function(X) {
        test <- CCdevs[CCdevs[, 1] == X & CCdevs[, 12] > 5 & CCdevs[, 13] > 3, ]
        sum(!is.na(test[6]))
    })

    EarlyCalls <- sapply(Cells, function(X) {
        test <- CCdevs[CCdevs[, 1] == X & CCdevs[, 12] < (-5) & CCdevs[, 13] < (-3), ]
        sum(!is.na(test[6]))
    })

    MissedCalls <- sapply(Cells, function(X) {
        test <- CCdevs[CCdevs[, 1] == X & CCdevs[, 15] > 5 & CCdevs[, 16] > 3, ]
        sum(!is.na(test[7]))
    }) # Check 7 (WT Div) or 4 (WT CC)? Originally 4 -> WT Div Time SD. Wait.

    EctopicCalls <- sapply(Cells, function(X) {
        test <- CCdevs[CCdevs[, 1] == X & is.na(CCdevs[, 7]) & !is.na(CCdevs[, 6]), ]
        sum(!is.na(test[6]))
    })

    # Lost / TrajectoryEnds Logic:
    # Defined as not dividing (CC is NA) AND track ends significantly earlier than Sister ( > 5min difference).
    LostCalls <- sapply(Cells, function(X) {
        # Subset for this cell
        myRows <- CCdevs[CCdevs[, 1] == X, , drop = FALSE]
        if (nrow(myRows) == 0) {
            return(0)
        }

        Sister <- GetSister(X)
        count <- 0

        for (i in 1:nrow(myRows)) {
            emb <- myRows[i, 2] # Embryo
            # Check if Terminal (CC is NA)
            if (is.na(myRows[i, 9])) { # Col 9 is "CC" (Mutant Cell Cycle)
                # Get my end time
                # My Start = DivTime - CC (if CC exist) OR DivTime - TerminalLen (if Term exist)?
                # Wait. Div Time in table is generally "Time cell divided".
                # If it didn't divide, MutantDivTimes is NA.
                # But MutantTerminalLengths (Col 14) implies length observed.
                # Birth Time = MyDivTime - MyCC.
                # But if I didn't divide, how do I know my BirthTime?
                # Ans: Parent's Div Time.
                # Easier: Compare EndTime directly.
                # MyEndTime = BirthTime + TerminalLength.
                # SisterEndTime = SisterBirthTime + SisterLength.
                # Since BirthTime == SisterBirthTime (same parent div),
                # Lost condition simplifies to: MyTerminalLength < SisterLength - 5.

                myTermLen <- myRows[i, 14] # Col 14: Terminal Cell Length

                # Get Sister Data for same embryo
                # We need to look up Sister in CCdevs for this embryo
                # Sister might divide (use CC, Col 9) or not (use TerminalLen, Col 14).
                # Also Sister might be dividing, so use MAX(CC, TerminalLen)?
                # Actually if sister divides, her "Length" is CC.

                # Fetch sister row from full CCdevs is slow inside loop.
                # Better: Pre-fetch or lookup.
                # Given limited N, lookup is okay.

                sisRow <- CCdevs[CCdevs[, 1] == Sister & CCdevs[, 2] == emb, ]
                if (nrow(sisRow) > 0) {
                    sisLen <- NA
                    if (!is.na(sisRow[1, 9])) { # Sister Divided
                        sisLen <- sisRow[1, 9]
                    } else if (!is.na(sisRow[1, 14])) { # Sister Terminal
                        sisLen <- sisRow[1, 14]
                    }

                    if (!is.na(myTermLen) && !is.na(sisLen)) {
                        if (myTermLen < (sisLen - 5)) {
                            count <- count + 1
                        }
                    }
                }
            }
        }
        return(count)
    })

    TotalDefects <- LateCalls + EarlyCalls + MissedCalls + EctopicCalls
    DeathCandidates <- MissedCalls + LostCalls

    write.csv(data.frame(Cells, ExpVals, WTDivMeans[Cells], LateCalls, EarlyCalls, MissedCalls, EctopicCalls, LostCalls, DeathCandidates, TotalDefects, CC_Dif_pValues, WTCCMeans[Cells], MutantCCMeans[Cells]), file = file.path(output_dir, paste0(Name, "CellDefectSummary.csv")), quote = F, row.names = F)

    ## quantitative deviations
    Deviations <- MutantCellCycles[Cells, ] - WTCCMeans[Cells]

    ## Ectopic divisions - enter CC and set to zero if not observed (conservative on inclusion but may be too high level
    # set divisions not observed in WT as zero
    #    Deviations[is.na(WTCCMeans[Cells]),]=0
    # set subset of those observed in mutants to mutant CC (consider instead making this a small # based on scale of CC devs ~10 min or so?)
    #    Deviations[!is.na(MutantCellCycles[Cells,]) & is.na(WTCCMeans[Cells])] = MutantCellCycles[!is.na(MutantCellCycles[Cells,])& is.na(WTCCMeans[Cells])]
    Deviations[!is.na(MutantCellCycles[Cells, ]) & is.na(WTCCMeans[Cells])] <- -50

    ## non-observed divisions (conservative evidence for late)
    NonObservedMinDevs <- MutantTerminalLengths[Cells, 2:length(MutantTerminalLengths)] - WTCCMeans[Cells]
    Deviations[NonObservedMinDevs > 0 & !is.na(NonObservedMinDevs)] <- NonObservedMinDevs[NonObservedMinDevs > 0 & !is.na(NonObservedMinDevs)]


    # UMAP / PCA Visualizations (skip if dim_reduction is FALSE)
    if (dim_reduction) {
        metadata <- NULL
        if (file.exists(embryo_metadata_file)) metadata <- read.csv(embryo_metadata_file, stringsAsFactors = FALSE)

        # Define metrics to analyze
        vis_list <- list(
            list(
                title = "Cell Cycle Lengths",
                data = load_and_merge_data(
                    file.path(mutant_input_dir, paste0(Name, "CCLengthNorm.tsv")),
                    file.path(wt_ref_dir, paste0(wt_prefix, "CCLengthNorm.tsv"))
                )
            ),
            list(
                title = "Division Times",
                data = load_and_merge_data(
                    file.path(mutant_input_dir, paste0(Name, "DivTimeNorm.tsv")),
                    file.path(wt_ref_dir, paste0(wt_prefix, "DivTimeNorm.tsv"))
                )
            ),
            list(
                title = "CC Deviations",
                data = read.table(file.path(output_dir, paste0(Name, "CCdev.txt")), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)[, grepl("^[0-9]|X[0-9]", colnames(read.table(file.path(output_dir, paste0(Name, "CCdev.txt")), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))), drop = FALSE]
            )
        )

        # Filter CC Deviations specifically to remove Sulston (if not already handled by grepl)
        if (!is.null(vis_list[[3]]$data)) {
            vis_list[[3]]$data <- vis_list[[3]]$data[, !colnames(vis_list[[3]]$data) %in% c("20081128_sulston", "X20081128_sulston"), drop = FALSE]
        }

        run_dimensionality_reduction_report(vis_list, metadata,
            output_pdf = file.path(output_dir, paste0(Name, "_Lineage_Visualizations_Grid.pdf")),
            report_title = "Lineage Kinetics",
            mutant_ids = MutantEmbryos,
            mutant_label = Name,
            Expression = ExpVals
        )
    } else {
        message("AnalyzeDivTimes: dim_reduction=FALSE; skipping UMAP/PCA report.")
    }

    return(list(cc_dev = Deviations, Cells = Cells))
}
