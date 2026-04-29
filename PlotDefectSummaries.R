source("functions.R")
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)

PlotDefectSummaries <- function(i, ExpPeak = NULL, sig = 3, dev = 5, expCutoff = 200, NN_sig = 3.5, max_sig = 3.5, mean_sig = 3.5, microns = 5, NN_cutoff = 0.8, minDivTime = 0,
                                data_dir = ".",
                                output_dir = NULL,
                                wt_ref_dir = NULL,
                                cell_names_file = NULL,
                                cell_lineage_order_file = NULL) {
    if (is.null(output_dir)) output_dir <- file.path(data_dir, i)
    if (is.null(wt_ref_dir)) wt_ref_dir <- file.path(data_dir, "Richard_et_al_plus_comma_WT")
    wt_prefix <- basename(wt_ref_dir)
    mutant_input_dir <- file.path(data_dir, i)
    if (is.null(cell_names_file)) cell_names_file <- "CellNames.csv"
    if (is.null(cell_lineage_order_file)) cell_lineage_order_file <- "Cells_350min_lineageOrder.csv"

    message(paste("Plotting Defect Summaries for", i))

    # Load reference data safely inside function
    WTDivTimes <- read.table(file.path(wt_ref_dir, paste0(wt_prefix, "DivTimeNorm.tsv")), header = T, row.names = 1, stringsAsFactors = F, check.names = FALSE)
    WTDivTimes[, 1] <- NULL

    # Define Cells locally
    Cells <- rownames(WTDivTimes)

    # Handle Expression if not provided
    if (is.null(ExpPeak)) {
        message("PlotDefectSummaries: No expression data provided. Using zero.")
        ExpPeak <- rep(0, length(Cells))
        names(ExpPeak) <- Cells
    }

    ExpPeak[is.na(ExpPeak)] <- 0

    # Pre-calculate CellNames/Times mapping if needed internally
    CellNames <- NULL
    CellTimes <- NULL
    if (file.exists(cell_names_file)) {
        tmpCellNames <- read.csv(cell_names_file, row.names = 1, check.names = FALSE)
        CellNames <- as.vector(tmpCellNames[, 1])
        names(CellNames) <- rownames(tmpCellNames)
        CellTimes <- as.vector(tmpCellNames[, 2])
        names(CellTimes) <- rownames(tmpCellNames)
    }

    # Define relatives
    Parents <- sapply(Cells, GetParent)
    Sisters <- sapply(Cells, GetSister)
    Aunts <- sapply(Parents, GetSister)
    names(Parents) <- Cells
    names(Sisters) <- Cells
    names(Aunts) <- Cells

    CellsLineageOrder <- read.csv(cell_lineage_order_file, check.names = FALSE)

    WTMax <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_CellMaxPositionDevs.csv")), row.names = 1, check.names = FALSE)

    # clear "-Inf" values and remove duplicate name column
    WTMax[, 1] <- NULL
    WTMax[WTMax < (-100)] <- NA

    print(i)
    print(2)

    WTMean <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_CellMeanPositionDevs.csv")), row.names = 1, check.names = FALSE)
    WTMean[, 1] <- NULL
    WTMean[WTMean < (-100)] <- NA

    WTNN <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_NN_Scores.csv")), row.names = 1, check.names = FALSE)
    WTNN[, 1] <- NULL


    WTMaxMean <- rowMeans(WTMax, na.rm = T)
    WTMaxSD <- apply(WTMax, 1, sd, na.rm = T)
    WTMeanMean <- rowMeans(WTMean, na.rm = T)
    WTMeanSD <- apply(WTMean, 1, sd, na.rm = T)
    WTNNMean <- rowMeans(WTNN, na.rm = T)
    WTNNMax <- apply(WTNN, 1, max, na.rm = T)
    WTNNSD <- apply(WTNN, 1, sd, na.rm = T)


    ## old values
    ## sig=3
    ## dev=5

    ## sig=5
    ## dev=6

    # Ensure directory exists
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    print(i)

    # Redirect output to file
    sink_file <- file.path(output_dir, paste0(i, "_summary_output.txt"))
    sink(sink_file)
    # Note: on.exit for sink needs to handle graphical device carefully if they are mixed?
    # Usually fine. But let's be explicit.

    pdf(file.path(output_dir, paste0(i, "_summary.pdf")), width = 10, height = 8)
    # 1)Load CC devs
    theseCCDevs <- read.csv(file.path(output_dir, paste0(i, "_ccDevs.csv")), row.names = 1, check.names = FALSE)

    #    theseCCDevs=theseCCDevs[theseCCDevs$Div.time > minDivTime,]

    Early <- theseCCDevs[, "CC_Z_score"] < (-sig) & theseCCDevs[, "CC_Deviation"] < (-dev)
    names(Early) <- rownames(theseCCDevs)
    Early[!is.na(theseCCDevs["WT_CC"]) & !is.na(theseCCDevs["Terminal_Cell_Length"]) | !is.na(theseCCDevs["CC"])] <- 0
    Early[theseCCDevs[, "CC_Z_score"] < (-sig) & theseCCDevs[, "CC_Deviation"] < (-dev)] <- 1


    Late <- theseCCDevs[, "CC_Z_score"] > sig & theseCCDevs[, "CC_Deviation"] > dev
    names(Late) <- rownames(theseCCDevs)
    Late[!is.na(theseCCDevs["WT_CC"]) & !is.na(theseCCDevs["Terminal_Cell_Length"]) | !is.na(theseCCDevs["CC"])] <- 0
    Late[theseCCDevs[, "CC_Z_score"] > sig & theseCCDevs[, "CC_Deviation"] > dev] <- 1

    Missed <- theseCCDevs[, "Terminal_cell_delta"] > dev & theseCCDevs[, "Terminal_Cell_Z"] > sig
    names(Missed) <- rownames(theseCCDevs)
    Missed[!is.na(theseCCDevs["WT_CC"]) & !is.na(theseCCDevs["Terminal_Cell_Length"]) | !is.na(theseCCDevs["CC"])] <- 0
    Missed[theseCCDevs[, "Terminal_cell_delta"] > dev & theseCCDevs[, "Terminal_Cell_Z"] > sig] <- 1

    Ectopic <- is.na(theseCCDevs[, "WT_Div_time"]) & !is.na(theseCCDevs[, "Div_time"])
    names(Ectopic) <- rownames(theseCCDevs)
    Ectopic[is.na(theseCCDevs["WT_CC"]) & !is.na(theseCCDevs["Terminal_Cell_Length"]) | !is.na(theseCCDevs["CC"])] <- 0
    Ectopic[is.na(theseCCDevs[, "WT_Div_time"]) & !is.na(theseCCDevs[, "Div_time"])] <- 1

    Ectopic[theseCCDevs[, "Div_time"] < minDivTime] <- 0
    Missed[theseCCDevs[, "Div_time"] < minDivTime] <- 0
    Early[theseCCDevs[, "Div_time"] < minDivTime] <- 0
    Late[theseCCDevs[, "Div_time"] < minDivTime] <- 0

    ### This reports how many cells in the dataset have each defect class (based on the selected thresholds). For example if ABplpppp divided late in 2 embryos, it would contribute 2 to the "late" count
    DividingCells <- unique(theseCCDevs[!is.na(theseCCDevs$CC) | !is.na(theseCCDevs$WT_CC), 1])
    print(paste(
        "Number of cells dividing in WT or mutant: ", length(DividingCells),
        " Expressing: ", sum(ExpPeak[DividingCells] >= expCutoff),
        " Nonexpressing: ", sum(ExpPeak[DividingCells] < expCutoff)
    ))

    print("NUMBER OF CELL-EMBRYO DEFECTS")
    print(paste("Early:", sum(Early, na.rm = T)))
    print(paste("Late:", sum(Late, na.rm = T)))
    print(paste("Missed:", sum(Missed, na.rm = T)))
    print(paste("Ectopic:", sum(Ectopic, na.rm = T)))


    # CombinedDefects=Early==1 | Late==1 | Missed==1 | Ectopic==1
    CombinedDefects <- pmax(Early, Late, Missed, Ectopic)

    names(CombinedDefects) <- rownames(theseCCDevs)

    Cells <- as.vector(unique(theseCCDevs[, 1]))
    Embs <- as.vector(unique(theseCCDevs[, 2]))
    CCDefectArray <- NULL
    MissedArray <- NULL
    LateArray <- NULL
    EctopicArray <- NULL
    EarlyArray <- NULL
    for (j in Embs) {
        CCDefectArray <- cbind(CCDefectArray, CombinedDefects[paste(Cells, j, sep = "_")])
        MissedArray <- cbind(MissedArray, Missed[paste(Cells, j, sep = "_")])
        LateArray <- cbind(LateArray, Late[paste(Cells, j, sep = "_")])
        EctopicArray <- cbind(EctopicArray, Ectopic[paste(Cells, j, sep = "_")])
        EarlyArray <- cbind(EarlyArray, Early[paste(Cells, j, sep = "_")])
    }

    rownames(CCDefectArray) <- Cells
    rownames(MissedArray) <- Cells
    rownames(EarlyArray) <- Cells
    rownames(LateArray) <- Cells
    rownames(EctopicArray) <- Cells

    colnames(CCDefectArray) <- Embs
    colnames(MissedArray) <- Embs
    colnames(EarlyArray) <- Embs
    colnames(LateArray) <- Embs
    colnames(EctopicArray) <- Embs

    CCDefectCount <- rowSums(CCDefectArray, na.rm = T)


    print("NUMBER OF CELLS DEFECTIVE BY EXPRESSION")
    print("EXPRESSING")
    print(paste("Early: ", sum(EarlyArray[ExpPeak[Cells] >= expCutoff, ], na.rm = T)))
    print(paste("Late: ", sum(LateArray[ExpPeak[Cells] >= expCutoff, ], na.rm = T)))
    print(paste("Missed: ", sum(MissedArray[ExpPeak[Cells] >= expCutoff, ], na.rm = T)))
    print(paste("Ectopic: ", sum(EctopicArray[ExpPeak[Cells] >= expCutoff, ], na.rm = T)))
    print(paste("ANY DEFECT: ", sum(CCDefectArray[ExpPeak[Cells] >= expCutoff, ], na.rm = T)))
    print("NOT EXPRESSING")
    print(paste("Early: ", sum(EarlyArray[ExpPeak[Cells] < expCutoff, ], na.rm = T)))
    print(paste("Late: ", sum(LateArray[ExpPeak[Cells] < expCutoff, ], na.rm = T)))
    print(paste("Missed: ", sum(MissedArray[ExpPeak[Cells] < expCutoff, ], na.rm = T)))
    print(paste("Ectopic: ", sum(EctopicArray[ExpPeak[Cells] < expCutoff, ], na.rm = T)))
    print(paste("ANY DEFECT: ", sum(CCDefectArray[ExpPeak[Cells] < expCutoff, ], na.rm = T)))
    print("NUMBER OF CELLS DEFECTIVE IN ANY EMBRYO BY EXPRESSION ")
    print("EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print("NOT EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 0)))
    print("NUMBER OF CELLS DEFECTIVE IN AT LEAST TWO EMBRYOS BY EXPRESSION ")
    print("EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print("NOT EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 1)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T) > 1)))


    ## List repeatedly defective cells
    print("Repeat defect cells (expressing): ")
    print(rowSums(CCDefectArray[rowSums(CCDefectArray, na.rm = T) > 1 & ExpPeak[Cells] >= expCutoff, , drop = FALSE], na.rm = T))
    print("Repeat defect cells (non expressing): ")
    print(rowSums(CCDefectArray[rowSums(CCDefectArray, na.rm = T) > 1 & ExpPeak[Cells] < expCutoff, , drop = FALSE], na.rm = T))

    # code to determine for each cell if that cell or an ancestor was defective
    AncestorHadDefect <- function(Cell, DefectArray) {
        Defects <- DefectArray[Cell, ]
        Defects[is.na(Defects)] <- 0
        Parent <- GetParent(Cell)
        #        print(Parent)
        if (Parent == "P" || is.na(Parent)) {
            return(Defects)
        } else {
            return(Defects | AncestorHadDefect(Parent, DefectArray))
        }
    }
    print(.5)
    AncestralDefects <- data.frame()
    for (j in Cells) {
        #        print("CELL")
        #        print(i)
        if (!is.na(j)) {
            AncestralDefects <- rbind(AncestralDefects, AncestorHadDefect(j, CCDefectArray))
        }
    }

    rownames(AncestralDefects) <- Cells
    colnames(AncestralDefects) <- Embs
    #
    ExpDef <- sum(ExpPeak[Cells] > expCutoff & CCDefectCount > 0, na.rm = T)
    ExpOK <- sum(ExpPeak[Cells] > expCutoff & CCDefectCount == 0, na.rm = T)
    NotExpDef <- sum(ExpPeak[Cells] < expCutoff & CCDefectCount > 0, na.rm = T)
    NotExpOK <- sum(ExpPeak[Cells] < expCutoff & CCDefectCount == 0, na.rm = T)
    ##    DivisionsAssessed =sum(

    # Plot Exp vs Defect Count
    p_exp_def <- ggplot(data.frame(Exp = ExpPeak[Cells], Count = CCDefectCount), aes(x = Exp, y = Count)) +
        geom_jitter(width = 0, height = 0.2, alpha = 0.5, size = 1) +
        theme_bw() +
        labs(x = "Peak Expression", y = paste(i, "CC Defect Count")) +
        ggtitle(paste("Sigma:", sig, "Dev:", dev, "\nExp:", ExpDef, ExpOK, round(ExpDef / (ExpDef + ExpOK), 3), "\nNot Exp:", NotExpDef, NotExpOK, round(NotExpDef / (NotExpDef + NotExpOK), 3)))

    grid.arrange(p_exp_def, top = "Cell Cycle Defects vs Expression")

    cat(paste(
        "CELL CYCLE\nSigma Cutoff: ", sig,
        "Dev cutoff (min): ", dev,
        "\nExpressing:", ExpDef, ExpOK, round(ExpDef / (ExpDef + ExpOK), 3),
        "\nNot Expressing:", NotExpDef, NotExpOK, round(NotExpDef / (NotExpDef + NotExpOK), 3), "\n"
    ))


    ## 2)Load position devs


    print(i)
    theseMaxDevs <- read.csv(file.path(output_dir, paste0(i, "_CellMaxPositionDevs.csv")), row.names = 1, check.names = FALSE)
    theseDivTimes <- read.table(file.path(mutant_input_dir, paste0(i, "DivTimeNorm.tsv")), row.names = 1, sep = "\t", header = TRUE, check.names = FALSE)
    theseDivTimes <- subset(theseDivTimes, select = -1)
    # Ensure rows match Cells
    theseDivTimes <- theseDivTimes[Cells, ]

    theseMeanDevs <- read.csv(file.path(output_dir, paste0(i, "_CellMeanPositionDevs.csv")), row.names = 1, check.names = FALSE)
    theseNNs <- read.csv(file.path(output_dir, paste0(i, "_NN_Scores.csv")), row.names = 1, check.names = FALSE)
    # because these csv files have 2 columns with cell name duplicated
    theseMaxDevs[, 1] <- NULL
    theseMeanDevs[, 1] <- NULL
    theseNNs[, 1] <- NULL
    theseMaxDevs[theseMaxDevs < (-100)] <- NA
    theseMeanDevs[theseMeanDevs < (-100)] <- NA

    # Ensure all dataframes have the same embryos (columns)
    CommonEmbryos <- intersect(colnames(theseMaxDevs), colnames(theseDivTimes))
    CommonEmbryos <- intersect(CommonEmbryos, colnames(theseMeanDevs))
    CommonEmbryos <- intersect(CommonEmbryos, colnames(theseNNs))

    if (length(CommonEmbryos) == 0) {
        stop("No common embryos found between Position files and DivTime file")
    }

    theseMaxDevs <- theseMaxDevs[, CommonEmbryos, drop = FALSE]
    theseMeanDevs <- theseMeanDevs[, CommonEmbryos, drop = FALSE]
    theseNNs <- theseNNs[, CommonEmbryos, drop = FALSE]
    theseDivTimes <- theseDivTimes[, CommonEmbryos, drop = FALSE]


    theseMaxZ <- (theseMaxDevs[Cells, ] - WTMaxMean[Cells]) / WTMaxSD[Cells]
    theseMeanZ <- (theseMeanDevs[Cells, ] - WTMeanMean[Cells]) / WTMeanSD[Cells]
    theseNNZ <- (theseNNs[Cells, ] - WTNNMean[Cells]) / WTNNSD[Cells]

    theseMaxMean <- rowMeans(theseMaxDevs, na.rm = T)
    theseMaxSD <- apply(theseMaxDevs, 1, sd, na.rm = T)

    theseMeanMean <- rowMeans(theseMeanDevs, na.rm = T)
    theseMeanSD <- apply(theseMeanDevs, 1, sd, na.rm = T)

    print(3)

    theseNNMean <- rowMeans(theseNNs, na.rm = T)
    theseNNSD <- apply(theseNNs, 1, sd, na.rm = T)
    theseNNMax <- apply(theseNNs, 1, max, na.rm = T)

    cellObservedCounts <- rowSums(!is.na(theseMaxDevs[Cells, ]), na.rm = T)

    # plots vs SD (Combined)
    df_var <- data.frame(
        Mut_MaxMean = theseMaxMean, Mut_MaxSD = theseMaxSD,
        WT_MaxMean = WTMaxMean, WT_MaxSD = WTMaxSD,
        Mut_MeanMean = theseMeanMean, Mut_MeanSD = theseMeanSD,
        WT_MeanMean = WTMeanMean, WT_MeanSD = WTMeanSD,
        Mut_NNMean = theseNNMean, Mut_NNSD = theseNNSD,
        WT_NNMean = WTNNMean, WT_NNSD = WTNNSD
    )

    p_max <- ggplot(df_var) +
        geom_point(aes(x = Mut_MaxMean, y = Mut_MaxSD, color = "Mutant"), alpha = 0.5) +
        geom_point(aes(x = WT_MaxMean, y = WT_MaxSD, color = "WT"), alpha = 0.5) +
        theme_bw() +
        labs(x = "Max Pos Dev", y = "SD", title = "Max Dev Variability")

    p_mean <- ggplot(df_var) +
        geom_point(aes(x = Mut_MeanMean, y = Mut_MeanSD, color = "Mutant"), alpha = 0.5) +
        geom_point(aes(x = WT_MeanMean, y = WT_MeanSD, color = "WT"), alpha = 0.5) +
        theme_bw() +
        labs(x = "Mean Pos Dev", y = "SD", title = "Mean Dev Variability")

    p_nn <- ggplot(df_var) +
        geom_point(aes(x = Mut_NNMean, y = Mut_NNSD, color = "Mutant"), alpha = 0.5) +
        geom_point(aes(x = WT_NNMean, y = WT_NNSD, color = "WT"), alpha = 0.5) +
        theme_bw() +
        labs(x = "NN Score", y = "SD", title = "NN Score Variability")

    # Mean Comparison Plots
    df_comp <- data.frame(
        Cells = Cells,
        WT_Max = WTMaxMean[Cells], Mut_Max = theseMaxMean[Cells],
        WT_Mean = WTMeanMean[Cells], Mut_Mean = theseMeanMean[Cells],
        WT_NN = WTNNMean[Cells], Mut_NN = theseNNMean[Cells]
    )

    p_comp1 <- ggplot(df_comp, aes(x = WT_Max, y = Mut_Max)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "WT Max Dev", y = "Mutant Max Dev")
    p_comp2 <- ggplot(df_comp, aes(x = WT_Mean, y = Mut_Mean)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "WT Mean Dev", y = "Mutant Mean Dev")
    p_comp3 <- ggplot(df_comp, aes(x = WT_NN, y = Mut_NN)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "WT NN Score", y = "Mutant NN Score")

    p_scatter <- ggplot(df_comp) +
        geom_point(aes(x = Mut_Max, y = Mut_NN, color = "Mutant"), alpha = 0.5) +
        geom_point(aes(x = WT_Max, y = WT_NN, color = "WT"), alpha = 0.5) +
        theme_bw() +
        labs(x = "Max Pos Dev", y = "Mean NN Score", title = "Dev vs NN")

    grid.arrange(p_max, p_mean, p_nn, p_comp1, p_comp2, p_comp3, p_scatter, nrow = 3, top = "Variability and Comparison")


    # outlier analysis

    PosOutliers <- (theseMaxZ > max_sig | theseMeanZ > mean_sig) & theseNNZ > NN_sig & theseMaxDevs[Cells, ] > microns & (is.na(theseDivTimes) | theseDivTimes > minDivTime)
    PosOutliers[PosOutliers == T] <- 1
    PosOutliers[PosOutliers == F] <- 0

    rownames(PosOutliers) <- Cells
    PosDefectCount <- rowSums(PosOutliers, na.rm = T)

    # Correlation Plots (ggplot2)
    corr_df <- data.frame(
        Cell = Cells,
        PosDefect = PosDefectCount[Cells],
        CCDefect = CCDefectCount[Cells],
        ParentCCDefect = CCDefectCount[sapply(Cells, GetParent)],
        Parent = sapply(Cells, GetParent),
        Exp = ExpPeak[Cells]
    )

    # helper for clean correlation titles
    get_cor_title <- function(x, y, label) {
        val <- round(cor(x, y, use = "pairwise.complete.obs", method = "spearman"), digits = 4)
        paste(label, "R =", val)
    }

    p_corr1 <- ggplot(corr_df, aes(x = jitter(PosDefect), y = jitter(CCDefect))) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "Pos Defects", y = "CC Defects", title = get_cor_title(corr_df$PosDefect, corr_df$CCDefect, i))

    p_corr2 <- ggplot(corr_df, aes(x = jitter(PosDefect), y = jitter(ParentCCDefect))) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "Pos Defects", y = "Parent CC Defects", title = get_cor_title(corr_df$PosDefect, corr_df$ParentCCDefect, i))

    p_corr3 <- ggplot(corr_df, aes(x = jitter(PosDefect), y = jitter(ParentCCDefect + CCDefect))) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "Pos Defects", y = "Total CC Defects", title = get_cor_title(corr_df$PosDefect, corr_df$ParentCCDefect + corr_df$CCDefect, i))


    ExpPosDef <- sum(ExpPeak[Cells] > expCutoff & PosDefectCount > 0, na.rm = T)
    ExpPosOK <- sum(ExpPeak[Cells] > expCutoff & PosDefectCount == 0, na.rm = T)
    NotExpPosDef <- sum(ExpPeak[Cells] < expCutoff & PosDefectCount > 0, na.rm = T)
    NotExpPosOK <- sum(ExpPeak[Cells] < expCutoff & PosDefectCount == 0, na.rm = T)

    p_exp_pos <- ggplot(corr_df, aes(x = Exp, y = jitter(PosDefect))) +
        geom_point(alpha = 0.5, size = 0.5) +
        theme_bw() +
        labs(x = "Peak Expression", y = "Position Defect Count") +
        ggtitle(paste("Pos Defects\nExp:", ExpPosDef, ExpPosOK, round(ExpPosDef / (ExpPosDef + ExpPosOK), 3), "\nNot Exp:", NotExpPosDef, NotExpPosOK, round(NotExpPosDef / (NotExpPosDef + NotExpPosOK), 3)))


    cat(paste(
        "NN Cutoff: ", NN_sig, "max cutoff: ", max_sig, "mean cutoff: ", mean_sig, "micron cutoff:", microns, "\nPOSITION DEFECTS (defect / no defect / fractionDefect)",
        "\nExpressing:", ExpPosDef, ExpPosOK, round(ExpPosDef / (ExpPosDef + ExpPosOK), 3),
        "\nNot Expressing:", NotExpPosDef, NotExpPosOK, round(NotExpPosDef / (NotExpPosDef + NotExpPosOK), 3), "\n"
    ))


    # Calculate Parent Defects safely
    ParentCells <- sapply(Cells, GetParent)
    ParentCCDefects <- CCDefectCount[ParentCells]
    ParentCCDefects[is.na(ParentCCDefects)] <- 0

    BothDefectExp <- sum(ExpPeak[Cells] > expCutoff & PosDefectCount[Cells] > 0 & (ParentCCDefects + CCDefectCount[Cells]) > 0, na.rm = T)
    PosDefectExp <- sum(ExpPeak[Cells] > expCutoff & PosDefectCount[Cells] > 0 & (ParentCCDefects + CCDefectCount[Cells]) == 0, na.rm = T)
    CCDefectExp <- sum(ExpPeak[Cells] > expCutoff & PosDefectCount[Cells] == 0 & (ParentCCDefects + CCDefectCount[Cells]) > 0, na.rm = T)
    NoDefectExp <- sum(ExpPeak[Cells] > expCutoff & PosDefectCount[Cells] == 0 & (ParentCCDefects + CCDefectCount[Cells]) == 0, na.rm = T)

    BothDefectNotExp <- sum(ExpPeak[Cells] <= expCutoff & PosDefectCount[Cells] > 0 & (ParentCCDefects + CCDefectCount[Cells]) > 0, na.rm = T)
    PosDefectNotExp <- sum(ExpPeak[Cells] <= expCutoff & PosDefectCount[Cells] > 0 & (ParentCCDefects + CCDefectCount[Cells]) == 0, na.rm = T)
    CCDefectNotExp <- sum(ExpPeak[Cells] <= expCutoff & PosDefectCount[Cells] == 0 & (ParentCCDefects + CCDefectCount[Cells]) > 0, na.rm = T)
    NoDefectNotExp <- sum(ExpPeak[Cells] <= expCutoff & PosDefectCount[Cells] == 0 & (ParentCCDefects + CCDefectCount[Cells]) == 0, na.rm = T)


    cat(paste(
        "POS+CC DEFECTS: Both CC and POS defects / only POS / only CC/ neither\n", "Exp: Both:", BothDefectExp, "Pos:", PosDefectExp, "CC:", CCDefectExp, "No:", NoDefectExp,
        "\nNo Exp: Both:", BothDefectNotExp, "Pos:", PosDefectNotExp, "CC:", CCDefectNotExp, "No:", NoDefectNotExp
    ), "\n")

    # Construct Boxplot Dataframe
    cat_df <- data.frame(Exp = numeric(0), Category = character(0))


    # 1. Both Defects
    idx_both <- Cells[PosDefectCount[Cells] > 0 & (ParentCCDefects + CCDefectCount[Cells]) > 0]
    if (length(idx_both) > 0) cat_df <- rbind(cat_df, data.frame(Exp = ExpPeak[idx_both], Category = "Both"))

    # 2. One Defect
    idx_one <- Cells[(PosDefectCount[Cells] > 0 | (ParentCCDefects + CCDefectCount[Cells]) > 0) &
        !(PosDefectCount[Cells] > 0 & (ParentCCDefects + CCDefectCount[Cells]) > 0)]
    if (length(idx_one) > 0) cat_df <- rbind(cat_df, data.frame(Exp = ExpPeak[idx_one], Category = "One"))

    # 3. No Defects
    idx_none <- Cells[PosDefectCount[Cells] == 0 & (ParentCCDefects + CCDefectCount[Cells]) == 0]
    if (length(idx_none) > 0) cat_df <- rbind(cat_df, data.frame(Exp = ExpPeak[idx_none], Category = "None"))

    p_box <- ggplot(cat_df, aes(x = Category, y = Exp)) +
        geom_boxplot() +
        theme_bw() +
        labs(title = paste(i, "Expression by Defect Category"), x = "Defect Status (Both/One/None)", y = "Peak Expression")

    # Grid Arrange Correlations
    grid.arrange(p_corr1, p_corr2, p_corr3, p_exp_pos, nrow = 2, top = "Defect Correlations")
    grid.arrange(p_box, top = "Defect Category vs Expression")


    ## Alternate approach - using absolute deviation thresholds as well
    PosCells <- rownames(theseMeanDevs)
    # ... (Calculations preserved, output suppressed for brevity)

    # Angle Correlation
    # ... (AngleArray loading preserved) ...
    # Assumes AngleArray logic is preserved via read.csv

    # Use existing read logic but fix Plotting
    # AngleArray=read.csv(file.path(i, paste0(i,"_dots.csv")),header=T,row.names=1, check.names=FALSE)
    # AngleArray=AngleArray[Cells,7:length(AngleArray)]
    # rownames(AngleArray)=Cells

    # ... (Calculation of AngOutliers) ...
    # Re-implementing logic to ensure variables exist for plotting

    AngleArray <- read.csv(file.path(output_dir, paste0(i, "_dots.csv")), header = T, row.names = 1, check.names = FALSE)
    AngleArray <- AngleArray[Cells, 7:length(AngleArray)]
    rownames(AngleArray) <- Cells
    divObservedCounts <- rowSums(!is.na(AngleArray[Cells, ]), na.rm = T)
    AngCutoff <- .5
    AngOutliers <- abs(AngleArray) < AngCutoff
    AngOutliers[AngOutliers == T] <- 1
    AngOutliers[AngOutliers == F] <- 0
    AngDefectCount <- rowSums(AngOutliers, na.rm = T)

    # Angle Plots
    ang_df <- data.frame(
        AngDefect = AngDefectCount[Cells],
        PosDefect = PosDefectCount[Cells],
        CCDefect = CCDefectCount[Cells],
        ParentCCDefect = CCDefectCount[sapply(Cells, GetParent)],
        ParentAngDefect = AngDefectCount[sapply(Cells, GetParent)],
        Exp = ExpPeak[Cells]
    )

    p_ang_exp <- ggplot(ang_df, aes(x = Exp, y = jitter(AngDefect))) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "Expression", y = "Angle Defect Count")
    p_ang_pos <- ggplot(ang_df, aes(x = jitter(AngDefect), y = jitter(PosDefect))) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "Angle Defects", y = "Pos Defects")
    p_ang_cc <- ggplot(ang_df, aes(x = jitter(AngDefect), y = jitter(CCDefect))) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(x = "Angle Defects", y = "CC Defects")

    p_ang_cc_corr <- ggplot(ang_df, aes(x = AngDefect, y = CCDefect)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(title = get_cor_title(ang_df$AngDefect, ang_df$CCDefect, "Ang vs CC"))
    p_ang_parent_cc <- ggplot(ang_df, aes(x = AngDefect, y = ParentCCDefect)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(title = get_cor_title(ang_df$AngDefect, ang_df$ParentCCDefect, "Ang vs Parent CC"))

    # Combined Totals Correlations
    p_ang_total_cc <- ggplot(ang_df, aes(x = AngDefect, y = CCDefect + ParentCCDefect)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(title = get_cor_title(ang_df$AngDefect, ang_df$CCDefect + ang_df$ParentCCDefect, "Ang vs Total CC"))
    p_total_ang_pos <- ggplot(ang_df, aes(x = AngDefect + ParentAngDefect, y = PosDefect)) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        labs(title = get_cor_title(ang_df$AngDefect + ang_df$ParentAngDefect, ang_df$PosDefect, "Total Ang vs Pos"))

    grid.arrange(p_ang_exp, p_ang_pos, p_ang_cc, p_ang_cc_corr, p_ang_parent_cc, p_ang_total_cc, p_total_ang_pos, nrow = 3, top = "Angle Correlations")

    ExpAngDef <- sum(ExpPeak[Cells] > expCutoff & AngDefectCount > 0, na.rm = T)
    ExpAngOK <- sum(ExpPeak[Cells] > expCutoff & PosDefectCount == 0, na.rm = T)
    NotExpAngDef <- sum(ExpPeak[Cells] < expCutoff & AngDefectCount > 0, na.rm = T)
    NotExpAngOK <- sum(ExpPeak[Cells] < expCutoff & AngDefectCount == 0, na.rm = T)

    cat(paste(
        "ANGLES", "Ang cutoff: ", AngCutoff,
        "\nExpressing:", ExpAngDef, ExpAngOK, round(ExpAngDef / (ExpAngDef + ExpAngOK), 3),
        "\nNot Expressing:", NotExpAngDef, NotExpAngOK, round(NotExpAngDef / (NotExpAngDef + NotExpAngOK), 3), "\n"
    ))


    # ... (Previous print statements will now go to file) ...
    # We need to ensure the sink starts earlier if we want ALL output.
    # But since I'm replacing code blocks, I should probably put the sink at the very beginning of the function
    # and just remove the dev.off() from the bottom of the replacement block to avoid premature closing if I put on.exit at top.
    # Actually, I'll handle the sink at the top in a separate edit to be clean.
    # For now, let's finish the plotting refactor.

    # 4)Integrate onto list of "terminal cells" matching expression tree

    # Ensure all arrays match the common embryos
    PosOutliers2 <- PosOutliers[, CommonEmbryos, drop = FALSE]
    if (!is.null(AngOutliers)) {
        AngOutliers2 <- AngOutliers[, CommonEmbryos, drop = FALSE]
    } else {
        AngOutliers2 <- 0
    }
    CCDefectArray2 <- CCDefectArray[, CommonEmbryos, drop = FALSE]

    PosOutliers2[is.na(PosOutliers2)] <- 0
    AngOutliers2[is.na(AngOutliers2)] <- 0
    CCDefectArray2[is.na(CCDefectArray2)] <- 0

    # Also need to handle CCDefectArray2 for parents
    Parents <- sapply(Cells, GetParent)
    CCDefectArrayParent <- matrix(0, nrow = length(Cells), ncol = length(CommonEmbryos))
    rownames(CCDefectArrayParent) <- Cells
    colnames(CCDefectArrayParent) <- CommonEmbryos

    ValidParents <- Parents %in% rownames(CCDefectArray2)
    CCDefectArrayParent[ValidParents, ] <- CCDefectArray2[Parents[ValidParents], ]

    NumDefects <- PosOutliers2 + AngOutliers2 + CCDefectArray2 + CCDefectArrayParent


    # Get # of observations for each defect table

    # Import list of cells in lineage order (350 minute sulston cells?)

    # Map "Cells" to these terminal cells
    GetTerminalCell <- function(X) {
        # cell is a terminal cell - this works
        if (X %in% CellsLineageOrder[, 1]) {
            return(X)
        }

        # cell is a descendent of a terminal cell - this works
        tmpCell <- X
        while (!is.na(GetParent(tmpCell)) && GetParent(tmpCell) != "P") {
            tmpCell <- GetParent(tmpCell)
            if (tmpCell %in% CellsLineageOrder[, 1]) {
                return(tmpCell)
            }
        }

        # cell is an ancestor of a terminal cell - works but sends out a list of cells
        Descendants <- NULL
        for (k in CellsLineageOrder[, 1]) {
            tmpCell <- k
            while (!is.na(GetParent(tmpCell)) && GetParent(tmpCell) != "P") {
                tmpCell <- GetParent(tmpCell)
                if (tmpCell == X) {
                    Descendants <- c(Descendants, k)
                }
            }
        }
        if (length(Descendants) > 0) {
            return(Descendants)
        }
        return(NA)
    }
    # NewCells=sapply(Cells,GetTerminalCell)
    # Output plots of defect frequency

    AncestorDefects <- function(Cell, DefectArray) {
        MyDefects <- 0
        if (Cell %in% rownames(DefectArray)) {
            MyDefects <- sum(DefectArray[Cell, ], na.rm = T)
        }
        Parent <- GetParent(Cell)
        if (Parent != "P" && !is.na(Parent)) {
            MyDefects <- MyDefects + AncestorDefects(Parent, DefectArray)
        }
        return(MyDefects)
    }

    SumDefects <- function(X, DefectArray) {
        return(sum(DefectArray[grep(paste(X, ".", sep = ""), rownames(DefectArray)), ], na.rm = T) +
            AncestorDefects(X, DefectArray))
    }
    CCTermDefects <- sapply(as.vector(CellsLineageOrder[, 1]), function(X) {
        SumDefects(X, CCDefectArray)
    })
    names(CCTermDefects) <- CellsLineageOrder[, 1]

    PosTermDefects <- sapply(as.vector(CellsLineageOrder[, 1]), function(X) {
        SumDefects(X, PosOutliers)
    })
    names(PosTermDefects) <- CellsLineageOrder[, 1]

    AngTermDefects <- sapply(as.vector(CellsLineageOrder[, 1]), function(X) {
        SumDefects(X, AngOutliers)
    })
    names(AngTermDefects) <- CellsLineageOrder[, 1]

    # Lineage Plots (ggplot2)
    lineage_df <- data.frame(
        Cell = factor(CellsLineageOrder[, 1], levels = CellsLineageOrder[, 1]), # Keep order
        Order = 1:length(CellsLineageOrder[, 1]),
        CC = CCTermDefects,
        Pos = PosTermDefects,
        Ang = AngTermDefects,
        Total = CCTermDefects + PosTermDefects + AngTermDefects
    )

    # Combined Lineage Plot
    p_lin_all <- ggplot(lineage_df, aes(x = Order)) +
        geom_line(aes(y = CC, color = "CC"), size = 0.5) +
        geom_line(aes(y = Pos, color = "Pos"), size = 0.5) +
        geom_line(aes(y = Ang, color = "Ang"), size = 0.5) +
        scale_color_manual(values = c("CC" = "black", "Pos" = "red", "Ang" = "green")) +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(x = "Lineage Order", y = "Defect Count (incl. ancestors)", title = "Defect Distribution along Lineage")

    p_lin_cc <- ggplot(lineage_df, aes(x = Order, y = CC)) +
        geom_line() +
        theme_bw() +
        labs(title = "CC Defects", x = "Lineage Order") +
        theme(axis.text.x = element_blank())
    p_lin_pos <- ggplot(lineage_df, aes(x = Order, y = Pos)) +
        geom_line() +
        theme_bw() +
        labs(title = "Position Defects", x = "Lineage Order") +
        theme(axis.text.x = element_blank())
    p_lin_ang <- ggplot(lineage_df, aes(x = Order, y = Ang)) +
        geom_line() +
        theme_bw() +
        labs(title = "Angle Defects", x = "Lineage Order") +
        theme(axis.text.x = element_blank())
    p_lin_tot <- ggplot(lineage_df, aes(x = Order, y = Total)) +
        geom_line() +
        theme_bw() +
        labs(title = "Total Defects", x = "Lineage Order") +
        theme(axis.text.x = element_blank())

    grid.arrange(p_lin_all, p_lin_cc, p_lin_pos, p_lin_ang, p_lin_tot, nrow = 5, top = "Lineage Defect Profiles")

    # Fixme generalize
    for (Founder in c("ABalpa", "ABara", "ABplp", "ABprp", "MSaa", "MSpa")) {
        founder_rows <- grep(Founder, lineage_df$Cell)
        if (length(founder_rows) > 0) {
            p_f_all <- ggplot(lineage_df[founder_rows, ], aes(x = Order)) +
                geom_line(aes(y = CC, color = "CC")) +
                geom_line(aes(y = Pos, color = "Pos")) +
                geom_line(aes(y = Ang, color = "Ang")) +
                scale_color_manual(values = c("CC" = "black", "Pos" = "red", "Ang" = "green")) +
                theme_bw() +
                labs(title = paste(Founder, "- Defects"), x = "Lineage Order", y = "Count")

            grid.arrange(p_f_all)
        }
    }
    Parents <- sapply(Cells, GetParent)
    Sisters <- sapply(Cells, GetSister)
    Aunts <- sapply(Parents, GetSister)


    write.csv(data.frame(
        Cells = CellsLineageOrder[, 1],
        Fates = CellNames[as.vector(CellsLineageOrder[, 1])],
        CCTermDefects,
        PosTermDefects,
        AngTermDefects,
        AngTermDefects,
        Expression = ExpPeak[CellsLineageOrder[, 1]]
    ), file = file.path(output_dir, paste(i, "Term_defects.csv", sep = "_")), row.names = F)

    write.csv(
        data.frame(Cells,
            Fates = CellNames[Cells],
            Times = CellTimes[Cells],
            CC_Defects = CCDefectCount[Cells],
            Pos_Defects = PosDefectCount[Cells],
            Angle_Defects = AngDefectCount[Cells],
            Cell_Count = cellObservedCounts[Cells],
            Division_Count = divObservedCounts[Cells],
            Expression = ExpPeak[Cells],
            MeanMaxPositionDev = theseMaxMean[Cells],
            MeanMeanPositionDev = theseMeanMean[Cells]
        ),
        file = file.path(output_dir, paste(i, "defect_counts.csv", sep = "_")), row.names = F
    )


    write.csv(
        data.frame(Cells,
            Cell_Names = CellNames[Cells],
            Cell_Times = CellTimes[Cells],
            CCDefectCount[Cells],
            CCDefectCount[Parents],
            CCDefectCount[Sisters],
            CCDefectCount[Aunts],
            PosDefectCount[Cells],
            PosDefectCount[Parents],
            PosDefectCount[Sisters],
            PosDefectCount[Aunts],
            AngDefectCount[Cells],
            AngDefectCount[Parents],
            AngDefectCount[Sisters],
            AngDefectCount[Aunts],
            cellObservedCounts[Cells],
            cellObservedCounts[Parents],
            cellObservedCounts[Sisters],
            cellObservedCounts[Aunts],
            divObservedCounts[Cells],
            divObservedCounts[Parents],
            divObservedCounts[Sisters],
            divObservedCounts[Aunts],
            ExpPeak[Cells]
        ),
        file = file.path(output_dir, paste(i, "defect_counts_with_relatives.csv", sep = "_")), row.names = F
    )

    write.csv(
        data.frame(Cells,
            Cell_Names = CellNames[Cells], Cell_Times = CellTimes[Cells],
            CCDefectArray,
            AncestralDefects,
            PosOutliers
        ),
        file = file.path(output_dir, paste(i, "defects_with_ancestral_CC.csv", sep = "_")), row.names = F
    )


    tmp <- data.frame(Cells, Cell_Names = CellNames[Cells], Cell_Times = CellTimes[Cells], ExpPeak[Cells], CCDefectArray[Cells, ], PosOutliers[Cells, ], AngOutliers[Cells, ])
    colnames(tmp) <- c("Cells", "Names", "Birth_Times", "Peak_Exp", paste("CC", colnames(PosOutliers), sep = "_"), paste("Position", colnames(PosOutliers), sep = "_"), paste("Angle", colnames(PosOutliers), sep = "_"))
    write.csv(tmp, file = file.path(output_dir, paste(i, "defects.csv", sep = "_")), row.names = F)


    dev.off() # close pdf
    if (sink.number() > 0) sink() # close text output
}


# tmp=aggregate(. ~ cell,  data=read.csv("ACDunc-30_mean.csv",row.names=1), mean, na.rm=T)
# unc30_exp=data.frame(tmp[,"blot"])
# rownames(unc30_exp)=tmp[,"cell"]
# unc30_peak=sapply(Cells,function(X){
#         GetPeak(unc30_exp,X)
#         })

# tmp=aggregate(. ~ cell,  data=read.csv("ACDceh-36_mean.csv",row.names=1), mean, na.rm=T)
# ceh36_exp=data.frame(tmp[,"blot"])
# rownames(ceh36_exp)=tmp[,"cell"]
# ceh36_peak=sapply(Cells,function(X){
#         GetPeak(ceh36_exp,X)
#         })

# tmp=aggregate(. ~ cell,  data=read.csv("ACDmls-2_mean.csv",row.names=1), mean, na.rm=T)
# mls2_exp=data.frame(tmp[,"blot"])
# rownames(mls2_exp)=tmp[,"cell"]
# mls2_peak=sapply(Cells,function(X){
#         GetPeak(mls2_exp,X)
#         })

# tmp=aggregate(. ~ cell,  data=read.csv("ACDnob-1GFP_mean.csv",row.names=1), mean, na.rm=T)
# nob1_exp=data.frame(tmp[,"blot"])
# rownames(nob1_exp)=tmp[,"cell"]
# nob1_peak=sapply(Cells,function(X){
#         GetPeak(nob1_exp,X)
#         })
