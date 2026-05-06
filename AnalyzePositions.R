source("functions.R")
source("DimensionalityReductionHelpers.R")
library(rgl)
library(gdata)
library(plotrix)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)

source("GetPeak.R")

# Note: ceh36_peak global calculation moved inside AnalyzePositions or handled by LineagePhenotyping.R
# CellNames loading can stay if it's a static mapping, but let's functionalize it if possible.
LoadCellNames <- function(cell_names_file = "CellNames.csv") {
    if (file.exists(cell_names_file)) {
        tmpCellNames <- read.csv(cell_names_file)
        CellNames <- as.vector(tmpCellNames[, 2])
        names(CellNames) <- tmpCellNames[, 1]
        return(CellNames)
    }
    return(NULL)
}

PlotExpVsDev <- function(Name, outfile, exp = NULL, type = "Mean", eGain = 1000, ylim = c(0, 2), xlim = c(0, 10),
                         data_dir = ".", output_dir = NULL) {
    if (is.null(output_dir)) output_dir <- file.path(data_dir, Name)
    if (is.null(exp)) {
        message("PlotExpVsDev: No expression data provided. Skipping color coding.")
    }
    pdf(file.path(output_dir, outfile))

    devs <- read.csv(file.path(output_dir, paste0(Name, "_Cell", type, "PositionDevs.csv")), row.names = 1, check.names = FALSE)
    devs[, 1] <- NULL
    NNs <- read.csv(file.path(output_dir, paste0(Name, "_NN_Scores.csv")), row.names = 1, check.names = FALSE)
    NNs[, 1] <- NULL

    meanDevs <- rowMeans(devs, na.rm = T)
    meanNNs <- rowMeans(NNs, na.rm = T)

    # Helper function for plotting
    draw_plot <- function(data, title, add_labels = FALSE) {
        # Ensure row names are columns for labeling
        data$Cell <- rownames(data)

        p <- ggplot(data, aes(x = dev, y = NN)) +
            geom_point(aes(color = e), size = 1) +
            theme_classic() +
            scale_color_gradient(name = "Expression", low = "blue", high = "red") +
            ggtitle(title) +
            labs(x = "Mean Deviation (microns)", y = "Mean Neighbor Deviation Score") +
            xlim(xlim) +
            ylim(ylim)

        if (add_labels) {
            # Filter for outliers: NN > 0.8 OR Dev > 3
            outliers <- subset(data, NN > 0.8 | dev > 3)
            p <- p + geom_text_repel(data = outliers, aes(label = Cell), size = 2, max.overlaps = Inf)
        }
        return(p)
    }

    # Prepare Data
    mean_data <- data.frame(dev = meanDevs, NN = meanNNs[names(meanDevs)], e = pmax(0, pmin(exp[names(meanDevs)], eGain)))
    rownames(mean_data) <- names(meanDevs)

    # 1. Generate Standard PDF
    pdf(file.path(output_dir, outfile))
    plot(draw_plot(mean_data, paste(Name, type, "Mean", sep = "_"), add_labels = FALSE))

    for (i in colnames(devs)) {
        embryo_data <- data.frame(dev = devs[names(meanDevs), i], NN = NNs[names(meanDevs), i], e = pmax(0, pmin(exp[names(meanDevs)], eGain)))
        rownames(embryo_data) <- names(meanDevs)
        plot(draw_plot(embryo_data, i, add_labels = FALSE))
    }
    dev.off()

    # 2. Generate Labeled PDF
    labeled_outfile <- sub("\\.pdf$", "_labeled.pdf", outfile)
    pdf(file.path(output_dir, labeled_outfile))
    plot(draw_plot(mean_data, paste(Name, type, "Mean Labeled", sep = "_"), add_labels = TRUE))

    for (i in colnames(devs)) {
        embryo_data <- data.frame(dev = devs[names(meanDevs), i], NN = NNs[names(meanDevs), i], e = pmax(0, pmin(exp[names(meanDevs)], eGain)))
        rownames(embryo_data) <- names(meanDevs)
        plot(draw_plot(embryo_data, i, add_labels = TRUE))
    }
    dev.off()

    # 3. Export Data to CSV
    csv_outfile <- sub("\\.pdf$", ".csv", outfile)
    export_df <- devs
    export_df$Mean_Dev <- meanDevs
    export_df$Mean_NN <- meanNNs
    export_df$Expression <- exp[rownames(devs)]
    write.csv(export_df, file = file.path(output_dir, csv_outfile))
}


# --- Helper Functions ---

LoadPositionData <- function(FilePath) {
    if (!file.exists(FilePath)) stop(paste("File not found:", FilePath))
    Positions <- read.table(FilePath, header = T, stringsAsFactors = F, sep = "\t")
    rownames(Positions) <- paste(Positions[, 1], Positions[, 2], sep = ":")
    return(Positions)
}

CalculateWTStats <- function(WTPositions, Name, output_dir = Name) {
    wtX <- WTPositions[, substr(colnames(WTPositions), start = 1, stop = 1) == "X"]
    wtEmbryos <- sub("X_", "", colnames(wtX))
    wtY <- WTPositions[, substr(colnames(WTPositions), start = 1, stop = 1) == "Y"]
    wtZ <- WTPositions[, substr(colnames(WTPositions), start = 1, stop = 1) == "Z"]
    colnames(wtX) <- wtEmbryos
    colnames(wtY) <- wtEmbryos
    colnames(wtZ) <- wtEmbryos

    wtXMeans <- rowMeans(wtX, na.rm = T)
    wtYMeans <- rowMeans(wtY, na.rm = T)
    wtZMeans <- rowMeans(wtZ, na.rm = T)

    wtXCellMeans <- tapply(wtXMeans, WTPositions[, 1], mean, na.rm = T)
    wtYCellMeans <- tapply(wtYMeans, WTPositions[, 1], mean, na.rm = T)
    wtZCellMeans <- tapply(wtZMeans, WTPositions[, 1], mean, na.rm = T)

    WTXdevs <- wtX - wtXMeans
    wtX_SD <- apply(WTXdevs, 1, sd, na.rm = T)
    WTYdevs <- wtY - wtYMeans
    wtY_SD <- apply(WTYdevs, 1, sd, na.rm = T)
    WTZdevs <- wtZ - wtZMeans
    wtZ_SD <- apply(WTZdevs, 1, sd, na.rm = T)

    WTDdevs <- sqrt(WTXdevs^2 + WTYdevs^2 + WTZdevs^2)
    WTDmeans <- apply(WTDdevs, 1, mean, na.rm = T)
    WTDsds <- apply(WTDdevs, 1, sd, na.rm = T)
    WTDcounts <- apply(WTDdevs, 1, function(X) {
        sum(!is.na(X))
    })

    # Plot WT Stats (ggplot2)
    pdf(file.path(output_dir, paste0(Name, "_WT_stats.pdf")), width = 10, height = 8)

    plot_df <- data.frame(
        MeanDev = WTDmeans,
        SDDev = WTDsds,
        Counts = WTDcounts
    )

    p1 <- ggplot(plot_df, aes(x = MeanDev, y = SDDev)) +
        geom_point(alpha = 0.5, size = 0.5) +
        theme_bw() +
        labs(x = "Mean WT Deviation (µm)", y = "SD of WT Deviations (µm)", title = "WT Deviation Variability")

    p2 <- ggplot(plot_df, aes(x = SDDev, y = Counts)) +
        geom_point(position = position_jitter(height = 0.2), alpha = 0.5, size = 0.5) +
        theme_bw() +
        labs(x = "SD of WT Deviations (µm)", y = "WT Observed Counts", title = "SD vs Observation Count")

    p3 <- ggplot(plot_df, aes(x = factor(Counts), y = SDDev)) +
        geom_boxplot(outlier.size = 0.5) +
        theme_bw() +
        labs(x = "WT Observed Counts", y = "SD of WT Deviations (µm)", title = "SD by Count")

    p4 <- ggplot(plot_df, aes(x = factor(Counts), y = MeanDev)) +
        geom_boxplot(outlier.size = 0.5) +
        theme_bw() +
        labs(x = "WT Observed Counts", y = "Mean WT Deviation (µm)", title = "Mean Deviation by Count")

    grid.arrange(p1, p2, p3, p4, nrow = 2, top = textGrob(paste("WT Reference Statistics -", Name), gp = gpar(fontsize = 14)))

    dev.off()

    # Filter SDs
    WTDsds <- (WTDsds[WTDcounts > 2])[names(WTDcounts)]

    return(list(
        wtXMeans = wtXMeans, wtYMeans = wtYMeans, wtZMeans = wtZMeans,
        wtX_SD = wtX_SD, wtY_SD = wtY_SD, wtZ_SD = wtZ_SD,
        wtXCellMeans = wtXCellMeans, wtYCellMeans = wtYCellMeans, wtZCellMeans = wtZCellMeans,
        WTDmeans = WTDmeans, WTDsds = WTDsds, WTDcounts = WTDcounts,
        wtEmbryos = wtEmbryos
    ))
}

CalculateNearestNeighbors <- function(MutantPositions, mutantX, mutantY, mutantZ, wtXMeans, wtYMeans, wtZMeans, Cells, embryo_name, Name) {
    # Optimized Vectorized NN Implementation

    message(paste("Processing embryo (NN):", embryo_name))

    mX <- mutantX[, embryo_name]
    mY <- mutantY[, embryo_name]
    mZ <- mutantZ[, embryo_name]
    names(mX) <- rownames(mutantX)
    names(mY) <- rownames(mutantY)
    names(mZ) <- rownames(mutantZ)

    # Pre-calculate Neighbor Lookups per Time
    Neighbors_Cache <- split(rownames(MutantPositions), MutantPositions[, 2])

    CellDeviationScores <- sapply(Cells, function(X) {
        myTimes <- MutantPositions[MutantPositions[, 1] == X, 2]
        if (length(myTimes) <= 1) {
            return(NA)
        }
        # sample 4 time points
        testTimes <- c(myTimes[0.1 * length(myTimes)], myTimes[0.35 * length(myTimes)], myTimes[0.7 * length(myTimes)], myTimes[0.9 * length(myTimes)])

        scores <- sapply(testTimes, function(Y) {
            key_X <- paste(X, Y, sep = ":")
            if (is.na(mX[key_X])) {
                return(NA)
            }

            Neighbors_IDs <- Neighbors_Cache[[as.character(Y)]]
            if (is.null(Neighbors_IDs)) {
                return(NA)
            }

            # Mutant Distances
            dX <- mX[Neighbors_IDs] - mX[key_X]
            dY <- mY[Neighbors_IDs] - mY[key_X]
            dZ <- mZ[Neighbors_IDs] - mZ[key_X]
            theseNeighborDistances <- sqrt(dX^2 + dY^2 + dZ^2)

            # WT Distances
            dX_WT <- wtXMeans[Neighbors_IDs] - wtXMeans[key_X]
            dY_WT <- wtYMeans[Neighbors_IDs] - wtYMeans[key_X]
            dZ_WT <- wtZMeans[Neighbors_IDs] - wtZMeans[key_X]
            theseWTNeighborDistances <- sqrt(dX_WT^2 + dY_WT^2 + dZ_WT^2)

            mean(tail(sort(abs(log2(theseNeighborDistances / theseWTNeighborDistances))), 10), na.rm = T)
        })
        return(median(scores, na.rm = TRUE))
    })
    return(CellDeviationScores)
}


OptimizeRotations <- function(i, Name, Times, mutantX, mutantY, mutantZ, wtXMeans, wtYMeans, wtZMeans, MutantPositions, WTPositions, PercentCellsRequired, wtXCellMeans, wtYCellMeans, wtZCellMeans, WTDmeans, WTDsds) {
    pdf(paste(Name, i, "PositionPlots.pdf", sep = "_"), width = 12, height = 10)

    Thetas <- NULL
    Distances <- NULL
    UnrotatedDistances <- NULL
    RotatedMutant <- NULL

    for (j in sort(Times)) {
        Pos <- data.frame(mutantX[MutantPositions[2] == j, i], mutantY[MutantPositions[2] == j, i], mutantZ[MutantPositions[2] == j, i], row.names = MutantPositions[MutantPositions[2] == j, 1])
        colnames(Pos) <- c("X", "Y", "Z")

        if (nrow(Pos) == 0) {
            Thetas <- rbind(Thetas, c(j, 0))
            rownames(Thetas) <- Thetas[, 1]
            Distances <- rbind(Distances, c(j, NA))
            rownames(Distances) <- Distances[, 1]
            UnrotatedDistances <- rbind(UnrotatedDistances, c(j, NA))
            rownames(UnrotatedDistances) <- UnrotatedDistances[, 1]
            next
        }
        WTMeanPos <- data.frame(wtXMeans[WTPositions[2] == j], wtYMeans[WTPositions[2] == j], wtZMeans[WTPositions[2] == j], row.names = WTPositions[WTPositions[2] == j, 1])
        colnames(WTMeanPos) <- c("X", "Y", "Z")
        WTMeanPos <- WTMeanPos[rownames(Pos), ]
        WTMutantDevs <- WTMeanPos - Pos
        WTMutantDists <- sqrt(rowSums(WTMutantDevs^2))
        MeanDist <- mean(WTMutantDists, na.rm = T)
        UnrotatedMeanDist <- MeanDist
        BestTheta <- 0
        if (j > 2) BestTheta <- Thetas[as.character(j - 1), 2]
        if (is.na(BestTheta)) BestTheta <- 0

        testRange <- pi * (-60:60) / 60
        bestPos <- Pos

        for (k in testRange) {
            thesePos <- rbind(rotate3d(as.matrix(Pos[, 1:3]), k, 1, 0, 0))
            colnames(thesePos) <- c("X", "Y", "Z")
            rownames(thesePos) <- rownames(Pos)
            theseWTMutantDevs <- WTMeanPos - thesePos
            thisMeanDist <- mean(sqrt(rowSums(theseWTMutantDevs^2)), na.rm = T)

            if (!is.na(thisMeanDist) && !is.na(MeanDist)) {
                if (thisMeanDist < MeanDist) {
                    BestTheta <- k
                    bestPos <- thesePos
                    MeanDist <- thisMeanDist
                }
            }
        }
        Thetas <- rbind(Thetas, c(j, BestTheta))
        rownames(Thetas) <- Thetas[, 1]
        Distances <- rbind(Distances, c(j, MeanDist))
        rownames(Distances) <- Distances[, 1]
        UnrotatedDistances <- rbind(UnrotatedDistances, c(j, UnrotatedMeanDist))
        rownames(UnrotatedDistances) <- UnrotatedDistances[, 1]

        rownames(bestPos) <- paste(rownames(bestPos), j, sep = ":")
        RotatedMutant <- rbind(RotatedMutant, bestPos)
    }

    # --- Plot 1: Rotation and Distance Metrics ---
    rot_metrics_df <- data.frame(
        Time = as.numeric(Thetas[, 1]),
        Theta = Thetas[, 2],
        Dist_Mean = Distances[, 2],
        Dist_Unrot = UnrotatedDistances[, 2]
    )

    p_theta <- ggplot(rot_metrics_df, aes(x = Time, y = Theta)) +
        geom_line() +
        geom_point() +
        theme_bw() +
        labs(x = "Time (min)", y = "Rotation Theta (rad)", title = "Optimal Rotation vs Time")

    p_dist <- ggplot(rot_metrics_df, aes(x = Time)) +
        geom_line(aes(y = Dist_Unrot, color = "Unrotated"), linetype = "dashed") +
        geom_line(aes(y = Dist_Mean, color = "Rotated"), size = 1) +
        scale_color_manual(values = c("Unrotated" = "black", "Rotated" = "red")) +
        theme_bw() +
        theme(legend.title = element_blank(), legend.position = "bottom") +
        labs(x = "Time (min)", y = "Mean Deviation (µm)", title = "Alignment Quality")


    RotatedMutant <- RotatedMutant[rownames(mutantX), ]
    rownames(RotatedMutant) <- rownames(mutantX)

    # --- Plot 2: Cell Counts ---
    WTCellCounts <- sapply(Times, function(X) {
        sum(!is.na(wtXMeans[WTPositions[, 2] == X]))
    })
    MutantCellCounts <- sapply(Times, function(X) {
        sum(!is.na(RotatedMutant[MutantPositions[, 2] == X, 1]))
    })

    counts_df <- data.frame(
        Time = Times,
        WT_Count = WTCellCounts,
        Mutant_Count = MutantCellCounts
    )

    p_counts <- ggplot(counts_df, aes(x = Time)) +
        geom_line(aes(y = WT_Count, color = "WT Reference")) +
        geom_line(aes(y = Mutant_Count, color = "This Embryo")) +
        geom_point(aes(y = WT_Count, color = "WT Reference"), size = 1) +
        geom_point(aes(y = Mutant_Count, color = "This Embryo"), size = 1) +
        scale_color_manual(values = c("WT Reference" = "black", "This Embryo" = "red")) +
        theme_bw() +
        theme(legend.title = element_blank(), legend.position = "bottom") +
        labs(x = "Time (min)", y = "Cell Count", title = "Cell Count Verification")

    # Filtering
    Keep <- MutantCellCounts / WTCellCounts > PercentCellsRequired
    Keep2 <- Keep[MutantPositions[rownames(RotatedMutant), 2]]
    names(Keep2) <- rownames(RotatedMutant)
    to_remove <- which(Keep2 == FALSE & !is.na(RotatedMutant[, 2]))
    if (length(to_remove) > 0) RotatedMutant[to_remove, ] <- NA

    # --- Plot 3: Global Position Correlations ---
    # Dataframe for correlation plots
    corr_df <- data.frame(
        WT_X = wtXMeans[rownames(mutantX)],
        WT_Y = wtYMeans[rownames(mutantX)],
        WT_Z = wtZMeans[rownames(mutantX)],
        Mut_X = mutantX[, i],
        Mut_Y = mutantY[, i],
        Mut_Z = mutantZ[, i],
        Rot_Y = RotatedMutant[rownames(mutantX), 2],
        Rot_Z = RotatedMutant[rownames(mutantX), 3]
    )

    p_x <- ggplot(corr_df, aes(x = WT_X, y = Mut_X)) +
        geom_point(alpha = 0.2, size = 0.5) +
        theme_bw() +
        labs(x = "WT Mean X", y = "Mutant X")

    p_y <- ggplot(corr_df) +
        geom_point(aes(x = WT_Y, y = Mut_Y, color = "Original"), alpha = 0.2, size = 0.5) +
        geom_point(aes(x = WT_Y, y = Rot_Y, color = "Rotated"), alpha = 0.2, size = 0.5) +
        scale_color_manual(values = c("Original" = "black", "Rotated" = "red")) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(x = "WT Mean Y", y = "Mutant Y")

    p_z <- ggplot(corr_df) +
        geom_point(aes(x = WT_Z, y = Mut_Z, color = "Original"), alpha = 0.2, size = 0.5) +
        geom_point(aes(x = WT_Z, y = Rot_Z, color = "Rotated"), alpha = 0.2, size = 0.5) +
        scale_color_manual(values = c("Original" = "black", "Rotated" = "red")) +
        theme_bw() +
        theme(legend.position = "none") +
        labs(x = "WT Mean Z", y = "Mutant Z")

    # Cell Means Correlation
    xCellMeans <- tapply(mutantX[, i], MutantPositions[, 1], mean, na.rm = T)
    yCellMeans <- tapply(mutantY[, i], MutantPositions[, 1], mean, na.rm = T)
    zCellMeans <- tapply(mutantZ[, i], MutantPositions[, 1], mean, na.rm = T)

    # Calculate R-squared for title
    r_sq <- tryCatch(
        {
            mod <- lm(wtXCellMeans[names(xCellMeans)] ~ xCellMeans)
            round(summary(mod)$adj.r.squared, 3)
        },
        error = function(e) NA
    )

    cell_corr_df <- data.frame(WT_X_Mean = wtXCellMeans[names(xCellMeans)], Mut_X_Mean = xCellMeans)
    p_cell_corr <- ggplot(cell_corr_df, aes(x = WT_X_Mean, y = Mut_X_Mean)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", se = FALSE, color = "blue", size = 0.5) +
        theme_bw() +
        labs(x = "WT Cell Mean X", y = "Mutant Cell Mean X", title = paste("X Correlation (R2 =", r_sq, ")"))

    # Calculate Deviations (Rotated)
    xDevs <- RotatedMutant[rownames(MutantPositions), 1] - wtXMeans[rownames(MutantPositions)]
    yDevs <- RotatedMutant[rownames(MutantPositions), 2] - wtYMeans[rownames(MutantPositions)]
    zDevs <- RotatedMutant[rownames(MutantPositions), 3] - wtZMeans[rownames(MutantPositions)]
    dDevs <- sqrt(xDevs^2 + yDevs^2 + zDevs^2)

    # Aggregates
    xDevMeans <- tapply(xDevs, MutantPositions[, 1], mean, na.rm = T)
    yDevMeans <- tapply(yDevs, MutantPositions[, 1], mean, na.rm = T)
    zDevMeans <- tapply(zDevs, MutantPositions[, 1], mean, na.rm = T)
    dDevMeans <- tapply(dDevs, MutantPositions[, 1], mean, na.rm = T)
    dDevMax <- tapply(dDevs, MutantPositions[, 1], max, na.rm = T)
    dZ <- (dDevs - WTDmeans[rownames(MutantPositions)]) / WTDsds[rownames(MutantPositions)]
    dZMeans <- tapply(dZ, MutantPositions[, 1], mean, na.rm = T)

    # --- Plot 4: Deviation Summaries ---
    dev_summ_df <- data.frame(x = xDevMeans, y = yDevMeans, z = zDevMeans, d = dDevMeans, dMax = dDevMax, dZ = dZMeans)

    p_xy <- ggplot(dev_summ_df, aes(x = x, y = y)) +
        geom_point(size = 0.5) +
        theme_bw() +
        labs(x = "Mean X Dev", y = "Mean Y Dev")
    p_xz <- ggplot(dev_summ_df, aes(x = x, y = z)) +
        geom_point(size = 0.5) +
        theme_bw() +
        labs(x = "Mean X Dev", y = "Mean Z Dev")
    p_yz <- ggplot(dev_summ_df, aes(x = y, y = z)) +
        geom_point(size = 0.5) +
        theme_bw() +
        labs(x = "Mean Y Dev", y = "Mean Z Dev")
    p_d_dmax <- ggplot(dev_summ_df, aes(x = d, y = dMax)) +
        geom_point(size = 0.5) +
        theme_bw() +
        labs(x = "Mean Total Dev", y = "Max Total Dev")
    p_d_dz <- ggplot(dev_summ_df, aes(x = d, y = dZ)) +
        geom_point(size = 0.5) +
        theme_bw() +
        labs(x = "Mean Total Dev", y = "Mean Z-Score")

    # Layout Pages
    grid.arrange(p_theta, p_dist, p_counts, p_cell_corr, nrow = 2, top = textGrob(paste("Alignment -", i), gp = gpar(fontsize = 14)))
    grid.arrange(p_x, p_y, p_z, p_xy, p_xz, p_yz, nrow = 2)
    grid.arrange(p_d_dmax, p_d_dz, nrow = 2, top = textGrob("Deviation Distributions", gp = gpar(fontsize = 14)))

    dev.off()

    return(list(
        RotatedMutant = RotatedMutant,
        xDevs = xDevs, yDevs = yDevs, zDevs = zDevs, dDevs = dDevs
    ))
}


# --- Main Function ---

AnalyzePositions <- function(Name, CalculateNeighbors = FALSE, Expression = NULL, PercentCellsRequired = 0.25, expCutoff = 500,
                             WT_Ref_File = NULL,
                             data_dir = ".",
                             output_dir = NULL,
                             wt_ref_dir = NULL,
                             embryo_metadata_file = NULL,
                             default_expression_file = NULL,
                             dim_reduction = TRUE) {
    if (is.null(output_dir)) output_dir <- file.path(data_dir, Name)
    if (is.null(wt_ref_dir)) wt_ref_dir <- file.path(data_dir, "Richard_et_al_plus_comma_WT")
    wt_prefix <- basename(wt_ref_dir)
    mutant_input_dir <- file.path(data_dir, Name)
    if (is.null(WT_Ref_File)) WT_Ref_File <- file.path(wt_ref_dir, paste0(wt_prefix, "positions.txt"))
    if (is.null(embryo_metadata_file)) embryo_metadata_file <- "embryo_metadata.csv"
    if (is.null(default_expression_file)) default_expression_file <- file.path(data_dir, "data/CA20120714_JIM136_L3.csv")

    message(paste("Analyzing Positions -", Name, sep = " "))

    # 1. Load WT Data
    WTPositions <- LoadPositionData(WT_Ref_File)

    # Determine Cells (union) early to ensure downstream consistency
    mutant_cell_names <- read.table(file.path(mutant_input_dir, paste0(Name, "positions.txt")), header = T, stringsAsFactors = F, sep = "\t", check.names = FALSE)[, 1]
    Cells <- sort(union(unique(WTPositions[, 1]), unique(mutant_cell_names)))

    # Handle Expression if not provided
    if (is.null(Expression)) {
        if (file.exists(default_expression_file)) {
            message("Loading default ceh-36 expression...")
            ceh36_exp <- read.csv(default_expression_file, row.names = 2, check.names = FALSE)["blot"]
            Expression <- sapply(Cells, function(X) GetPeak(ceh36_exp, X))
        } else {
            Expression <- rep(NA, length(Cells))
            names(Expression) <- Cells
        }
    }

    # Ensure output directory exists
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    mutant_cell_names <- read.table(file.path(mutant_input_dir, paste0(Name, "positions.txt")), header = T, stringsAsFactors = F, sep = "\t", check.names = FALSE)[, 1]
    Cells <- sort(union(unique(WTPositions[, 1]), unique(mutant_cell_names)))
    Times <- unique(WTPositions[, 2])

    WTStats <- CalculateWTStats(WTPositions, Name, output_dir = output_dir)
    wtXMeans <- WTStats$wtXMeans
    wtYMeans <- WTStats$wtYMeans
    wtZMeans <- WTStats$wtZMeans
    wtEmbryos <- WTStats$wtEmbryos
    WTDmeans <- WTStats$WTDmeans
    WTDsds <- WTStats$WTDsds
    WTDcounts <- WTStats$WTDcounts

    WTDir <- dirname(WT_Ref_File)

    # 2. Load Mutant Data
    MutantPositions <- read.table(file.path(mutant_input_dir, paste0(Name, "positions.txt")), header = T, stringsAsFactors = F, sep = "\t", check.names = FALSE)
    rownames(MutantPositions) <- paste(MutantPositions[, 1], MutantPositions[, 2], sep = ":")

    mutantX <- MutantPositions[, substr(colnames(MutantPositions), start = 1, stop = 1) == "X"]
    mutantEmbryos <- sub("X_", "", colnames(mutantX))
    mutantY <- MutantPositions[, substr(colnames(MutantPositions), start = 1, stop = 1) == "Y"]
    mutantZ <- MutantPositions[, substr(colnames(MutantPositions), start = 1, stop = 1) == "Z"]
    colnames(mutantX) <- mutantEmbryos
    colnames(mutantY) <- mutantEmbryos
    colnames(mutantZ) <- mutantEmbryos

    # 3. Calculate Neighbors (Optimized)
    if (CalculateNeighbors == TRUE) {
        message("Calculating neighbors (Optimized)")
        NN_distance_summaries <- data.frame(Cells)

        for (i in mutantEmbryos) {
            scores <- CalculateNearestNeighbors(MutantPositions, mutantX, mutantY, mutantZ, wtXMeans, wtYMeans, wtZMeans, Cells, i, Name)
            NN_distance_summaries <- cbind(NN_distance_summaries, scores)
            colnames(NN_distance_summaries)[ncol(NN_distance_summaries)] <- i
        }
        write.csv(NN_distance_summaries, file = file.path(output_dir, paste0(Name, "_NN_Scores.csv")))
    }


    # 4. Rotation and Deviation Analysis
    MutantDevs <- MutantPositions[, 1:2]
    MutantDevs <- cbind(MutantDevs, WTDmeans[rownames(MutantPositions)], WTDsds[rownames(MutantPositions)], WTDcounts[rownames(MutantPositions)])
    MutantCellDevMaxs <- data.frame(Cells)
    MutantCellDevMeans <- data.frame(Cells)
    colnames(MutantDevs) <- c("Cell", "Time", "WT D mean", "WT D SD", "WT D count")

    RotatedX <- NULL
    RotatedY <- NULL
    RotatedZ <- NULL

    message("Optimizing rotation and calculating deviations")
    for (i in mutantEmbryos) {
        print(i)

        pdf(file.path(output_dir, paste0(Name, "_", i, "_PositionPlots.pdf")))

        Thetas <- NULL
        Distances <- NULL
        UnrotatedDistances <- NULL
        RotatedMutant <- NULL

        # We MUST loop over all times present in MutantPositions to ensure RotatedMutant has all rows
        MutantTimes <- sort(unique(MutantPositions[, 2]))

        for (j in MutantTimes) {
            Pos <- data.frame(mutantX[MutantPositions[, 2] == j, i], mutantY[MutantPositions[, 2] == j, i], mutantZ[MutantPositions[, 2] == j, i], row.names = MutantPositions[MutantPositions[, 2] == j, 1])
            colnames(Pos) <- c("X", "Y", "Z")
            if (nrow(Pos) == 0) {
                Thetas <- rbind(Thetas, c(j, 0))
                rownames(Thetas) <- Thetas[, 1]
                Distances <- rbind(Distances, c(j, NA))
                rownames(Distances) <- Distances[, 1]
                UnrotatedDistances <- rbind(UnrotatedDistances, c(j, NA))
                rownames(UnrotatedDistances) <- UnrotatedDistances[, 1]
                next
            }

            # WT reference might be missing for this time point
            WT_Time_Indices <- WTPositions[, 2] == j
            if (sum(WT_Time_Indices) > 0) {
                WTMeanPos <- data.frame(wtXMeans[WT_Time_Indices], wtYMeans[WT_Time_Indices], wtZMeans[WT_Time_Indices], row.names = WTPositions[WT_Time_Indices, 1])
                colnames(WTMeanPos) <- c("X", "Y", "Z")
                common_cells <- intersect(rownames(Pos), rownames(WTMeanPos))

                if (length(common_cells) > 0) {
                    WTMeanPos_Sub <- WTMeanPos[common_cells, , drop = FALSE]
                    Pos_Sub <- Pos[common_cells, , drop = FALSE]

                    WTMutantDevs <- WTMeanPos_Sub - Pos_Sub
                    WTMutantDists <- sqrt(rowSums(WTMutantDevs^2))
                    MeanDist <- mean(WTMutantDists, na.rm = T)
                    UnrotatedMeanDist <- MeanDist
                    BestTheta <- 0
                    if (j > 2 && as.character(j - 1) %in% rownames(Thetas)) BestTheta <- Thetas[as.character(j - 1), 2]
                    if (is.na(BestTheta)) BestTheta <- 0

                    testRange <- pi * (-60:60) / 60

                    # Track best across multiple starting orientations
                    BestDist <- MeanDist
                    BestThetaFinal <- BestTheta
                    # Orientation codes: 0=Normal, 1=180-around-Y (X/Z flip), 2=180-around-Z (X/Y flip)
                    BestOrientation <- 0

                    orientations <- list(
                        list(code = 0, x_mul = 1, y_mul = 1, z_mul = 1),
                        list(code = 1, x_mul = -1, y_mul = 1, z_mul = -1), # 180 around Y
                        list(code = 2, x_mul = -1, y_mul = -1, z_mul = 1) # 180 around Z
                    )

                    for (ort in orientations) {
                        P_test <- Pos_Sub
                        P_test$X <- P_test$X * ort$x_mul
                        P_test$Y <- P_test$Y * ort$y_mul
                        P_test$Z <- P_test$Z * ort$z_mul

                        for (k in testRange) {
                            thesePos <- rbind(rotate3d(as.matrix(P_test[, 1:3]), k, 1, 0, 0))
                            colnames(thesePos) <- c("X", "Y", "Z")
                            rownames(thesePos) <- rownames(P_test)
                            theseWTMutantDevs <- WTMeanPos_Sub - thesePos
                            thisMeanDist <- mean(sqrt(rowSums(theseWTMutantDevs^2)), na.rm = T)
                            if (!is.na(thisMeanDist) && thisMeanDist < BestDist) {
                                BestThetaFinal <- k
                                BestDist <- thisMeanDist
                                BestOrientation <- ort$code
                            }
                        }
                    }

                    # Apply best found orientation and rotation
                    BestTheta <- BestThetaFinal
                    MeanDist <- BestDist
                    if (BestOrientation == 1) {
                        Pos$X <- -Pos$X
                        Pos$Z <- -Pos$Z
                    } else if (BestOrientation == 2) {
                        Pos$X <- -Pos$X
                        Pos$Y <- -Pos$Y
                    }

                    bestPos <- rbind(rotate3d(as.matrix(Pos[, 1:3]), BestTheta, 1, 0, 0))
                    colnames(bestPos) <- c("X", "Y", "Z")
                    rownames(bestPos) <- rownames(Pos)

                    Thetas <- rbind(Thetas, c(j, BestTheta))
                    rownames(Thetas) <- Thetas[, 1]
                    Distances <- rbind(Distances, c(j, MeanDist))
                    rownames(Distances) <- Distances[, 1]
                    UnrotatedDistances <- rbind(UnrotatedDistances, c(j, UnrotatedMeanDist))
                    rownames(UnrotatedDistances) <- UnrotatedDistances[, 1]
                } else {
                    BestTheta <- 0
                    if (j > 2 && as.character(j - 1) %in% rownames(Thetas)) BestTheta <- Thetas[as.character(j - 1), 2]
                    bestPos <- Pos
                    Thetas <- rbind(Thetas, c(j, 0))
                    rownames(Thetas) <- Thetas[, 1]
                    Distances <- rbind(Distances, c(j, NA))
                    rownames(Distances) <- Distances[, 1]
                    UnrotatedDistances <- rbind(UnrotatedDistances, c(j, NA))
                    rownames(UnrotatedDistances) <- UnrotatedDistances[, 1]
                }
            } else {
                # No WT data for this time point
                BestTheta <- 0
                if (j > 2 && as.character(j - 1) %in% rownames(Thetas)) BestTheta <- Thetas[as.character(j - 1), 2]
                bestPos <- Pos
                Thetas <- rbind(Thetas, c(j, 0))
                rownames(Thetas) <- Thetas[, 1]
                Distances <- rbind(Distances, c(j, NA))
                rownames(Distances) <- Distances[, 1]
                UnrotatedDistances <- rbind(UnrotatedDistances, c(j, NA))
                rownames(UnrotatedDistances) <- UnrotatedDistances[, 1]
            }

            rownames(bestPos) <- paste(rownames(bestPos), j, sep = ":")
            RotatedMutant <- rbind(RotatedMutant, bestPos)
        }

        plot(Thetas, main = i, xlab = "Time", ylab = "Least square distance Theta (rotation)", cex = 0.5)
        plot(UnrotatedDistances, main = i, xlab = "Time", ylab = "Mean deviations (Std black, rotated red)", cex = 0.5)
        points(Distances, col = 2, cex = 0.5)

        # Now RotatedMutant MUST have all rows from mutantX
        # But we use intersect just in case of weirdness
        common_rows <- intersect(rownames(mutantX), rownames(RotatedMutant))
        RotatedMutant <- RotatedMutant[common_rows, , drop = FALSE]

        # Ensure it has exactly the same rows as mutantX for downstream steps
        # If any are missing, fill with NA
        missing_rows <- setdiff(rownames(mutantX), rownames(RotatedMutant))
        if (length(missing_rows) > 0) {
            NA_matrix <- matrix(NA, nrow = length(missing_rows), ncol = 3, dimnames = list(missing_rows, c("X", "Y", "Z")))
            RotatedMutant <- rbind(RotatedMutant, NA_matrix)
        }
        RotatedMutant <- RotatedMutant[rownames(mutantX), , drop = FALSE]

        # Plot # of WT and mutant cells
        WTCellCounts <- sapply(unique(WTPositions[, 2]), function(X) {
            sum(!is.na(wtXMeans[WTPositions[, 2] == X]))
        })
        plot(unique(WTPositions[, 2]), WTCellCounts, xlab = "Time", ylab = "Cell counts - WT (blk), Mutant (red)")
        MutantCellCounts <- sapply(Times, function(X) {
            sum(!is.na(RotatedMutant[MutantPositions[, 2] == X, 1]))
        })
        points(Times, MutantCellCounts, col = 2)

        # Filter
        Keep <- MutantCellCounts / WTCellCounts > PercentCellsRequired
        Keep2 <- Keep[as.character(MutantPositions[rownames(RotatedMutant), 2])]
        names(Keep2) <- rownames(RotatedMutant)
        to_remove <- which(Keep2 == FALSE & !is.na(RotatedMutant[, 2]))
        if (length(to_remove) > 0) RotatedMutant[to_remove, ] <- NA

        plot(wtXMeans[rownames(mutantX)], mutantX[, i], cex = 0.1, xlab = "WT mean X", ylab = "Mutant X")
        plot(wtYMeans[rownames(mutantX)], mutantY[, i], cex = 0.1, xlab = "WT mean Y", ylab = "Mutant Y (red=rotated)")
        points(wtYMeans[rownames(mutantX)], RotatedMutant[rownames(mutantX), 2], col = 2, cex = 0.1)
        plot(wtZMeans[rownames(mutantX)], mutantZ[, i], cex = 0.1, xlab = "WT mean Z", ylab = "Mutant Z (red=rotated)")
        points(wtZMeans[rownames(mutantX)], RotatedMutant[rownames(mutantX), 3], col = 2, cex = 0.1)

        xCellMeans <- tapply(mutantX[, i], MutantPositions[, 1], mean, na.rm = T)
        yCellMeans <- tapply(mutantY[, i], MutantPositions[, 1], mean, na.rm = T)
        zCellMeans <- tapply(mutantZ[, i], MutantPositions[, 1], mean, na.rm = T)
        yRCellMeans <- tapply(RotatedMutant[rownames(MutantPositions), 2], MutantPositions[, 1], mean, na.rm = T)
        zRCellMeans <- tapply(RotatedMutant[rownames(MutantPositions), 3], MutantPositions[, 1], mean, na.rm = T)

        # Accumulate Rotated Vectors
        RotatedX <- cbind(RotatedX, RotatedMutant[rownames(MutantPositions), 1])
        RotatedY <- cbind(RotatedY, RotatedMutant[rownames(MutantPositions), 2])
        RotatedZ <- cbind(RotatedZ, RotatedMutant[rownames(MutantPositions), 3])

        print("Calculating deviations")
        xDevs <- RotatedMutant[rownames(MutantPositions), 1] - wtXMeans[rownames(MutantPositions)]
        xDevMeans <- tapply(xDevs, MutantPositions[, 1], mean, na.rm = T)
        xDevMax <- tapply(xDevs, MutantPositions[, 1], max, na.rm = T)
        xZ <- abs(xDevs) / WTStats$wtX_SD[rownames(MutantPositions)]
        xZMax <- tapply(xZ, MutantPositions[, 1], max, na.rm = T)

        yDevs <- RotatedMutant[rownames(MutantPositions), 2] - wtYMeans[rownames(MutantPositions)]
        yDevMeans <- tapply(yDevs, MutantPositions[, 1], mean, na.rm = T)
        yDevMax <- tapply(yDevs, MutantPositions[, 1], max, na.rm = T)
        yZ <- abs(yDevs) / WTStats$wtY_SD[rownames(MutantPositions)]
        yZMax <- tapply(yZ, MutantPositions[, 1], max, na.rm = T)

        zDevs <- RotatedMutant[rownames(MutantPositions), 3] - wtZMeans[rownames(MutantPositions)]
        zDevMeans <- tapply(zDevs, MutantPositions[, 1], mean, na.rm = T)
        zDevMax <- tapply(zDevs, MutantPositions[, 1], max, na.rm = T)
        zZ <- abs(zDevs) / WTStats$wtZ_SD[rownames(MutantPositions)]
        zZMax <- tapply(zZ, MutantPositions[, 1], max, na.rm = T)

        dDevs <- sqrt(xDevs^2 + yDevs^2 + zDevs^2)
        dDevMeans <- tapply(dDevs, MutantPositions[, 1], mean, na.rm = T)
        dDevMax <- tapply(dDevs, MutantPositions[, 1], max, na.rm = T)
        dZ <- (dDevs - WTDmeans[rownames(MutantPositions)]) / WTDsds[rownames(MutantPositions)]
        dZMeans <- tapply(dZ, MutantPositions[, 1], mean, na.rm = T)
        dZMax <- tapply(dZ, MutantPositions[, 1], max, na.rm = T)

        # Plot summaries
        plot(xDevMeans, yDevMeans, main = i)
        plot(xDevMeans, zDevMeans, main = i)
        plot(yDevMeans, zDevMeans, main = i)
        plot(dDevMeans, dDevMax, main = i)
        plot(dDevMeans, dZMeans, main = i)

        MutantDevs <- cbind(MutantDevs, dDevs)
        colnames(MutantDevs)[ncol(MutantDevs)] <- i

        test <- (dDevMax[Cells])
        names(test) <- Cells
        MutantCellDevMaxs <- cbind(MutantCellDevMaxs, test)
        test <- (dDevMeans[Cells])
        names(test) <- Cells
        MutantCellDevMeans <- cbind(MutantCellDevMeans, test)
        colnames(MutantCellDevMaxs)[ncol(MutantCellDevMaxs)] <- i
        colnames(MutantCellDevMeans)[ncol(MutantCellDevMeans)] <- i

        dev.off()
    }

    write.csv(MutantCellDevMaxs, file = file.path(output_dir, paste(Name, "CellMaxPositionDevs.csv", sep = "_")))
    write.csv(MutantCellDevMeans, file = file.path(output_dir, paste(Name, "CellMeanPositionDevs.csv", sep = "_")))
    write.csv(MutantDevs, file = file.path(output_dir, paste(Name, "PositionDevs.csv", sep = "_")))

    # 5. Defect Summaries
    print("Reading NN Scores and Defect Summary")
    MC_NN_scores <- read.csv(file.path(output_dir, paste(Name, "_NN_Scores.csv", sep = "")), row.names = 1, check.names = FALSE)
    MC_NN_scores[, 1] <- NULL

    WT_NN_scores <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_NN_Scores.csv")), row.names = 1, check.names = FALSE)
    WT_NN_scores[, 1] <- NULL

    MC_Devs <- read.csv(file.path(output_dir, paste(Name, "_CellMeanPositionDevs.csv", sep = "")), row.names = 1, check.names = FALSE)
    MC_Devs[, 1] <- NULL

    WT_Devs <- read.csv(file.path(wt_ref_dir, paste0(wt_prefix, "_CellMeanPositionDevs.csv")), row.names = 1, check.names = FALSE)
    WT_Devs[, 1] <- NULL

    CommonEmbryos <- intersect(colnames(MC_NN_scores), colnames(MC_Devs))
    if (length(CommonEmbryos) == 0) stop("No common embryos found between NN scores and Position Devs")
    MC_NN_scores <- MC_NN_scores[, CommonEmbryos, drop = FALSE]
    MC_Devs <- MC_Devs[, CommonEmbryos, drop = FALSE]
    mutantEmbryos <- CommonEmbryos
    colnames(WT_NN_scores) <- wtEmbryos
    colnames(WT_Devs) <- wtEmbryos

    pdf(file.path(output_dir, paste(Name, "positionDefects.pdf", sep = "_")), width = 8, height = 8)
    Cells <- rownames(MC_NN_scores)
    # P-Values
    NN_wilcox_p <- sapply(Cells, function(X) {
        try(return(wilcox.test(as.numeric(MC_NN_scores[X, ]), as.numeric(WT_NN_scores[X, ]))$p.value), silent = T)
        return(NA)
    })
    Dev_wilcox_p <- sapply(Cells, function(X) {
        try(return(wilcox.test(as.numeric(MC_Devs[X, ]), as.numeric(WT_Devs[X, ]))$p.value), silent = T)
        return(NA)
    })

    # Create Plotting DataFrame
    defects_df <- data.frame(
        Cell = Cells,
        LogP_NN = log10(NN_wilcox_p),
        LogP_Dev = log10(Dev_wilcox_p),
        Expression = Expression[Cells]
    )

    # Handle potentially missing expression or P-values
    defects_df <- defects_df %>% filter(!is.na(LogP_NN) & !is.na(LogP_Dev))

    p_pvals <- ggplot(defects_df, aes(x = LogP_Dev, y = LogP_NN)) +
        geom_point(aes(color = Expression), size = 2, alpha = 0.7) +
        scale_color_gradient(low = "blue", high = "red", na.value = "grey50") +
        geom_text_repel(aes(label = Cell), size = 2, max.overlaps = 20) +
        geom_vline(xintercept = -2, linetype = "dashed", color = "red") +
        geom_hline(yintercept = -2, linetype = "dashed", color = "red") +
        theme_bw() +
        labs(x = "Log10(P-value) Position Deviation", y = "Log10(P-value) Neighbor Deviation", title = "Significant Defects")

    # Z-Scores
    NN_Z <- (MC_NN_scores[Cells, ] - apply(WT_NN_scores[Cells, ], 1, mean, na.rm = T)) / apply(WT_NN_scores[Cells, ], 1, sd, na.rm = T)
    Dev_Z <- (MC_Devs[Cells, ] - apply(WT_Devs[Cells, ], 1, mean, na.rm = T)) / apply(WT_Devs[Cells, ], 1, sd, na.rm = T)
    max_nn_Z <- apply(NN_Z, 1, max, na.rm = T)
    max_dev_Z <- apply(Dev_Z, 1, max, na.rm = T)

    z_df <- data.frame(
        Cell = Cells,
        MaxZ_NN = max_nn_Z,
        MaxZ_Dev = max_dev_Z
    )

    p_z <- ggplot(z_df, aes(x = MaxZ_NN, y = MaxZ_Dev)) +
        geom_point(alpha = 0.6) +
        geom_text_repel(aes(label = Cell), size = 2, max.overlaps = 10) +
        theme_bw() +
        labs(x = "Max Neighbor Deviation Z-Score", y = "Max Position Deviation Z-Score", title = "Z-Score Severity")

    grid.arrange(p_pvals, p_z, nrow = 2)
    dev.off()

    # 6. Trajectory Plots
    rownames(RotatedX) <- rownames(MutantPositions)
    rownames(RotatedY) <- rownames(MutantPositions)
    rownames(RotatedZ) <- rownames(MutantPositions)
    colnames(RotatedX) <- mutantEmbryos
    colnames(RotatedY) <- mutantEmbryos
    colnames(RotatedZ) <- mutantEmbryos
    write.csv(RotatedX, file = file.path(output_dir, paste(Name, "_rotatedX.csv", sep = "")))
    write.csv(RotatedY, file = file.path(output_dir, paste(Name, "_rotatedY.csv", sep = "")))
    write.csv(RotatedZ, file = file.path(output_dir, paste(Name, "_rotatedZ.csv", sep = "")))

    # (Simplified plotting for final trajectory pdfs)


    # 7. UMAP / PCA Visualizations (Spatial Phenotyping) - skip if dim_reduction is FALSE
    if (!dim_reduction) {
        message("AnalyzePositions: dim_reduction=FALSE; skipping UMAP/PCA report.")
        return(MutantDevs)
    }

    metadata <- NULL
    if (file.exists(embryo_metadata_file)) metadata <- read.csv(embryo_metadata_file, stringsAsFactors = FALSE)

    # Metrics list
    vis_list <- list(
        list(
            title = "Mean Position Deviations",
            data = load_and_merge_data(
                file.path(output_dir, paste0(Name, "_CellMeanPositionDevs.csv")),
                file.path(wt_ref_dir, paste0(wt_prefix, "_CellMeanPositionDevs.csv"))
            )
        ),
        list(
            title = "Neighbor Deviation (NN) Scores",
            data = load_and_merge_data(
                file.path(output_dir, paste0(Name, "_NN_Scores.csv")),
                file.path(wt_ref_dir, paste0(wt_prefix, "_NN_Scores.csv"))
            )
        )
    )

    # Bonus: Terminal Position Analysis (Embryo Shape)
    # We reconstruct a [Cells x Embryos] matrix using the Euclidean distance
    # of each cell's terminal position from the WT mean terminal position.
    # Or more simply, just UMAP the terminal coordinates themselves (3 coords per cell).
    # To keep it simple and parallel, we'll use Mean Position Deviations as the primary 'Shape' metric.

    run_dimensionality_reduction_report(vis_list, metadata,
        output_pdf = file.path(output_dir, paste0(Name, "_Spatial_Visualizations_Grid.pdf")),
        report_title = "Spatial Phenotyping",
        Expression = Expression
    )

    # 7. UMAP / PCA Visualizations (Spatial Phenotyping)
    message("Generating Spatial UMAP/PCA Reports...")

    # Load metadata locally as this script might be run standalone
    metadata <- NULL
    if (file.exists(embryo_metadata_file)) metadata <- read.csv(embryo_metadata_file, stringsAsFactors = FALSE)

    # Define metrics to analyze
    # Using cell mean position deviations and NN scores
    vis_list <- list(
        list(
            title = "Mean Position Deviations",
            data = load_and_merge_data(
                file.path(output_dir, paste0(Name, "_CellMeanPositionDevs.csv")),
                file.path(wt_ref_dir, paste0(wt_prefix, "_CellMeanPositionDevs.csv"))
            )
        ),
        list(
            title = "Neighbor Deviation (NN) Scores",
            data = load_and_merge_data(
                file.path(output_dir, paste0(Name, "_NN_Scores.csv")),
                file.path(wt_ref_dir, paste0(wt_prefix, "_NN_Scores.csv"))
            )
        )
    )

    output_pdf_umap <- file.path(output_dir, paste0(Name, "_Spatial_Visualizations_Grid.pdf"))

    # Wrap in tryCatch to prevent analysis failure if UMAP fails (e.g. not enough data)
    tryCatch(
        {
            # Extract mutant IDs from one of the datasets (they are columns, merged with WT)
            # But load_and_merge_data merges them. The returned DF has headers.
            # We need the LIST of mutant IDs.
            # We can infer them as: All columns EXCEPT those in WT?
            # Or safely, we read the mutant file directly to get names.
            mut_for_names <- read.csv(file.path(output_dir, paste0(Name, "_CellMeanPositionDevs.csv")), check.names = FALSE, nrows = 1)
            # Assuming rownames is col 1, verify?
            # read.csv check.names=FALSE might keep X prefix if in file?
            # Usually file has Embryo names in header.
            # If row.names=1 used in load_and_merge, column names are embryos.
            # Let's read header.
            mut_ids <- colnames(mut_for_names) # includes "Cells" maybe?
            mut_ids <- mut_ids[mut_ids != "Cells" & mut_ids != ""]

            run_dimensionality_reduction_report(vis_list, metadata,
                output_pdf = output_pdf_umap,
                report_title = "Spatial Phenotyping",
                mutant_ids = mut_ids,
                mutant_label = Name,
                Expression = Expression
            )
        },
        error = function(e) {
            message("Error generating Spatial UMAPs:")
            print(e)
        }
    )

    return(MutantDevs)
}
