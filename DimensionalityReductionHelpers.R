# DimensionalityReductionHelpers.R
#
# UMAP / PCA visualization helpers used by AnalyzeDivTimes and AnalyzePositions
# to summarize per-embryo phenotypes. The two main entry points are:
#
#   .dim_red_prepare_meta(embryos, metadata, mutant_ids, mutant_label)
#       Build a tidy embryo-level metadata table from raw column names plus
#       optional user-provided embryo_metadata.csv. Handles "X20120714..."
#       style R-mangled names and the "WT (Reference)" / Mutant grouping.
#
#   run_dimensionality_reduction_report(vis_list, metadata, output_pdf,
#                                       report_title, mutant_ids,
#                                       mutant_label, Expression)
#       Take a list of metric data frames (cells x embryos), KNN-impute
#       missing values, run UMAP and PCA, and write a multi-page PDF plus
#       per-metric *_PCA_EigenCells.csv loadings tables to output_pdf's
#       directory.
#
# Inputs the report consumes are written by AnalyzeDivTimes (CCLengthNorm,
# DivTimeNorm, CCdev) and AnalyzePositions (CellMeanPositionDevs, NN_Scores)
# in the per-dataset output directory; the WT counterparts come from
# wt_ref_dir.
#
# Dependencies: umap, ggplot2, tidyr, dplyr, gridExtra, ggrepel, grid.

library(umap)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(grid)

# Internal helper: Prepare plotting metadata
.dim_red_prepare_meta <- function(embryos, metadata, mutant_ids = NULL, mutant_label = "Mutant") {
    if (is.null(embryos) || length(embryos) == 0) {
        return(NULL)
    }

    all_meta <- data.frame(Embryo_Raw = embryos, stringsAsFactors = FALSE)
    all_meta <- all_meta %>%
        mutate(Embryo_Clean = ifelse(grepl("^X[0-9]", .data$Embryo_Raw), substring(.data$Embryo_Raw, 2), .data$Embryo_Raw)) %>%
        mutate(Embryo_Match = gsub("\\.", "-", .data$Embryo_Clean))

    # Default Grouping
    all_meta$Group <- "WT (Reference)" # Default for WT control library

    # Apply Mutant Label if IDs provided
    if (!is.null(mutant_ids)) {
        # Match against raw or cleaned names. Usually raw names from colnames match raw names in stored list.
        # To be robust, we check both.
        all_meta$Group[all_meta$Embryo_Raw %in% mutant_ids] <- mutant_label
    }

    # Join with user metadata
    if (!is.null(metadata)) {
        all_meta <- left_join(all_meta, metadata, by = c("Embryo_Match" = "Embryo"))

        # Handle potential Group collision if metadata also has Group column
        if ("Group.y" %in% colnames(all_meta) && "Group.x" %in% colnames(all_meta)) {
            # Prefer metadata Group (Group.y) if present, fallback to our default (Group.x)
            all_meta$Group <- ifelse(!is.na(all_meta$Group.y), as.character(all_meta$Group.y), as.character(all_meta$Group.x))
            all_meta$Group.x <- NULL
            all_meta$Group.y <- NULL
        }
    } else {
        all_meta$Temperature <- NA
    }

    # Final Label Logic:
    # If Temperature is available (and valid 15/25), use it.
    # Otherwise use Group.

    # Check if Temperature exists (it might not if metadata passed but lacked it)
    if (!"Temperature" %in% colnames(all_meta)) all_meta$Temperature <- NA

    # If Temperature is NA, fill with Group (Mutant or WT)
    all_meta$Label <- as.character(all_meta$Temperature)

    # Logic: If Label is NA => Use Group.
    # If Label was "15" or "25", keep it.
    all_meta$Label[is.na(all_meta$Label)] <- all_meta$Group[is.na(all_meta$Label)]

    # Map "WT (20C)" from old metadata default?
    # If metadata had NAs, we want those to be Group (Mutant/WT) not forced to WT(20C) unless it IS WT.
    # My previous code logic: all_meta$Temperature[is.na(all_meta$Temperature)] <- "WT (20C)"
    # New logic: Only fallback to WT (20C) if we think it is one.
    # But with 'Group', we distinguish.

    # For backward compatibility with tbx-35 (where mutants have metadata "15"/"25" and WTs have NA -> WT 20C):
    # If Group is "Mutant" and Temperature is NA -> "Mutant" (Good)
    # If Group is "WT" and Temperature is NA -> "WT (20C)" (or just keep WT?)
    # Let's standardize WT label.
    all_meta$Label[all_meta$Label == "WT (Reference)"] <- "WT (20C)"

    # Ensure Label is Factor
    # Priority levels: 15, 25, Mutant, WT (20C)
    all_meta$Label <- factor(all_meta$Label, levels = unique(c("15", "25", mutant_label, "WT (20C)", all_meta$Label)))

    return(all_meta %>% distinct(Embryo_Raw, .keep_all = TRUE))
}

# Internal helper: Get legend grob
.dim_red_get_legend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) == 0) {
        return(NULL)
    }
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# --- New Imputation Logic ---
impute_knn_pearson <- function(df, k = 5) {
    # df: embryos as columns, cells as rows
    # We want to find similar embryos.

    mat <- as.matrix(df)

    # Check for empty columns (embryos with NO data)
    valid_cols <- colSums(!is.na(mat)) > 0
    mat <- mat[, valid_cols, drop = FALSE]

    # Calculate Correlation Matrix (pairwise.complete.obs)
    # We use this to find neighbors
    corr_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

    # Diagnostic: Check for NAs in correlation matrix (if no shared observations)
    corr_mat[is.na(corr_mat)] <- 0 # Treat as no correlation

    mat_imputed <- mat

    # Loop through each embryo (column) that has missing values
    na_indices <- which(is.na(mat), arr.ind = TRUE)

    if (nrow(na_indices) > 0) {
        # Optimization: Iterate by column (embryo)
        cols_with_na <- unique(na_indices[, 2])

        for (j in cols_with_na) {
            # Find K nearest neighbors based on correlation
            # Get correlations for this embryo
            corrs <- corr_mat[, j]

            # Exclude self
            corrs[j] <- -Inf

            # Sort decreasing
            top_k_indices <- order(corrs, decreasing = TRUE)[1:min(k, length(corrs) - 1)]

            # Identify missing rows (cells) for this embryo
            missing_rows <- which(is.na(mat[, j]))

            for (i in missing_rows) {
                # Get values of these neighbors for this cell
                neighbor_vals <- mat[i, top_k_indices]

                # Calculate mean (ignoring NAs in neighbors)
                replacement <- mean(neighbor_vals, na.rm = TRUE)

                # Fallback: If neighbors also missing, use global row mean
                if (is.nan(replacement) || is.na(replacement)) {
                    replacement <- mean(mat[i, ], na.rm = TRUE)
                }

                # Second Fallback: If row is entirely empty, use 0
                if (is.nan(replacement) || is.na(replacement)) {
                    replacement <- 0
                }

                mat_imputed[i, j] <- replacement
            }
        }
    }

    return(as.data.frame(mat_imputed, check.names = FALSE))
}

# Main Visualization Function
generate_dim_red_plots <- function(df, metadata, title, custom_colors = NULL, max_overlaps = 10) {
    # df: embryos as columns, cells as rows
    df_t <- t(df)
    embryos <- rownames(df_t)
    all_meta <- .dim_red_prepare_meta(embryos, metadata, mutant_ids = attr(df, "mutant_ids"), mutant_label = attr(df, "mutant_label"))

    # Filter cells (at least 20% observed and not constant)
    # df_t is Embryos x Cells
    valid_cells <- apply(df_t, 2, function(x) sum(!is.na(x)) > (nrow(df_t) * 0.2) && length(unique(x[!is.na(x)])) > 1)
    df_filtered <- df_t[, valid_cells, drop = FALSE]

    if (ncol(df_filtered) < 2) {
        message(paste("Skipping", title, "- too few valid cells"))
        return(NULL)
    }

    # Imputation using KNN (Pearson)
    # impute_knn_pearson expects Embryos as Columns (to find neighbors among them)
    # df_filtered is Embryos x Cells (Rows x Cols)
    # So we transpose it to Cells x Embryos
    df_for_imp <- t(df_filtered)

    message(paste("Imputing", title, "using KNN (Pearson)..."))
    df_imputed_t <- impute_knn_pearson(df_for_imp, k = 5)

    # Transpose back to Embryos x Cells for UMAP
    df_imputed <- t(df_imputed_t)

    # Handle zero-variance columns after imputation to prevent prcomp scale. = TRUE error
    valid_imp_cols <- apply(df_imputed, 2, function(x) var(x, na.rm = TRUE) > 1e-8)
    valid_imp_cols[is.na(valid_imp_cols)] <- FALSE
    df_imputed <- df_imputed[, valid_imp_cols, drop = FALSE]

    if (ncol(df_imputed) < 2) {
        message(paste("Skipping", title, "- too few valid cells after imputation"))
        return(NULL)
    }

    # 1. Run UMAP
    custom.config <- umap.defaults
    custom.config$n_neighbors <- min(15, nrow(df_imputed) - 1)
    um_res <- try(umap(df_imputed, config = custom.config), silent = TRUE)

    umap_df <- NULL
    if (!is(um_res, "try-error")) {
        umap_df <- as.data.frame(um_res$layout)
        colnames(umap_df) <- c("Dim1", "Dim2")
        umap_df$Embryo_Raw <- rownames(df_imputed)
        umap_df <- left_join(umap_df, all_meta, by = "Embryo_Raw")
    }

    # 2. Run PCA
    pca_res <- prcomp(df_imputed, scale. = TRUE)
    pca_df <- as.data.frame(pca_res$x[, 1:min(3, ncol(pca_res$x))])
    pca_df$Embryo_Raw <- rownames(df_imputed)
    pca_df <- left_join(pca_df, all_meta, by = "Embryo_Raw")

    plot_func <- function(d, x_col, y_col, x_lab, y_lab, method_name) {
        p <- ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]], color = Label)) +
            geom_point(size = 2.5, alpha = 0.7) +
            geom_text_repel(aes(label = Embryo_Match), size = 2, show.legend = FALSE, max.overlaps = max_overlaps) +
            theme_bw() +
            labs(title = paste(method_name, "-", title), x = x_lab, y = y_lab, color = "Group") +
            theme(plot.title = element_text(size = 10))

        # Dynamic Colors: Apply manual scale ONLY if standard labels are present
        current_labels <- unique(as.character(d$Label))
        # Known colors + Custom
        color_map <- c(
            "15" = "blue", "25" = "red", "WT (20C)" = "black", "WT" = "black", "WT (Reference)" = "black",
            "ceh-51" = "darkorange",
            "tbx-35 (15C)" = "blue",
            "tbx-35 (25C)" = "red",
            "Mutant" = "purple"
        )

        if (!is.null(custom_colors)) {
            # Merge custom colors, overwriting defaults if conflict
            color_map <- c(color_map, custom_colors)
            # Unique names, keeping last (custom)
            color_map <- color_map[unique(names(color_map))]
        }

        # Apply logic to pick needed colors
        if (any(current_labels %in% names(color_map))) {
            my_colors <- color_map
            others <- setdiff(current_labels, names(color_map))
            if (length(others) > 0) {
                extra_cols <- c("darkorange", "purple", "brown", "cyan4")
                new_map <- structure(extra_cols[1:length(others)], names = others)
                my_colors <- c(my_colors, new_map)
            }
            p <- p + scale_color_manual(values = my_colors)
        }

        return(p)
    }

    p_umap <- NULL
    if (!is.null(umap_df)) p_umap <- plot_func(umap_df, "Dim1", "Dim2", "UMAP1", "UMAP2", "UMAP")

    p_pca12 <- plot_func(pca_df, "PC1", "PC2", "PC1", "PC2", "PCA")
    p_pca13 <- plot_func(pca_df, "PC1", "PC3", "PC1", "PC3", "PCA")
    p_pca23 <- NULL
    if ("PC3" %in% colnames(pca_df)) p_pca23 <- plot_func(pca_df, "PC2", "PC3", "PC2", "PC3", "PCA")

    # Return results including data for export
    return(list(
        plots = list(p_umap, p_pca12, p_pca13, p_pca23),
        loadings = pca_res$rotation,
        scores = pca_df
    ))
}

# Logic to merge mutant run files with a WT control file (filtering Sulston)
load_and_merge_data <- function(mutant_path, wt_path) {
    if (!file.exists(mutant_path)) {
        message(paste("File missing:", mutant_path))
        return(NULL)
    }
    if (!file.exists(wt_path)) {
        message(paste("File missing:", wt_path))
        return(NULL)
    }

    # Helper to read based on extension
    read_auto <- function(p) {
        if (grepl("\\.csv$", p, ignore.case = TRUE)) {
            d <- read.csv(p, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            # Remove 'Cells' column if present (common in position files)
            if ("Cells" %in% colnames(d)) d$Cells <- NULL
            return(d)
        } else {
            return(read.table(p, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
        }
    }

    mut <- read_auto(mutant_path)
    wt <- read_auto(wt_path)

    # Exclude 20081128_sulston
    wt <- wt[, !colnames(wt) %in% c("20081128_sulston", "X20081128_sulston"), drop = FALSE]
    mut <- mut[, !colnames(mut) %in% c("20081128_sulston", "X20081128_sulston"), drop = FALSE]

    all_cells <- unique(c(rownames(wt), rownames(mut)))
    cols <- unique(c(colnames(wt), colnames(mut)))

    if (length(cols) == 0) {
        return(NULL)
    }

    combined <- as.data.frame(matrix(NA, nrow = length(all_cells), ncol = length(cols)))
    rownames(combined) <- all_cells
    colnames(combined) <- cols
    for (cn in colnames(wt)) combined[rownames(wt), cn] <- wt[, cn]
    for (cn in colnames(mut)) combined[rownames(mut), cn] <- mut[, cn]
    return(combined)
}

# High-level convenience function
run_dimensionality_reduction_report <- function(data_list, metadata, output_pdf, report_title = "Analysis", mutant_ids = NULL, mutant_label = "Mutant", Expression = NULL, custom_colors = NULL, max_overlaps = 10) {
    # data_list: list(list(data=df, title="Metric 1"), ...)

    for (i in seq_along(data_list)) {
        if (!is.null(data_list[[i]]$data)) {
            attr(data_list[[i]]$data, "mutant_ids") <- mutant_ids
            attr(data_list[[i]]$data, "mutant_label") <- mutant_label
        }
    }

    message(paste("Creating PDF report:", output_pdf))
    pdf(output_pdf, width = 14, height = 14)

    output_dir <- dirname(output_pdf)
    output_base <- sub("\\.pdf$", "", basename(output_pdf))

    for (item in data_list) {
        if (is.null(item$data)) next

        # generate_dim_red_plots now returns list(plots=..., loadings=..., scores=...)
        res <- generate_dim_red_plots(item$data, metadata, item$title, custom_colors, max_overlaps)
        if (is.null(res) || is.null(res$plots)) next

        plots <- res$plots

        # ------------------------
        # Export Loadings CSV
        # ------------------------
        sanitized_title <- gsub(" ", "_", item$title)
        export_loadings <- as.data.frame(res$loadings)
        export_loadings$Cell <- rownames(export_loadings)

        # 1. Lineage metadata
        get_lineage <- function(x) {
            # Explicitly categorize specific early cells as "Other"
            if (x %in% c("P1", "P2", "P3", "P4", "Z2", "Z3", "EMS")) {
                return("Other")
            }
            if (grepl("^AB", x)) {
                return("AB")
            }
            if (grepl("^MS", x)) {
                return("MS")
            }
            if (grepl("^E", x)) {
                return("E")
            }
            if (grepl("^C", x)) {
                return("C")
            }
            if (grepl("^D", x)) {
                return("D")
            }
            if (grepl("^P", x)) {
                return("P")
            }
            if (grepl("^Z", x)) {
                return("Z")
            }
            return("Other")
        }
        export_loadings$Lineage <- sapply(export_loadings$Cell, get_lineage)

        # 2. Expression (if provided)
        if (!is.null(Expression)) {
            # Lookup from named vector
            export_loadings$Expression <- Expression[export_loadings$Cell]
        }

        # 3. WT Birth Time
        # Assume file is in current or parent dir?
        # Standard location: "SupplementalTable2_DivisionTimes.txt" in root of workspace?
        # We try a few locations.
        wt_locs <- c("SupplementalTable2_DivisionTimes.txt", "../SupplementalTable2_DivisionTimes.txt")
        wt_file_found <- NULL
        for (f in wt_locs) if (file.exists(f)) wt_file_found <- f

        if (!is.null(wt_file_found)) {
            wt_data <- try(read.delim(wt_file_found, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
            if (!inherits(wt_data, "try-error")) {
                # Find Birth col
                birth_col <- grep("Birth", colnames(wt_data), value = TRUE, ignore.case = TRUE)[1]
                if (!is.na(birth_col) && "Cell" %in% colnames(wt_data)) {
                    birth_map <- wt_data[, c("Cell", birth_col)]
                    colnames(birth_map) <- c("Cell", "WT_Birth_Time")
                    export_loadings <- left_join(export_loadings, birth_map, by = "Cell")
                }
            }
        }

        # Reorder
        first_cols <- c("Cell", "Lineage")
        if ("WT_Birth_Time" %in% colnames(export_loadings)) first_cols <- c(first_cols, "WT_Birth_Time")
        if ("Expression" %in% colnames(export_loadings)) first_cols <- c(first_cols, "Expression")

        other_cols <- setdiff(colnames(export_loadings), first_cols)
        export_loadings <- export_loadings[, c(first_cols, other_cols)]

        csv_name <- paste0(output_base, "_", sanitized_title, "_PCA_EigenCells.csv")
        csv_path <- file.path(output_dir, csv_name)
        write.csv(export_loadings, csv_path, row.names = FALSE)
        message(paste("Exported loadings:", csv_path))

        # ------------------------
        # Plotting
        # ------------------------
        legend <- .dim_red_get_legend(plots[[2]] + theme(legend.position = "bottom"))
        plots_clean <- lapply(plots, function(p) if (!is.null(p)) p + theme(legend.position = "none") else NULL)

        grid.arrange(
            plots_clean[[1]], plots_clean[[2]],
            plots_clean[[3]], plots_clean[[4]],
            ncol = 2, nrow = 2,
            top = textGrob(paste(report_title, "-", item$title), gp = gpar(fontsize = 16, font = 2)),
            bottom = legend
        )
    }

    dev.off()
    message("Done.")
}
