#!/usr/bin/env Rscript

# PlotComparisonBoxplots.R
#
# Comparative boxplots of Cell Cycle deviations, grouped by Major Lineage,
# Depth, and per-Cell. Adds Wilcoxon p-values (Bonferroni-corrected) for
# every non-WT group vs the WT reference.
#
# Two ways to invoke:
#
#   1. As a function (preferred, called by run_pipeline.R):
#        source("PlotComparisonBoxplots.R")
#        PlotComparisonBoxplots(name, data_dir = ..., output_dir = ...,
#                               wt_ref_dir = ..., embryo_metadata_file = ...)
#
#   2. As a standalone Rscript (legacy):
#        Rscript PlotComparisonBoxplots.R <analysis_dir>
#      where <analysis_dir> is both the input/output directory and the dataset
#      name (current-working-directory-relative). All other paths fall back
#      to the legacy "Richard_et_al_plus_comma_WT/" / "embryo_metadata.csv".

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(grid)

# Source visualization helpers for metadata logic
if (file.exists("DimensionalityReductionHelpers.R")) {
    source("DimensionalityReductionHelpers.R")
} else {
    stop("DimensionalityReductionHelpers.R not found in current directory.")
}


# ------------------------------------------------------------------
# HELPER FUNCTIONS (top-level so they can be tested independently)
# ------------------------------------------------------------------

.pcb_get_major_lineage <- function(cell_name) {
    if (grepl("^AB", cell_name)) return("AB")
    if (grepl("^MS", cell_name)) return("MS")
    if (grepl("^E", cell_name)) return("E")
    if (grepl("^C", cell_name)) return("C")
    if (grepl("^D", cell_name)) return("D")
    if (grepl("^P", cell_name)) return("P")
    return("Other")
}

.pcb_has_data_for_test <- function(vals1, vals2) {
    length(na.omit(vals1)) >= 3 && length(na.omit(vals2)) >= 3
}

.pcb_calculate_all_stats <- function(data, group_cols) {
    unique_groups_in_data <- unique(as.character(data$Group))
    ref_group <- "WT (Reference)"

    comparisons <- list()
    mutant_groups <- setdiff(unique_groups_in_data, ref_group)

    if (ref_group %in% unique_groups_in_data) {
        for (m in mutant_groups) {
            comparisons[[length(comparisons) + 1]] <- c(m, ref_group)
        }
    }

    if ("15" %in% mutant_groups && "25" %in% mutant_groups) {
        comparisons[[length(comparisons) + 1]] <- c("15", "25")
    }

    if (length(comparisons) == 0) return(NULL)

    data$GroupID <- apply(data[, group_cols, drop = FALSE], 1, paste, collapse = ":")
    unique_ids <- unique(data$GroupID)

    results <- list()

    for (grp in unique_ids) {
        sub_data <- data[data$GroupID == grp, ]

        for (pair in comparisons) {
            g1 <- pair[1]
            g2 <- pair[2]

            p_val <- NA
            if (all(pair %in% unique(sub_data$Group))) {
                vals1 <- sub_data$CC_Deviation[sub_data$Group == g1]
                vals2 <- sub_data$CC_Deviation[sub_data$Group == g2]

                if (.pcb_has_data_for_test(vals1, vals2)) {
                    res <- try(wilcox.test(vals1, vals2), silent = TRUE)
                    if (!inherits(res, "try-error")) {
                        p_val <- res$p.value
                    }
                }
            }

            comparison_name <- paste(g1, "vs", g2)

            results[[length(results) + 1]] <- data.frame(
                GroupID = grp,
                Comparison = comparison_name,
                Group1 = g1,
                Group2 = g2,
                P_raw = p_val,
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(results) == 0) return(NULL)

    results_df <- do.call(rbind, results)

    results_df <- results_df %>%
        group_by(Comparison) %>%
        mutate(
            N_tests = sum(!is.na(P_raw)),
            P_adj = p.adjust(P_raw, method = "bonferroni")
        ) %>%
        ungroup()

    return(results_df)
}

.pcb_get_sig_brackets <- function(data, x_var, stats_df, group_cols) {
    if (is.null(stats_df)) return(data.frame())

    brackets <- data.frame()

    if (!is.factor(data[[x_var]])) data[[x_var]] <- factor(data[[x_var]])
    x_levels <- levels(data[[x_var]])
    group_levels <- levels(data$Group)

    y_max_global <- max(data$CC_Deviation, na.rm = TRUE)
    y_min_global <- min(data$CC_Deviation, na.rm = TRUE)
    step <- (y_max_global - y_min_global) * 0.08
    if (is.na(step) || step == 0) step <- 1

    for (lvl in x_levels) {
        sub_data <- data[data[[x_var]] == lvl, ]
        if (nrow(sub_data) == 0) next

        grp_id <- unique(sub_data$GroupID)
        if (length(grp_id) > 1) next

        my_stats <- stats_df %>% filter(GroupID == grp_id, !is.na(P_adj), P_adj < 0.05)

        if (nrow(my_stats) == 0) next

        current_y <- max(sub_data$CC_Deviation, na.rm = TRUE)

        for (k in 1:nrow(my_stats)) {
            st <- my_stats[k, ]
            g1 <- st$Group1
            g2 <- st$Group2

            if (all(c(g1, g2) %in% unique(sub_data$Group))) {

                p_val <- st$P_adj
                if (p_val < 0.001) {
                    label <- sprintf("p=%.1e", p_val)
                } else {
                    label <- sprintf("p=%.3f", p_val)
                }

                current_y <- current_y + step

                x_int <- match(lvl, x_levels)

                get_offset <- function(t, all_grps) {
                    n_grps <- length(all_grps)
                    idx <- match(t, all_grps)
                    if (n_grps == 1) return(0)
                    seq(-0.3, 0.3, length.out = n_grps)[idx]
                }

                x1 <- x_int + get_offset(g1, group_levels)
                x2 <- x_int + get_offset(g2, group_levels)

                brackets <- rbind(brackets, data.frame(
                    x = x1, xend = x2,
                    y = current_y,
                    label = label
                ))
            }
        }
    }
    return(brackets)
}

.pcb_plot_lineage_boxplot <- function(data, lineage_name, title_suffix, stats_df, x_col = "Depth", group_cols) {

    brackets <- .pcb_get_sig_brackets(data, x_col, stats_df, group_cols)

    known_colors <- c("15" = "blue", "25" = "red", "WT (20C)" = "grey50", "WT (Reference)" = "grey30")
    unique_grps <- levels(data$Group)

    final_colors <- known_colors[intersect(names(known_colors), unique_grps)]

    others <- setdiff(unique_grps, names(known_colors))
    if (length(others) > 0) {
        new_cols <- scales::hue_pal()(length(others))
        names(new_cols) <- others
        final_colors <- c(final_colors, new_cols)
    }

    p <- ggplot(data, aes(x = factor(.data[[x_col]]), y = CC_Deviation, fill = Group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.5, position = position_dodge(width = 0.8)) +
        geom_point(aes(color = Group),
                   position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                   size = 0.8, alpha = 0.6) +
        scale_fill_manual(values = final_colors) +
        scale_color_manual(values = final_colors) +
        theme_bw() +
        labs(
            title = paste("Lineage:", lineage_name, "-", title_suffix),
            x = ifelse(x_col == "Depth", "Lineage Depth (Generation)", "Cell Identity"),
            y = "Deviation (min)",
            caption = "Significance: Wilcoxon Rank Sum with Bonferroni Correction"
        ) +
        theme(legend.position = "bottom")

    if (nrow(brackets) > 0) {
        max_y_data <- max(data$CC_Deviation, na.rm = TRUE)
        max_y_bracket <- max(brackets$y, na.rm = TRUE)
        upper_lim <- max(max_y_data, max_y_bracket) + (max(data$CC_Deviation, na.rm = TRUE) - min(data$CC_Deviation, na.rm = TRUE)) * 0.1

        p <- p +
            geom_segment(data = brackets, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
            geom_text(data = brackets, aes(x = (x + xend) / 2, y = y, label = label), inherit.aes = FALSE, vjust = 0, size = 3) +
            coord_cartesian(ylim = c(NA, upper_lim))
    }

    return(p)
}


# ------------------------------------------------------------------
# MAIN ENTRY POINT
# ------------------------------------------------------------------

PlotComparisonBoxplots <- function(name,
                                   data_dir = ".",
                                   output_dir = NULL,
                                   wt_ref_dir = NULL,
                                   embryo_metadata_file = NULL) {
    if (is.null(output_dir)) output_dir <- file.path(data_dir, name)
    if (is.null(wt_ref_dir)) wt_ref_dir <- file.path(data_dir, "Richard_et_al_plus_comma_WT")
    wt_prefix <- basename(wt_ref_dir)
    if (is.null(embryo_metadata_file)) embryo_metadata_file <- "embryo_metadata.csv"

    message(paste("Generating Consolidated Comparison Boxplots for:", name))

    # 1. Load Data
    cc_devs_file <- file.path(output_dir, paste0(name, "_ccDevs.csv"))
    if (!file.exists(cc_devs_file)) {
        stop(paste("Data file not found:", cc_devs_file))
    }

    df_cc <- read.csv(cc_devs_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (any(colnames(df_cc) == "")) {
        colnames(df_cc)[colnames(df_cc) == ""] <- "RowID"
    }
    colnames(df_cc) <- gsub(" ", "_", colnames(df_cc))

    # 2. Metadata
    embryo_list <- unique(df_cc$Embryo)

    meta_file_specific <- file.path(output_dir, paste0(name, "_embryo_metadata.csv"))
    meta_file_generic <- embryo_metadata_file

    metadata_raw <- NULL
    if (file.exists(meta_file_specific)) {
        message(paste("Loading metadata from:", meta_file_specific))
        metadata_raw <- read.csv(meta_file_specific, stringsAsFactors = FALSE)
    } else if (file.exists(meta_file_generic)) {
        message(paste("Loading metadata from:", meta_file_generic))
        metadata_raw <- read.csv(meta_file_generic, stringsAsFactors = FALSE)
    }

    if (!is.null(metadata_raw)) {
        if ("Temperature" %in% colnames(metadata_raw)) {
            metadata_raw$Group <- as.character(metadata_raw$Temperature)
            message("Grouping by 'Temperature' column found in metadata.")
        } else if ("Group" %in% colnames(metadata_raw)) {
            metadata_raw$Group <- as.character(metadata_raw$Group)
            message("Grouping by 'Group' column found in metadata.")
        } else {
            metadata_raw$Group <- "Mutant"
            message("No standard grouping column (Temperature/Group) found. Treating all as 'Mutant'.")
        }
    }

    meta_clean <- .dim_red_prepare_meta(embryo_list, metadata_raw, mutant_ids = embryo_list, mutant_label = "Mutant")

    if ("Temperature" %in% colnames(meta_clean)) {
        if (is.null(metadata_raw)) {
            meta_clean$Group <- "Mutant"
        } else {
            if (!"Group" %in% colnames(meta_clean)) meta_clean$Group <- meta_clean$Temperature
        }
    }

    df_merged <- left_join(df_cc, meta_clean, by = c("Embryo" = "Embryo_Raw"))

    # 2b. Load WT reference data
    wt_file <- file.path(wt_ref_dir, paste0(wt_prefix, "_ccDevs.csv"))

    if (file.exists(wt_file)) {
        message("Loading WT Reference data...")
        df_wt <- read.csv(wt_file, stringsAsFactors = FALSE, check.names = FALSE)
        if (any(colnames(df_wt) == "")) colnames(df_wt)[colnames(df_wt) == ""] <- "RowID"
        colnames(df_wt) <- gsub(" ", "_", colnames(df_wt))

        df_wt$Group <- "WT (Reference)"
        if (!"Embryo" %in% colnames(df_wt)) df_wt$Embryo <- df_wt$RowID

        cols_to_keep <- c("Cell", "Embryo", "CC_Deviation", "Group")

        if (!"Group" %in% colnames(df_merged)) df_merged$Group <- "Mutant"

        df_merged$Group <- as.character(df_merged$Group)
        df_wt$Group <- as.character(df_wt$Group)

        df_merged <- bind_rows(
            df_merged[, cols_to_keep[cols_to_keep %in% colnames(df_merged)]],
            df_wt[, cols_to_keep[cols_to_keep %in% colnames(df_wt)]]
        )

        message("Merged Data Group distribution:")
        print(table(df_merged$Group, useNA = "always"))
    }

    # 3. Feature engineering
    df_merged <- df_merged %>%
        mutate(
            MajorLineage = sapply(Cell, .pcb_get_major_lineage),
            Depth = str_count(Cell, "[a-z]")
        )

    target_lineages <- c("AB", "MS", "C", "D", "E")
    df_plot <- df_merged %>% filter(MajorLineage %in% target_lineages)

    all_groups <- unique(df_plot$Group)
    wt_ref <- "WT (Reference)"
    others <- sort(setdiff(all_groups, wt_ref))

    if ("WT (20C)" %in% others) {
        others <- c(setdiff(others, "WT (20C)"), "WT (20C)")
    }

    final_levels <- c(others, wt_ref)
    df_plot$Group <- factor(df_plot$Group, levels = final_levels)

    message("Plot Data Factor Levels distribution:")
    print(table(df_plot$Group, useNA = "always"))

    df_plot <- df_plot %>% filter(!is.na(CC_Deviation))

    message("Calculating Stats for Lineage/Depth Grouping...")
    stats_lineage_depth <- .pcb_calculate_all_stats(df_plot, c("MajorLineage", "Depth"))

    message("Calculating Stats for Individual Cells Grouping...")
    stats_cells <- .pcb_calculate_all_stats(df_plot, c("Cell"))

    output_pdf <- file.path(output_dir, paste0(name, "_comparative_boxplots.pdf"))
    message(paste("Creating Combined PDF:", output_pdf))

    pdf(output_pdf, width = 12, height = 8)

    grid.newpage()
    grid.text("Comparative Lineage Analysis", y = 0.6, gp = gpar(fontsize = 24, fontface = "bold"))
    grid.text(paste("Dataset:", name), y = 0.5, gp = gpar(fontsize = 16))
    grid.text("1. All Data Points (By Lineage)", y = 0.35, gp = gpar(fontsize = 12, col = "blue"))
    grid.text("2. Cell Means (By Lineage)", y = 0.30, gp = gpar(fontsize = 12, col = "blue"))
    grid.text("3. Individual Cells (By Lineage & Depth)", y = 0.25, gp = gpar(fontsize = 12, col = "blue"))
    grid.text("Note: Significance is Bonferroni corrected.", y = 0.15, gp = gpar(fontsize = 10, fontface = "italic"))

    # SET 1: ALL DATA POINTS
    message("Plotting Set 1: All Points...")
    df_plot$GroupID <- paste(df_plot$MajorLineage, df_plot$Depth, sep = ":")

    for (lineage in target_lineages) {
        df_sub <- df_plot %>% filter(MajorLineage == lineage)
        if (!is.null(df_sub) && nrow(df_sub) > 0) {
            print(.pcb_plot_lineage_boxplot(df_sub, lineage, "All Data Points", stats_lineage_depth, "Depth", c("MajorLineage", "Depth")))
        }
    }

    # SET 2: CELL MEANS
    message("Plotting Set 2: Cell Means...")
    df_means <- df_plot %>%
        group_by(Cell, Group, MajorLineage, Depth) %>%
        summarise(CC_Deviation = mean(CC_Deviation, na.rm = TRUE), Count = n(), .groups = "drop")

    df_means$GroupID <- paste(df_means$MajorLineage, df_means$Depth, sep = ":")

    message("Calculating Stats for Cell Means Grouping...")
    stats_cell_means <- .pcb_calculate_all_stats(df_means, c("MajorLineage", "Depth"))

    for (lineage in target_lineages) {
        df_sub <- df_means %>% filter(MajorLineage == lineage)
        if (!is.null(df_sub) && nrow(df_sub) > 0) {
            print(.pcb_plot_lineage_boxplot(df_sub, lineage, "Cell Means", stats_cell_means, "Depth", c("MajorLineage", "Depth")))
        }
    }

    # SET 3: INDIVIDUAL CELLS
    message("Plotting Set 3: Individual Cells...")
    df_plot$GroupID <- df_plot$Cell

    for (lineage in target_lineages) {
        depths <- sort(unique(df_plot$Depth[df_plot$MajorLineage == lineage]))

        for (d in depths) {
            df_sub <- df_plot %>% filter(MajorLineage == lineage, Depth == d)
            if (!is.null(df_sub) && nrow(df_sub) > 0) {
                df_sub$Cell <- factor(df_sub$Cell, levels = sort(unique(df_sub$Cell)))
                p <- .pcb_plot_lineage_boxplot(df_sub, lineage, paste("Cells - Depth", d), stats_cells, "Cell", "Cell")
                p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
                print(p)
            }
        }
    }

    dev.off()
    message("PlotComparisonBoxplots: Done.")
    invisible(output_pdf)
}


# ------------------------------------------------------------------
# Standalone-script entry point (legacy CLI)
# ------------------------------------------------------------------
# Run only when sourced via Rscript / R CMD BATCH, not when source()'d
# from another script (sys.nframe() == 0 means we are at the top frame).

if (!interactive() && sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    analysis_dir <- if (length(args) > 0) args[1] else "tbx-35_all"
    PlotComparisonBoxplots(name = analysis_dir)
}
