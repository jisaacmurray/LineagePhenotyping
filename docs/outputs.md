# Output reference

A successful pipeline run writes ~60 files into `output_dir` (defaults
to `<data_dir>/<name>/`). The most important ones, in roughly the
order you'll consult them:

## Cell-cycle defect tables

### `<name>_ccDevs.csv`
Long-format table: one row per (cell, embryo) pair. Key columns:

| Column                | Meaning                                                                |
| --------------------- | ---------------------------------------------------------------------- |
| `Cell`                | Cell name                                                              |
| `Embryo`              | Embryo name (one of the movies in your list)                           |
| `CC_Deviation`        | Mutant cell-cycle length minus WT mean (minutes)                       |
| `CC_Z_score`          | Z-score relative to WT cell-cycle distribution                         |
| `Terminal_Cell_Length`| For non-dividing cells, length of observed terminal cycle              |
| `Terminal_Cell_Z`     | Z-score for the terminal-length test                                   |
| `Terminal_cell_delta` | Terminal length minus WT mean (positive = held longer than expected)   |
| `WT_CC`, `WT_Div_time`| WT mean cell-cycle length and division time for the cell               |
| `Div_time`            | Mutant division time                                                   |

Filtered downstream into the per-cell summary and the 3-sigma table.

### `<name>CellDefectSummary.csv`
One row per cell. Sort by `TotalDefects` to triage.

| Column          | Meaning                                                            |
| --------------- | ------------------------------------------------------------------ |
| `Cells`         | Cell name                                                          |
| `ExpVals`       | Peak expression value (from your expression file)                  |
| `LateCalls`     | Count of embryos in which the cell divided late                    |
| `EarlyCalls`    | Count of embryos in which the cell divided early                   |
| `MissedCalls`   | Cell observed but not seen to divide (held terminal)               |
| `EctopicCalls`  | Cell observed dividing in mutant but absent in WT reference        |
| `LostCalls`     | Sister-pair end-time deviation > 5 min (likely lost)               |
| `DeathCandidates` | Plausible death calls                                            |
| `TotalDefects`  | Sum of the above                                                   |
| `CC_Dif_pValues`| t-test p-value for CC length difference                            |

### `<name>5min_3sd_Devs.csv`
Pre-filtered subset of `<name>_ccDevs.csv` showing only the cell-embryo
combinations with both `CC_Z_score` > 3 and absolute `CC_Deviation` > 5.
The "manual QC" target list — open each in AceTree and verify.

### `<name>CC.txt` / `<name>CCdev.txt`
Tab-separated wide-format dumps of per-cell CC lengths and deviations
across embryos. `<name>CCdev.txt` is the file to open and sort by
embryo when iterating with AceTree.

## Position / spatial defect tables

| File                                    | Contents                                                          |
| --------------------------------------- | ----------------------------------------------------------------- |
| `<name>_PositionDevs.csv`               | Per (cell, time) deviation magnitude for each embryo              |
| `<name>_CellMeanPositionDevs.csv`       | Per-cell mean across times                                        |
| `<name>_CellMaxPositionDevs.csv`        | Per-cell max across times                                         |
| `<name>_NN_Scores.csv`                  | Neighbor-deviation (log2 ratio of mutant/WT pairwise distances)   |
| `<name>_rotated{X,Y,Z}.csv`             | Mutant positions rotated to align with the WT axes                |
| `<name>_dots.csv`                       | Per-cell sister-vector dot products and t-test p-values vs WT     |

## Defect-classification dumps

| File                                          | Contents                                                |
| --------------------------------------------- | ------------------------------------------------------- |
| `<name>_defects.csv`                          | Wide table of per-cell defect classifications           |
| `<name>_defect_counts.csv`                    | Aggregated counts per cell                              |
| `<name>_defect_counts_with_relatives.csv`     | Counts including sister/parent/aunt context             |
| `<name>_defects_with_ancestral_CC.csv`        | Adds ancestral cell-cycle context                       |
| `<name>_Term_defects.csv`                     | Terminal-cell defect detail                             |

## Plots

### Multi-page summaries

| File                              | What it shows                                                       |
| --------------------------------- | ------------------------------------------------------------------- |
| `<name>_summary.pdf`              | Spatial defect maps, expression-vs-defect heatmaps, frequency grids |
| `<name>CC_plots.pdf`              | Mutant vs WT cell-cycle and division-time scatter (one per embryo) |
| `<name>_WT_stats.pdf`             | WT reference variability (sanity check)                             |
| `<name>_positionDefects.pdf`      | Significance and severity scatterplots for spatial defects          |
| `<name>_posDirectionPlots.pdf`    | Time-resolved per-axis position bias                                |
| `<name>_comparative_boxplots.pdf` | Mutant-vs-WT boxplots by lineage and depth (Wilcoxon, Bonferroni)   |

### Arrow / projection plots

| File                                       | What it shows                                                           |
| ------------------------------------------ | ----------------------------------------------------------------------- |
| `<name>_mean_arrows.pdf`                   | Mean WT→mutant displacement vectors, colored by deviation magnitude     |
| `<name>_mean_arrows_lineage.pdf`           | Same arrows, colored by lineage founder (AB[apl], MS, C, D, E)          |
| `<name>_mean_arrows_exp.pdf`               | Same arrows, colored by peak expression                                 |
| `<datestamp>_..._arrows.pdf`               | One per individual mutant embryo                                        |
| `<name>_<embryo>_PositionPlots.pdf`        | Per-embryo, per-time-point projection grids                             |

### UMAP / PCA reports (only if `dim_reduction: true`)

| File                                                          | What it shows                                  |
| ------------------------------------------------------------- | ---------------------------------------------- |
| `<name>_Lineage_Visualizations_Grid.pdf`                      | Embryo embeddings over kinetics metrics        |
| `<name>_Spatial_Visualizations_Grid.pdf`                      | Embryo embeddings over spatial metrics         |
| `<name>_*_PCA_EigenCells.csv`                                 | Cell loadings for the PCA components           |

### Expression vs deviation

| File                              | What it shows                                                       |
| --------------------------------- | ------------------------------------------------------------------- |
| `<name>.ExpVsDev.pdf`             | Mean position deviation vs peak expression, color-coded             |
| `<name>.ExpVsDev_labeled.pdf`     | Same, with text labels on outlier cells                             |
| `<name>.ExpVsDev.csv`             | Underlying data                                                     |

## Run metadata

| File              | Contents                                                    |
| ----------------- | ----------------------------------------------------------- |
| `run_log.txt`     | Resolved config + parameter values + start/end timestamps   |
| `<name>_summary_output.txt` | Captured stdout / printed stats from `PlotDefectSummaries` |
