# LineagePhenotyping

Identify cell-cycle, position, and division-orientation defects in lineaged
*C. elegans* embryos by comparing them against a wild-type reference.

The toolchain has two halves:

1. **Per-embryo data extraction (Perl).** Process AceTree-edited movies into
   per-dataset summary files (positions, division times, cell-cycle lengths,
   angles). This part is unchanged from the original protocol — see
   [`README.txt`](README.txt) and [`docs/protocol.md`](docs/protocol.md).

2. **Statistical analysis and plotting (R).** Compare a mutant dataset
   against the bundled WT reference and produce per-cell defect calls plus
   summary plots. This is the part this README focuses on.

## Quick start

### Install dependencies

```bash
Rscript dependencies.R
```

This installs all required CRAN packages (ggplot2, dplyr, tidyr, stringr,
ggrepel, gridExtra, umap, yaml, plotly, ggtree, patchwork, scatterplot3d,
plotrix, rgl, gdata, reshape2). It is safe to re-run; already-installed
packages are skipped.

### Run the pipeline on one dataset

The recommended entry point is config-driven:

```bash
Rscript run_pipeline.R configs/example_jim721.yaml
```

Copy [`configs/template.yaml`](configs/template.yaml) to a new file in
`configs/`, fill in the paths for your dataset, and run it. The full schema
is documented in the template; the minimum is:

```yaml
name: my_dataset
data_dir: /absolute/path/to/data/root
expression_file: data/CA<my_marker>.csv
params:
  expCutoff:   200
  sig:         3
  microns:     5
  posDevTime:  250
```

Outputs land in `<data_dir>/<name>/` by default; override with `output_dir`
in the config to keep generated files separate from inputs.

### Or, the legacy CLI (still supported)

```bash
Rscript LineagePhenotyping.R <name> <expression_csv> <expCutoff> <sig> <microns> <posDevTime>
```

This expects the working directory to contain the analysis scripts,
the WT reference (`Richard_et_al_plus_comma_WT/`), and a subdirectory
`<name>/` holding the dataset's perl-pipeline outputs (positions.txt,
DivTimeNorm.tsv, etc.). Outputs are written into `<name>/`.

## Data layout

```
<data_dir>/
├── Richard_et_al_plus_comma_WT/    # WT reference (bundled)
│   ├── Richard_et_al_plus_comma_WTpositions.txt
│   ├── Richard_et_al_plus_comma_WTDivTimeNorm.tsv
│   ├── Richard_et_al_plus_comma_WTCCLengthNorm.tsv
│   ├── Richard_et_al_plus_comma_WT_rotated{X,Y,Z}.csv
│   ├── Richard_et_al_plus_comma_WT_NN_Scores.csv
│   ├── Richard_et_al_plus_comma_WT_CellMean/MaxPositionDevs.csv
│   └── Richard_et_al_plus_comma_WT_ccDevs.csv
├── data/
│   └── CA<marker>.csv               # one per expression marker (per-cell peaks)
├── <name>/                          # one per mutant dataset
│   ├── <name>positions.txt          # produced by Perl pipeline
│   ├── <name>DivTimeNorm.tsv
│   ├── <name>CCLengthNorm.tsv
│   └── <name>CCLengthMinTerminal.tsv
└── embryo_metadata.csv              # optional, used for boxplot grouping
```

`output_dir` (defaults to `<data_dir>/<name>/`, but overridable) collects
the analysis-stage outputs:

```
<output_dir>/
├── <name>_ccDevs.csv                # main cell-cycle defects table
├── <name>CellDefectSummary.csv      # per-cell defect counts (sort by TotalDefects)
├── <name>5min_3sd_Devs.csv          # filtered outliers
├── <name>_PositionDevs.csv          # per-cell-time deviations
├── <name>_CellMean/MaxPositionDevs.csv
├── <name>_NN_Scores.csv             # neighbor-deviation scores
├── <name>_dots.csv                  # division-orientation t-tests
├── <name>_rotated{X,Y,Z}.csv        # mutant positions rotated to WT frame
├── <name>_summary.pdf               # multi-page defect summary
├── <name>CC_plots.pdf               # cell-cycle scatter plots
├── <name>_mean_arrows*.pdf          # 3D arrow projections
├── <name>_Lineage_Visualizations_Grid.pdf   # UMAP/PCA over kinetics
├── <name>_Spatial_Visualizations_Grid.pdf   # UMAP/PCA over positions
├── <name>_comparative_boxplots.pdf  # mutant vs WT lineage/depth boxplots
└── run_log.txt                      # resolved config + timestamps
```

A more detailed catalogue is in [`docs/outputs.md`](docs/outputs.md).

## Configuration parameters

The five numeric knobs that affect calls:

| Field          | Default | What it does                                                            |
| -------------- | ------- | ----------------------------------------------------------------------- |
| `expCutoff`    | 200     | Minimum peak expression to call a cell "expressing"                     |
| `sig`          | 3       | Z-score threshold for division-time defects                             |
| `microns`      | 5       | Position-deviation threshold (µm)                                       |
| `posDevTime`   | 250     | Time point used for arrow / projection plots                            |
| `minDivTime`   | 70      | Cells dividing earlier than this are ignored (skip P0/AB-level noise)   |

Plus two flags:

| Field           | Default | What it does                                                         |
| --------------- | ------- | -------------------------------------------------------------------- |
| `peak_recalc`   | true    | Re-run `ReadPeakExpression` after the full mutant cell list is known |
| `dim_reduction` | true    | Run UMAP / PCA reports (~3-5 min extra; needs the `umap` package)    |

See [`docs/parameters.md`](docs/parameters.md) for tuning notes.

## Recommended downstream workflow

After a successful run, the prioritization sketch in the original
[`README.txt`](README.txt) still applies:

1. **Rate of development.** Look at `<name>CC_plots.pdf` and the rate-by-lineage
   summaries.
2. **Defect candidates.** `<name>CellDefectSummary.csv` sorted by `TotalDefects`,
   then `<name>5min_3sd_Devs.csv` for outliers, then open suspicious cells in
   AceTree to manually check for editing errors.
3. **Position/orientation phenotypes.** `<name>_summary.pdf`,
   `<name>_mean_arrows_lineage.pdf`, `<name>_dots.csv`.
4. **Embryo-level grouping.** `<name>_comparative_boxplots.pdf` and the
   UMAP/PCA grids if `dim_reduction: true`.

The historical workflow protocol lives at the
[Google Doc](https://docs.google.com/document/d/1XM8kchZvxdCX81YEYu5rknO62mZlo48bTNOYDQB5e5U/edit);
[`docs/protocol.md`](docs/protocol.md) ports the parts that are still
accurate against the new pipeline.

## Repository layout

| Path                          | What it is                                                |
| ----------------------------- | --------------------------------------------------------- |
| `run_pipeline.R`              | Single config-driven entry point (recommended)            |
| `LineagePhenotyping.R`        | Legacy positional-CLI entry point (still supported)       |
| `AnalyzeDivTimes.R`           | Step 1: division times and cell-cycle deviations          |
| `AnalyzePositions.R`          | Step 2: per-time-point spatial deviations + neighbor analysis |
| `AnalyzeRotation.R`           | Step 3: division-orientation t-tests vs WT                |
| `PlotDefectSummaries.R`       | Step 4: `<name>_summary.pdf` and the per-defect CSVs      |
| `PlotPositionDevs.R`          | Step 5: arrow projections, exp-vs-deviation               |
| `PlotComparisonBoxplots.R`    | Step 6: mutant-vs-WT boxplots (also a function)           |
| `DimensionalityReductionHelpers.R` | UMAP / PCA helpers used by steps 1 and 2             |
| `functions.R`, `GetPeak.R`, `FindNeighborDeviation.R` | shared lineage and analysis utilities |
| `Mutant_Time_Analysis.R`      | Legacy mega-script; kept as historical reference          |
| `*.pl`, `MakeDB.pm`           | Perl extraction toolchain (data preparation)              |
| `configs/`                    | YAML pipeline configs                                     |
| `docs/`                       | Protocol, parameters, outputs                             |
| `data/`, `die-1/`, `ceh-32_mutant/` | Example datasets shipped with the repo              |
| `Richard_et_al_plus_comma_WT/` | Wild-type reference (Richard et al.)                     |

## Citing

If you use this toolchain please cite the original lineage-phenotyping
work that established the WT reference and the protocol:

> Richard et al., *(citation TBD)*

## Limitations

- The lineage-tree helpers (`GetSister`, `GetParent`, `GetFounder`) and
  cell-name parsing are **specific to *C. elegans***. Generalizing to a
  different organism would require an adapter.
- The three large core functions (`AnalyzeDivTimes`, `AnalyzePositions`,
  `PlotDefectSummaries`) currently combine compute, IO, and plotting. A
  future cleanup will split these into pure compute + thin IO/plot layers
  to enable proper unit tests and `R CMD check`-clean packaging.
