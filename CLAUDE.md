# LineagePhenotyping — agent orientation

R-based pipeline for *C. elegans* embryo lineage phenotyping. Compares
mutant datasets (lineaged from AceTree movies via the **Perl** extraction
half) against the Richards et al. 2013 wild-type reference and produces
per-cell defect calls for cell cycle, position, and division orientation,
plus a battery of trees, arrows, UMAP/PCA grids, boxplots, and per-defect
CSVs.

Repo: <https://github.com/jisaacmurray/LineagePhenotyping>. Owner: John
Murray (`jisaacmurray`). `gh` is authenticated as `jisaacmurray`; the
remote URL is HTTPS with a PAT, so pushes work without prompts.

This document is for AI assistants picking up the project. The pipeline
is the **R analysis half** of a two-half toolchain. The extraction half
(historically Perl, in the Dropbox tree) is now partly owned by this repo:
`build_inputs.R` (a standalone R port of `CompareDivTime.pl` +
`ComparePositions.pl`) builds the four input tables `run_pipeline.R`
consumes, from a freeze directory of `dats/*.csv`. See
"Extraction half — `build_inputs.R`" below. `GetACD.pl` (coordinate
alignment) is being rewritten in R separately and is **not** in scope here;
`GetAngles_revRotate.pl` is **retired** — `AnalyzeRotation.R` recomputes
division angles internally, so it is not consumed downstream.

---

## Quickstart

```bash
# from the repo root:
Rscript run_pipeline.R configs/example_jim721.yaml
```

That single command runs the full 8-step pipeline against the JIM721
dataset and writes outputs to `<output_dir>/`. Total runtime ~15-20 min
on a typical Mac. See `configs/template.yaml` for the full YAML schema.

To make a new run, copy `configs/template.yaml` to `configs/<your>.yaml`,
edit `name` + `data_dir` (and `expression_file` if you have one), then
run as above. Configs are portable (relative paths resolve against the
config file's directory).

---

## Extraction half — `build_inputs.R`

Modern replacement for the Perl extraction step that feeds
`run_pipeline.R`. Standalone Rscript (does **not** source the analysis
library), so a collaborator with only a freeze of `dats/*.csv` and no
embryoDB install can build inputs.

```bash
Rscript build_inputs.R <freeze_dir> <name> [output_dir] [list_file]
# defaults: output_dir = <freeze_dir>/<name>
#           list_file  = <freeze_dir>/<name>.list  (one series_name per line)
```

The freeze dir holds per-series CSVs copied flat: `CD<series>.csv`,
`ACD<series>.csv`, `TIME<series>.csv` (and the reference series
`20081128_sulston`). `embryodb phenotyping freeze <dataset>` produces
exactly this layout plus a `configs/<name>.yaml` and a list file.

Outputs (in `output_dir`):

| File | Ported from | Consumed by `run_pipeline.R`? |
|---|---|---|
| `<name>DivTimeNorm.tsv` | `CompareDivTime.pl` | yes |
| `<name>CCLengthNorm.tsv` | `CompareDivTime.pl` | yes |
| `<name>CCLengthMinTerminal.tsv` | `CompareDivTime.pl` | yes |
| `<name>positions.txt` | `ComparePositions.pl` | yes |
| `<name>{CCLength,DivTime,CCLinNorm,STATS,CellsVsTime}.tsv`, `<name>radialpositions.txt` | both | no (parity extras the Perl scripts also emitted) |

**TIME handling:** if `TIME<series>.csv` is absent, frame index → minutes
via `minutes_per_timepoint` (read from a `minutes_per_timepoint:` line in
`<freeze_dir>/configs/<name>.yaml`, or env `BUILD_INPUTS_MPT`). Series
with neither TIME nor a fallback are skipped, matching `CompareDivTime.pl`.

**Byte-diff validation:** `die-1/` and `ceh-32_mutant/` ship the
Perl-generated golden outputs (the four `.tsv`/`.txt` above). The port
produces **byte-identical** output for both (verified). Mind R's
`check.names` sanitization (see "R column-name sanitization" below) —
`build_inputs.R` builds header strings by hand to avoid it.

**List order = column order (division tables only).** `DivTime*` and
`CCLength*` order their per-series columns by the *list-file order*
(after the leading `sulston` reference), exactly as the Perl did;
`positions.txt` sorts series in byte order and is therefore
order-independent. So a freeze whose `.list` is in a different order than
a historical Perl run yields identical column *values* but shifted column
*order* — to reproduce an old golden byte-for-byte, feed the same list
order.

**Performance:** parsing is row-by-row over the large `CD*`/`ACD*` files
(~80–100k lines each), so a 4–5-series dataset takes a few minutes.
Per-cell accumulators are split/collapsed once (not grown with `c()` in
the loop) to avoid O(n²) blowup on the ~100k-row `positions.txt`.

The Perl `int()` semantics are mirrored exactly: `sig_figs()` truncates
toward zero, the least-squares fit replays `Statistics::Descriptive`'s
plain-double left-to-right sums, and sorts use byte/ASCII order
(`sort(method="radix")`) to match Perl's default sort.

---

## Project layout

Top-level R sources are the analysis library. They `source()` each other
relatively; `run_pipeline.R` handles loading them in dependency order.

| File | Role |
|---|---|
| `build_inputs.R` | **Extraction-half entry point.** Standalone Rscript port of `CompareDivTime.pl` + `ComparePositions.pl`. Reads a freeze dir of `dats/*.csv` (per-series `CD*.csv`/`ACD*.csv`/`TIME*.csv`) + a series list, emits the four tables `run_pipeline.R` reads (`<name>DivTimeNorm.tsv`, `<name>CCLengthNorm.tsv`, `<name>CCLengthMinTerminal.tsv`, `<name>positions.txt`) plus parity extras. Does not source the analysis library. See section below. |
| `run_pipeline.R` | **Analysis-half entry point.** Parses YAML, sources helpers, runs the 8-step pipeline. |
| `functions.R` | Shared utilities (`GetParent`, `GetSister`, `GetFounder`, name parsing). Sourced by every other helper. |
| `plot_paths.R` | `.plots_dir()` helper that routes display outputs into `plots/<kind>/` under `subdir_layout: by_kind` (default). |
| `AnalyzeDivTimes.R` | Step 1: cell-cycle deviation. Writes `_ccDevs.csv`, `_CCdev_matrix.txt`, `CC_plots.pdf`, calls dim-reduction for lineage kinetics. **~390 LOC**, mixes compute + IO + plotting. |
| `AnalyzePositions.R` | Step 2: position deviation + per-mutant rotation alignment. Writes `_PositionDevs.csv`, `_Cell{Mean,Max}PositionDevs.csv`, `_NN_Scores.csv`, `_rotated{X,Y,Z}.csv`, calls dim-reduction for spatial phenotyping, plots WT stats + per-embryo PositionPlots + positionDefects. **~960 LOC**. Contains `OptimizeRotations` (per-timepoint theta optimization) and `PlotExpVsDev`. |
| `AnalyzeRotation.R` | Step 3: rotated division-orientation analysis. Writes `_dots.csv`. |
| `PlotDefectSummaries.R` | Step 4: integrates CC + position + angle defects into per-cell summary table + summary PDF. **~840 LOC**. |
| `PlotPositionDevs.R` | Step 5: 3D arrow plots (`PlotDeviationsSingle` does the rendering). |
| `PlotComparisonBoxplots.R` | Step 7: mutant-vs-WT boxplots by lineage/depth/cell with Wilcoxon stats. |
| `DimensionalityReductionHelpers.R` | UMAP/PCA reports. Called from `AnalyzeDivTimes` (lineage kinetics grid) and `AnalyzePositions` (spatial grid). |
| `DefectScoreFrames.R` | Phase 5 (new). Adapters converting LineagePhenotyping outputs into the canonical `(cell, time, value[, embryo])` frame for `LIVEtools::plot_lineage_tree()`. Functions: `cc_dev_long`, `position_dev_long_from_csv`, `position_dev_per_t_long`, `position_dev_per_t_components_long` (AP-signed + radial-signed), `dot_dev_long`, `arbitrary_long`. |
| `PlotDefectTrees.R` | Phase 5 (new). Orchestrator: `PlotDefectTrees()` renders the full defect-tree set; `defect_scheme()` builds color schemes per defect kind; `replot_defect_subset()` re-renders a single sublineage without a pipeline rerun. |
| `CellCountDiagnostic.R` | Phase 5.3 (new). `analyze_cell_counts_per_t()` + `plot_cell_counts_per_t()`. Auto-runs in every pipeline; useful for spotting late-timepoint tracking cliffs. |
| `dependencies.R` | One-shot installer. CRAN + Bioconductor + GitHub (LIVEtools). |
| `configs/template.yaml` | YAML reference — every knob is documented here. |
| `configs/example_jim721.yaml` | Canonical end-to-end test config. |
| `tests/test_defect_trees.R` | Synthetic-data smoke test (16+ assertions). |
| `docs/outputs.md` | What each output file contains. |
| `docs/protocol.md`, `docs/parameters.md` | High-level docs. |

Reference data shipped with the code:

| File | Purpose |
|---|---|
| `CellNames.csv` | Sulston cell name table. |
| `Cells_350min_lineageOrder.csv` | Canonical cell ordering for plots. |
| `SupplementalTable2_DivisionTimes.txt` | WT division times (for cross-referencing in PCA loadings). |
| `data/CA*.csv` | Real lab CA (one-row-per-cell `blot`) files. When the YAML has no `expression_file`, `run_pipeline.R` falls back to the first `data/CA*.csv` (sorted) as a **placeholder** — reporter-specific, so set `expression_file` for a meaningful peak-expression annotation. |
| ~~`ceh36_peak.csv`~~ | **Retired as the fallback** — it was malformed (header `"","x"`, all-NA), crashing `ReadPeakExpression` with `duplicate 'row.names'`. The bundled `data/CA*.csv` files replace it. |

---

## Pipeline steps

`run_pipeline.R` runs these in order. The defect-tree section (5/5.1+)
is new; the rest comes from the earlier phase-2 refactor.

1. **`AnalyzeDivTimes(...)`** → CC deviation tables + lineage UMAP/PCA grid
2. **`AnalyzePositions(...)`** → position deviation tables + rotated CSVs + spatial UMAP/PCA grid
3. **`AnalyzeRotation(...)`** → division-orientation deviation (uses rotated XYZ from step 2)
4. **`PlotDefectSummaries(...)`** → integrated defect summary PDF + counts CSVs
5. **`PlotDeviationsList(...)`** (from `PlotPositionDevs.R`) → 3D arrows
6. **`PlotExpVsDev(...)`** (from `AnalyzePositions.R`) → expression-vs-deviation scatter
7. **`PlotComparisonBoxplots(...)`** → mutant-vs-WT boxplots with stats
7b. **Cell-count cliff diagnostic** (`CellCountDiagnostic.R`)
8. **`PlotDefectTrees(...)`** → defect-colored lineage trees (CC, position, angle, AP-signed, radial-signed)

Each step takes explicit `data_dir`, `output_dir`, `wt_ref_dir` args
(plus a few per-step ones). The runner threads YAML values through.
Pipeline gates `dim_reduction:` (skip UMAP/PCA reports) and
`defect_trees:` (skip Phase 5 trees) for faster iteration.

---

## Output structure (Phase 5.2 default: `subdir_layout: by_kind`)

```
<output_dir>/
  *.csv, *.tsv, *.txt, run_log.txt            ← data files (read by downstream code)
  plots/
    cc/                                       ← CC dev plots
      <name>CC_plots.pdf
      <name>_summary.pdf
      <name>_DefectTrees_CCDev_{mean,max}.pdf
      per_embryo/
    position/
      cell/                                   ← per-cell aggregates
        <name>_DefectTrees_PositionDev{Mean,Max}.pdf
        <name>_positionDefects.pdf
        <name>_posDirectionPlots.pdf
        per_embryo/<name>_<emb>_PositionPlots.pdf
      per_t/                                  ← per-timepoint trees
        <name>_DefectTrees_PositionDev_per_t.pdf
        <name>_DefectTrees_PositionDev_AP_per_t.pdf
        <name>_DefectTrees_PositionDev_radial_per_t.pdf
        per_embryo/
      arrows/                                 ← 3D arrow plots
        <name>_mean_arrows*.pdf
        per_embryo/<name>_<emb>_arrows.pdf
    angle/
      <name>_DefectTrees_DotDev_mean.pdf
      per_embryo/
    expression/
      <name>_ExpVsDev{,_labeled}.pdf
    boxplots/
      <name>_comparative_boxplots.pdf
    summary/                                  ← multi-panel grids
      <name>_Spatial_Visualizations_Grid.pdf
      <name>_Lineage_Visualizations_Grid.pdf
      <name>_WT_stats.pdf
    _diagnostics/
      <name>_cell_counts_vs_t.pdf
```

For byte-identical pre-Phase-5.2 output layout (everything flat under
`<output_dir>`), set `subdir_layout: flat` in the YAML.

---

## Key YAML knobs

See `configs/template.yaml` for the complete schema. Most-used:

```yaml
name: my_dataset
data_dir: /path/to/data
# output_dir, wt_ref_dir, script_dir, expression_file, embryo_metadata_file all optional

params:
  expCutoff:   200          # min expression to call a cell "expressing"
  sig:         3            # Z-score threshold for division-time defects
  microns:     5            # position deviation threshold (µm)
  posDevTime:  250          # time point for arrow projections
  minDivTime:  70           # ignore cells dividing before this time
  peak_recalc: true         # recompute peak after AnalyzeDivTimes finds full cell set
  dim_reduction: true       # run UMAP/PCA reports (skip with false to save 3-5 min)

  # Phase 5: defect-colored lineage trees
  defect_trees:        true
  defect_tree_root:    P0
  defect_tree_split_roots: [ABa, ABp, P1]
  # defect_tree_milestone: { lineage: MS, n_cells: 32 }   # optional crop

  # Phase 5.2: output organization
  subdir_layout:   by_kind  # or "flat"
  keep_arrow_jpgs: false    # delete temp arrow jpg dirs after PDFs sealed

  # Phase 5.3: cliff filter (drops per-t values where embryo has < N cells alive)
  min_cells_per_timepoint: 200

  # Phase 5.4: grey-band sensitivity overrides
  # defect_scheme_bands:
  #   position_dev_AP:            1.5
  #   position_dev_radial_signed: 1.5
  #   cc_dev:                     5
```

---

## Development workflow

### Branches and PRs

- `main` is the merged branch.
- Active development branch (Phase 5+): `phase5-defect-trees` — PR #7.
- Use `gh pr view --json url`, `gh pr list`, etc.

### Run the smoke test

```bash
Rscript tests/test_defect_trees.R
```

Synthetic-data assertions for `LIVEtools::plot_lineage_tree()`, the
adapters in `DefectScoreFrames.R`, and `.plots_dir()` routing. Runs
in ~10s. Fast feedback loop for any change to these files.

### Run the canonical end-to-end test

```bash
Rscript run_pipeline.R configs/example_jim721.yaml
```

JIM721 is the canonical dataset. After Phase 5.4/5.5 the expected
output is ~72 defect-tree PDFs + the rest of the standard outputs.
Total runtime ~15-20 min.

Phase-2 baseline outputs (pre-refactor reference for byte-diff
verification) live at:

```
~/Library/CloudStorage/Dropbox/Reverse_Genetics/LineagePhenotyping/JIM721_baseline_2026-04-27/
```

Data-file outputs (`.csv`, `.tsv`, `.txt`) should be byte-identical to
this baseline for any pure refactor; display files (PDF/PNG/JPG) may
differ trivially due to embedded paths or layout reorg.

### Verifying a change locally

1. Syntax check first: `Rscript -e "parse(file='X.R')"`.
2. Smoke test if you touched `DefectScoreFrames.R` / `PlotDefectTrees.R`
   / `plot_paths.R`.
3. JIM721 end-to-end if you touched a step in the pipeline.

---

## Dependencies

`dependencies.R` is the one-shot installer:

- CRAN: `ggplot2`, `ggrepel`, `dplyr`, `tidyr`, `stringr`, `reshape2`,
  `scales`, `yaml`, `plotly`, `patchwork`, `scatterplot3d`, `plotrix`,
  `rgl`, `umap`, `gdata`, `ape`, `remotes`.
- Bioconductor: `ggtree`.
- GitHub: `LIVEtools` from `johnmurraylab/LIVE_tools` (provides
  `plot_lineage_tree`, `paginate_tree_plots`, etc. — required by
  Phase 5 defect-tree plotting).

Run once with `Rscript dependencies.R`.

---

## Known quirks / gotchas

### Sparse checkout — the working tree may look smaller than HEAD

The side-copy at `/Users/jmurr/MacTools/LineagePhenotyping_update/`
uses git sparse-checkout: only top-level files plus `configs/`,
`docs/`, `tests/` are materialized. Most of the legacy Perl scripts
and subdirectories are hidden.

**Pattern at `.git/info/sparse-checkout`:**

```
/*
!/*/
/configs/
/docs/
/tests/
```

**DO NOT** run `git read-tree -m -u HEAD` after editing this file —
it overwrites the working tree from HEAD and **deletes untracked
files**. To extend sparse-checkout safely, either:

- Use `git sparse-checkout add <path>` (modern command), or
- Append to the file and use `git checkout` instead of `read-tree -u`.

### Dropbox sibling vs side-copy

Two separate copies of this codebase exist on the user's machine:

| Path | Purpose |
|---|---|
| `~/Library/CloudStorage/Dropbox/Reverse_Genetics/LineagePhenotyping/` | **Active analysis tree.** Where data + outputs live. Has ~890 working-tree changes pre-update. **Treat as untouchable** — Dropbox version history is the rollback. Don't `git pull` between these two. |
| `/Users/jmurr/MacTools/LineagePhenotyping_update/` (this side-copy) | **Development copy.** Sparse-checkout from the GitHub remote. All branch / push activity happens here. |

### Naming conventions to respect

- `Richard_et_al_plus_comma_WT/` (no 's' in the directory name) — the
  WT reference bundle. Don't rename: it's the file-name prefix for
  every WT file (e.g. `Richard_et_al_plus_comma_WTDivTimeNorm.tsv`).
  The paper's first author is "Richards" (with 's') — known typo,
  but the directory name is load-bearing.
- `<name>_ccDevs.csv` is the rich per-row deviation table.
- `<name>_CCdev_matrix.txt` is the cells × embryos matrix (Phase 5.2
  rename from `<name>CCdev.txt`).
- `_CCdev_matrix.txt` is **only** consumed by `AnalyzeDivTimes.R`
  itself; nothing else reads it.
- `embryo_metadata.csv` is per-dataset, lives in `data_dir`. Optional.
- `<name>_dots.csv` is the division-orientation dot-product table.

### Latent issues flagged for future cleanup

- **AnalyzePositions.R has a duplicate UMAP/PCA call.** The Spatial
  Visualizations Grid is generated twice (lines ~926 and ~977); both
  write to the same PDF, so the second clobbers the first. Doubles
  runtime of that step but is functionally harmless. Memory file
  notes this as "two duplicate UMAP/PCA blocks". Folds into the
  eventual decomposition phase.

- **The three giants** (`AnalyzeDivTimes` ~390 LOC, `AnalyzePositions`
  ~960 LOC, `PlotDefectSummaries` ~840 LOC) each mix compute + IO +
  plotting in a single function. Splitting them is the prerequisite
  to an R-package conversion. Not yet started; flagged as Phase 7.

- **Phase 5.3H2(b) — alignment fix is deferred.** `OptimizeRotations`
  re-optimizes per-frame theta independently; when n_cells drops at
  late timepoints, the optimization is dominated by survivors and
  produces artefactual deviations. The Phase 5.3H2(a) viz filter
  (`min_cells_per_timepoint`) masks the issue at visualization time,
  but the proper fix (freeze theta when n_cells < threshold) is not
  done — needs a separate PR with byte-diff data verification.

### R column-name sanitization

This burned us in Phase 5.7. `data.frame(...)` with default
`check.names = TRUE` rewrites `tab-1_L4` → `tab.1_L4` and adds an
`X` prefix to names starting with a digit. `read.table(...,
check.names = FALSE)` preserves whatever's in the file. The CC
Deviations panel of the Lineage Visualizations Grid was rendering
every point black for months because the embryo-ID match in
`DimensionalityReductionHelpers.R` failed silently against the
sanitized names. The fix at line 49-65 of that file matches against
both raw and normalized forms. If you add new dim-reduction inputs,
be aware of this.

### LIVEtools is a separate repo

`LIVEtools::plot_lineage_tree()` is consumed by Phase 5. It lives at
`~/Library/CloudStorage/Dropbox/Papers/YG_3D_micro/LIVE_tools/` (a
git clone of `johnmurraylab/LIVE_tools`). The active LIVEtools branch
is `continuous-tree-plots` (PR #9). When debugging defect-tree
rendering you may need to read/edit code there too. The tree-plot
function defaults: `value_col = "blot"`, `resolution = "cell"`.

---

## Output decision tree

If you want to answer "where is X produced?":

| Output | Producer | Subdir |
|---|---|---|
| `_ccDevs.csv` | `AnalyzeDivTimes` | top-level |
| `_CCdev_matrix.txt` | `AnalyzeDivTimes` | top-level |
| `CC_plots.pdf` | `AnalyzeDivTimes` | `plots/cc/` |
| `_PositionDevs.csv`, `_Cell{Mean,Max}PositionDevs.csv`, `_rotated{X,Y,Z}.csv`, `_NN_Scores.csv` | `AnalyzePositions` | top-level |
| `_WT_stats.pdf`, `_Spatial_Visualizations_Grid.pdf`, `_Lineage_Visualizations_Grid.pdf` | dim-reduction helper | `plots/summary/` |
| `_positionDefects.pdf`, `_posDirectionPlots.pdf`, per-embryo `PositionPlots.pdf` | `AnalyzePositions` + `PlotPositionDevs` | `plots/position/cell/` |
| `_mean_arrows*.pdf`, per-embryo `_arrows.pdf` | `PlotPositionDevs::PlotDeviationsSingle` | `plots/position/arrows/` |
| `_dots.csv` | `AnalyzeRotation` | top-level |
| `_DefectTrees_*.pdf` | `PlotDefectTrees` | `plots/<cc\|position\|angle>/[per_embryo/]` |
| `_summary.pdf` | `PlotDefectSummaries` | `plots/summary/` |
| `_comparative_boxplots.pdf` | `PlotComparisonBoxplots` | `plots/boxplots/` |
| `_ExpVsDev{,_labeled}.{pdf,csv}` | `AnalyzePositions::PlotExpVsDev` | `plots/expression/` |
| `_cell_counts_vs_t.{pdf,csv}` | `CellCountDiagnostic` | `plots/_diagnostics/` (PDF) + top-level (CSV) |

`docs/outputs.md` has full per-file column descriptions.

---

## Where to look first when…

- **A defect-tree color looks wrong** → `PlotDefectTrees.R::defect_scheme()` and `LIVEtools::plot_lineage_tree`'s `color_values` / `value_min` / `value_max` args. YAML override: `defect_scheme_bands`.
- **PCA/UMAP plot has wrong/missing coloring** → `DimensionalityReductionHelpers.R::.dim_red_prepare_meta()` line 49-65. Check whether the data frame's column names match `mutant_ids` after any sanitization.
- **A new output file is missing from `plots/<kind>/`** → search for the file's basename in the analysis scripts. Any `pdf()` / `ggsave()` call site should be wrapped in `file.path(.plots_dir(output_dir, "<kind>"), ...)`. If not, add it.
- **Per-timepoint defect tree shows artefactual late spikes** → set `min_cells_per_timepoint: 200` (or higher) in the YAML. Or run `CellCountDiagnostic.R` standalone to inspect the cliff.
- **A tree's branches end too long with grey NA tails** → `plot_lineage_tree(...)` defaults to `truncate_to_last_data = TRUE`; check that's not been overridden.
- **An empty leaf cell shows as grey in a CC plot** → CC scheme uses `drop_empty_cells = TRUE` by default in Phase 5.1+. If something regressed, check `defect_scheme("cc_dev")`'s output.
- **You need to re-render trees for a sublineage without rerunning the pipeline** → `PlotDefectTrees.R::replot_defect_subset(name, output_dir, sublineage = "MSpap")`.
- **Pipeline is too slow during dev** → set `dim_reduction: false` and `defect_trees: false` in your YAML.
- **A path resolution looks wrong** → `script_dir` (Phase 5.5) controls where reference files are located; defaults to where `run_pipeline.R` lives. If your scripts and data are in separate trees, set `script_dir:` explicitly.

---

## Project history (for context)

The Phase-2 refactor (PRs #2-#5, merged) modernized the original
per-dataset wrapper scripts (`LineagePhenotyping_<dataset>.R`) into
a single YAML-driven runner with explicit-path args, `dependencies.R`,
and docs. Phases 5/5.1-5.7 (PR #7, in flight) add the continuous-
resolution defect-tree plotting (via the new LIVEtools function),
per-defect color schemes, output reorg, output-rename + jpg cleanup,
a sublineage replot helper, fixed-axes arrows, the cell-count
diagnostic, the late-timepoint viz filter, the `script_dir` config
from a fork, narrower default sensitivity bands, and the CC-Dev
color fix.

Future roadmap (planned but not started):

- **Phase 5.3H2(b)**: root-cause alignment fix (per-t theta freezing).
- **Phase 6**: Spatial + trajectory views via more LIVEtools wraps.
- **Phase 7**: Decompose the three giant functions into compute / IO /
  plot units.
- **Phase 8**: R-package conversion of LineagePhenotyping.
- **Phases 9-13**: organism generalization, testthat harness,
  cross-machine validation, UMAP cleanup, full docs port.
- **Long horizon**: EmbryoDB.jar integration calling the package-ified
  plotters.
