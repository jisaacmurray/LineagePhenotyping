# Protocol

Full pipeline from raw AceTree-edited movies to a defect-call summary
PDF. The Perl half (steps 1–5) is unchanged from the historical workflow
documented in [`../README.txt`](../README.txt) and the
[Google Doc](https://docs.google.com/document/d/1XM8kchZvxdCX81YEYu5rknO62mZlo48bTNOYDQB5e5U/edit).
Steps 6–7 use the new config-driven runner.

## 1. Edit the movies

Edit, check, extract, reopen in AceTree, resave, and re-extract each
movie. Make sure your partial-editing codes are correctly specified in
embryoDB.

## 2. Make a movie list

Plain text file with one movie name per line. Conventionally lives in
`<HPC_root>/lists/<your_dataset>`. The list filename becomes the
dataset's `name` for the rest of the pipeline.

Optionally re-extract all movies at once via `sageExtract.pl <list>`
or `sageNotextract.pl <list>` (run from the `tools3` directory).

## 3. Validate trees

```bash
PrintTrees.pl /path/to/lists/<name> -200 500 rainbow 5
```

(Run from the `tools3` directory.) Visually confirm the trees and
editing look right before continuing.

## 4. Extract per-cell summary data

```bash
GetACD.pl /path/to/lists/<name>
```

(Also from `tools3`.) Produces the per-movie ACD outputs that
`RunAllAnalysis.pl` will aggregate.

## 5. Aggregate to a per-dataset summary

```bash
RunAllAnalysis.pl /path/to/lists/<name>
```

Run this from the directory containing the LineagePhenotyping R scripts.
It writes a `<name>/` subdirectory containing:

- `<name>positions.txt`
- `<name>DivTimeNorm.tsv`
- `<name>CCLengthNorm.tsv`
- `<name>CCLengthMinTerminal.tsv`
- (plus per-embryo angles files used by future analyses)

## 6. (Optional) Provide an expression file

If you care about whether defects are enriched in expressing cells,
copy a `CA<gene>.csv` file for the marker of interest into
`<data_dir>/data/`. Otherwise the pipeline runs without expression
coloring (or you can fall back to `ceh36_peak.csv` already in the repo).

## 7. Run the analysis

### Recommended: config-driven

```bash
cp configs/template.yaml configs/<name>.yaml
# Edit <name>.yaml: name, data_dir, expression_file, params
Rscript run_pipeline.R configs/<name>.yaml
```

### Legacy: positional CLI

Still supported. Equivalent to invoking step 7 with all defaults plus
the listed positional args:

```bash
Rscript LineagePhenotyping.R <name> <expression_csv> <expCutoff> <sig> <microns> <posDevTime>
```

For example, the historical command-line for the ceh-32 dataset was:

```bash
Rscript LineagePhenotyping.R ceh-32_mutant data/CA_CEH-32_SYS85.csv 200 3 3 250 \
    > ceh-32_mutant/ceh-32_mutant.stats.txt
```

## 8. Triage the outputs

The output catalogue is in [outputs.md](outputs.md). The recommended
prioritization is:

1. **Rate of development.** `<name>CC_plots.pdf` and the rate-by-lineage
   summaries.
2. **Defect candidates.** Sort `<name>CellDefectSummary.csv` by the
   `TotalDefects` column. For deep dives use `<name>5min_3sd_Devs.csv`.
3. **QC the top hits in AceTree.** Open `<name>CCdev.txt`, sort by
   embryo, and manually verify each suspect cell — many high-Z hits are
   editing artefacts that need to be fixed at the AceTree stage.
4. **Position / orientation phenotypes.** `<name>_summary.pdf`,
   `<name>_mean_arrows_lineage.pdf`, `<name>_dots.csv`.
5. **Embryo grouping / treatment effects.** `<name>_comparative_boxplots.pdf`,
   plus the UMAP/PCA grids if `dim_reduction: true`.

## Iterating

The first pass usually surfaces editing errors in addition to real
phenotypes. Fix the AceTree edits → re-run from step 4 (`GetACD.pl`) →
re-run step 7. Most defect calls converge after one or two iterations.
