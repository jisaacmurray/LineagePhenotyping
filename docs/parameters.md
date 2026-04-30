# Parameter tuning

The five numeric parameters that drive defect calls, plus the two
boolean flags. All live under `params:` in a config file.

## `expCutoff` (default 200)

Minimum peak expression to call a cell "expressing" the marker. Used
for filtering and color-coding. The right value is roughly **5–10% of
the maximum peak intensity in your expression file**, but is highly
marker-dependent.

- *ceh-36* / *blot* peak: 200 works well.
- *ceh-27* / *OP135* peak: ~1000 (raw values up to ~10800).
- *tab-1* / *SYS674* peak: 500.

If `<name>CellDefectSummary.csv` shows large numbers of "expressing
defects" in cells that look quiet by eye, raise this cutoff. If you
seem to be missing real expressing cells, lower it.

## `sig` (default 3)

Z-score threshold for division-time defect calls. A cell is flagged if
its mutant-vs-WT Z-score for any of (CC length, division time,
NN-deviation, mean position deviation) exceeds `sig`. Default 3 picks
up most real defects without too many false positives in well-edited
data.

Raise to 4 or 5 to focus on the strongest hits when first triaging a
new mutant.

## `microns` (default 5)

Magnitude cutoff (µm) for position deviation defects. Used together
with `sig` — a cell needs both Z > sig AND deviation > microns to be
called. This avoids flagging cells whose Z-score is high simply
because the WT distribution is very tight.

A reasonable range is 3–10. Raise it for noisy data, lower it for
mutants with subtle phenotypes.

## `posDevTime` (default 250)

Time point (minutes after first cleavage) used for the arrow and
projection plots. 250 catches most cells before terminal differentiation;
use earlier (e.g. 200) for early-acting genes, later (e.g. 350) for
later-acting ones.

## `minDivTime` (default 70)

Cells dividing earlier than this time are excluded. Set high enough to
skip P0/AB-level noise where small absolute deviations dominate the
relative metrics.

## `peak_recalc` (default true)

After `AnalyzeDivTimes` discovers the full union of cells present in
the dataset, re-run `ReadPeakExpression` over that cell list so cells
not directly observed in the expression file get imputed from
ancestors. Set to **false** when you need to reproduce legacy runs
that did not include this step (e.g. for diffing against the
JIM721 baseline).

## `dim_reduction` (default true)

Run the UMAP/PCA reports
(`<name>_Lineage_Visualizations_Grid.pdf` and
`<name>_Spatial_Visualizations_Grid.pdf`). Adds ~3–5 minutes of
runtime and requires the `umap` package. Turn off if you are
iterating on the kinetics or position calls and don't need the
embryo-level summaries each pass.

## Worked example: tuning a new dataset

A reasonable first pass on an unfamiliar mutant:

```yaml
params:
  expCutoff:     <peak_max * 0.05>   # 5% of marker max
  sig:           4                    # high-confidence calls only
  microns:       5
  posDevTime:    250
  minDivTime:    70
  peak_recalc:   true
  dim_reduction: false                # skip on first pass for speed
```

After triage:
- Lower `sig` to 3 to widen the net.
- Adjust `expCutoff` based on what looks "expressing" in your hands.
- Re-run with `dim_reduction: true` for the summary plots.
