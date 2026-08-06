# scExplorer - User Guide

A walkthrough of every tab and every option in the app: what each control does,
what information it gives you, and what you'd use it for. scExplorer is an
upload-based browser app for exploring any Seurat object - no coding needed.

---

## Loading data (the "Get started" screen)

- Upload your own data → drag-and-drop or Browse a Seurat object saved as `.rds`.
- Use example data → loads a small bundled demo object.

Once a valid Seurat object is read, the tabs appear. The uploader is only available
*before* a dataset is loaded - to switch datasets, reload the page.

### Summary cards
Three cards summarise the object:

| Card | What it shows |
|---|---|
| Cells | Number of cells (`ncol` of the object) |
| Genes | Number of genes/features (`nrow`) |
| Samples | Number of distinct samples/batches, detected from a sample/batch column (shows `-` if none found) |

### Assays/layers strip
Below the cards, a strip lists every assay and the layers it contains, e.g.:

> Assays/layers: RNA: `counts`  SCT: `counts` `data` `scale.data`  integrated: `data` `scale.data`

This tells you at a glance which assays exist and which data layers are available in
each - useful because the Assay/Data-layer menus throughout the app only offer what's
actually present. (E.g. an RNA assay with only `counts` has never been normalized.)

---

## Downloading plots and tables (works on every tab)

- Plots → PNG. Every plot has a ⬇ PNG button in its top-right corner. It saves
  a static PNG named after the plot's on-screen title (e.g. `Malat1_·_Normalized_log.png`),
  so files are self-describing.
- Download all as one PNG. When a tab shows more than one panel (Gene Expression,
  Split by Group), a ⬇ Download all as one PNG button appears above the grid. It
  composites every panel (plus the shared colorbar) into a single image.
- Tables → Excel. Every table (Cluster Markers, Contrasts, Aggregated) has a
  Download Excel button that saves the entire table (all rows, respecting any
  search filter) as an `.xlsx` - opens in Excel, LibreOffice, or Google Sheets.

---

## Two concepts used across many tabs

### Assay
Seurat objects can hold several assays (e.g. `RNA`, `SCT`, `integrated`, `ADT`).
Where an Assay dropdown appears, it selects which assay's values to use. Notes:

- It only appears when the object has more than one assay.
- On Gene Expression it includes *all* assays (viewing any assay is harmless).
- On Cluster Markers and Contrasts it excludes `integrated` - differential
  expression on batch-corrected values is statistically invalid - and an explanatory
  note is shown.

### Data layer
Each assay can store the same genes at different processing stages, as layers:

| Layer | Label in app | What it is | Typical use |
|---|---|---|---|
| `data` | Normalized (log) | Log-normalized expression | The default for viewing and for most marker tests |
| `counts` | Raw counts | Raw UMI counts | Count-based models (Poisson); seeing raw depth |
| `scale.data` | Scaled (z-score) | Z-scored (mean 0, unit variance), variable features only | Comparing relative pattern; note many genes are absent here |

Wherever a Data layer menu appears, the app shows the raw values of that layer
(no extra z-scoring) - so you literally see the numbers at that stage.

---

## 1. Clusters / Metadata

The map of your cells: each point is one cell, positioned by the UMAP embedding, and
coloured by a metadata column you choose.

| Option | What it does | Information / use |
|---|---|---|
| Resolution / clustering *(appears only if the object has clustering columns)* | Quick-pick a clustering column (e.g. `SCT_snn_res.0.1`), listed as `assay · res X (N clusters)` | Recolour the map by that clustering. Higher resolution = more, finer clusters; the points don't move, only their colours. Stays in sync with "Color UMAP by" (picking a clustering in either box updates the other). |
| Color UMAP by | Colour cells by any metadata column (clusters, cell type, sample, phase, condition, QC scores…) | See structure and how any variable is distributed across the map. |
| Show categories (checkboxes) *(categorical columns)* | Toggle individual categories on/off; unticked ones turn grey | Isolate one or a few groups - e.g. grey out everything except "T cells". Toggling is instant (no re-render). |
| Min / Max boxes *(continuous columns)* | Grey out cells whose value falls outside the typed range | Focus on a value window - e.g. only cells with a QC score above a threshold. |

Categorical vs continuous is auto-detected. Discrete labels (clusters, cell types)
get distinct colours + the "Show categories" checkboxes. Numeric columns (counts, QC
metrics, prediction scores - even numbers stored as text) get a viridis gradient +
colorbar, exact hover value, and the Min/Max range filter. Clustering columns are
always treated as categorical, however many clusters they have (a 60-cluster leiden
result is 60 categories, not a scale).

Category order. Wherever categories are listed (the legend, the "Show categories"
checkboxes, and the group pickers on the other tabs), they are sorted naturally -
numerically when the labels are numbers (so `0, 1, 2, …, 10, 11`, not `0, 1, 10, 11, 2`),
otherwise alphabetically.

Used for: understanding the overall structure - what clusters/cell types exist, how
they're arranged, and how any metadata variable maps onto them.

---

## 2. Gene Expression

One UMAP per gene, coloured by that gene's expression on a diverging
blue → white → red scale. Values are the raw values of the chosen data layer
(no app-side z-scoring), so:

- For Normalized (log) / Raw counts (non-negative): the scale runs 0 → max
  (low = blue, high = red).
- For Scaled (z-score): the scale is symmetric about 0 (blue = below average,
  white = average, red = above average).

| Option | What it does | Information / use |
|---|---|---|
| Compare | Two modes (below) | Choose whether to compare *genes* or *layers*. |
| &nbsp;&nbsp;• *Genes (in one data layer)* | One UMAP per selected gene, all in one layer | Compare several genes' expression patterns. |
| &nbsp;&nbsp;• *One gene across data layers* | The first selected gene shown in each ticked layer, side by side | See the same gene as raw counts vs normalized vs scaled - a teaching tool for what normalization does. |
| Assay *(if >1 assay)* | Which assay to read from | Compare RNA vs SCT vs integrated expression. |
| Sort genes by | Orders the gene menu: Alphabetical, or Highest / Lowest average expression | Find the most/least expressed genes without knowing names. |
| Gene(s) | Pick one or more genes | Each gene gets its own UMAP (in *genes* mode). |
| Data layer *(genes mode)* | Which layer to show (Normalized / Raw counts / Scaled) | See values before/after normalization. |
| Layers to show *(layers mode)* | Tick which layers to display for the first gene | Starts empty - tick the ones you want. (Rendering all at once is heavy; opt in to avoid slow-downs.) |
| Plots per row | Grid width (1-4) | Layout only. |
| Color scale → Same scale for all genes | All plots share one colorbar | Compare magnitude across genes (which is higher). |
| Color scale → Scale adjusted per gene | Each plot uses its own range | See the pattern even for weakly-expressed genes. |

> Genes missing from `scale.data` (only variable features are scaled) show a short
> message in that panel - expected, not an error.

Used for: seeing *where* a gene is expressed - confirming a cluster's identity from
a marker, or checking a gene of interest's activity across the map.

---

## 3. Cluster Markers

Genes that define each group - markers of one group vs *all other cells*. Two
versions depending on the object:

A) Precomputed (scdrake objects with `@misc$cluster_markers`):

| Option | Meaning |
|---|---|
| Clustering resolution | Which clustering the markers were computed on |
| Test type | Which precomputed statistical test |
| Cluster | Which cluster's markers to show |
| Search in / Search | Filter the table by a chosen column (e.g. a gene symbol) |

B) Computed live (any object):

| Option | Meaning / use |
|---|---|
| Grouping column | The metadata column that defines the groups - pick your cell-type column to get markers per cell type, or a clustering column for per-cluster markers |
| Assay *(if >1)* | Which assay to test on (`integrated` excluded - invalid for DE) |
| Data layer | Which layer the test reads (only the layers present in the assay are offered) |
| Group / cluster | Which group to find markers for (that group vs everyone else) |
| Test type | Statistical test (see the list below) |
| Search | Filter the results table |

Markers are auto-computed on selection. Table columns: `avg_log2FC` (how much higher
in this group), `pct.1` / `pct.2` (fraction of cells expressing it, in the group vs the
rest), `p_val` and `p_val_adj` (significance; very significant genes can display as `0`
when they underflow double precision).

Used for: annotating cell types / clusters - "what genes make this group
different?" → look them up to name it.

---

## 4. Contrasts (DEA)

Differential expression between two groups you choose (A vs B) - a head-to-head
comparison ("DEA" = Differential Expression Analysis). Precomputed or live:

A) Precomputed (scdrake `@misc$contrasts`): Contrast (the A vs B pair),
Test type, Search in / Search.

B) Computed live (any object):

| Option | Meaning / use |
|---|---|
| Grouping column | The metadata column to compare within |
| Assay *(if >1)* | Assay to test on (`integrated` excluded) |
| Data layer | Layer the test reads |
| Target (A) | The first group |
| Reference (B) | The group to compare against |
| Test type | Statistical test |
| Search | Filter results |

Difference from Cluster Markers: Markers = one group vs *everyone else*;
Contrasts = one group vs *one specific other group*.

Used for: condition comparisons (treated vs control), or one cell type vs another
directly.

---

## The "Test type" options (Markers & Contrasts)

All answer "which genes differ," with different statistics:

- Wilcoxon - the robust standard. Start here.
- Student's t - assumes roughly normal data; faster, less robust.
- ROC / AUC - ranks genes by how well they *classify* the group (reports an AUC
  "power", no p-value). Good for clean markers.
- Logistic regression (LR) - can account for covariates.
- Likelihood ratio (bimod) - a single-cell-specific test.
- Poisson (counts) - for raw UMI counts.

Standard tests (Wilcoxon/t/…) run on the normalized (`data`) layer; count models
(Poisson) run on counts. Installing [presto](https://github.com/immunogenomics/presto)
makes these substantially faster (Seurat uses it automatically). If an assay lacks a
`data` layer, the app normalizes on the fly so the tests still run.

---

## 5. Aggregated Expression

A summary heatmap + table: for each selected gene, one value per group (a
pseudobulk-style summary).

| Option | What it does | Information / use |
|---|---|---|
| Aggregate by | Metadata column defining the groups (heatmap columns) | e.g. cluster, cell type, condition |
| Sort genes by | Orders the gene menu | convenience |
| Gene(s) | Genes = rows of the heatmap | your genes of interest |
| Assay *(if >1)* | Which assay to read | choose the expression source |
| Data layer | Which layer to summarise (raw values) | Normalized / Raw counts / Scaled |
| Summary per group | Mean / Sum / Median per group | Mean = average level; Sum = total (pseudobulk); Median = robust middle |

- Heatmap: genes (rows) × groups (columns) on the blue → white → red scale, with the
  colour range pinned to the values (0→max for counts/data; symmetric for scale.data).
  Its title states the summary, gene count, grouping column and layer.
- Table: the same numbers (downloadable to Excel).

Used for: comparing a panel of genes across conditions/clusters at a glance -
"is this marker set higher in condition X?" - without scanning many separate UMAPs.

---

## 6. Split by Group

One gene shown as a separate expression UMAP per group - e.g. one panel per cell
type, condition, or sample. Each panel shows only that group's cells (coloured by the
gene's expression) over a grey backdrop of all cells, on a shared colour scale so
panels are directly comparable. Values are the raw values of the chosen data layer
(same colour meaning as Gene Expression).

| Option | What it does | Information / use |
|---|---|---|
| Split by | Metadata column whose groups become the panels | e.g. cell type, condition, sample |
| Then by (optional) | A second column to split again → a grid (e.g. cell type × condition) | leave "(none)" for a single split |
| Sort genes by | Orders the gene menu | find genes quickly |
| Gene | The single gene to show, split across the groups | one gene per view |
| Panels per row | Grid width | layout only |
| Assay *(if >1)* | Which assay to read | choose the expression source |
| Data layer | Which layer's raw values to show | Normalized / Raw counts / Scaled |
| Groups to show | Pick which "Split by" groups get a panel | starts empty - you choose which to compare |
| Groups to show (secondary) | Appears when "Then by" is set - pick which of *its* groups to include | only the combinations you pick are drawn |

- Panels are titled by the gene and group with the cell count, e.g.
  `Malat1 · Immature B cell (n=2,213)`. A combination with no cells (e.g. a KO×celltype
  that's unannotated) shows just the grey backdrop with `n=0`.
- Both pickers start empty - nothing renders until you add groups.
- Panel limit: at most 12 panels at once (each is a WebGL plot and browsers limit
  how many can be live). Pick the few groups you actually want to compare.

Used for: seeing how one gene's expression differs across groups - "is it higher
in KO than WT?", "which cell types express it?" - same gene each time, panels differ only
by which cells are shown.

---

## Quick reference - which tab for which question

| Question | Tab |
|---|---|
| What clusters / cell types are there? | Clusters / Metadata |
| Where is gene X expressed (and in which layer)? | Gene Expression |
| What defines this cluster / cell type? | Cluster Markers |
| How do groups A and B differ (gene table)? | Contrasts (DEA) |
| How does a gene panel compare across conditions? | Aggregated Expression |
| How does one gene differ across groups/cell types? | Split by Group |

---

## A note on colour scales

- Expression plots (Gene Expression, Aggregated, Split) use a diverging
  blue → white → red scale. For non-negative layers (counts/normalized) it runs
  0 → max; for scaled (z-score) data it's symmetric about 0 (white = average).
- Continuous metadata in Clusters / Metadata uses a viridis gradient.
- Categorical metadata uses distinct discrete colours with a toggleable legend.
