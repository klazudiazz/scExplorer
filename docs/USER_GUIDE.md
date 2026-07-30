# scExplorer — User Guide

A walkthrough of every tab and option in the app: what it does, what information
it gives you, and what you'd use it for.

---

## Loading data (the "Get started" screen)

- **Upload your own data** → drop or Browse a Seurat object `.rds`.
- **Use example data** → loads the bundled demo (~400 cells).

Once a valid object is loaded, the tabs appear. Which tabs show depends on the
object (see **Cluster Markers** / **Contrasts** below). The uploader is only
available *before* a dataset is loaded; to switch datasets, reload the page.

**Summary cards.** Just below the "Loaded…" line, four cards summarise the
object: **Cells**, **Genes**, **Clusters**, and **Samples**. Cells and Genes come
straight from the object; Clusters and Samples are detected from the object's
metadata (the clustering column and a sample/batch column) and show "-" when the
object has no such column.

---

## 1. Clusters / Metadata

The map of your cells. Each point is one cell, positioned by the UMAP embedding.

| Option | What it does | Info it gives / use |
|---|---|---|
| **Color UMAP by** | Colors cells by a chosen metadata column (clusters, cell type, batch, phase, condition…) | See the cluster structure and how any category is distributed across the map. Identify clusters, spot batch effects, see where a cell type or condition sits. |
| **Show categories** (checkboxes) | Toggle individual categories on/off (unchecked ones go grey) | Isolate one or a few groups — e.g. grey out everything except "T cells" to see only where they are. |

**Used for:** understanding the overall structure of the dataset — what
clusters/cell types exist and how they're arranged.

---

## 2. Gene Expression

One UMAP **per gene**, colored by that gene's **z-scored** expression on a
diverging blue → white → red scale: **blue = below the gene's average,
white = average, red = above average**. (Z-scoring centres each gene at 0, so
the colours are comparable in meaning across genes.)

| Option | What it does | Info it gives / use |
|---|---|---|
| **Sort genes by** | Orders the gene dropdown: Alphabetical, or by Highest / Lowest average expression | Quickly find the most (or least) expressed genes without knowing names. |
| **Gene(s)** | Pick one or more genes to plot | Each selected gene gets its own UMAP. |
| **Plots per row** | Grid width (1–4) | Layout only. |
| **Color scale → Same scale for all genes** | All plots share one colorbar (same ± range) | Compare **magnitude** between genes (which gene deviates more from its mean). |
| **Color scale → Scale adjusted per gene** | Each plot scales to its own ± max | See the **pattern** even for weakly-varying genes (each clearly visible). |

**Used for:** seeing *where* a gene is expressed relative to its average — e.g.
confirming a cluster is T cells because CD3 is high there, or checking your
gene of interest's activity.

---

## 3. Cluster Markers

Genes that **define each cluster** (markers of one group vs all other cells).
Two versions, depending on your object:

**A) If the object has scdrake's precomputed tables:**

| Option | Meaning |
|---|---|
| **Clustering resolution** | Which clustering to use |
| **Test type** | Which precomputed statistical test |
| **Cluster** | Which cluster's markers to show |
| **Search in / Search** | Filter the table by a column (e.g. type a gene symbol) |

**B) If not (computed live from any object):**

| Option | Meaning |
|---|---|
| **Grouping column** | Which metadata column defines the groups (clusters) |
| **Group / cluster** | Which group to find markers for |
| **Test type** | Statistical test (see the test list below) |
| **Search** | Filter the results table |

**Table columns (computed mode):** `avg_log2FC` (how much higher in this group),
`pct.1` / `pct.2` (fraction of cells expressing it in the group vs the rest),
`p_val_adj` (significance).

**Used for:** **annotating cell types** — "what genes make this cluster
different?" → look them up to name the cluster.

---

## 4. Contrasts

Genes that differ between **two groups you choose** (A vs B) — a head-to-head
comparison. Also precomputed (scdrake) or live:

| Option | Meaning |
|---|---|
| **Grouping column** | Which metadata column to compare within |
| **Target (A)** / **Reference (B)** | The two groups to compare |
| **Test type** | Statistical test |
| **Search** | Filter results |

**Difference from Cluster Markers:** Markers = one group vs *everyone else*.
Contrasts = one group vs *one specific other group*.

**Used for:** condition comparisons — treated vs control, or cluster 1 vs
cluster 2 directly.

---

## The "Test type" options (Markers & Contrasts)

All answer "which genes differ," but with different statistics:

- **Wilcoxon** — the standard, robust default. Start here.
- **Student's t** — assumes roughly normal data; faster, less robust.
- **ROC / AUC** — ranks genes by how well they *classify* the group (gives an AUC
  "power", no p-value). Good for finding clean markers.
- **Logistic regression (LR)** — can account for covariates.
- **Likelihood ratio (bimod)** — single-cell-specific test.
- **Poisson (counts)** — for raw UMI counts.

Installing [presto](https://github.com/immunogenomics/presto) makes these tests
substantially faster (Seurat uses it automatically).

---

## 5. Aggregated Expression

A **summary heatmap + table**: for each selected gene, one value per group
(pseudobulk-style).

| Option | What it does | Info / use |
|---|---|---|
| **Aggregate by** | Metadata column defining the groups (heatmap columns) | e.g. cluster, condition |
| **Sort genes by** | Orders the gene dropdown | convenience |
| **Gene(s)** | Genes = rows of the heatmap | pick your genes of interest |
| **Summary per group** | Mean / Sum / Median per group | Mean for average level; Sum for total (pseudobulk); Median for a robust middle |

- **Heatmap:** genes (rows) × groups (columns), coloured by the aggregated value
  on the same blue → white → red scale (here it runs low → high across the actual
  values — these are raw aggregates, not z-scored).
- **Table:** the same numbers.

**Used for:** comparing a set of genes across conditions/clusters **at a glance**
— "is this marker panel higher in condition X?" — without scanning many separate
UMAPs.

---

## 6. Split by Group

One gene shown as a **separate expression UMAP per group** — e.g. one panel per
cell type, condition, or sample. Each panel shows only that group's cells
(coloured by the gene's expression) over a grey backdrop of all cells, on a
**shared colour scale** so the panels are directly comparable.

The gene's expression is shown **z-scored** (centred at 0), with a diverging
blue → white → red scale: **blue = below average, white = average, red = above
average** for that gene.

| Option | What it does | Info / use |
|---|---|---|
| **Split by** | Metadata column whose groups become the panels | e.g. cell type, condition, sample |
| **Then by (optional)** | A second column to split again — turns the panels into a grid (e.g. cell type × condition) | leave as "(none)" for a single split |
| **Sort genes by** | Orders the gene dropdown (Alphabetical / Highest / Lowest average expression) | find genes quickly |
| **Gene** | The single gene to show, split across the groups | one gene per view |
| **Panels per row** | Grid width | layout only |
| **Groups to show** | Pick which of the "Split by" groups get a panel | starts empty — you choose the groups to compare |
| **Groups to show** (secondary) | Appears when "Then by" is set — pick which of *its* groups to include | only the combinations you pick are drawn |

- Panels are **titled by the group name** (both groups when "Then by" is set) with
  the cell count, e.g. `Immature B cell - WT (n=2,213)`. A combination with no
  cells (e.g. `KO - Plasma Cells (n=0)`) simply shows the grey backdrop.
- Both pickers start **empty**; nothing renders until you add groups.
- **Panel limit:** the tab shows **at most 12 panels at once** (each panel is a
  WebGL plot and browsers limit how many can be live) — the pickers won't let you
  exceed it. Pick the few groups you actually want to compare rather than all. To
  change the limit, edit `split_max_panels` in `inst/app/app.R`.

**Used for:** seeing how one gene's expression **differs across groups** —
"is this gene higher in the KO than the WT?", "which cell types express it?".
It's the same gene each time, so the panels differ only by which cells are shown.

---

## Quick reference — which tab for which question

| Question | Tab |
|---|---|
| What clusters / cell types are there? | Clusters / Metadata |
| Where is gene X expressed? | Gene Expression |
| What defines this cluster? | Cluster Markers |
| How do groups A and B differ (gene table)? | Contrasts |
| How does this gene panel compare across conditions? | Aggregated Expression |
| How does one gene's expression differ across groups/cell types? | Split by Group |
