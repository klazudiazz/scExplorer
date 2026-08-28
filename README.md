# scExplorer

**Interactive explorer for Seurat objects.**

A Shiny app (`scExplorer`) that lets you upload any Seurat object `.rds` (or load
a bundled example) and explore it in the browser.

## Documentation

See the **[User Guide](docs/USER_GUIDE.md)** for a full walkthrough of every tab and
option: what each control does, the information it gives, and what you'd use it for.

## Features

The analysis is fixed to the **RNA assay**, **Normalized (log)** values and the
**Wilcoxon** test (falling back to SCT / the object's default assay if RNA has no
normalized layer).

- **Clusters / Metadata** - a chosen reduction (`umap`, `pca`, ...) coloured by any
  metadata column (clusters, cell types, sample/batch, phase, QC scores...), with
  per-category show/hide checkboxes, a Min/Max filter for numeric columns, and an
  optional **Split by group** view (one small-multiple panel per group).
- **Gene Expression** - two sub-tabs: **Per-gene maps** (one UMAP per gene, on a
  shared or per-gene colour scale) and **Aggregated (pseudobulk)** (a heatmap + table
  summarising genes Mean / Sum / Median per group).
- **Cluster Markers** - marker genes for one group vs the rest, computed live behind a
  Compute button; searchable, exportable table.
- **Marker Plot** - a dot plot (dot size = % of cells expressing, colour = z-scored
  mean expression) or violins of chosen genes across a grouping column.
- **Differential expression** - two chosen groups (target A vs reference B): a
  **pseudobulk DESeq2** analysis (sums raw counts per sample - the recommended way,
  and the default, for comparing conditions across samples) or a per-cell test.
- **Downloads** - every plot has a toolbar (the camera icon saves a PNG), multi-panel
  grids export all panels into one PNG, and every table exports to Excel (`.xlsx`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("klazudiazz/scExplorer")
```

If the repository is **private**, pass a token:
`remotes::install_github("klazudiazz/scExplorer", auth_token = "ghp_...")`.

## Usage

```r
library(scexplorer)
run_app()
```

`run_app()` starts the app and **opens it in your web browser** - it runs
locally, no server or internet connection needed. Then either upload a Seurat
`.rds` or click **"Use example data"**.

**If the app doesn't open on its own** (e.g. running R from a plain terminal or a
headless / conda environment with no default browser), look for the line

```
Listening on http://127.0.0.1:XXXX
```

in the R console and **copy that URL into your browser manually** - the app is
already running there. You can also disable the auto-open attempt with
`run_app(launch.browser = FALSE)`, or point R at a browser first, e.g.
`options(browser = "firefox")`.

## Dependencies

Core packages used by the app: **shiny**, **Seurat**, **plotly**, **DT**,
**shinycssloaders**. (Excel export is done in the browser by DT, so no extra Excel
package is needed.)

The **Differential expression** tab uses **pseudobulk DESeq2** by default, so
**DESeq2** (Bioconductor) is required to run it:

```r
# install.packages("BiocManager")
BiocManager::install("DESeq2")
```

Optional:

- **Cluster Markers / Marker Plot / Differential expression (per-cell)** run
  `Seurat::FindMarkers`. Installing
  **[presto](https://github.com/immunogenomics/presto)** makes those tests
  substantially faster (Seurat picks it up automatically):

  ```r
  remotes::install_github("immunogenomics/presto")
  ```

## Preparing a Seurat object (e.g. from CellRanger)

If you're starting from 10x CellRanger output (a `filtered_feature_bc_matrix`
folder with `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`), read it in
with Seurat, run the standard pipeline so the object has a UMAP, and save it as
`.rds`:

```r
library(Seurat)

counts <- Read10X(data.dir = "path/to/filtered_feature_bc_matrix")
seu <- CreateSeuratObject(counts = counts, project = "my_sample")

seu <- NormalizeData(seu)                 # or SCTransform(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)   # adds a clustering column
seu <- RunUMAP(seu, dims = 1:30)             # adds the required "umap" reduction

saveRDS(seu, "my_object.rds")
```

Then upload `my_object.rds` in the app.

## What your Seurat object needs

- At least one **dimensionality reduction** (`umap`, `tsne`, `pca`, ...). The
  map-based views plot on it; if several are present a "Reduction" dropdown lets you
  switch. A fresh object built straight from counts has none - run `RunPCA()` then
  `RunUMAP()` first.
- One or more **`meta.data`** columns to colour / group by. The app automatically
  offers categorical columns (clusters, cell types, sample/batch...) and numeric
  columns (shown as a gradient); very high-cardinality text columns (cell barcodes)
  are excluded.
- For gene expression and markers the app uses the **RNA** assay's **normalized
  (`data`)** layer. If RNA has only `counts`, run `NormalizeData(obj, assay = "RNA")`;
  otherwise the app falls back to the SCT / default assay.
- *(Optional)* `SYMBOL` / `ENSEMBL` columns in the assay's feature metadata for nicer
  gene labels (falls back to the gene rowname otherwise).

## Customising

The app is a single, heavily-commented file: `inst/app/app.R`. Search it for the
tag **`# TWEAK:`** to find the values you're most likely to change - theme
colour, point sizes, which metadata columns appear, plot sizing, the expression
colour gradient, table columns, and the example-data path.

## Notes

- The bundled example lives at `inst/app/example_data/seurat_example.rds`
  (~90 MB - GitHub accepts it but shows a large-file warning, and it is pulled on
  every `install_github`). All other `.rds` files are git-ignored.
- To regenerate the package help page after editing roxygen comments, run
  `devtools::document()`.

## Citation

scExplorer was inspired by **scdrake**, a single-cell RNA-seq analysis pipeline, and
it can read the cluster-marker and contrast tables that scdrake stores in a Seurat
object's `@misc` slot. 

- Kubovčiak J, Kolář M, Novotný J (2023). Scdrake: a reproducible and scalable
  pipeline for scRNA-seq data analysis. *Bioinformatics Advances* 3(1):vbad089.
  https://doi.org/10.1093/bioadv/vbad089

Tools and data scExplorer relies on:

**Software** (run `citation("<pkg>")` in R for the exact, version-correct reference):

- **Seurat** (v5) - the single-cell framework the app is built on. Hao Y, et al.
  (2024). Dictionary learning for integrative, multimodal and scalable single-cell
  analysis. *Nature Biotechnology* 42:293-304.
  https://doi.org/10.1038/s41587-023-01767-y
- **DESeq2** - Love MI, Huber W, Anders S (2014). Moderated estimation of fold change
  and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
  https://doi.org/10.1186/s13059-014-0550-8 (used for pseudobulk differential
  expression).
- **presto** - Korsunsky I et al. Fast Wilcoxon rank-sum test, which speeds up marker
  detection when installed. https://github.com/immunogenomics/presto
- **plotly** - Sievert C (2020). Interactive Web-Based Data Visualization with R,
  plotly, and shiny. Chapman and Hall/CRC. https://plotly-r.com
- **shiny** - Chang W, et al. (2026). shiny: Web Application Framework for R.
  R package version 1.14.0. https://shiny.posit.co/
- **DT** - Xie Y, Cheng J, Tan X, Aden-Buie G (2025). DT: A Wrapper of the JavaScript
  Library 'DataTables'. R package version 0.34.0. https://github.com/rstudio/DT

**Methods**

- Squair JW et al. (2021). Confronting false discoveries in single-cell differential
  expression. *Nature Communications* 12:5692.

**Example data**

The bundled example is the human pancreas dataset from:

- Baron M et al. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse
  Pancreas Reveals Inter- and Intra-cell Population Structure. *Cell Systems*
  3(4):346-360.

obtained via CZ CELLxGENE Discover (Chan Zuckerberg Initiative,
https://cellxgene.cziscience.com).
