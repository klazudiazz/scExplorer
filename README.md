# scExplorer

**Interactive explorer for Seurat objects.**

A Shiny app (`scExplorer`) that lets you upload any Seurat object `.rds` (or load
a bundled example) and explore it in the browser - no coding required. Built to
work with plain Seurat objects.

## Documentation

See the **[User Guide](docs/USER_GUIDE.md)** for a full walkthrough of every tab and
option: what each control does, the information it gives, and what you'd use it for.

## Features

- **Clusters / Metadata** - UMAP coloured by any metadata column (clusters,
  cell types, sample/batch, phase, QC scores...), with a resolution quick-pick,
  per-category show/hide checkboxes, and a Min/Max filter for numeric columns.
- **Gene Expression** - one UMAP per selected gene, with an **Assay** and
  **Data layer** (raw counts / normalized / scaled) selector, a shared or
  per-gene colour scale, and a "one gene across data layers" comparison mode.
- **Cluster Markers** - searchable marker-gene tables for one group vs the rest;
  precomputed (scdrake `@misc$cluster_markers`) or computed live from any object.
- **Contrasts (DEA)** - searchable differential-expression tables between two
  chosen groups; precomputed (`@misc$contrasts`) or computed live.
- **Aggregated Expression** - a heatmap + table summarising selected genes
  (Mean / Sum / Median) per level of a metadata column.
- **Split by Group** - one gene shown as a separate expression UMAP per group
  (optionally a two-way grid), on a shared scale.
- **Downloads** - every plot exports to **PNG** (or all panels into one PNG), and
  every table exports to **Excel** (`.xlsx`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("klazudiazz/scexplorer")
```

If the repository is **private**, pass a token:
`remotes::install_github("klazudiazz/scexplorer", auth_token = "ghp_...")`.

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
**shinycssloaders**.

For the **Cluster Markers** / **Contrasts** tabs (differential expression via
`Seurat::FindMarkers`), the app uses
**[presto](https://github.com/immunogenomics/presto)** - Seurat picks it up
automatically to make those tests substantially faster. Install it with:

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

- A dimensionality reduction named **`umap`** (`Embeddings(obj, "umap")`).
- One or more **`meta.data`** columns to colour/group by. The app automatically
  offers suitable categorical columns (clusters, cell types, sample/batch...) and
  numeric columns (shown as a gradient); very high-cardinality text columns
  (e.g. cell barcodes) are excluded.
- *(Optional)* `obj@misc$cluster_markers` / `obj@misc$contrasts` data frames
  (scdrake format) enable the precomputed marker/contrast tables - otherwise those
  tabs compute results live with `Seurat::FindMarkers`.
- *(Optional)* `SYMBOL` / `ENSEMBL` columns in the assay's feature metadata for
  nicer gene labels (falls back to the gene rowname otherwise).

## Customising

The app is a single, heavily-commented file: `inst/app/app.R`. Search it for the
tag **`# TWEAK:`** to find the values you're most likely to change - theme
colour, point sizes, which metadata columns appear, plot sizing, the expression
colour gradient, table columns, and the example-data path.

## Notes

- The bundled example lives at `inst/app/example_data/seurat_example.rds`
  (small enough to commit). All other `.rds` files are git-ignored.
- To regenerate the package help page after editing roxygen comments, run
  `devtools::document()`.

## License

MIT (c) 2026 Klaudia Zieba
