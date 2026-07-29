# scexplorer

**Interactive explorer for Seurat objects.**

A Shiny app (`scExplorer`) that lets you upload any Seurat object `.rds` (or load
a bundled example) and explore it in the browser - no coding required. Built to
work with plain Seurat objects and with objects produced by the
[scdrake](https://github.com/bioinfocz/scdrake) pipeline.

## Features

- **Clusters / Metadata** - UMAP coloured by a metadata column (`batch`,
  `phase`, cell-type annotations, ...), with per-category show/hide filtering.
- **Gene Expression** - one UMAP per selected gene, with a shared or per-gene
  colour scale.
- **Cluster Markers** - searchable per-cluster marker-gene tables
  *(shown only if the object carries `@misc$cluster_markers`)*.
- **Contrasts** - searchable differential-expression tables
  *(shown only if the object carries `@misc$contrasts`)*.

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

## What your Seurat object needs

- A dimensionality reduction named **`umap`** (`Embeddings(obj, "umap")`).
- One or more **`meta.data`** columns to colour by. By default the app offers
  columns matching `^batch$`, `^phase$`, or `*_labels_pruned`; edit
  `useful_patterns` in `inst/app/app.R` to expose others.
- *(Optional)* `obj@misc$cluster_markers` / `obj@misc$contrasts` data frames
  (scdrake format) to enable the marker and contrast tabs - they auto-hide when
  absent.
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
