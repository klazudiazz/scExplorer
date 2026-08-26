# =============================================================================
# scExplorer - interactive Shiny app for exploring a Seurat object
# =============================================================================
#
# WHAT IT DOES
#   Load a Seurat .rds (upload or bundled example) and explore it in the browser:
#     * Clusters / Metadata : UMAP coloured by a metadata column, with a
#                             per-category show/hide filter.
#     * Gene Expression     : one UMAP per chosen gene, coloured by expression.
#     * Cluster Markers     : per-cluster marker-gene tables  (shown if present).
#     * Contrasts           : differential-expression tables  (shown if present).
#
# WHAT THE OBJECT NEEDS
#   * A reduction named "umap"  (accessed via Embeddings(obj, "umap")).
#   * meta.data columns to colour by. Every categorical column (factor /
#     character / logical) plus low-cardinality numerics (e.g. cluster numbers)
#     is offered; see `grouping_cols` in output$app_body to change the rule.
#   * Optional obj@misc$cluster_markers / obj@misc$contrasts data frames
#     (scdrake format) enable the last two tabs; they auto-hide when absent.
#   * Optional rowData SYMBOL/ENSEMBL columns give nicer gene labels.
#
# HOW A SHINY APP IS STRUCTURED (quick primer)
#   `ui`     describes the page (inputs + output placeholders).
#   `server` fills the outputs. An output$x <- render*({...}) block re-runs
#            automatically whenever an input$* or reactive() it reads changes.
#   reactive({...}) is a cached value shared by several outputs.
#   More: https://shiny.posit.co
#
# WHERE TO CUSTOMISE
#   Search this file for the tag  "# TWEAK:"  to jump to the things you're most
#   likely to want to change (theme colour, point sizes, which columns appear,
#   plot sizing, expression colour gradient, default genes-per-row, ...).
# =============================================================================

# ---- Packages ----
library(shiny)             # web-app framework
library(Seurat)            # single-cell object (Embeddings, FetchData, ...)
library(DT)                # interactive marker / contrast tables
library(shinycssloaders)   # the "Rendering plot..." loading spinners
library(plotly)            # interactive UMAPs + colorbars

# Raise the upload size limit (default is 5MB; Seurat objects might be much bigger)
options(shiny.maxRequestSize = 5 * 1024^3)  # 5 GB
# The marker/contrast/aggregated tables render CLIENT-side (renderDT server = FALSE)
# on purpose, so the DataTables "Download Excel" button exports ALL rows (not just
# the current page). DT warns when a client-side table is large -- silence that
# expected warning here.
options(DT.warn.size = FALSE)

# ---- Example data: read once per R session, reuse on every app (re)start ----
# The example .rds is large, so instead of readRDS()-ing it every time the app
# starts we cache the parsed object in the global environment. The first use in
# an R session reads it; any later app restart in the SAME session reuses the
# in-memory object (no re-read). A full R restart clears memory and re-reads once.
# WHERE "Use example data" LOADS FROM. The app tries these paths in order and
# uses the first that exists.
#
# >>> OTHER USERS: change the FIRST path below to the absolute location of the
# >>> example .rds on YOUR machine. Find it (in the bash terminal) with:
# >>>     find ~ -name "seurat_example*.rds"
# >>> If you leave it as-is it simply won't match on your system, and the app
# >>> falls back to the copy bundled in the package (2nd path) -- which is
# >>> usually exactly what you want.
.example_candidate_paths <- c(
  "/home/rstudio/inst/app/example_data/seurat_example.rds",                     # <-- CHANGE THIS per machine
  system.file("app/example_data/seurat_example.rds", package = "scexplorer"),   # bundled in the package (portable)
  "example_data/seurat_example.rds",                                            # relative to the app's working dir
  "/home/rstudio/scdrake_shiny_app/example_data/seurat_example.rds"             # legacy fallback
)

get_example_seu <- function() {
  if (exists(".seurat_example_obj", envir = globalenv(), inherits = FALSE)) {
    return(get(".seurat_example_obj", envir = globalenv()))
  }
  example_path <- .example_candidate_paths[file.exists(.example_candidate_paths)][1]
  if (is.na(example_path)) return(NULL)
  obj <- readRDS(example_path)
  assign(".seurat_example_obj", obj, envir = globalenv())
  obj
}

# Turn a meta.data column name into a friendly label for the UI: underscores ->
# spaces, first letter capitalised (e.g. "immgen_labels_pruned" -> "Immgen labels
# pruned"). Only the DISPLAY changes -- the selected value stays the real name.
pretty_label <- function(x) {
  x <- gsub("_", " ", x)
  sub("^(\\w)", "\\U\\1", x, perl = TRUE)
}

# Sort category levels for display: numerically when every level is a number
# (so 2 comes before 10, not 10 before 2), otherwise alphabetically. Used
# everywhere categorical levels are listed (legend, checkboxes, group pickers).
sort_levels <- function(x) {
  x <- unique(as.character(x))
  xn <- suppressWarnings(as.numeric(x))
  if (length(x) > 0 && all(!is.na(xn))) x[order(xn)] else sort(x)
}

# Choices for every "Sort genes by" dropdown (Gene Expression, Aggregated, Split,
# Marker Plot). ordered_gene_choices() interprets these values.
GENE_SORT_CHOICES <- c("Alphabetical", "Highest average expression", "Lowest average expression")

# Diverging blue -> white -> red colour scale used for expression on EVERY plot
# (Gene Expression, Aggregated Expression heatmap, Split by Group, and their
# colorbars). TWEAK: edit these stops to re-theme the expression gradient app-wide.
EXPR_COLORSCALE <- list(
  c(0.00, "#3B4CC0"),   # low  - blue
  c(0.25, "#8DB0FE"),
  c(0.50, "#F7F7F7"),   # mid  - white
  c(0.75, "#F49A7B"),
  c(1.00, "#B40426")    # high - red
)

# Colour limits for RAW layer values: non-negative layers (counts, data) start at
# 0; a layer with negatives (scale.data = z-scores) gets a symmetric range so 0
# maps to white on the blue-white-red scale.
expr_limits <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(c(0, 1))
  if (min(v) < 0) {
    lim <- max(abs(v)); if (lim <= 0) lim <- 1
    c(-lim, lim)
  } else {
    hi <- max(v); if (hi <= 0) hi <- 1
    c(0, hi)
  }
}

# Readable name for an assay layer
layer_label <- function(L) switch(L,
  data = "Normalized (log)", counts = "Raw counts", scale.data = "Scaled (z)", L)

# Labels for a "Data layer" dropdown
LAYER_LABELS <- c(data = "Normalized (log)", counts = "Raw counts", scale.data = "Scaled (z-score)")

# Named choices for a "Data layer" selector: the layers present in `assay_use`
# of `seu_obj`, in preferred order (data, counts, scale.data), labelled by
# LAYER_LABELS. Shared by every Data-layer dropdown in the app.
layer_choices <- function(seu_obj, assay_use) {
  lyrs <- tryCatch(SeuratObject::Layers(seu_obj[[assay_use]]), error = function(e) "data")
  ord <- intersect(c("data", "counts", "scale.data"), lyrs)
  if (length(ord) == 0) ord <- lyrs[1]
  setNames(ord, ifelse(ord %in% names(LAYER_LABELS), LAYER_LABELS[ord], ord))
}

# ---- Download helpers -------------------------------------------------------
# Tables: build a DT with a Excel download button. Uses the
# DataTables Buttons extension; the .xlsx is generated in the browser via JSZip
# (bundled with the DT package) -> no server-side Excel writer / no Java needed,
# and the file opens in Excel / LibreOffice / Google Sheets alike. `filename` is
# the .xlsx base name; `page = "all"` exports every row, not just the current page.
dt_with_excel <- function(df, filename) {
  DT::datatable(
    df, rownames = FALSE, extensions = "Buttons",
    options = list(
      dom = "Bfrtip",
      buttons = list(list(extend = "excel", text = "Download Excel",
                          filename = filename, title = NULL,
                          exportOptions = list(modifier = list(page = "all")))),
      pageLength = 15, searching = TRUE, scrollX = TRUE,
      processing = TRUE, language = list(processing = "Loading table...")
    )
  )
}

# Plots: a small "PNG" button that exports a plotly plot as a static PNG entirely
# in the browser (Plotly.downloadImage -> see the dlPlot() JS in tags$head).`id` = the plotlyOutput id, `fname` = PNG name.
safe_fname <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
plot_dl_btn <- function(id, fname) {
  tags$button(class = "png-dl-btn", type = "button", title = "Download PNG",
              onclick = sprintf("dlPlot('%s','%s')", id, safe_fname(fname)),
              HTML("&#8615; PNG"))
}
# "Download all as one PNG" button for a multi-panel grid: composites every
# listed panel (+ the shared colorbar, if given) into a single image laid out in
# `ncol` columns. `ids` = character vector of plot ids.
plot_dl_all_btn <- function(ids, ncol, colorbar_id = NULL, filename = "plots",
                            label = "Download all as one PNG") {
  ids_json <- jsonlite::toJSON(as.character(ids))   # always a JSON array
  cb <- if (is.null(colorbar_id)) "null" else sprintf("'%s'", colorbar_id)
  tags$button(class = "png-dl-all-btn", type = "button",
              onclick = sprintf("dlPlotsCombined(%s, %d, %s, '%s')",
                                ids_json, as.integer(ncol), cb, safe_fname(filename)),
              HTML(paste0("&#8615; ", label)))
}

# ---- In-app help ------------------------------------------------------------
# Link to the full user guide (edit if you fork / rename the repo).
GUIDE_URL <- "https://github.com/klazudiazz/scExplorer/blob/main/docs/USER_GUIDE.md"

# Plain-language tooltips for the terms users most often stumble on, shown as a
# small circled-i next to the control label (native browser tooltip -- no deps).
INFO_ASSAY <- "Assay = WHICH data matrix. Assays can hold different genes and mean different things: RNA (raw counts), SCT (normalized), integrated (batch-corrected; not valid for DE), ADT (protein), etc. Pick the assay first, then a data layer within it."
INFO_LAYER <- "Data layer = HOW the assay's values are transformed (same genes and cells): Normalized (log) is the usual expression, Raw counts are the integers, Scaled (z-score) is centered/scaled (variable genes only, has negatives)."
INFO_TEST  <- "Statistical test for differential expression. Wilcoxon (default) and t / LR / bimod run on normalized values; Poisson uses raw counts; ROC reports an AUC 'power' instead of a p-value."
INFO_RES   <- "A clustering column. The resolution sets HOW MANY clusters (higher = more, finer); the prefix (e.g. SCT vs integrated) is which assay's neighbor graph it was computed on. Colours the map by that clustering."
INFO_RED_CMP <- "The two embeddings (reductions) to place side by side on the SAME cells, e.g. the standard UMAP vs a foundation-model embedding (scGPT/Geneformer). Both panels share colours; highlighting a group lights it up in both."

# A small circled-i carrying a hover tooltip (native `title`, no JS).
info_icon <- function(text) tags$span(class = "info-ic", title = text, "ⓘ")
# A control label with the info icon appended (pass to a *Input's `label`).
lbl <- function(text, info) tagList(text, info_icon(info))

# A collapsible "How to use this tab" panel (native <details>, no JS), with a
# link to the full guide appended after the body.
help_details <- function(..., .summary = "How to use this tab", .guide = TRUE) {
  tags$details(class = "tab-help",
    tags$summary(.summary),
    tags$div(class = "tab-help-body", ...,
             if (.guide) tags$p(tags$a(href = GUIDE_URL, target = "_blank", "Open the full user guide ↗"))))
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Hide all browser scrollbars (deletes the gray bar on the right) */
      ::-webkit-scrollbar { width: 0px; background: transparent; }
      body { -ms-overflow-style: none; scrollbar-width: none; }

      /* Force hide Plotly's modebar (the gray line with options) */
      .plotly .modebar { display: none !important; }

      /* Restyle DataTables' own 'processing' indicator (the default blue pulsing
         dots) to match the app's green spinner. */
      @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
      .dataTables_processing, div.dt-processing {
        color: #333 !important;
        font-size: 14px !important;
        border: none !important;
        box-shadow: none !important;
        background: rgba(255, 255, 255, 0.9) !important;
        /* Sit below the header row (in the body area) instead of centred over
           the still-empty table, so it doesn't overlap the column names. */
        top: 130px !important;
        bottom: auto !important;
        margin-top: 0 !important;
      }
      /* hide DataTables' default animated bars/dots */
      .dataTables_processing > div, div.dt-processing > div { display: none !important; }
      /* green spinner in front of the text, matching the other loaders */
      .dataTables_processing::before, div.dt-processing::before {
        content: '';
        display: inline-block; width: 22px; height: 22px;
        margin-right: 10px; vertical-align: middle;
        border: 3px solid #eee; border-top: 3px solid #39FF14; border-radius: 50%;
        animation: spin 0.8s linear infinite;
      }

      .selectize-control { min-width: 220px; }
      .selectize-dropdown { min-width: 220px; }
      @media (max-width: 768px) {
        .sidebar { margin-bottom: 20px; }
      }

      /* Neon green accent theme.
         # TWEAK: to re-theme the app, replace the accent colour #39FF14 (and the
         darker hover #2ecc0f) everywhere in this <style> block. */
      .progress-bar {
        background-color: #39FF14 !important;
        font-weight: 600;
      }

      .shiny-file-input-progress { position: relative; }
      .shiny-file-input-progress .progress-bar { color: transparent !important; }
      .shiny-file-input-progress::after {
        content: 'Uploading data...';
        position: absolute;
        left: 0; right: 0; top: 0; bottom: 0;
        display: flex;
        align-items: center;
        justify-content: center;
        color: #1a1a1a;
        font-weight: 600;
        pointer-events: none;
      }
      .shiny-file-input-progress.upload-done::after {
        content: 'Uploaded the data. Reading...';
      }
      .shiny-file-input-progress.app-ready {
        display: none;
      }
      .nav-tabs > li.active > a, .nav-tabs > li.active > a:hover, .nav-tabs > li.active > a:focus {
        border-top: 3px solid #39FF14 !important;
      }
      .nav-tabs > li > a { color: #1a1a1a; }
      .nav-tabs > li > a:hover { color: #2ecc0f; }
      .btn-default:hover, .btn-default:focus { border-color: #39FF14; color: #2ecc0f; }
      .selectize-control.single .selectize-input.focus { border-color: #39FF14 !important; box-shadow: 0 0 0 1px #39FF14; }

      /* Per-plot 'Download PNG' button (top-right corner, always visible). */
      .png-dl-btn {
        position: absolute; top: 4px; right: 4px; z-index: 10;
        font-size: 11px; line-height: 1; padding: 3px 8px; cursor: pointer;
        background: rgba(255,255,255,0.9); color: #1a1a1a;
        border: 1px solid #39FF14; border-radius: 4px;
        opacity: 0.85; transition: opacity 0.15s, background 0.15s;
      }
      .png-dl-btn:hover { opacity: 1; background: #39FF14; }

      /* 'Download all PNGs' button shown above a multi-panel grid. */
      .png-dl-all-btn {
        font-size: 13px; padding: 5px 14px; cursor: pointer; margin-bottom: 10px;
        background: rgba(255,255,255,0.95); color: #1a1a1a;
        border: 1px solid #39FF14; border-radius: 4px;
      }
      .png-dl-all-btn:hover { background: #39FF14; }

      /* DataTables 'Download Excel' button, themed to match. */
      .dt-buttons .dt-button {
        background: rgba(255,255,255,0.9) !important; color: #1a1a1a !important;
        border: 1px solid #39FF14 !important; border-radius: 4px !important;
        padding: 4px 12px !important; margin-bottom: 8px !important; font-size: 13px !important;
      }
      .dt-buttons .dt-button:hover { background: #39FF14 !important; }

      /* Summary stat cards (Cells / Genes / Clusters / Samples) */
      .stat-cards { display: flex; gap: 12px; flex-wrap: wrap; margin: 6px 0 14px; }
      .stat-card {
        flex: 1 1 120px; min-width: 120px;
        background: #f7f7f7; border: 1px solid #e6e6e6; border-top: 3px solid #39FF14;
        border-radius: 8px; padding: 10px 16px;
      }
      .stat-card .stat-label {
        font-size: 12px; color: #8a8a8a; text-transform: uppercase; letter-spacing: 0.04em;
      }
      .stat-card .stat-value { font-size: 26px; font-weight: 700; color: #1a1a1a; line-height: 1.15; }

      /* Nested sub-tabs (scFM Explorer, Gene Expression): small rounded 'pill'
         chips, so they read as a compact secondary selector for the active tab,
         clearly smaller and different from the main tab bar. */
      .sub-tabs { margin: 4px 0 8px; }
      .sub-tabs .nav-pills { margin-bottom: 12px; }
      .sub-tabs .nav-pills > li > a {
        font-size: 14px; padding: 7px 18px; margin-right: 8px; color: #333 !important;
        background: #f2f2f2 !important; border: 1px solid #39FF14; border-radius: 4px;
      }
      .sub-tabs .nav-pills > li > a:hover {
        color: #2ecc0f !important; background: #eaeaea !important; border-color: #2ecc0f;
      }
      .sub-tabs .nav-pills > li.active > a,
      .sub-tabs .nav-pills > li.active > a:hover,
      .sub-tabs .nav-pills > li.active > a:focus {
        color: #1a1a1a !important; background: #fafafa !important;
        border: 2px solid #39FF14; font-weight: 600;
      }

      /* In-app help: circled-i tooltip icon + 'How to use this tab' panel */
      .info-ic { cursor: help; color: #2ecc0f; font-weight: 700; margin-left: 4px; font-size: 12px; }
      .tab-help { margin: 8px 0 12px; background: #fafafa; border: 1px solid #e6e6e6;
        border-left: 3px solid #39FF14; border-radius: 6px; padding: 6px 12px; }
      .tab-help > summary { cursor: pointer; font-weight: 600; color: #1a1a1a; outline: none; }
      .tab-help > summary:hover { color: #2ecc0f; }
      .tab-help-body { margin-top: 8px; color: #444; font-size: 13px; line-height: 1.5; }
      .tab-help-body p { margin: 4px 0; }
      .tab-help-body ol, .tab-help-body ul { margin: 6px 0; padding-left: 20px; }
      .tab-help-body ol li, .tab-help-body ul li { margin: 3px 0; }
      .tab-help-body a { color: #2ecc0f; font-weight: 600; }

      /* Assay / layer info strip */
      .assay-layer-strip { margin: 0 0 10px; padding: 8px 14px; background: #fafafa;
        border: 1px solid #e6e6e6; border-radius: 6px; font-size: 13px; color: #555; }
      .assay-layer-strip .assay-name { font-weight: 600; color: #1a1a1a; margin-right: 2px; }
      .assay-layer-strip .layer-badge { display: inline-block; padding: 1px 7px; margin: 0 2px;
        background: #e8f5e9; border: 1px solid #a5d6a7; border-radius: 4px;
        font-size: 11px; color: #2e7d32; font-family: monospace; }
      .assay-layer-strip .assay-block { display: inline-block; margin-right: 18px; }

      /* Make each spinner-wrapped feature plot fill its (square) grid cell */
      #dynamic_feature_plots .shiny-spinner-output-container { width: 100% !important; height: 100% !important; }
      #dynamic_feature_plots .html-widget { width: 100% !important; height: 100% !important; }
      #umap_plot_box .shiny-spinner-output-container { width: 100% !important; height: 100% !important; }
      #umap_plot_box .html-widget { width: 100% !important; height: 100% !important; }
      /* Split-by-group panels: make the plotly widget fill its fixed square box. */
      .umap-split-cell .shiny-spinner-output-container { width: 100% !important; height: 100% !important; }
      .umap-split-cell .html-widget { width: 100% !important; height: 100% !important; }
    ")),
    tags$script(HTML("
      document.addEventListener('DOMContentLoaded', function() {
        var observer = new MutationObserver(function() {
          document.querySelectorAll('.shiny-file-input-progress').forEach(function(container) {
            var bar = container.querySelector('.progress-bar');
            if (bar && bar.innerText && bar.innerText.trim() === 'Upload complete') {
              container.classList.add('upload-done');
            }
          });
          if (document.querySelector('.nav-tabs')) {
            document.querySelectorAll('.shiny-file-input-progress').forEach(function(container) {
              container.classList.add('app-ready');
            });
          }
        });
        observer.observe(document.body, { childList: true, subtree: true, characterData: true });

        // Once 'Use example data' is chosen, disable the 'Upload your own data'
        // option for the rest of the session (switching back is not supported;
        // reload the page to upload instead).
        document.addEventListener('change', function(e) {
          if (e.target && e.target.name === 'data_source' && e.target.value === 'example') {
            var up = document.querySelector('input[name=data_source][value=upload]');
            if (up) {
              up.disabled = true;
              var lbl = up.closest('label') || up.parentElement;
              if (lbl) { lbl.style.opacity = '0.5'; lbl.style.cursor = 'not-allowed'; }
            }
          }
        });
      });

      // Turn any string into a safe file base name.
      function pngSanitize(s) {
        return String(s).replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '') || 'plot';
      }
      // Read a plotly plot's on-screen title so the PNG name matches what the
      // user sees (falls back to the R-provided name if there's no title).
      function pngName(gd, fallback) {
        try {
          var t = gd && gd.layout && gd.layout.title;
          var txt = t ? (typeof t === 'string' ? t : t.text) : null;
          if (txt) { return pngSanitize(txt); }
        } catch (e) {}
        return pngSanitize(fallback);
      }
      // Indices of gradient traces (a numeric-coloured scattergl trace) whose
      // colorbar is currently hidden -- these are the panels that, in shared-scale
      // mode, rely on the ONE colorbar drawn to the side.
      function hiddenScaleTraces(gd) {
        var idx = [];
        (gd.data || []).forEach(function(tr, i) {
          if (tr.marker && tr.marker.colorscale && Array.isArray(tr.marker.color)
              && tr.marker.showscale === false) { idx.push(i); }
        });
        return idx;
      }
      // Export ONE plotly plot as a static PNG.
      // If the panel's own colorbar is hidden (shared-scale mode), turn it on just
      // for the capture and restore it after, so every downloaded PNG has a scale
      // without changing what's on screen.
      window.dlPlot = function(id, fallback) {
        var el = document.getElementById(id);
        if (!el) { return; }
        var gd = el.classList.contains('js-plotly-plot') ? el : el.querySelector('.js-plotly-plot');
        if (!gd || !window.Plotly) { return; }
        var opts = { format: 'png', filename: pngName(gd, fallback), scale: 2 };
        var hidden = hiddenScaleTraces(gd);
        if (hidden.length === 0) {                         // already has a colorbar
          window.Plotly.downloadImage(gd, opts);
          return;
        }
        Promise.resolve(window.Plotly.restyle(gd, { 'marker.showscale': true }, hidden))
          .then(function() { return window.Plotly.downloadImage(gd, opts); })
          .then(function() { window.Plotly.restyle(gd, { 'marker.showscale': false }, hidden); })
          .catch(function() { window.Plotly.restyle(gd, { 'marker.showscale': false }, hidden); });
      };
      // Export ALL displayed panels as ONE combined PNG: render each panel to an
      // image via Plotly (correct WebGL readback), then paint them onto a single
      // canvas in the grid layout (+ the shared colorbar on the right, if any).
      // ids = panel plot ids; ncol = panels per row; colorbarId = shared colorbar
      // id or null; filename = output base name.
      window.dlPlotsCombined = function(ids, ncol, colorbarId, filename) {
        if (!ids || !ids.length || !window.Plotly) { return; }
        var SCALE = 2, GAP = 20;
        function gdOf(id) {
          var el = document.getElementById(id); if (!el) { return null; }
          return el.classList.contains('js-plotly-plot') ? el : el.querySelector('.js-plotly-plot');
        }
        function snap(id) {
          var gd = gdOf(id); if (!gd) { return Promise.resolve(null); }
          var w = gd.clientWidth || 500, h = gd.clientHeight || 500;
          return window.Plotly.toImage(gd, { format: 'png', width: w, height: h, scale: SCALE });
        }
        function toImg(url) {
          return new Promise(function(res) {
            if (!url) { res(null); return; }
            var im = new Image();
            im.onload = function() { res(im); };
            im.onerror = function() { res(null); };
            im.src = url;
          });
        }
        Promise.all(ids.map(snap)).then(function(urls) { return Promise.all(urls.map(toImg)); })
          .then(function(panelImgs) {
            panelImgs = panelImgs.filter(Boolean);
            if (!panelImgs.length) { return; }
            var proceed = function(cbImg) {
              var n = panelImgs.length;
              var cols = Math.max(1, Math.min(ncol || 1, n));
              var rows = Math.ceil(n / cols);
              var cw = panelImgs[0].naturalWidth, ch = panelImgs[0].naturalHeight;
              var gridW = cols * cw + (cols - 1) * GAP, gridH = rows * ch + (rows - 1) * GAP;
              var cbW = cbImg ? (cbImg.naturalWidth + GAP) : 0;
              var canvas = document.createElement('canvas');
              canvas.width = gridW + cbW; canvas.height = gridH;
              var ctx = canvas.getContext('2d');
              ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, canvas.width, canvas.height);
              panelImgs.forEach(function(im, i) {
                var r = Math.floor(i / cols), c = i % cols;
                ctx.drawImage(im, c * (cw + GAP), r * (ch + GAP), cw, ch);
              });
              if (cbImg) { ctx.drawImage(cbImg, gridW + GAP, Math.max(0, (gridH - cbImg.naturalHeight) / 2)); }
              canvas.toBlob(function(blob) {
                if (!blob) { return; }
                var a = document.createElement('a');
                a.href = URL.createObjectURL(blob);
                a.download = (filename || 'plots') + '.png';
                document.body.appendChild(a); a.click();
                setTimeout(function() { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
              }, 'image/png');
            };
            if (colorbarId) { snap(colorbarId).then(toImg).then(proceed); }
            else { proceed(null); }
          });
      };
    "))
  ),

  titlePanel("scExplorer"),

  # "Get started" controls -- hidden once a dataset is successfully loaded
  # (output.data_ready is TRUE). conditionalPanel only hides them (CSS), so the
  # inputs keep their values.
  conditionalPanel(
    condition = "!output.data_ready",
    radioButtons("data_source", "Get started:",
                 choices = c("Upload your own data" = "upload",
                             "Use example data" = "example"),
                 selected = character(0)),
    conditionalPanel(
      condition = "input.data_source == 'upload'",
      fileInput("seurat_file", "Upload a Seurat object (.rds)", accept = ".rds")
    )
  ),

  # Once data is loaded, show a short status line in place of the controls.
  conditionalPanel(
    condition = "output.data_ready",
    uiOutput("data_status"),
    uiOutput("summary_cards")
  ),

  # The whole app body (tabs + plots) is built on the server side once data is
  # loaded, so nothing shows until a Seurat object is available.
  withSpinner(uiOutput("app_body"), caption = "Reading the data...", color = "#39FF14")
)

# =============================================================================
# SERVER: all reactive logic
# =============================================================================
server <- function(input, output, session) {

  # seu(): the loaded Seurat object. Everything downstream depends on this.
  # Reads either the uploaded file or the cached example, and validates it.
  seu <- reactive({
    req(input$data_source)
    if (input$data_source == "example") {
      obj <- tryCatch(get_example_seu(), error = function(e) {
        validate(need(FALSE, paste("Failed to read example dataset:", conditionMessage(e))))
      })
      validate(need(!is.null(obj), "Example dataset not found."))
    } else {
      req(input$seurat_file)
      obj <- tryCatch(readRDS(input$seurat_file$datapath), error = function(e) {
        validate(need(FALSE, paste("Failed to read uploaded file:", conditionMessage(e))))
      })
    }
    validate(need(inherits(obj, "Seurat"), "This file is not a Seurat object."))
    obj
  })

  # TRUE only once a valid Seurat object is loaded. Drives hiding the
  # "Get started" controls (and showing the status line) in the UI.
  output$data_ready <- reactive({
    !is.null(tryCatch(seu(), error = function(e) NULL))
  })
  outputOptions(output, "data_ready", suspendWhenHidden = FALSE)

  # Short status line shown in place of the controls once data is loaded.
  output$data_status <- renderUI({
    if (identical(input$data_source, "example")) {
      msg <- "Loaded example data"
    } else {
      req(input$seurat_file)
      msg <- paste0("Loaded file: ", input$seurat_file$name)
    }
    tags$p(tags$strong(msg))
  })

  # Summary stat cards shown under the status line once data is loaded. Cells and
  # Genes come straight from the object; Clusters / Samples are detected from the
  # first matching meta.data column (common names), and show "-" when absent so
  # the cards stay sensible for any Seurat object.
  # TWEAK: extend these candidate lists to match your own column names.
  output$summary_cards <- renderUI({
    req(seu())
    md <- seu()@meta.data
    distinct_of <- function(col) {
      if (is.na(col) || !col %in% colnames(md)) return(NA_integer_)
      length(unique(md[[col]][!is.na(md[[col]])]))
    }
    # Pick a column by exact name, else by a regex over the column names.
    pick_col <- function(exact, pattern) {
      col <- intersect(exact, colnames(md))[1]
      if (is.na(col)) {
        hit <- grep(pattern, colnames(md), ignore.case = TRUE, value = TRUE)
        col <- if (length(hit)) hit[1] else NA_character_
      }
      col
    }
    n_samples  <- distinct_of(pick_col(c("sample", "orig.ident", "samples", "batch", "donor", "patient"),
                                       "sample|batch|orig\\.ident"))

    fmt  <- function(x) if (is.na(x)) "-" else format(x, big.mark = ",")
    card <- function(label, value) {
      div(class = "stat-card",
          div(class = "stat-label", label),
          div(class = "stat-value", value))
    }
    # Assay / layer inventory
    all_assays <- tryCatch(SeuratObject::Assays(seu()), error = function(e) names(seu()@assays))
    assay_layer_info <- lapply(all_assays, function(a) {
      lyrs <- tryCatch(SeuratObject::Layers(seu()[[a]]), error = function(e) "?")
      tags$span(class = "assay-block",
                tags$span(class = "assay-name", paste0(a, ":")),
                lapply(lyrs, function(l) tags$span(class = "layer-badge", l)))
    })

    tagList(
      div(class = "stat-cards",
          card("Cells",   fmt(ncol(seu()))),
          card("Genes",   fmt(nrow(seu()))),
          card("Samples", fmt(n_samples))),
      div(class = "assay-layer-strip",
          tags$strong("Assays/layers: "), assay_layer_info)
    )
  })

  # gene_lookup(): maps each gene rowname to a nice display label
  # "SYMBOL (ENSEMBL)" when those rowData columns exist, else just the rowname.
  # Used to populate the searchable gene dropdown.
  gene_lookup <- reactive({
    req(seu())
    genes <- rownames(seu())
    assay_name <- DefaultAssay(seu())
    mf <- tryCatch(seu()[[assay_name]][[]], error = function(e) NULL)
    if (!is.null(mf) && all(c("SYMBOL", "ENSEMBL") %in% colnames(mf))) {
      sym <- as.character(mf$SYMBOL)
      ens <- as.character(mf$ENSEMBL)
      # Fall back to the gene's rowname whenever SYMBOL/ENSEMBL is missing (NA),
      # so a label never reads "NA (NA)".
      label <- ifelse(!is.na(sym) & !is.na(ens), paste0(sym, " (", ens, ")"),
                      ifelse(!is.na(sym), sym,
                             ifelse(!is.na(ens), ens, genes)))
      data.frame(value = genes, label = label, stringsAsFactors = FALSE)
    } else {
      data.frame(value = genes, label = genes, stringsAsFactors = FALSE)
    }
  })

  # gene_avg_expr(): mean expression per gene, used to offer the "sort genes by
  # highest/lowest average expression" options. Prefers the normalized "data"
  # layer of the default assay, but falls back to "counts" when "data" is missing
  # or empty (e.g. RNA in an SCT-only pipeline). 
  gene_avg_expr <- reactive({
    req(seu())
    try_layer <- function(L) {
      m <- suppressWarnings(tryCatch(GetAssayData(seu(), layer = L),
                                     error = function(e) NULL))
      if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) NULL else m
    }
    m <- try_layer("data")
    if (is.null(m)) m <- try_layer("counts")
    if (is.null(m)) return(setNames(rep(0, nrow(seu())), rownames(seu())))  # last resort
    Matrix::rowMeans(m)
  })

  # found_genes(): the subset of the user-selected genes that actually exist in
  # the object. shared_limits(): the pooled value range across those genes (of the
  # chosen data layer), used by BOTH the individual plots and the single shared
  # colorbar so their scales match in "Same scale for all genes" mode.
  found_genes <- reactive({
    req(seu(), input$genes)
    input$genes[input$genes %in% rownames(seu())]
  })

  # ---- Per-tab assay / layer selectors (shared factories) -------------------
  # Several tabs each have an "Assay" + "Data layer" pair with identical logic;
  # these factories build the backing reactives from just the input id, so each
  # tab is a one-liner instead of a repeated block.

  # Preferred default assay for viewing: RNA (log-normalized) when it has a `data`
  # layer -- the standard, interpretable choice for DE and cross-sample comparison;
  # else SCT (which usually has `data`); else the object's DefaultAssay. Requiring a
  # `data` layer keeps the app generic (an RNA assay with only counts falls through).
  preferred_view_assay <- function(seu_obj) {
    all_a <- tryCatch(SeuratObject::Assays(seu_obj), error = function(e) names(seu_obj@assays))
    has_data <- function(a) isTRUE("data" %in% tryCatch(SeuratObject::Layers(seu_obj[[a]]), error = function(e) character(0)))
    if ("RNA" %in% all_a && has_data("RNA")) return("RNA")
    if ("SCT" %in% all_a && has_data("SCT")) return("SCT")
    tryCatch(Seurat::DefaultAssay(seu_obj), error = function(e) all_a[1])
  }

  # Selected assay, defaulting to preferred_view_assay() until the input exists
  # (used by tabs that may view ANY assay: Gene Expression, Aggregated, Marker Plot).
  make_assay_val <- function(id) reactive({
    a <- input[[id]]
    all_a <- tryCatch(SeuratObject::Assays(seu()), error = function(e) names(seu()@assays))
    if (!is.null(a) && a %in% all_a) a else preferred_view_assay(seu())
  })

  # Selected assay for the DE tabs (Cluster Markers, Contrasts): restricted to
  # de_assay_choices() (which excludes `integrated`), defaulting to default_de_assay().
  make_de_assay_val <- function(id) reactive({
    a <- input[[id]]
    if (!is.null(a) && a %in% de_assay_choices()) a else default_de_assay()
  })

  # Selected data layer, defaulting to "data" until the input renders.
  make_layer_val <- function(id) reactive({
    l <- input[[id]]
    if (is.null(l) || !nzchar(l)) "data" else l
  })

  feature_assay_val <- make_assay_val("feature_assay")

  # The data layer for the Gene Expression plots is fixed to "data" (Normalized log);
  # the selector was removed, so make_layer_val falls back to "data".
  feature_layer <- make_layer_val("feature_layer")

  # Shared colour limits for "Same scale for all genes" mode: the pooled range of
  # the RAW values (of the chosen layer) across the selected genes.
  shared_limits <- reactive({
    fg <- found_genes()
    req(length(fg) > 0)
    assay_use <- feature_assay_val()
    lmat <- GetAssayData(seu(), layer = feature_layer(), assay = assay_use)   # fetch once
    m <- as.matrix(lmat[intersect(fg, rownames(lmat)), , drop = FALSE])
    expr_limits(as.vector(m))
  })

  # ---- Embedding (reduction) selection ----------------------------------------
  # The 2D+ reductions in the object (umap, tsne, pca, and foundation-model
  # embeddings like scgpt / geneformer if they were saved as reductions). Every
  # embedding plot in the app draws from the ONE selected here.
  available_reductions <- reactive({
    req(seu())
    reds <- tryCatch(SeuratObject::Reductions(seu()), error = function(e) names(seu()@reductions))
    reds[vapply(reds, function(r)
      tryCatch(ncol(Embeddings(seu(), r)) >= 2, error = function(e) FALSE), logical(1))]
  })

  # The chosen reduction, defaulting to "umap" (then the first available).
  reduction_val <- reactive({
    reds <- available_reductions()
    r <- input$reduction
    if (!is.null(r) && r %in% reds) r
    else if ("umap" %in% reds) "umap"
    else if (length(reds)) reds[1] else "umap"
  })

  # Memoised embedding matrix of the selected reduction -- every plot pulls from
  # this one reactive instead of re-extracting Embeddings() on each render.
  umap_emb <- reactive({
    req(seu())
    Embeddings(seu(), reduction_val())
  })

  # Axis titles for the current embedding (its own column names, e.g. UMAP_1 /
  # scgptumap_1), so plots are labelled for whatever reduction is shown.
  axis_labels <- reactive({
    cn <- colnames(umap_emb())
    if (length(cn) >= 2) cn[1:2] else c("Dim 1", "Dim 2")
  })

  # Display name of the current embedding for plot titles (UMAP/TSNE/PCA upper-cased,
  # anything else as-is, e.g. a foundation-model reduction name).
  embedding_label <- reactive({
    r <- reduction_val()
    if (tolower(r) %in% c("umap", "tsne", "pca")) toupper(r) else r
  })

  # One common square axis window (padded full extent of the UMAP embedding),
  # applied to BOTH axes of EVERY plot so all UMAPs are framed identically with
  # the same x/y range and scale.
  umap_range <- reactive({
    emb <- umap_emb()[, 1:2]
    r <- range(emb, na.rm = TRUE)
    pad <- diff(r) * 0.03
    c(r[1] - pad, r[2] + pad)
  })

  # A fixed, "nice" tick spacing for the UMAP axes (aiming for ~6 gridlines).
  # Pinning it makes every colouring draw gridlines at the SAME values, instead
  # of plotly auto-choosing different spacing per plot from the pixel width.
  umap_dtick <- reactive({
    raw <- diff(umap_range()) / 6
    mag <- 10^floor(log10(raw))
    nice <- c(1, 2, 2.5, 5, 10)
    mag * nice[which.min(abs(nice - raw / mag))]
  })

  # ordered_gene_choices(): the gene dropdown's choices in a given sort order --
  # "Alphabetical" (fast; no expression computation) or by average expression
  # (uses gene_avg_expr(), the slow rowMeans, only when actually requested).
  # Returns a named vector: names = display labels, values = gene rownames.
  ordered_gene_choices <- function(sort_mode) {
    lookup <- gene_lookup()
    if (identical(sort_mode, "Alphabetical")) {
      ordered_values <- sort(lookup$value)
    } else {
      avg_expr <- gene_avg_expr()   # only computed when sorting by expression
      ascending_order <- order(avg_expr, names(avg_expr), na.last = TRUE, method = "radix")
      genes_low_to_high <- names(avg_expr)[ascending_order]
      ordered_values <- if (identical(sort_mode, "Highest average expression")) {
        rev(genes_low_to_high)
      } else {
        genes_low_to_high
      }
    }
    lookup_ordered <- lookup[match(ordered_values, lookup$value), ]
    setNames(lookup_ordered$value, lookup_ordered$label)
  }

  # Keep each tab's gene dropdown in sync with its own sort selector. (Server-side
  # selectize so it stays fast even with tens of thousands of genes.)
  observeEvent(list(seu(), input$gene_sort), {
    req(seu(), input$gene_sort)
    updateSelectizeInput(session, "genes",
                         choices = ordered_gene_choices(input$gene_sort), server = TRUE)
  })
  observeEvent(list(seu(), input$agg_gene_sort), {
    req(seu(), input$agg_gene_sort)
    updateSelectizeInput(session, "agg_genes",
                         choices = ordered_gene_choices(input$agg_gene_sort), server = TRUE)
  })

  # Assays available for differential expression (FindMarkers).
  # Excludes 'integrated' (and any assay starting with 'integrated') because DE on
  # batch-integrated data is statistically invalid (it's batch-reconstructed
  # and loses true single-cell variance structure; Seurat warns against it).
  de_assay_choices <- reactive({
    req(seu())
    all_a <- tryCatch(SeuratObject::Assays(seu()), error = function(e) names(seu()@assays))
    valid_a <- all_a[!grepl("^integrated", tolower(all_a))]
    if (length(valid_a) == 0) valid_a <- all_a
    valid_a
  })

  default_de_assay <- reactive({
    req(seu())
    choices <- de_assay_choices()
    # Prefer RNA (log-normalized) as the default DE assay, else SCT; excludes
    # integrated already. prepare_obj_for_de auto-normalizes RNA if it lacks `data`.
    if ("RNA" %in% choices) return("RNA")
    if ("SCT" %in% choices) return("SCT")
    def <- tryCatch(Seurat::DefaultAssay(seu()), error = function(e) NULL)
    if (!is.null(def) && def %in% choices) def else choices[1]
  })

  gmarker_assay_val <- make_de_assay_val("gmarker_assay")
  gcontrast_assay_val <- make_de_assay_val("gcontrast_assay")

  # Clustering / resolution columns in the metadata (e.g. SCT_snn_res.0.1,
  # seurat_clusters, final_clusters). Shared by the "Resolution / clustering"
  # quick-pick and the two-way sync with "Color UMAP by" below.
  clust_cols_r <- reactive({
    req(seu())
    grep("snn_res|_res\\.[0-9]|cluster", colnames(seu()@meta.data), ignore.case = TRUE, value = TRUE)
  })

  # app_body(): builds the whole tabbed UI once data is loaded. The
  # tabs Cluster Markers and Contrasts are only added when the object carries the
  # matching @misc data, so plain Seurat objects simply won't show them.
  output$app_body <- renderUI({
    req(input$data_source)
    # If "upload" is selected but no file has been chosen yet, show nothing.
    # This clears any previously-loaded example view (so a fresh upload isn't
    # masked by stale tabs) and removes the tabs, which also lets the file-input
    # progress bar reappear for the new upload.
    if (input$data_source == "upload" && is.null(input$seurat_file)) return(NULL)
    req(seu())
    all_cols <- colnames(seu()@meta.data)
    # Which meta.data columns appear in "Color UMAP by": every column that makes
    # sense to colour by categorically - factor / character / logical columns,
    # plus low-cardinality numerics (e.g. cluster numbers). Continuous numerics
    # with many distinct values (nCount, nFeature, doublet_score, percent.mt, ...)
    # are excluded, because the UMAP assigns one colour per value and those would
    # produce thousands of categories.
    # Grouping columns (Split / Aggregate / Markers / Contrasts) = REAL categories:
    # factor/character/logical or low-cardinality numeric, and NOT a continuous
    # column. Capped at `max_categories` distinct values so barcode / cell-name
    # columns (thousands of levels) are excluded - they'd make thousands of
    # panels/heatmap columns and trigger Shiny's "large number of options" warning.
    # TWEAK: raise/lower max_categories to change what counts as a grouping column.
    max_categories <- 100
    grouping_cols <- all_cols[vapply(all_cols, function(col) {
      if (col %in% clust_cols_r()) return(TRUE)                     # clusterings: always a category, any count
      x <- seu()@meta.data[[col]]
      if (is_continuous_col(col)) return(FALSE)                     # a score, not a category
      (is.factor(x) || is.character(x) || is.logical(x) || is.numeric(x)) &&
        length(unique(x[!is.na(x)])) <= max_categories
    }, logical(1))]
    grouping_choices <- setNames(grouping_cols, pretty_label(grouping_cols))
    grouping_choices <- grouping_choices[order(names(grouping_choices))]

    # "Color UMAP by" = the grouping columns PLUS continuous columns (numeric, or
    # numeric-like scores) which are drawn as a gradient (see is_continuous_col).
    # Still excludes huge text columns (barcodes) so the UMAP can't get thousands
    # of discrete colours.
    color_cols <- all_cols[vapply(all_cols, function(col) {
      if (col %in% clust_cols_r()) return(TRUE)                     # clusterings are categorical colourings (offered here now)
      x <- seu()@meta.data[[col]]
      if (is_continuous_col(col)) return(TRUE)                      # gradient
      (is.factor(x) || is.character(x) || is.logical(x) || is.numeric(x)) &&
        length(unique(x[!is.na(x)])) <= max_categories
    }, logical(1))]
    color_choices <- setNames(color_cols, pretty_label(color_cols))
    color_choices <- color_choices[order(names(color_choices))]

    # Embeddings (reductions) to offer -- only shown when the object has >1, so
    # the picker appears exactly when there's a choice (e.g. umap + a
    # foundation-model embedding). Applies to EVERY embedding plot in the app.
    red_choices <- available_reductions()

    # Assay layers present (counts / data / scale.data) for the Gene Expression
    # "Data layer" toggle -- view raw values before/after normalization.
    # All assays (for the Gene Expression assay selector).
    all_assays_ge <- tryCatch(SeuratObject::Assays(seu()), error = function(e) names(seu()@assays))
    def_assay_ge <- preferred_view_assay(seu())

    has_markers <- !is.null(seu()@misc$cluster_markers)
    has_contrasts <- !is.null(seu()@misc$contrasts)

    de_assays <- de_assay_choices()
    def_de_assay <- default_de_assay()

    # DE tests offered by the generic Markers/Contrasts tabs.
    # All built into Seurat and fast enough to auto-compute. (poisson uses raw
    # counts; roc reports an AUC power instead of a p-value.)
    # NOTE: the slow per-gene tests (negbinom, MAST, DESeq2) were removed for now
    # -- they take minutes and these tabs auto-compute. Re-add here if wanted
    # (the counts-slot logic below already handles negbinom/DESeq2).
    de_test_choices <- c("Wilcoxon" = "wilcox", "Student's t" = "t", "ROC / AUC" = "roc",
                         "Logistic regression" = "LR", "Likelihood ratio (bimod)" = "bimod",
                         "Poisson (counts)" = "poisson")

    tabs <- list(
      tabPanel(
        "Clusters / Metadata",
        help_details(
          tags$p("Colours the UMAP (one point per cell) by any metadata column, so you can see how cell types, conditions, clusters, or scores are laid out."),
          tags$ol(
            tags$li("Pick a column in \"Color UMAP by\" to colour the map (clusterings are in this list too)."),
            tags$li("Choose which embedding to plot on in \"Reduction (embedding)\" (for example umap or pca), when the object has more than one."),
            tags$li("For a categorical column, tick or untick the categories in the sidebar to show or hide them."),
            tags$li("For a numeric column (scores, QC), type a Min or Max to grey out cells outside that range."),
            tags$li("To break the map out per group, tick \"Split by group\", pick a column, then choose which groups in \"Groups to show\" (starts empty; one panel per chosen group over a grey backdrop, up to 12)."),
            tags$li("Hover a point to read its value; use the PNG button to save the plot."))),
        # Controls stacked vertically: Reduction, then Color UMAP by + its category
        # filter, then the Split-by option (revealed by a checkbox). Plot renders below.
        wellPanel(
          if (length(red_choices) > 1)
            fluidRow(column(4, selectInput("reduction", "Reduction:", choices = red_choices,
                                           selected = if ("umap" %in% red_choices) "umap" else red_choices[1], width = "100%"))),
          fluidRow(column(4, selectInput("group_var", "Color UMAP by:", choices = color_choices, width = "100%"))),
          uiOutput("category_filter_ui"),
          fluidRow(column(6, checkboxInput("split_on", "Split by group", value = FALSE))),
          conditionalPanel(
            condition = "input.split_on",
            fluidRow(column(4, selectInput("split_by", "Split by:", choices = grouping_choices, width = "100%"))),
            fluidRow(column(6, uiOutput("umap_split_levels_ui")))
          )
        ),
        # Not splitting: the single interactive UMAP (category toggles restyle it in place).
        conditionalPanel(
          condition = "!input.split_on",
          div(style = "display: flex; justify-content: center; width: 100%; overflow: hidden;",
              div(id = "umap_plot_box", class = "plot-dl-wrap", style = "position: relative; width: calc(min(100vh - 340px, 60vw) + 260px); height: min(100vh - 340px, 60vw);",
                  plot_dl_btn("umap_plot", "UMAP"),
                  withSpinner(plotlyOutput("umap_plot", width = "100%", height = "100%"), caption = "Rendering plot...", color = "#39FF14")))
        ),
        # Splitting: one UMAP panel per chosen group (up to 12).
        conditionalPanel(
          condition = "input.split_on",
          uiOutput("umap_split_grid")
        )
      ),
      tabPanel(
        "Gene Expression",
        help_details(
          tags$p("Explore gene expression two ways, using the sub tabs below (the smaller green buttons). Both use the RNA assay and Normalized (log) values. Each sub tab has its own \"How to use\" panel with steps."),
          tags$ul(
            tags$li("Per-gene maps: one UMAP per gene, coloured by expression."),
            tags$li("Aggregated (pseudobulk): the selected genes summarised per group as a heatmap and table."))),
        div(class = "sub-tabs", tabsetPanel(
          type = "pills",
          tabPanel(
            "Per-gene maps",
            help_details(.summary = "How to use: Per-gene maps", .guide = FALSE,
              tags$p("One UMAP per gene, coloured by expression on a blue to white to red scale. Uses the RNA assay and Normalized (log) values."),
              tags$ol(
                tags$li("Type gene symbols in the \"Gene(s):\" box."),
                tags$li("Set \"Plots per row\" and \"Color scale\" (same scale for all genes to compare, or adjusted per gene)."),
                tags$li("Save a panel with its PNG button, or all panels with \"Download all as one PNG\"."))),
        wellPanel(
          fluidRow(
            column(4, selectInput("gene_sort", "Sort genes by:", width = "100%",
                                  choices = GENE_SORT_CHOICES)),
            column(8, selectizeInput("genes", "Gene(s):", choices = NULL, multiple = TRUE, width = "100%",
                                     options = list(placeholder = "Type one or more gene symbols...", sortField = list(list(field = "$order")))))
          ),
          fluidRow(
            # TWEAK: default / allowed plots-per-row for the Gene Expression grid.
            column(2, numericInput("feature_ncol", "Plots per row:", value = 2, min = 1, max = 4, width = "100%")),
            column(7, radioButtons("scale_mode", "Color scale:",
                                   choices = c("Same scale for all genes (comparable)" = "shared",
                                               "Scale adjusted per gene (all clearly visible)" = "per_gene"),
                                   selected = "shared"))
          )
        ),
        uiOutput("dynamic_feature_plots")   # (centres its own grid internally)
          ),
          # ---- Aggregated (pseudobulk): heatmap + table, nested under Gene Expression ----
          tabPanel(
            "Aggregated (pseudobulk)",
            help_details(.summary = "How to use: Aggregated (pseudobulk)", .guide = FALSE,
              tags$p("Summarises the selected genes per group as a heatmap and a table, a compact pseudobulk overview. Uses the RNA assay and Normalized (log) values."),
              tags$ol(
                tags$li("Pick the \"Aggregate by\" column (the groups)."),
                tags$li("Type genes in the \"Gene(s):\" box."),
                tags$li("Choose the \"Summary per group\" (mean, sum, or median)."),
                tags$li("Save the heatmap with the PNG button, or the table with \"Download Excel\"."))),
            wellPanel(
              fluidRow(
                column(3, selectInput("agg_group_var", "Aggregate by:", choices = grouping_choices, width = "100%")),
                column(3, selectInput("agg_gene_sort", "Sort genes by:", width = "100%", choices = GENE_SORT_CHOICES)),
                column(6, selectizeInput("agg_genes", "Gene(s):", choices = NULL, multiple = TRUE, width = "100%",
                                         options = list(placeholder = "Type one or more gene symbols...", sortField = list(list(field = "$order")))))
              ),
              fluidRow(
                column(6, radioButtons("agg_fun", "Summary per group:",
                                       choices = c("Mean" = "mean", "Sum" = "sum", "Median" = "median"),
                                       selected = "mean", inline = TRUE))
              )
            ),
            div(class = "plot-dl-wrap", style = "position: relative;",
                plot_dl_btn("agg_heatmap", "aggregated_expression"),
                withSpinner(plotlyOutput("agg_heatmap", height = "600px"), caption = "Rendering plot...", color = "#39FF14")),
            tags$hr(),
            withSpinner(DTOutput("agg_table"), caption = "Loading table...", color = "#39FF14")
          )
        ))
      )
    )

    # Cluster Markers tab. If the object carries scdrake's precomputed
    # @misc$cluster_markers we show those (instant). Otherwise we compute markers
    # on demand from any Seurat object (FindMarkers, behind a Run button).
    if (has_markers) {
      tabs <- c(tabs, list(tabPanel(
        "Cluster Markers",
        help_details(
          tags$p("Precomputed marker tables: genes that mark each cluster versus all other cells (Wilcoxon test)."),
          tags$ol(
            tags$li("Pick a \"Clustering resolution\" (which clustering to use)."),
            tags$li("Pick a cluster to view its markers."),
            tags$li("Type in the table's search box to filter; \"Download Excel\" exports every row."))),
        wellPanel(
          fluidRow(
            column(3, selectInput("marker_resolution", "Clustering resolution:", width = "100%",
                                  choices = unique(seu()@misc$cluster_markers$source_column))),
            column(2, uiOutput("marker_cluster_ui"))
          )
        ),
        withSpinner(DTOutput("marker_table"), caption = "Loading table...", color = "#39FF14")
      )))
    } else {
      tabs <- c(tabs, list(tabPanel(
        "Cluster Markers",
        help_details(
          tags$p("Finds genes that mark one group versus all other cells, computed live from the object. Uses the RNA assay, Normalized (log) values, and the Wilcoxon test."),
          tags$ol(
            tags$li("Pick the \"Grouping column\" and the group whose markers you want."),
            tags$li("Click \"Compute markers\". It is heavy on large objects, so it does not run automatically."),
            tags$li("Type in the table's search box to filter; \"Download Excel\" exports every row."))),
        wellPanel(
          fluidRow(
            column(2, selectInput("gmarker_group", "Grouping column:", choices = grouping_choices, width = "100%")),
            column(2, uiOutput("gmarker_cluster_ui"))
          ),
          fluidRow(
            column(3, div(style = "margin-top: 5px;",
                          actionButton("gmarker_compute", "Compute markers",
                                       style = "background-color: #fafafa; color: #000000; border: 3px solid #39FF14; font-weight: 600;")))
          )
        ),
        withSpinner(DTOutput("gmarker_table"), caption = "Computing markers...", color = "#39FF14")
      )))
    }

    # Marker Plot tab: a dot plot of chosen genes across a grouping column -- dot
    # size = % of cells expressing, colour = average expression (z-scored per gene).
    # The classic "markers per cell type" visualisation; pair it with the marker
    # tables (copy top genes from Cluster Markers into the gene picker here).
    tabs <- c(tabs, list(tabPanel(
      "Marker Plot",
      help_details(
        tags$p("Compares chosen genes across groups, as a dot plot or violins. Dot plot: dot size is the percent of cells in the group expressing the gene, colour is mean expression z scored per gene. Violin: the expression distribution per group, one row per gene. Uses the RNA assay and Normalized (log) values."),
        tags$ol(
          tags$li("Choose a \"Plot type\": dot plot or violin."),
          tags$li("Pick the \"Group by\" column."),
          tags$li("Dot plot: choose genes manually, or pick \"Top markers per group (auto)\" and click \"Compute markers\" (Wilcoxon). Violin: type the genes to show."),
          tags$li("Use the PNG button to save the plot."))),
      wellPanel(
        # Plot type is the FIRST thing to choose; everything below adapts to it.
        fluidRow(
          column(12, radioButtons("marker_plot_type", "Plot type:", inline = TRUE,
                                  choices = c("Dot plot" = "dot", "Violin" = "violin"), selected = "dot"))
        ),
        fluidRow(
          column(4, selectInput("dot_group_var", "Group by:", choices = grouping_choices, width = "100%"))
        ),
        # ---- DOT PLOT options: auto top markers, or a manual gene set ----
        conditionalPanel(
          condition = "input.marker_plot_type == 'dot'",
          fluidRow(
            column(12, radioButtons("dot_gene_mode", "Genes:", inline = TRUE,
                                    choices = c("Top markers per group (auto)" = "auto",
                                                "Choose manually" = "manual"), selected = "auto"))
          ),
          # AUTO: find each group's top markers (gated behind a button -- it computes
          # markers for EVERY group).
          conditionalPanel(
            condition = "input.dot_gene_mode == 'auto'",
            fluidRow(
              column(3, numericInput("dot_top_n", "Top markers / group:", value = 3, min = 1, max = 15, width = "100%")),
              column(3, div(style = "margin-top: 25px;",
                            actionButton("dot_compute", "Compute markers",
                                         style = "background-color: #fafafa; color: #000000; border: 3px solid #39FF14; font-weight: 600;")))
            )
          ),
          # MANUAL: pick genes yourself (e.g. from the Cluster Markers table).
          conditionalPanel(
            condition = "input.dot_gene_mode == 'manual'",
            fluidRow(
              column(3, selectInput("dot_gene_sort", "Sort genes by:", width = "100%",
                                    choices = GENE_SORT_CHOICES)),
              column(9, selectizeInput("dot_genes", "Gene(s):", choices = NULL, multiple = TRUE, width = "100%",
                                       options = list(placeholder = "Type marker genes (e.g. from Cluster Markers)...", sortField = list(list(field = "$order")))))
            )
          )
        ),
        # ---- VIOLIN options: just pick the genes to show ----
        conditionalPanel(
          condition = "input.marker_plot_type == 'violin'",
          fluidRow(
            column(3, selectInput("violin_gene_sort", "Sort genes by:", width = "100%",
                                  choices = GENE_SORT_CHOICES)),
            column(9, selectizeInput("violin_genes", "Genes to show (violin):",
                                     choices = NULL, multiple = TRUE, width = "100%",
                                     options = list(placeholder = "Type or pick genes to show as violins...", sortField = list(list(field = "$order")))))
          )
        )
      ),
      conditionalPanel(
        condition = "input.marker_plot_type == 'dot'",
        div(class = "plot-dl-wrap", style = "position: relative;",
            plot_dl_btn("dot_plot", "marker_dotplot"),
            withSpinner(plotlyOutput("dot_plot", height = "600px"),
                        caption = "Computing markers / rendering...", color = "#39FF14"))
      ),
      conditionalPanel(
        condition = "input.marker_plot_type == 'violin'",
        div(class = "plot-dl-wrap", style = "position: relative; width: 100%;",
            plot_dl_btn("violin_plot", "marker_violin"),
            withSpinner(uiOutput("violin_plot_ui"),
                        caption = "Computing markers / rendering...", color = "#39FF14"))
      )
    )))

    # Contrasts tab. Precomputed scdrake @misc$contrasts if present, otherwise
    # differential expression between two chosen groups computed on demand.
    if (has_contrasts) {
      tabs <- c(tabs, list(tabPanel(
        "Differential expression",
        help_details(
          tags$p("Precomputed differential expression tables between two groups (target versus reference). Unlike cluster markers (a group versus everyone else), a contrast is one group versus one specific other group."),
          tags$ol(
            tags$li("Pick a \"Contrast\" (the A vs B pair)."),
            tags$li("Type in the table's search box to filter; \"Download Excel\" exports every row (Wilcoxon test)."))),
        wellPanel(
          fluidRow(
            column(4, selectInput("contrast_choice", "Contrast:", width = "100%",
                                  choices = unique(paste0(seu()@misc$contrasts$target, " vs ", seu()@misc$contrasts$reference))))
          )
        ),
        withSpinner(DTOutput("contrast_table"), caption = "Loading table...", color = "#39FF14")
      )))
    } else {
      tabs <- c(tabs, list(tabPanel(
        "Differential expression",
        help_details(
          tags$p("Finds genes that differ between two groups you choose (target A versus reference B), for example one condition versus another. Unlike cluster markers (a group versus everyone else), this is one group versus one specific other group."),
          tags$ol(
            tags$li("Choose a \"Method\": \"Per cell test\" is fast and exploratory; \"Pseudobulk DESeq2\" is recommended for comparing conditions across samples, because per cell tests inflate false positives there."),
            tags$li("Pick the \"Grouping column\" that holds your two groups, then the target group (A) and the reference group (B)."),
            tags$li("Per cell test uses the RNA assay, Normalized (log) values, and the Wilcoxon test. Pseudobulk: set the \"Sample / replicate column\" (it needs at least two samples per group), and optionally \"Restrict to\" one cell type; it sums raw RNA counts and runs DESeq2."),
            tags$li("Click \"Compute contrast\" (or \"Compute pseudobulk DE\"). It does not run automatically."),
            tags$li("Type in the table's search box to filter; \"Download Excel\" exports every row."))),
        wellPanel(
          fluidRow(column(12, radioButtons("gcontrast_method", "Method:", inline = TRUE,
                                  choices = c("Per-cell test (exploratory)" = "percell",
                                              "Pseudobulk DESeq2 (recommended for conditions)" = "pseudobulk"),
                                  selected = "percell"))),
          fluidRow(
            column(4, selectInput("gcontrast_group", "Grouping column:", choices = grouping_choices, width = "100%")),
            column(4, uiOutput("gcontrast_target_ui")),
            column(4, uiOutput("gcontrast_reference_ui"))
          ),
          # ---- Per-cell test options ----
          conditionalPanel(
            condition = "input.gcontrast_method == 'percell'",
            fluidRow(
              column(3, div(style = "margin-top: 5px;",
                            actionButton("gcontrast_compute", "Compute contrast",
                                         style = "background-color: #fafafa; color: #000000; border: 3px solid #39FF14; font-weight: 600;")))
            )
          ),
          # ---- Pseudobulk (DESeq2) options ----
          conditionalPanel(
            condition = "input.gcontrast_method == 'pseudobulk'",
            fluidRow(
              column(3, selectInput("gcontrast_pb_sample", "Sample / replicate column:", choices = grouping_choices, width = "100%")),
              column(3, selectInput("gcontrast_pb_subset_col", "Restrict to (optional):", choices = c("(all cells)" = "", grouping_choices), width = "100%")),
              column(3, uiOutput("gcontrast_pb_subset_level_ui")),
              column(3, div(style = "margin-top: 25px;",
                            actionButton("gcontrast_pb_compute", "Compute pseudobulk DE",
                                         style = "background-color: #fafafa; color: #000000; border: 3px solid #39FF14; font-weight: 600;")))
            )
          )
        ),
        conditionalPanel(
          condition = "input.gcontrast_method == 'percell'",
          withSpinner(DTOutput("gcontrast_table"), caption = "Computing contrast...", color = "#39FF14")),
        conditionalPanel(
          condition = "input.gcontrast_method == 'pseudobulk'",
          uiOutput("gcontrast_pb_caption"),
          withSpinner(DTOutput("gcontrast_pb_table"), caption = "Running DESeq2...", color = "#39FF14"))
      )))
    }


    do.call(tabsetPanel, tabs)
  })

  # Which metadata column the Clusters/Metadata UMAP is coloured by. TWO controls
  # can drive it and they hold DIFFERENT columns: "Color UMAP by" (non-clustering
  # columns only) and the "Resolution / clustering" quick-pick (clustering columns
  # only). We colour by whichever the user changed most recently -- so clusterings
  # appear ONLY in the Resolution picker, never duplicated in "Color UMAP by".
  # The active colour column for the Clusters/Metadata UMAP (clusterings now live
  # in "Color UMAP by" too, so there's no separate Resolution picker to sync with).
  active_group <- reactive({
    req(input$group_var)
    input$group_var
  })

  # The per-category show/hide checkboxes under "Color UMAP by". Rebuilt every
  # time the colouring column changes; NA (unlabelled) cells become an "NA" box.
  output$category_filter_ui <- renderUI({
    gv <- active_group()
    req(gv)
    # Continuous columns: a gradient + a value-range slider. Cells outside the
    # range go grey; cells inside keep the gradient. Defaults span the full range.
    if (is_continuous_col(gv)) {
      x <- suppressWarnings(as.numeric(as.character(seu()@meta.data[[gv]])))
      rng <- range(x, na.rm = TRUE)
      step <- signif(diff(rng) / 100, 2)
      return(tagList(
        helpText("Type a value range to keep (others go grey); press Tab/Enter to apply:"),
        fluidRow(
          column(6, numericInput("cont_min", "Min:", value = rng[1], step = step)),
          column(6, numericInput("cont_max", "Max:", value = rng[2], step = step))
        )
      ))
    }
    vals <- as.character(seu()@meta.data[[gv]])
    vals[is.na(vals)] <- "NA"
    vals <- sort_levels(vals)
    checkboxGroupInput("categories_shown", "Show categories:", choices = vals, selected = vals, inline = TRUE)
  })

  # ---------------------------------------------------------------------------
  # Clusters / Metadata UMAP (one plotly scatter coloured by a metadata column)
  # ---------------------------------------------------------------------------

  # Turn a metadata column into categories for colouring. Categorical columns
  # (factor/character/logical, low-cardinality numeric) are used as-is; a
  # CONTINUOUS numeric (>50 distinct values) is split into quantile RANGES so it
  # can be shown categorically (one colour per range) instead of per value.
  # Returns the per-cell values, the ordered level list (present levels only,
  # NA last), and whether it was binned.
  # TWEAK: change `length.out = 7` for a different number of ranges (6 here).
  # A metadata column is treated as CONTINUOUS (coloured by a gradient) when its
  # values are numeric -- including numbers stored as character/factor, e.g.
  # prediction scores. Otherwise it's categorical (discrete colours + checkboxes).
  #   * clustering columns          -> always categorical
  #   * ">90% coerce to numbers"    -> text labels stay categorical.
  #   * has NON-INTEGER values      -> a real-valued score (e.g. 0.031) is
  #       continuous even with few distinct values (like a per-sample score).
  #   * OR >50 distinct integers    -> e.g. nCount_RNA is continuous, while small
  #       integer codes (cluster ids 0..17) stay categorical.
  # TWEAK: change the 50 to shift the integer categorical/continuous threshold.
  is_continuous_col <- function(col) {
    if (col %in% clust_cols_r()) return(FALSE)           # clusterings are categorical
    xn <- suppressWarnings(as.numeric(as.character(seu()@meta.data[[col]])))
    if (mean(!is.na(xn)) <= 0.9) return(FALSE)            # not numeric-like
    u <- unique(xn[!is.na(xn)])
    any(u != round(u)) || length(u) > 50                 # decimals, or many integers
  }

  # Base data for the UMAP, in a fixed (seeded) random point order. Depends ONLY
  # on the object + colouring column - NOT on which categories are ticked - so
  # that toggling categories doesn't rebuild it (the observer below just
  # recolours the existing plot in place). Branches on continuous vs categorical.
  umap_base <- reactive({
    gv <- active_group()
    req(gv)
    emb <- umap_emb()
    raw <- seu()@meta.data[[gv]]
    set.seed(1)
    ord <- sample(nrow(emb))             # fixed interleaved order (avoids stacking)
    df <- data.frame(UMAP_1 = emb[ord, 1], UMAP_2 = emb[ord, 2], stringsAsFactors = FALSE)

    if (is_continuous_col(gv)) {
      v <- suppressWarnings(as.numeric(as.character(raw)))[ord]   # handles numeric stored as text
      df$value <- v
      df$tooltip <- paste0(pretty_label(gv), ": ", round(v, 3))
      return(list(df = df, continuous = TRUE))
    }

    # Categorical: discrete colours (hex - plotly.js renders names like "grey85"
    # as black). NA cells become an explicit "NA" category.
    g <- as.character(raw)[ord]; g[is.na(g)] <- "NA"
    df$group <- g
    df$tooltip <- paste0("Group: ", g)
    cats <- sort_levels(g)
    non_na <- setdiff(cats, "NA")
    pal <- setNames(scales::hue_pal()(length(non_na)), non_na)
    if ("NA" %in% cats) pal[["NA"]] <- "#B3B3B3"   # unassigned cells (grey70)
    list(df = df, continuous = FALSE, cats = cats, palette = pal)
  })

  # Given a selection, returns the per-point colours (category colour if shown,
  # else grey #D9D9D9) plus the per-category legend-swatch colours. A NULL/empty/
  # stale selection falls back to "all shown" so the plot is never blank.
  umap_point_cols <- function(base, cats_shown) {
    if (is.null(cats_shown) || !any(base$cats %in% cats_shown)) cats_shown <- base$cats
    disp <- ifelse(base$cats %in% cats_shown, base$palette[base$cats], "#D9D9D9")
    names(disp) <- base$cats
    col <- unname(disp[base$df$group])
    col[is.na(col)] <- "#D9D9D9"        # safety net: never pass NA to plotly
    list(points = col, legend = disp)
  }

  # Axis range / tick step / labels for an arbitrary embedding matrix (used by the
  # compare view, which frames each reduction on its own scale).
  red_frame <- function(emb) {
    r <- range(emb[, 1:2], na.rm = TRUE); pad <- diff(r) * 0.03
    rng <- c(r[1] - pad, r[2] + pad)
    raw <- diff(rng) / 6; mag <- 10^floor(log10(raw)); nice <- c(1, 2, 2.5, 5, 10)
    cn <- colnames(emb)
    list(rng = rng, dtick = mag * nice[which.min(abs(nice - raw / mag))],
         xlab = if (length(cn) >= 2) cn[1] else "Dim 1",
         ylab = if (length(cn) >= 2) cn[2] else "Dim 2")
  }
  red_label <- function(r) if (tolower(r) %in% c("umap", "tsne", "pca")) toupper(r) else r

  # Build a full UMAP plotly from a base (colouring) + a frame -- used by the
  # side-by-side compare view. (The MAIN single plot has its own render with the
  # restyle optimisation; this one just re-renders, which is fine for a 2-up view.)
  build_umap_plot <- function(base, cats_shown, fr, title) {
    df <- base$df
    if (isTRUE(base$continuous)) {
      full <- range(df$value, na.rm = TRUE)
      lo <- if (is.null(input$cont_min) || is.na(input$cont_min)) full[1] else input$cont_min
      hi <- if (is.null(input$cont_max) || is.na(input$cont_max)) full[2] else input$cont_max
      in_rng <- !is.na(df$value) & df$value >= lo & df$value <= hi
      out_df <- df[!in_rng, ]; in_df <- df[in_rng, ]
      p <- plot_ly()
      if (nrow(out_df) > 0)
        p <- add_trace(p, data = out_df, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                       marker = list(size = 3, color = "#E8E8E8"), text = ~tooltip, hoverinfo = "text", showlegend = FALSE)
      if (nrow(in_df) > 0)
        p <- add_trace(p, data = in_df, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                       marker = list(size = 4, color = ~value, colorscale = "Viridis", cmin = full[1], cmax = full[2],
                                     colorbar = list(thickness = 12, len = 0.5, outlinewidth = 0,
                                                     x = 1.02, xanchor = "left", y = 0.5, yanchor = "middle")),
                       text = ~tooltip, hoverinfo = "text", showlegend = FALSE)
      showlg <- FALSE; lg <- NULL
    } else {
      cols <- umap_point_cols(base, cats_shown); df$col <- cols$points
      p <- plot_ly(df, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                   marker = list(size = 4, color = ~col), text = ~tooltip, hoverinfo = "text", showlegend = FALSE)
      for (cat in sort_levels(df$group)) {
        idx <- which(df$group == cat)[1]
        p <- add_trace(p, x = df$UMAP_1[idx], y = df$UMAP_2[idx], type = "scatter", mode = "markers",
                       marker = list(size = 0.1, color = cols$legend[[cat]]), name = cat,
                       showlegend = TRUE, hoverinfo = "skip", inherit = FALSE)
      }
      showlg <- TRUE; lg <- list(itemsizing = "constant", y = 0.5, yanchor = "middle")
    }
    p %>% plotly::layout(
      hovermode = "closest",
      title = list(text = title, x = 0, xref = "paper", xanchor = "left",
                   font = list(size = 14, family = "Open Sans, verdana, arial, sans-serif")),
      xaxis = list(title = fr$xlab, zeroline = FALSE, range = fr$rng, tick0 = 0, dtick = fr$dtick,
                   constrain = "domain", constraintoward = "center"),
      yaxis = list(title = fr$ylab, zeroline = FALSE, scaleanchor = "x", scaleratio = 1, tick0 = 0, dtick = fr$dtick,
                   constrain = "domain", constraintoward = "middle", range = fr$rng),
      showlegend = showlg, legend = lg, margin = list(l = 46, r = 10, t = 36, b = 40)
    ) %>% plotly::config(displayModeBar = FALSE)
  }

  output$umap_plot <- renderPlotly({
    base <- umap_base()
    df <- base$df
    gv <- active_group()

    # ---- Continuous column: viridis gradient + a typed Min/Max range. Cells in
    # range keep the gradient; cells outside go grey. Colour limits are the FULL
    # data range so a value's colour doesn't shift when you narrow the range. ----
    if (isTRUE(base$continuous)) {
      full <- range(df$value, na.rm = TRUE)
      lo <- if (is.null(input$cont_min) || is.na(input$cont_min)) full[1] else input$cont_min
      hi <- if (is.null(input$cont_max) || is.na(input$cont_max)) full[2] else input$cont_max
      in_rng <- !is.na(df$value) & df$value >= lo & df$value <= hi
      out_df <- df[!in_rng, ]; in_df <- df[in_rng, ]

      p <- plot_ly()
      if (nrow(out_df) > 0) {
        p <- add_trace(p, data = out_df, x = ~UMAP_1, y = ~UMAP_2,
                       type = "scattergl", mode = "markers",
                       marker = list(size = 3, color = "#E8E8E8"),
                       text = ~tooltip, hoverinfo = "text", showlegend = FALSE)
      }
      if (nrow(in_df) > 0) {
        p <- add_trace(p, data = in_df, x = ~UMAP_1, y = ~UMAP_2,
                       type = "scattergl", mode = "markers",
                       marker = list(size = 4, color = ~value, colorscale = "Viridis",
                                     cmin = full[1], cmax = full[2],
                                     # Title on the SIDE (vertical, along the bar)
                                     colorbar = list(title = list(text = pretty_label(gv), side = "right"),
                                                     thickness = 15, len = 0.5, outlinewidth = 0,
                                                     x = 1.02, xanchor = "left", y = 0.5, yanchor = "middle")),
                       text = ~tooltip, hoverinfo = "text", showlegend = FALSE)
      }
      p <- plotly::layout(
        p, hovermode = "closest",
        title = list(text = paste(embedding_label(), "colored by", pretty_label(gv)),
                     x = 0, xref = "paper", xanchor = "left",
                     font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
        xaxis = list(title = axis_labels()[1], zeroline = FALSE, range = umap_range(),
                     tick0 = 0, dtick = umap_dtick(),
                     constrain = "domain", constraintoward = "center"),
        yaxis = list(title = axis_labels()[2], zeroline = FALSE, scaleanchor = "x", scaleratio = 1,
                     tick0 = 0, dtick = umap_dtick(),
                     constrain = "domain", constraintoward = "middle", range = umap_range()),
        margin = list(l = 50, r = 10, t = 40, b = 40)
      )
      return(plotly::config(p, displayModeBar = FALSE))
    }

    # ---- Categorical column: discrete colours + toggleable legend (below). ----
    # Initial colours reflect the current selection but without a reactive
    # dependency on it (isolate) - the observer below updates them in place.
    cols <- umap_point_cols(base, isolate(input$categories_shown))
    df$col <- cols$points

    # A SINGLE marker trace (index 0) holds every point, so toggling a category
    # is just a colour restyle - the plot never reloads. (plotly's `color=`
    # would split into one trace per category and stack them, hiding one.)
    p <- plot_ly(
      df, x = ~UMAP_1, y = ~UMAP_2,
      type = "scattergl", mode = "markers",
      marker = list(size = 4, color = ~col),   # TWEAK: UMAP point size
      text = ~tooltip, hoverinfo = "text", showlegend = FALSE
    )

    # Legend on the right
    legend_order <- sort_levels(df$group)
    for (cat in legend_order) {
      idx <- which(df$group == cat)[1]
      p <- add_trace(
        p, x = df$UMAP_1[idx], y = df$UMAP_2[idx],
        type = "scatter", mode = "markers",
        marker = list(size = 0.1, color = cols$legend[[cat]]),
        name = cat, showlegend = TRUE, hoverinfo = "skip", inherit = FALSE
      )
    }

    p <- plotly::layout(
      p, hovermode = "closest",
      title = list(text = paste(embedding_label(), "colored by", pretty_label(gv)), x = 0, xref = "paper", xanchor = "left",
                   font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
      xaxis = list(title = axis_labels()[1], zeroline = FALSE, range = umap_range(),
                   tick0 = 0, dtick = umap_dtick(),
                   constrain = "domain", constraintoward = "center"),
      yaxis = list(title = axis_labels()[2], zeroline = FALSE, scaleanchor = "x", scaleratio = 1,
                   tick0 = 0, dtick = umap_dtick(),
                   constrain = "domain", constraintoward = "middle", range = umap_range()),
      showlegend = TRUE,
      legend = list(title = list(text = pretty_label(gv)), itemsizing = "constant",
                    y = 0.5, yanchor = "middle"),
      margin = list(l = 50, r = 10, t = 40, b = 40)
    )
    plotly::config(p, displayModeBar = FALSE)
  })

  # Toggling a category recolours the EXISTING plot via plotly restyle
  observeEvent(input$categories_shown, {
    req(active_group())
    base <- umap_base()
    if (isTRUE(base$continuous)) return()   # gradient mode has no category toggles
    cols <- umap_point_cols(base, input$categories_shown)
    legend_order <- sort_levels(base$df$group)
    proxy <- plotlyProxy("umap_plot", session)
    plotlyProxyInvoke(proxy, "restyle", list(marker.color = list(cols$points)), list(0))
    leg_cols <- unname(cols$legend[legend_order])
    plotlyProxyInvoke(proxy, "restyle",
                      list(marker.color = as.list(leg_cols)),
                      as.list(seq_along(legend_order)))
  }, ignoreInit = TRUE)

  # ---- Split by: one UMAP panel per level of a chosen column ----------------
  # Same idea as Seurat's DimPlot(split.by=): break the metadata-coloured map out
  # into small multiples, one per group, over a grey backdrop of all cells so the
  # panels share the same frame and are directly comparable.
  umap_split_max <- 12   # WebGL-context cap, like the Gene Expression grid.

  # All levels of the chosen split column (choices for the "groups to show" picker).
  umap_split_all_levels <- reactive({
    sb <- input$split_by
    req(!is.null(sb) && nzchar(sb))
    v <- as.character(seu()@meta.data[[sb]]); v[is.na(v)] <- "NA"
    sort_levels(unique(v))
  })

  # A "groups to show" picker (starts EMPTY): rendering every panel at once is a
  # WebGL overload that freezes the tab, so the user opts into the groups to draw.
  output$umap_split_levels_ui <- renderUI({
    req(input$split_by, nzchar(input$split_by))
    selectizeInput("umap_split_show", "Groups to show:", choices = umap_split_all_levels(),
                   selected = character(0), multiple = TRUE,
                   options = list(placeholder = "Pick groups to show (up to 12)..."))
  })

  # The selected groups, ignoring any stale picks left over from another column.
  umap_split_sel <- reactive({
    req(input$split_by, nzchar(input$split_by))
    intersect(input$umap_split_show, umap_split_all_levels())
  })

  # Build one split panel: grey backdrop (all cells) + this group's cells coloured
  # by the active "Color UMAP by" column, on the shared UMAP frame.
  build_umap_split_panel <- function(level) {
    base <- umap_base(); df <- base$df
    emb  <- umap_emb(); set.seed(1); ord <- sample(nrow(emb))   # same shuffle umap_base() used
    sb <- input$split_by
    sval <- as.character(seu()@meta.data[[sb]])[ord]; sval[is.na(sval)] <- "NA"
    ing <- sval == level
    title <- paste0(pretty_label(sb), ": ", level, " (n=", sum(ing), ")")
    p <- plot_ly()
    # Grey backdrop = all cells (context for sparse groups).
    p <- add_trace(p, data = df, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                   marker = list(size = 3, color = "#E8E8E8"), hoverinfo = "skip", showlegend = FALSE)
    if (isTRUE(base$continuous)) {
      full <- range(df$value, na.rm = TRUE)
      lo <- if (is.null(input$cont_min) || is.na(input$cont_min)) full[1] else input$cont_min
      hi <- if (is.null(input$cont_max) || is.na(input$cont_max)) full[2] else input$cont_max
      keep <- ing & !is.na(df$value) & df$value >= lo & df$value <= hi
      if (any(keep))
        p <- add_trace(p, data = df[keep, ], x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                       marker = list(size = 4, color = ~value, colorscale = "Viridis", cmin = full[1], cmax = full[2],
                                     colorbar = list(title = list(text = pretty_label(active_group()), side = "right"),
                                                     thickness = 12, len = 0.5, outlinewidth = 0,
                                                     x = 1.02, xanchor = "left", y = 0.5, yanchor = "middle")),
                       text = ~tooltip, hoverinfo = "text", showlegend = FALSE)
      showlg <- FALSE; lg <- NULL
    } else {
      # Category colours reflect the current show/hide selection (read once; split
      # panels don't live-restyle, so no reactive dependency on the checkboxes). No
      # per-panel legend -- the grid shows ONE shared legend so the map fills the square.
      cols <- umap_point_cols(base, isolate(input$categories_shown))
      sub <- df[ing, ]; sub$col <- cols$points[ing]
      if (nrow(sub) > 0)
        p <- add_trace(p, data = sub, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
                       marker = list(size = 4, color = ~col), text = ~tooltip, hoverinfo = "text", showlegend = FALSE)
      showlg <- FALSE; lg <- NULL
    }
    p %>% plotly::layout(
      hovermode = "closest",
      title = list(text = title, x = 0, xref = "paper", xanchor = "left",
                   font = list(size = 13, family = "Open Sans, verdana, arial, sans-serif")),
      xaxis = list(title = axis_labels()[1], zeroline = FALSE, range = umap_range(), tick0 = 0, dtick = umap_dtick(),
                   constrain = "domain", constraintoward = "center"),
      yaxis = list(title = axis_labels()[2], zeroline = FALSE, scaleanchor = "x", scaleratio = 1, tick0 = 0, dtick = umap_dtick(),
                   constrain = "domain", constraintoward = "middle", range = umap_range()),
      showlegend = showlg, legend = lg, margin = list(l = 44, r = 10, t = 34, b = 38)
    ) %>% plotly::config(displayModeBar = FALSE)
  }

  output$umap_split_grid <- renderUI({
    req(input$split_by, nzchar(input$split_by))
    levs <- umap_split_sel()
    validate(need(length(levs) >= 1, "Pick one or more groups in \"Groups to show\" (in the sidebar)."))
    validate(need(length(levs) <= umap_split_max,
                  paste0("Showing at most ", umap_split_max, " panels at once; pick fewer groups.")))
    # Every panel is the SAME fixed size as the single (non-split) plot. Panels sit
    # side by side and wrap to the next row when they don't fit, instead of shrinking.
    cells <- lapply(seq_along(levs), function(i) {
      id <- paste0("umap_split_", i)
      div(class = "plot-dl-wrap umap-split-cell", style = "position: relative; flex: 0 0 auto; width: min(100vh - 340px, 60vw); aspect-ratio: 1 / 1;",
          plot_dl_btn(id, paste0("UMAP_", levs[i])),
          withSpinner(plotlyOutput(id, width = "100%", height = "100%"), caption = "Rendering...", color = "#39FF14"))
    })
    # grid does NOT grow (flex:0 1 auto), so the shared legend sits snug on its right
    # instead of being pushed to the far edge. Panels wrap when they don't fit.
    grid <- div(style = "display:flex; flex-wrap:wrap; justify-content:center; gap:16px; flex:0 1 auto; min-width:0;", cells)
    # "Download all as one PNG" (composites the panels), shown when >1 panel.
    ids <- paste0("umap_split_", seq_along(levs))
    dl_all <- if (length(levs) > 1)
      div(style = "text-align:right; margin-bottom:6px;",
          plot_dl_all_btn(ids, min(length(levs), 2), NULL, "umap_split_panels"))
    # One shared legend for the whole grid (categorical only), so each panel's map
    # isn't squeezed by its own legend. Continuous panels carry a thin colorbar.
    base <- umap_base()
    if (isTRUE(base$continuous))
      return(tagList(dl_all, div(style = "display:flex; justify-content:center;", grid)))
    pal <- base$palette
    items <- lapply(sort_levels(base$cats), function(cc)
      tags$div(style = "display:flex; align-items:center; margin-bottom:4px; font-size:12px; white-space:nowrap;",
               tags$span(style = sprintf("width:11px; height:11px; border-radius:50%%; background:%s; display:inline-block; margin-right:6px; flex:0 0 auto;", pal[[cc]])),
               cc))
    legend <- div(style = "flex:0 0 auto; max-height:80vh; overflow-y:auto;",
                  tags$strong(style = "display:block; margin-bottom:6px;", paste0(pretty_label(active_group()), ":")), items)
    tagList(
      dl_all,
      # align-items:center -> legend sits vertically centred against the plots.
      div(style = "display:flex; align-items:center; justify-content:center; gap:18px;", grid, legend)
    )
  })

  lapply(seq_len(umap_split_max), function(i) {
    force(i)
    output[[paste0("umap_split_", i)]] <- renderPlotly({
      levs <- umap_split_sel()
      req(i <= length(levs))
      build_umap_split_panel(levs[i])
    })
  })

  # Max Gene Expression panels drawn at once. Each gene is a full-resolution WebGL
  # UMAP, and browsers cap live WebGL contexts (~8-16) TWEAK: raise cautiously.
  feature_max_panels <- 12

  # ---- Square grid for Gene Expression ----
  # Each plot is a square whose side is the SMALLER of the width available per
  # column and the height available per row (viewport-aware). Using the same
  # value for the column and row tracks makes every cell a true square (so
  # plotly's equal-scale axes render correctly, single plot included), and the
  # height term keeps the whole grid on screen - 4 plots fit without scrolling.
  output$dynamic_feature_plots <- renderUI({
    specs <- plot_specs()
    n <- length(specs)
    if (n == 0) {
      # In "layers" mode a gene may be chosen but no layer ticked yet -> nudge
      if (identical(input$compare_mode, "layers") && length(found_genes()) > 0)
        return(helpText("Tick one or more layers under \"Layers to show\" to display them."))
      return(NULL)
    }
    # Refuse to draw more than feature_max_panels WebGL UMAPs at once - past the
    # browser's live-context ceiling the tab freezes (same guard Split by Group uses).
    validate(need(n <= feature_max_panels,
                  sprintf("That would draw %d plots at once (max %d), which can freeze the browser. Select fewer genes.",
                          n, feature_max_panels)))

    # "layers" mode = one gene per layer, each with its own scale, so no single
    # shared colorbar. Otherwise the shared colorbar shows in "same scale" mode.
    layers_mode <- identical(input$compare_mode, "layers")
    shared <- !layers_mode && identical(input$scale_mode, "shared")
    ncol <- max(1, min(input$feature_ncol, n))
    nrow <- ceiling(n / ncol)

    # TWEAK: these offsets tune the (fixed) square plot size. If plots feel too
    # small/large, adjust side_offset (horizontal) or top_offset (the height cap).
    # (Controls now sit in a row ABOVE the plots, not in a left sidebar, so the
    # plots use nearly the full page width - small side_offset.)
    #
    # The size is driven by the width available per column and CAPPED to one
    # viewport height - it does NOT depend on the number of rows, so each plot
    # keeps the same size no matter how many genes are added; the extra rows
    # simply extend the page and it scrolls down.
    gap <- 10
    side_offset <- 90    # px taken up by page margins/gutters (no sidebar now)
    top_offset  <- 380   # px above the plots (title, status, tabs, controls row)
    cb_width    <- if (shared) 150 else 0   # room for the single shared colorbar

    w_term <- sprintf("calc((100vw - %dpx) / %d)", side_offset + cb_width + (ncol - 1) * gap, ncol)
    h_cap  <- sprintf("calc(100vh - %dpx)", top_offset)   # constant cap (independent of nrow)
    cell   <- sprintf("min(%s, %s)", w_term, h_cap)

    plot_list <- lapply(seq_len(n), function(i) {
      id <- paste0("feature_plot_individual_", i)
      div(class = "plot-dl-wrap", style = "position: relative; width: 100%; height: 100%;",
          plot_dl_btn(id, specs[[i]]$title),
          withSpinner(
            plotlyOutput(id, width = "100%", height = "100%"),
            caption = "Rendering plot...", color = "#39FF14"
          ))
    })

    grid_div <- div(
      style = paste0("display: grid; grid-template-columns: repeat(", ncol,
                     ", ", cell, "); grid-template-rows: repeat(", nrow,
                     ", ", cell, "); gap: ", gap, "px; justify-content: center;",
                     " align-content: start;"),
      plot_list)

    # In shared mode each plot's own colorbar is hidden; show one to the right.
    # Its column stretches to the full grid height and a one-cell-tall colorbar
    # is centred in it, so the scale sits in the vertical middle of ALL the
    # plots no matter how many genes are added.
    body <- if (shared) {
      div(style = "display: flex; align-items: stretch; justify-content: center; gap: 15px;",
          grid_div,
          # Column stretches to full grid height; inner box is one cell tall and
          # centred, so the colorbar sits in the vertical middle of all plots.
          # (height is set on this div's style, NOT on plotlyOutput, because the
          # output's height param rejects a min()/calc() expression.)
          div(style = paste0("width: ", cb_width, "px; display: flex; align-items: center; justify-content: center;"),
              div(style = paste0("width: 100%; height: ", cell, ";"),
                  plotlyOutput("shared_colorbar", width = "100%", height = "100%"))))
    } else {
      grid_div
    }

    # "Download all as one PNG" (all panels composited) - only for >1 panel.
    # In shared-scale mode the single side colorbar is included in the image.
    if (n > 1) {
      ids <- paste0("feature_plot_individual_", seq_len(n))
      cb <- if (shared) "shared_colorbar" else NULL
      tagList(div(style = "text-align: right; padding-right: 0;",
                  plot_dl_all_btn(ids, ncol, cb, "gene_expression_panels")), body)
    } else {
      body
    }
  })

  # Single shared colorbar, shown to the right of the grid in shared scale mode
  output$shared_colorbar <- renderPlotly({
    req(identical(input$scale_mode, "shared"))
    lims <- shared_limits()
    lab <- switch(feature_layer(), counts = "Counts", scale.data = "Scaled (z)", "Expression")
    plot_ly(
      x = c(0, 0), y = c(0, 1),
      type = "scatter", mode = "markers",
      marker = list(
        color = c(lims[1], lims[2]),
        # TWEAK: keep this gradient identical to the per-plot one above so the
        # shared colorbar matches the points.
        colorscale = EXPR_COLORSCALE,   # blue -> white -> red (see EXPR_COLORSCALE)
        cmin = lims[1], cmax = lims[2],
        opacity = 0,                                           # hide the dummy points
        showscale = TRUE,
        # x = 0 / xanchor left -> bar sits at the LEFT of its own image, so it
        # ends up right next to the plots (on screen and in the combined PNG).
        colorbar = list(title = list(text = lab), thickness = 15,
                        len = 0.5, outlinewidth = 0, x = 0, xanchor = "left")
      )
    ) %>%
      plotly::layout(
        xaxis = list(visible = FALSE, range = c(-1, 1)),
        yaxis = list(visible = FALSE, range = c(-1, 2)),
        margin = list(l = 0, r = 0, t = 10, b = 10),
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor = "rgba(0,0,0,0)"
      ) %>%
      plotly::config(displayModeBar = FALSE)
  })

  # ---------------------------------------------------------------------------
  # Gene Expression: one independent plot per gene (so each can have its own
  # colour scale in "per gene" mode). We pre-declare up to max_plots outputs;
  # dynamic_feature_plots (above) decides how many are actually shown.
  # TWEAK: raise max_plots if you ever need to show more than 50 genes at once.
  # ---------------------------------------------------------------------------
  max_plots <- 50

  # Assay layers present, in an order: data, counts, scale.data
  # Now respects the Gene Expression assay selector.
  available_layers <- reactive({
    assay_use <- feature_assay_val()
    al <- tryCatch(SeuratObject::Layers(seu()[[assay_use]]), error = function(e) NULL)
    if (is.null(al) || length(al) == 0) al <- "data"
    lo <- intersect(c("data", "counts", "scale.data"), al)
    if (length(lo) == 0) lo <- al[1]
    lo
  })

  # What each plot slot should show: a list of {gene, layer, title}.
  #  * "genes" mode  -> one plot per selected gene, all in the chosen layer.
  #  * "layers" mode -> one plot per CHOSEN layer (via "Layers to show"), of the
  #    FIRST gene. Defaults to all available layers until the user narrows it.
  plot_specs <- reactive({
    fg <- found_genes()
    req(length(fg) > 0)
    if (identical(input$compare_mode, "layers")) {
      g <- fg[[1]]
      avail <- available_layers()
      # Only the layers the user ticked, in preferred order (data, counts,
      # scale.data). Nothing ticked -> nothing rendered
      keep <- avail[avail %in% input$feature_layers_show]
      lapply(keep, function(L)
        list(gene = g, layer = L, title = paste0(g, " · ", layer_label(L))))
    } else {
      layer <- feature_layer()
      # Title shows the gene and the data layer, so a downloaded panel is self-
      # describing (e.g. "Malat1 · Normalized (log)")
      lapply(fg, function(g) list(gene = g, layer = layer,
                                  title = paste0(g, " · ", layer_label(layer))))
    }
  })

  # One value-stable reactive per plot slot (adding a gene only invalidates the
  # NEW slot - existing slots return the same spec, so their render is cached).
  slot_spec <- lapply(1:max_plots, function(i) {
    force(i)
    reactive({
      ps <- plot_specs()
      if (i <= length(ps)) ps[[i]] else NULL
    })
  })

  lapply(1:max_plots, function(i) {
    force(i)
    output[[paste0("feature_plot_individual_", i)]] <- renderPlotly({
      spec <- slot_spec[[i]]()
      req(!is.null(spec))
      g <- spec$gene; layer <- spec$layer

      umap_coords <- umap_emb()
      # RAW values of the slot's layer (counts / data / scale.data) - no z-scoring,
      # so you literally see the layer's values (e.g. before vs after normalization).
      assay_use <- feature_assay_val()
      lmat <- GetAssayData(seu(), layer = layer, assay = assay_use)
      validate(need(g %in% rownames(lmat),
                    paste0(g, " is not present in the '", layer, "' layer of assay '", assay_use, "'.")))
      plot_df <- data.frame(
        UMAP_1 = umap_coords[, 1],
        UMAP_2 = umap_coords[, 2],
        expression = as.numeric(lmat[g, ]),
        check.names = FALSE
      )

      # Colour limits. In "layers" mode each panel is a DIFFERENT layer (different
      # scale), so every panel gets its own scale + colorbar. Otherwise: shared
      # across genes (one colorbar) or per-gene. expr_limits() = 0..max, or a
      # symmetric range for scale.data (which has negatives).
      layers_mode <- identical(input$compare_mode, "layers")
      if (!layers_mode && identical(input$scale_mode, "shared")) {
        lims <- shared_limits(); show_bar <- FALSE
      } else {
        lims <- expr_limits(plot_df$expression); show_bar <- TRUE
      }

      # Draw low-value points first so high-value ones sit on top.
      plot_df <- plot_df[order(plot_df$expression), ]
      lab <- switch(layer, counts = "Counts", scale.data = "Scaled (z)", "Expression")

      # Hover reads the value straight from marker.color (the expression already
      # sent to plotly) via hovertemplate - so we don't ship a second ~19k-element
      # per-cell tooltip string array to the browser for every gene panel.
      p <- plot_ly(
        plot_df, x = ~UMAP_1, y = ~UMAP_2,
        type = "scattergl", mode = "markers",
        marker = list(
          size = 4,                                              # TWEAK: gene-expression point size
          color = ~expression,
          colorscale = EXPR_COLORSCALE,   # TWEAK: edit EXPR_COLORSCALE near the top
          cmin = lims[1], cmax = lims[2],
          showscale = show_bar,
          colorbar = list(title = list(text = lab), thickness = 15, len = 0.5, outlinewidth = 0)
        ),
        hovertemplate = paste0("Gene: ", g, "<br>", lab, ": %{marker.color:.3f}<extra></extra>")
      )

      # Same framing as the Clusters/Metadata UMAP: fixed square window, no tick
      # chrome, left-aligned title (matching font/size).
      p <- plotly::layout(
        p,
        title = list(text = spec$title, x = 0, xref = "paper", xanchor = "left",
                     font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
        xaxis = list(title = axis_labels()[1], zeroline = FALSE, showline = FALSE, ticks = "",
                     range = umap_range(), constrain = "range"),
        yaxis = list(title = axis_labels()[2], zeroline = FALSE, showline = FALSE, ticks = "",
                     scaleanchor = "x", scaleratio = 1, constrain = "range",
                     range = umap_range()),
        margin = list(l = 50, r = 10, t = 40, b = 30)
      )
      plotly::config(p, displayModeBar = FALSE)
    })
  })

  # ===========================================================================
  # Aggregated Expression (pseudobulk) tab.
  # For each selected gene, summarise its expression (from the "data" layer,
  # same source as the Gene Expression tab) across the cells of each level of
  # the chosen metadata column -> a genes x groups matrix shown as a heatmap
  # and a table.
  # ===========================================================================

  # Aggregated tab: fixed to RNA assay + "data" (Normalized log); selectors removed.
  agg_assay_val <- make_assay_val("agg_assay")
  agg_layer <- make_layer_val("agg_layer")

  # agg_matrix(): the genes x group-levels matrix of aggregated expression, using
  # the RAW values of the chosen assay/layer (no app-side z-scoring).
  agg_matrix <- reactive({
    req(input$agg_group_var)
    validate(need(length(input$agg_genes) > 0, "Select one or more genes to aggregate."))
    assay_use <- agg_assay_val()
    layer_use <- agg_layer()
    full <- GetAssayData(seu(), layer = layer_use, assay = assay_use)   # genes x cells
    # Filter against the chosen layer's genes: scale.data only holds variable
    # features, so some picked genes may legitimately be absent there.
    genes <- input$agg_genes[input$agg_genes %in% rownames(full)]
    validate(need(length(genes) > 0,
                  paste0("None of the selected genes are in the '", layer_use,
                         "' layer of assay '", assay_use, "'.")))

    fun <- switch(input$agg_fun, sum = sum, median = stats::median, mean)
    expr_mat <- as.matrix(full[genes, , drop = FALSE])     # genes x cells (chosen layer)
    grp <- as.character(seu()@meta.data[[input$agg_group_var]])
    grp[is.na(grp)] <- "NA"                                # un-annotated cells -> "NA" group
    groups <- sort_levels(grp)

    mat <- matrix(NA_real_, nrow = length(genes), ncol = length(groups),
                  dimnames = list(genes, groups))
    # tapply groups once per gene (vs subsetting per gene x group in a nested loop).
    for (g in genes) {
      agg_vals <- tapply(expr_mat[g, ], grp, fun)
      mat[g, names(agg_vals)] <- agg_vals
    }
    mat
  })

  # Short label for the summary function (colorbar / column headers).
  agg_fun_label <- reactive({
    switch(input$agg_fun, sum = "Sum", median = "Median", "Mean")
  })

  # Human-readable name of the chosen layer, for the heatmap title.
  agg_layer_label <- reactive({
    layer_label(agg_layer())
  })

  # Heatmap: genes (rows) x group levels (columns), blue -> white -> red like the
  # Gene Expression plots. Colour scale auto-fits the aggregated values.
  output$agg_heatmap <- renderPlotly({
    mat <- agg_matrix()
    # Fix the colour range to the RAW aggregated values: 0..max for non-negative
    # layers (counts / data), symmetric about 0 for scale.data so 0 reads white.
    lims <- expr_limits(as.vector(mat))
    plot_ly(
      x = colnames(mat), y = rownames(mat), z = mat,
      type = "heatmap",
      colorscale = EXPR_COLORSCALE,   # blue -> white -> red (see EXPR_COLORSCALE)
      zmin = lims[1], zmax = lims[2], zauto = FALSE,
      colorbar = list(title = list(text = agg_fun_label()), thickness = 15, outlinewidth = 0),
      hovertemplate = paste0("Gene: %{y}<br>", pretty_label(input$agg_group_var),
                             ": %{x}<br>", agg_fun_label(), ": %{z:.3f}<extra></extra>")
    ) %>%
      plotly::layout(
        title = list(text = paste0(agg_fun_label(), " expression of ", nrow(mat),
                                   " gene", if (nrow(mat) != 1) "s" else "", " per ",
                                   pretty_label(input$agg_group_var), "  ·  ",
                                   agg_layer_label()),
                     x = 0, xref = "paper", xanchor = "left",
                     font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
        xaxis = list(title = pretty_label(input$agg_group_var), tickangle = -45),
        yaxis = list(title = "Gene", autorange = "reversed"),   # first gene on top
        margin = list(l = 90, r = 10, t = 40, b = 100)
      ) %>%
      plotly::config(displayModeBar = FALSE)
  })

  # Table: same matrix as a data frame (Gene column + one column per group level).
  output$agg_table <- renderDT({
    mat <- agg_matrix()
    df <- data.frame(Gene = rownames(mat), round(mat, 3),
                     check.names = FALSE, row.names = NULL)
    dt_with_excel(df, "aggregated_expression")
  }, server = FALSE)   # Excel button exports ALL rows, not just the current page

  # ===========================================================================
  # Marker Plot tab: a dot plot of selected genes x a grouping column.
  #   dot size  = % of cells in the group that express the gene (raw layer > 0)
  #   dot colour = average expression, z-scored per gene across groups
  # (so the diverging blue-white-red scale reads as "low/high in this group").
  # ===========================================================================
  dot_assay_val <- make_assay_val("dot_assay")
  dot_layer <- make_layer_val("dot_layer")

  # Populate the (manual-mode) gene picker in the chosen sort order.
  observeEvent(list(seu(), input$dot_gene_sort), {
    req(seu(), input$dot_gene_sort)
    updateSelectizeInput(session, "dot_genes",
                         choices = ordered_gene_choices(input$dot_gene_sort), server = TRUE)
  })

  # AUTO mode: FindAllMarkers on the grouping column (each group vs the rest),
  # computed only when "Compute markers" is clicked (it's heavier than the
  # single-group Cluster Markers tab). Returns the raw results table (or an error :).
  dot_auto <- eventReactive(input$dot_compute, {
    req(input$dot_group_var)
    grp <- input$dot_group_var; assay_use <- dot_assay_val(); test <- "wilcox"
    obj <- prepare_obj_for_de(seu(), assay_use)
    Seurat::Idents(obj) <- grp                         # use the grouping column as identities
    res <- tryCatch(
      suppressWarnings(Seurat::FindAllMarkers(
        obj, only.pos = TRUE, test.use = test, assay = assay_use,
        min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)),   # drop sparse / low-FC genes
      error = function(e) e
    )
    # Remember what this was computed for, so we can detect stale results below.
    list(res = res, group = grp, assay = assay_use, test = test)
  })

  # The genes to plot: the manual picker, or the top-N markers per group from the
  # last auto computation (re-slicing top-N is instant; it doesn't recompute DE).
  dot_genes_final <- reactive({
    if (identical(input$dot_gene_mode, "auto")) {
      a <- dot_auto()
      # Stale-guard: the auto genes belong to the grouping/assay/test they were
      # computed with. If any changed since, prompt a refresh instead of showing
      # the previous grouping's genes.
      validate(need(!is.null(a) &&
                    identical(a$group, input$dot_group_var) &&
                    identical(a$assay, dot_assay_val()) &&
                    identical(a$test, "wilcox"),
                    "Click 'Compute markers' to find the top markers for the current grouping."))
      res <- a$res
      validate(need(!inherits(res, "error"),
                    paste("Could not compute markers:",
                          if (inherits(res, "error")) conditionMessage(res) else "")))
      validate(need(is.data.frame(res) && nrow(res) > 0, "No markers found for this grouping."))
      n <- if (is.null(input$dot_top_n) || is.na(input$dot_top_n)) 3 else input$dot_top_n
      # Keep significant markers (when the test gives a p-value; ROC doesn't), so
      # the top-N are real cluster markers, not high-FC-but-sparse junk genes.
      if ("p_val_adj" %in% names(res)) res <- res[!is.na(res$p_val_adj) & res$p_val_adj < 0.05, ]
      validate(need(nrow(res) > 0,
                    "No significant markers passed the filter. Try a different test or grouping."))
      res <- res[order(res$cluster, -res$avg_log2FC), ]
      unique(unlist(by(res, res$cluster, function(d) head(d$gene, n))))   # top N per group, block-ordered
    } else {
      input$dot_genes
    }
  })

  # dot_data(): per gene x group, the % expressing and the mean expression, plus
  # the mean z-scored per gene across groups (for colour).
  # Shared by the dot plot and violin: the (genes x cells) matrix for `genes_sel`
  # in the current assay/layer, plus the group vector and its sorted levels.
  marker_gene_matrix <- function(genes_sel) {
    assay_use <- dot_assay_val(); layer_use <- dot_layer()
    full <- GetAssayData(seu(), layer = layer_use, assay = assay_use)   # genes x cells
    genes <- genes_sel[genes_sel %in% rownames(full)]
    validate(need(length(genes) > 0,
                  paste0("None of the selected genes are in the '", layer_use,
                         "' layer of assay '", assay_use, "'.")))
    grp <- as.character(seu()@meta.data[[input$dot_group_var]]); grp[is.na(grp)] <- "NA"
    list(mat = full[genes, , drop = FALSE], genes = genes, grp = grp, groups = sort_levels(grp))
  }

  dot_data <- reactive({
    req(input$dot_group_var)
    gsel <- dot_genes_final()
    validate(need(length(gsel) > 0,
                  "Pick genes (manual mode) or click 'Compute markers' (auto mode)."))
    mm <- marker_gene_matrix(gsel)
    genes <- mm$genes; m <- mm$mat; grp <- mm$grp; groups <- mm$groups

    pct <- matrix(0, length(genes), length(groups), dimnames = list(genes, groups))
    avg <- matrix(0, length(genes), length(groups), dimnames = list(genes, groups))
    for (l in groups) {
      sub <- m[, grp == l, drop = FALSE]
      pct[, l] <- as.numeric(Matrix::rowMeans(sub > 0)) * 100
      avg[, l] <- as.numeric(Matrix::rowMeans(sub))
    }
    # z-score each gene's mean across groups (constant genes -> 0, i.e. white).
    scaled <- t(scale(t(avg)))
    scaled[!is.finite(scaled)] <- 0
    list(genes = genes, groups = groups, pct = pct, avg = avg, scaled = scaled)
  })

  output$dot_plot <- renderPlotly({
    d <- dot_data()
    df <- expand.grid(gene = d$genes, group = d$groups,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    idx <- cbind(match(df$gene, d$genes), match(df$group, d$groups))
    df$pct    <- d$pct[idx]
    df$avg    <- d$avg[idx]
    df$scaled <- d$scaled[idx]
    df$dotsize <- 2 + (df$pct / 100) * 16                       # 0% -> 2px, 100% -> 18px
    df$tt <- sprintf("Gene: %s<br>%s: %s<br>%% expressing: %.1f%%<br>Avg expr: %.3f",
                     df$gene, pretty_label(input$dot_group_var), df$group, df$pct, df$avg)
    lim <- max(abs(d$scaled), na.rm = TRUE); if (!is.finite(lim) || lim <= 0) lim <- 1

    plot_ly(df, x = ~gene, y = ~group, type = "scatter", mode = "markers",
            marker = list(size = ~dotsize, color = ~scaled, colorscale = EXPR_COLORSCALE,
                          cmin = -lim, cmax = lim, line = list(width = 0.5, color = "rgba(120,120,120,0.5)"),
                          colorbar = list(title = list(text = "Avg expr (z)"),
                                          thickness = 15, len = 0.6, outlinewidth = 0)),
            text = ~tt, hoverinfo = "text") %>%
      plotly::layout(
        title = list(text = paste0("Marker dot plot by ", pretty_label(input$dot_group_var)),
                     x = 0, xref = "paper", xanchor = "left",
                     font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
        xaxis = list(title = "", type = "category", tickangle = -45,
                     categoryorder = "array", categoryarray = d$genes),
        yaxis = list(title = "", type = "category",
                     categoryorder = "array", categoryarray = rev(d$groups)),
        margin = list(l = 130, r = 10, t = 50, b = 120)
      ) %>%
      plotly::config(displayModeBar = FALSE)
  })

  # ---- Marker Plot: Violin variant --------------------------------------------
  dot_genes_safe <- reactive({
    tryCatch(dot_genes_final(), error = function(e) character(0))
  })

  # Fill the violin gene picker with all genes (a searchable dropdown -- type to
  # find any gene). Starts EMPTY so the user chooses which genes to show; keeps the
  # user's own selection when it's already set (does NOT auto-fill any genes).
  observeEvent(list(seu(), input$marker_plot_type, input$violin_gene_sort), {
    req(seu(), identical(input$marker_plot_type, "violin"))
    sort_mode <- if (is.null(input$violin_gene_sort)) "Alphabetical" else input$violin_gene_sort
    cur <- isolate(input$violin_genes)
    sel <- cur[nchar(cur) > 0]   # keep the user's picks, else leave empty
    updateSelectizeInput(session, "violin_genes",
                         choices = ordered_gene_choices(sort_mode),
                         selected = sel, server = TRUE)
  })

  # One stacked row of violins per selected gene (x = group), sharing the group
  # axis. Cells are downsampled for responsiveness - a violin only needs the shape.
  violin_data <- reactive({
    req(input$dot_group_var)
    validate(need(length(input$violin_genes) > 0,
                  "Pick one or more genes above to show as violins."))
    mm <- marker_gene_matrix(head(input$violin_genes, 20))   # cap rows (stacked violins get tall)
    grp <- mm$grp
    keep <- seq_along(grp)
    if (length(keep) > 6000) { set.seed(1); keep <- sort(sample(keep, 6000)) }   # downsample (points are drawn)
    list(mat = as.matrix(mm$mat[, keep, drop = FALSE]),
         grp = grp[keep], genes = mm$genes, groups = mm$groups)
  })

  # Height grows with the number of genes so each violin row stays readable.
  output$violin_plot_ui <- renderUI({
    n <- length(violin_data()$genes)
    plotlyOutput("violin_plot", width = "100%", height = paste0(max(560, 340 * n), "px"))
  })

  output$violin_plot <- renderPlotly({
    d <- violin_data()
    pal <- setNames(scales::hue_pal()(length(d$groups)), d$groups)
    gf <- factor(d$grp, levels = d$groups)
    panels <- lapply(seq_along(d$genes), function(i) {
      g <- d$genes[i]
      df <- data.frame(group = gf, value = d$mat[i, ], stringsAsFactors = FALSE)
      plot_ly(df, x = ~group, y = ~value, type = "violin", color = ~group, colors = pal,
              # jittered individual cells + an inner box (median/quartiles), like a
              # standard single-cell violin. Light, small points so dense zero-heavy
              # genes don't turn into a black bar.
              points = "all", jitter = 0.35, pointpos = 0,
              # spanmode "hard" clips the violin outline to the real data min/max,
              # so the shape doesn't overshoot past the actual max (the default
              # "soft" KDE extends ~one bandwidth beyond the extreme points).
              spanmode = "hard",
              marker = list(size = 1.5, color = "rgba(50,50,50,0.18)"),
              box = list(visible = TRUE, width = 0.15), meanline = list(visible = TRUE),
              scalemode = "width",
              # hover on: points show a cell's exact value, the box shows the
              # quartile stats, the violin outline shows the local density value.
              showlegend = FALSE, hoveron = "violins+points+kde") %>%
        plotly::layout(
          yaxis = list(title = list(text = g, font = list(size = 12)), zeroline = FALSE),
          xaxis = list(title = "", tickangle = -45))
    })
    subplot(panels, nrows = length(d$genes), shareX = TRUE, titleY = TRUE, margin = 0.004) %>%
      plotly::layout(
        title = list(text = paste0("Marker violins by ", pretty_label(input$dot_group_var)),
                     x = 0, xref = "paper", xanchor = "left",
                     font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
        margin = list(l = 90, r = 10, t = 50, b = 120)) %>%
      plotly::config(displayModeBar = FALSE)
  })

  # ===========================================================================
  # Cluster Markers & Contrasts tabs (interactive DT tables).
  # Both read data frames from obj@misc (scdrake format): each row is one
  # cluster/contrast whose `markers` cell holds a per-gene stats table.
  # ===========================================================================

  # Cluster dropdown for the Markers tab, numerically sorted (1,2,...,10 not 1,10,2).
  output$marker_cluster_ui <- renderUI({
    req(input$marker_resolution)
    df <- seu()@misc$cluster_markers
    # Test type selector removed: default to wilcox if present, else the first test.
    mt <- if ("wilcox" %in% df$test_type) "wilcox" else df$test_type[1]
    df <- df[df$source_column == input$marker_resolution & df$test_type == mt, ]
    cluster_ids <- unique(df$group_level)
    cluster_ids <- cluster_ids[order(as.numeric(cluster_ids))]
    selectInput("marker_cluster", "Cluster:", choices = cluster_ids, width = "100%")
  })

  # Helper: put the identifier / key stat columns first, in a fixed order, then
  # append whatever else the table has. (Row names are the gene IDs, hence
  # rownames = FALSE in the datatable() calls below.)
  # TWEAK: change canonical_order to reorder / prioritise different columns.
  reorder_id_cols <- function(df) {
    canonical_order <- c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME", "self.average", "self.detected", "other.average", "other.detected", "Top", "logFC", "AUC", "p.value", "FDR", "lfc_summary", "auc_summary")
    present_in_order <- intersect(canonical_order, colnames(df))
    remaining <- setdiff(colnames(df), present_in_order)
    df[, c(present_in_order, remaining), drop = FALSE]
  }

  # Marker table: pick the row matching (resolution, test, cluster), pull its
  # per-gene stats table, sort by FDR, reorder columns. The DataTables built-in
  # search box (top-right of the table) handles filtering.
  # TWEAK: pageLength controls rows per page.
  output$marker_table <- renderDT({
    req(input$marker_resolution, input$marker_cluster)
    df <- seu()@misc$cluster_markers
    mt <- if ("wilcox" %in% df$test_type) "wilcox" else df$test_type[1]
    row <- df[df$source_column == input$marker_resolution & df$test_type == mt & df$group_level == input$marker_cluster, ]
    validate(need(nrow(row) == 1, "No data for this combination."))
    markers_df <- as.data.frame(row$markers[[1]])
    markers_df <- markers_df[order(markers_df$FDR), ]
    markers_df <- reorder_id_cols(markers_df)
    dt_with_excel(markers_df, paste0("cluster_markers_cluster_", input$marker_cluster))
  }, server = FALSE) 
  
  # Contrast table: same idea as the marker table, but the row is chosen by the
  # "target vs reference" pair (split back apart here) plus the test type.
  output$contrast_table <- renderDT({
    req(input$contrast_choice)
    parts <- strsplit(input$contrast_choice, " vs ")[[1]]
    df <- seu()@misc$contrasts
    ct <- if ("wilcox" %in% df$test_type) "wilcox" else df$test_type[1]
    row <- df[df$target == parts[1] & df$reference == parts[2] & df$test_type == ct, ]
    validate(need(nrow(row) == 1, "No data for this combination."))
    contrast_df <- as.data.frame(row$markers[[1]])
    contrast_df <- contrast_df[order(contrast_df$FDR), ]
    contrast_df <- reorder_id_cols(contrast_df)
    dt_with_excel(contrast_df, safe_fname(paste0("contrasts_DEA_", input$contrast_choice)))
  }, server = FALSE) 

  # ===========================================================================
  # Generic Cluster Markers & Contrasts.
  # Used when the object has NO precomputed @misc tables: differential
  # expression is computed on demand with Seurat::FindMarkers, gated by a Run
  # button (the computation is heavy) and cached until Run is pressed again.
  # ===========================================================================

  # Format a FindMarkers() result for display: sort by adjusted p-value, prepend
  # a gene id column (+ a friendly label when rowData gives one), round numbers.
  format_generic_de <- function(df) {
    # Sort by adjusted p-value when available; the ROC test has no p-value, so
    # fall back to its "power" (higher = better).
    if ("p_val_adj" %in% names(df)) {
      df <- df[order(df$p_val_adj), , drop = FALSE]
    } else if ("power" %in% names(df)) {
      df <- df[order(df$power, decreasing = TRUE), , drop = FALSE]
    }
    feat <- rownames(df)
    out <- data.frame(gene = feat, stringsAsFactors = FALSE, check.names = FALSE)
    lbl <- gene_lookup()$label[match(feat, gene_lookup()$value)]
    if (!all(is.na(lbl)) && !isTRUE(all(lbl == feat))) out$label <- lbl
    out <- cbind(out, df)
    num <- vapply(out, is.numeric, logical(1))
    out[num] <- lapply(out[num], function(x) signif(x, 3))
    rownames(out) <- NULL
    out
  }

  # ---- Generic Cluster Markers: chosen group vs all other cells ----
  output$gmarker_cluster_ui <- renderUI({
    req(input$gmarker_group)
    lvls <- sort_levels(seu()@meta.data[[input$gmarker_group]])
    selectInput("gmarker_cluster", "Group / cluster:", choices = lvls, width = "100%")
  })

  # Auto-computes whenever the grouping column or group changes (cached until
  # then, so editing the search box only re-filters - no recompute).
  # Helper: prepare a Seurat object for FindMarkers on a given assay.
  # Sets the default assay, joins split layers (Seurat v5 integration creates
  # per-sample layers like counts.1/data.1 that must be merged before DE),
  # normalizes if a 'data' layer is missing (e.g. RNA assay that was never
  # NormalizeData'd - common in SCT-only integration pipelines), and runs
  # PrepSCTFindMarkers when SCT is chosen.
  prepare_obj_for_de <- function(obj, assay_use) {
    Seurat::DefaultAssay(obj) <- assay_use
    # Seurat v5: join split layers (counts.1, data.1, ...) into single matrices.
    if (exists("JoinLayers", where = asNamespace("SeuratObject"))) {
      obj[[assay_use]] <- tryCatch(
        suppressWarnings(SeuratObject::JoinLayers(obj[[assay_use]])),
        error = function(e) obj[[assay_use]]
      )
    }
    # If the assay has no 'data' layer (only raw counts), log-normalize on the
    # fly so FindMarkers can use slot = "data".  This is fast (~seconds) and
    # won't touch SCT (which already has a data layer).
    avail <- tryCatch(SeuratObject::Layers(obj[[assay_use]]), error = function(e) character(0))
    if (!"data" %in% avail && "counts" %in% avail) {
      obj <- tryCatch(
        suppressWarnings(Seurat::NormalizeData(obj, assay = assay_use, verbose = FALSE)),
        error = function(e) obj
      )
    }
    if (identical(assay_use, "SCT")) {
      obj <- tryCatch(
        suppressWarnings(Seurat::PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)),
        error = function(e) obj
      )
    }
    obj
  }

  # The DE tabs are fixed to the RNA assay + "data" (Normalized log); their assay
  # and data-layer selectors were removed. input$gmarker_layer / input$gcontrast_layer
  # are absent, so the slot logic below falls back to "data".

  # Gated behind the "Compute markers" button. FindMarkers over the whole object
  # is heavy (seconds to minutes on large data) and Shiny is single-threaded, so
  # running it on every dropdown change would freeze the app.
  # eventReactive only recomputes on click; the result is stamped with the inputs
  # it was computed for so the table can flag a stale selection (see gmarker_table).
  gmarker_res <- eventReactive(input$gmarker_compute, {
    req(seu(), input$gmarker_group, input$gmarker_cluster)
    assay_use <- gmarker_assay_val()
    req(assay_use)

    obj <- prepare_obj_for_de(seu(), assay_use)

    # Use the user-chosen layer, falling back to auto-detection by test type.
    slot <- if (!is.null(input$gmarker_layer) && nzchar(input$gmarker_layer)) {
      input$gmarker_layer
    } else if ("wilcox" %in% c("negbinom", "poisson", "DESeq2")) {
      "counts"
    } else {
      "data"
    }
    res <- tryCatch(
      # suppressWarnings: some tests (e.g. negbinom) emit one benign per-gene
      # convergence warning each - hundreds of lines. Errors are still caught.
      suppressWarnings(
        Seurat::FindMarkers(obj, ident.1 = input$gmarker_cluster, group.by = input$gmarker_group,
                            test.use = "wilcox", slot = slot, assay = assay_use)
      ),
      error = function(e) e
    )
    # Stamp what this was computed for, so an input change since can be detected.
    list(res = res, group = input$gmarker_group, cluster = input$gmarker_cluster,
         assay = assay_use, layer = input$gmarker_layer, test = "wilcox")
  })

  output$gmarker_table <- renderDT({
    # Before the first compute, prompt rather than hanging on the spinner.
    validate(need(isTRUE(input$gmarker_compute > 0),
                  "Set the options above, then click 'Compute markers'."))
    a <- gmarker_res()
    # Stale-guard: the cached result belongs to the inputs it was computed with.
    # If any changed since, ask for a recompute instead of showing mislabeled data.
    validate(need(identical(a$group, input$gmarker_group) &&
                  identical(a$cluster, input$gmarker_cluster) &&
                  identical(a$assay, gmarker_assay_val()) &&
                  identical(a$layer, input$gmarker_layer) &&
                  identical(a$test, "wilcox"),
                  "Inputs changed. Click 'Compute markers' to update the table."))
    res <- a$res
    if (inherits(res, "error")) validate(need(FALSE, paste("Could not compute markers:", conditionMessage(res))))
    out <- format_generic_de(res)
    dt_with_excel(out, safe_fname(paste0("cluster_markers_", input$gmarker_cluster)))
  }, server = FALSE)

  # ---- Generic Contrasts: target group vs reference group ----
  output$gcontrast_target_ui <- renderUI({
    req(input$gcontrast_group)
    lvls <- sort_levels(seu()@meta.data[[input$gcontrast_group]])
    selectInput("gcontrast_target", "Target (A):", choices = lvls, width = "100%")
  })
  output$gcontrast_reference_ui <- renderUI({
    req(input$gcontrast_group)
    lvls <- sort_levels(seu()@meta.data[[input$gcontrast_group]])
    selectInput("gcontrast_reference", "Reference (B):", choices = lvls,
                selected = if (length(lvls) > 1) lvls[2] else lvls[1], width = "100%")
  })

  # Auto-computes whenever the groups or test change (cached until then, so the
  # search box only re-filters - no recompute).
  # Gated behind the "Compute contrast" button (same reason as gmarker_res above:
  # FindMarkers is heavy and Shiny is single-threaded). Only recomputes on click;
  # stamped with its inputs so the table can flag a stale selection.
  gcontrast_res <- eventReactive(input$gcontrast_compute, {
    req(seu(), input$gcontrast_group, input$gcontrast_target, input$gcontrast_reference)
    validate(need(input$gcontrast_target != input$gcontrast_reference, "Pick two different groups."))
    assay_use <- gcontrast_assay_val()
    req(assay_use)

    obj <- prepare_obj_for_de(seu(), assay_use)

    # Use the user-chosen layer, falling back to auto-detection by test type.
    slot <- if (!is.null(input$gcontrast_layer) && nzchar(input$gcontrast_layer)) {
      input$gcontrast_layer
    } else if ("wilcox" %in% c("negbinom", "poisson", "DESeq2")) {
      "counts"
    } else {
      "data"
    }
    res <- tryCatch(
      # suppressWarnings: silence the benign per-gene convergence notes (e.g.
      # negbinom's theta.ml "iteration limit reached"). Errors still caught.
      suppressWarnings(
        Seurat::FindMarkers(obj, ident.1 = input$gcontrast_target,
                            ident.2 = input$gcontrast_reference, group.by = input$gcontrast_group,
                            test.use = "wilcox", slot = slot, assay = assay_use)
      ),
      error = function(e) e
    )
    list(res = res, group = input$gcontrast_group, target = input$gcontrast_target,
         reference = input$gcontrast_reference, assay = assay_use,
         layer = input$gcontrast_layer, test = "wilcox")
  })

  output$gcontrast_table <- renderDT({
    validate(need(isTRUE(input$gcontrast_compute > 0),
                  "Set the options above, then click 'Compute contrast'."))
    a <- gcontrast_res()
    validate(need(identical(a$group, input$gcontrast_group) &&
                  identical(a$target, input$gcontrast_target) &&
                  identical(a$reference, input$gcontrast_reference) &&
                  identical(a$assay, gcontrast_assay_val()) &&
                  identical(a$layer, input$gcontrast_layer) &&
                  identical(a$test, "wilcox"),
                  "Inputs changed. Click 'Compute contrast' to update the table."))
    res <- a$res
    if (inherits(res, "error")) validate(need(FALSE, paste("Could not compute contrast:", conditionMessage(res))))
    out <- format_generic_de(res)
    dt_with_excel(out, safe_fname(paste0("contrasts_DEA_", input$gcontrast_target, "_vs_", input$gcontrast_reference)))
  }, server = FALSE)

  # ---- Pseudobulk differential expression (DESeq2) --------------------------
  # Sum raw counts per sample within each group, then DESeq2 on the group factor.
  # The recommended way to compare conditions across samples (avoids the pseudo-
  # replication false positives of per-cell tests; Squair et al. 2021).

  # Optional "restrict to a cell type / cluster" level picker.
  output$gcontrast_pb_subset_level_ui <- renderUI({
    col <- input$gcontrast_pb_subset_col
    if (is.null(col) || !nzchar(col)) return(helpText("(no restriction)"))
    selectInput("gcontrast_pb_subset_level", "Restrict to level:",
                choices = sort_levels(seu()@meta.data[[col]]), width = "100%")
  })

  gcontrast_pb_res <- eventReactive(input$gcontrast_pb_compute, {
    req(seu(), input$gcontrast_group, input$gcontrast_target, input$gcontrast_reference, input$gcontrast_pb_sample)
    validate(need(input$gcontrast_target != input$gcontrast_reference, "Pick two different groups."))
    validate(need(requireNamespace("DESeq2", quietly = TRUE),
                  "The DESeq2 package is required for pseudobulk mode (install it on the server)."))
    assay_use <- gcontrast_assay_val()
    cnt <- tryCatch(GetAssayData(seu(), layer = "counts", assay = assay_use), error = function(e) NULL)
    validate(need(!is.null(cnt) && nrow(cnt) > 0, paste0("Assay '", assay_use, "' has no raw 'counts' layer for pseudobulk.")))
    md <- seu()@meta.data
    cond <- as.character(md[[input$gcontrast_group]])
    samp <- as.character(md[[input$gcontrast_pb_sample]])
    keep <- cond %in% c(input$gcontrast_target, input$gcontrast_reference) & !is.na(samp)
    sub_col <- input$gcontrast_pb_subset_col
    if (!is.null(sub_col) && nzchar(sub_col) && !is.null(input$gcontrast_pb_subset_level))
      keep <- keep & !is.na(md[[sub_col]]) & as.character(md[[sub_col]]) == input$gcontrast_pb_subset_level
    validate(need(sum(keep) > 0, "No cells match the chosen groups / restriction."))
    cnt <- cnt[, keep, drop = FALSE]; cond <- cond[keep]; samp <- samp[keep]

    # Sum counts per sample (genes x samples) via an indicator matrix.
    sf <- factor(samp)
    ind <- Matrix::sparse.model.matrix(~ 0 + sf)
    pb <- as.matrix(cnt %*% ind)
    colnames(pb) <- levels(sf)
    # Each sample's group = its cells' (majority) condition.
    samp_group <- tapply(cond, samp, function(x) names(sort(table(x), decreasing = TRUE))[1])
    groups <- samp_group[colnames(pb)]
    ok <- groups %in% c(input$gcontrast_target, input$gcontrast_reference)
    pb <- pb[, ok, drop = FALSE]; groups <- groups[ok]
    n_t <- sum(groups == input$gcontrast_target); n_r <- sum(groups == input$gcontrast_reference)
    validate(need(n_t >= 2 && n_r >= 2, sprintf(
      "Pseudobulk needs at least 2 samples per group (found %d in '%s', %d in '%s'). Pick a sample column that actually has replicates.",
      n_t, input$gcontrast_target, n_r, input$gcontrast_reference)))

    pb <- round(pb)
    pb <- pb[rowSums(pb) >= 10, , drop = FALSE]   # drop all-but-empty genes
    coldata <- data.frame(group = factor(groups, levels = c(input$gcontrast_reference, input$gcontrast_target)))
    rownames(coldata) <- colnames(pb)
    res <- tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = pb, colData = coldata, design = ~ group)
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      as.data.frame(DESeq2::results(dds, contrast = c("group", input$gcontrast_target, input$gcontrast_reference)))
    }, error = function(e) e)
    list(res = res, n_t = n_t, n_r = n_r, group = input$gcontrast_group,
         target = input$gcontrast_target, reference = input$gcontrast_reference,
         sample = input$gcontrast_pb_sample, assay = assay_use,
         subset = if (!is.null(sub_col) && nzchar(sub_col)) paste0(sub_col, " = ", input$gcontrast_pb_subset_level) else "all cells")
  })

  output$gcontrast_pb_caption <- renderUI({
    req(isTRUE(input$gcontrast_pb_compute > 0))
    a <- gcontrast_pb_res()
    if (inherits(a$res, "error")) return(NULL)
    helpText(sprintf("Pseudobulk DESeq2: %d samples in '%s' vs %d in '%s', on %s (assay '%s', raw counts). log2FoldChange is %s vs %s.",
                     a$n_t, a$target, a$n_r, a$reference, a$subset, a$assay, a$target, a$reference))
  })

  output$gcontrast_pb_table <- renderDT({
    validate(need(isTRUE(input$gcontrast_pb_compute > 0),
                  "Set the options above, then click 'Compute pseudobulk DE'."))
    a <- gcontrast_pb_res()
    validate(need(identical(a$group, input$gcontrast_group) &&
                  identical(a$target, input$gcontrast_target) &&
                  identical(a$reference, input$gcontrast_reference) &&
                  identical(a$sample, input$gcontrast_pb_sample) &&
                  identical(a$assay, gcontrast_assay_val()),
                  "Inputs changed. Click 'Compute pseudobulk DE' to update the table."))
    res <- a$res
    if (inherits(res, "error")) validate(need(FALSE, paste("DESeq2 failed:", conditionMessage(res))))
    df <- data.frame(gene = rownames(res), res, row.names = NULL, check.names = FALSE)
    df <- df[order(df$padj, na.last = TRUE), ]
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], function(x) signif(x, 3))
    dt_with_excel(df, safe_fname(paste0("pseudobulk_DESeq2_", a$target, "_vs_", a$reference)))
  }, server = FALSE)

}

shinyApp(ui, server)
