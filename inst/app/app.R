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

# Raise the upload size limit (default is 5MB; Seurat objects are much bigger)
options(shiny.maxRequestSize = 5 * 1024^3)  # 5 GB

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

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Hide all browser scrollbars (deletes the gray bar on the right) */
      ::-webkit-scrollbar { width: 0px; background: transparent; }
      body { -ms-overflow-style: none; scrollbar-width: none; }

      /* Force hide Plotly's modebar (the gray line with options) */
      .plotly .modebar { display: none !important; }

      /* Restyle DataTables' own 'processing' indicator (the default blue pulsing
         dots) to match the app's green spinner. It shows while a large table
         draws on the client -- after the server render that the shinycssloaders
         spinner covers -- so this is the loading icon users actually see. */
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

      /* Make each spinner-wrapped feature plot fill its (square) grid cell */
      #dynamic_feature_plots .shiny-spinner-output-container { width: 100% !important; height: 100% !important; }
      #dynamic_feature_plots .html-widget { width: 100% !important; height: 100% !important; }
      #umap_plot_box .shiny-spinner-output-container { width: 100% !important; height: 100% !important; }
      #umap_plot_box .html-widget { width: 100% !important; height: 100% !important; }
    ")),
    tags$script(HTML("
      document.addEventListener('DOMContentLoaded', function() {
        // ---- Load Test render timer ----
        // When the 'Copies' slider changes, a full rebuild of 3*N plots follows.
        // Start a clock, then count the lt_ plot renders (shiny:value) until all
        // 3*N have drawn, and report the elapsed wall-clock time.
        var ltStart = null, ltExpected = 0, ltCount = 0;
        $(document).on('shiny:inputchanged', function(e) {
          if (e.name === 'lt_copies') {
            ltStart = performance.now();
            ltExpected = 3 * e.value;
            ltCount = 0;
            var t = document.getElementById('lt_timer');
            if (t) t.textContent = 'Rendering ' + e.value + ' copies (' + ltExpected + ' plots)...';
          }
        });
        $(document).on('shiny:value', function(e) {
          if (ltStart === null || !e.name) return;
          if (e.name.indexOf('lt_umap_') === 0 || e.name.indexOf('lt_gene_') === 0 ||
              e.name.indexOf('lt_heat_') === 0) {
            ltCount++;
            if (ltCount >= ltExpected) {
              var secs = ((performance.now() - ltStart) / 1000).toFixed(1);
              var t = document.getElementById('lt_timer');
              if (t) t.textContent = 'Rendered ' + (ltExpected / 3) + ' copies (' +
                ltExpected + ' plots) in ' + secs + ' s';
              ltStart = null;
            }
          }
        });

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
    uiOutput("data_status")
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

  # gene_avg_expr(): mean log-expression per gene, used to offer the
  # "sort genes by highest/lowest average expression" options.
  gene_avg_expr <- reactive({
    req(seu())
    Matrix::rowMeans(GetAssayData(seu(), layer = "data"))
  })

  # found_genes(): the subset of the user-selected genes that actually exist in
  # the object. shared_upper(): the pooled max expression across those genes,
  # used by BOTH the individual plots and the single shared colorbar so their
  # scales match in "Same scale for all genes" mode.
  found_genes <- reactive({
    req(seu(), input$genes)
    input$genes[input$genes %in% rownames(seu())]
  })

  shared_upper <- reactive({
    fg <- found_genes()
    req(length(fg) > 0)
    # Read straight from the data layer (faster than FetchData's data.frame path).
    all_expr <- as.vector(GetAssayData(seu(), layer = "data")[fg, , drop = FALSE])
    upper <- max(all_expr, na.rm = TRUE)              # true maximum across selected genes
    if (!is.finite(upper) || upper <= 0) upper <- 1   # guard: all-zero expression
    as.numeric(upper)
  })

  # One common square axis window (padded full extent of the UMAP embedding),
  # applied to BOTH axes of EVERY plot so all UMAPs are framed identically with
  # the same x/y range and scale.
  umap_range <- reactive({
    req(seu())
    emb <- Embeddings(seu(), "umap")[, 1:2]
    r <- range(emb, na.rm = TRUE)
    pad <- diff(r) * 0.03
    c(r[1] - pad, r[2] + pad)
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

  # app_body(): builds the whole tabbed UI once data is loaded. The last two
  # tabs (Cluster Markers, Contrasts) are only added when the object carries the
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
    # sense to colour by categorically -- factor / character / logical columns,
    # plus low-cardinality numerics (e.g. cluster numbers). Continuous numerics
    # with many distinct values (nCount, nFeature, doublet_score, percent.mt, ...)
    # are excluded, because the UMAP assigns one colour per value and those would
    # produce thousands of categories.
    # TWEAK: to offer literally EVERY column instead, use `grouping_cols <- all_cols`.
    # To raise/lower what counts as "categorical", change the 50 below.
    grouping_cols <- all_cols[vapply(seu()@meta.data, function(x) {
      is.factor(x) || is.character(x) || is.logical(x) ||
        (is.numeric(x) && length(unique(x[!is.na(x)])) <= 50)
    }, logical(1))]
    # Friendly, A-Z-sorted labels for the dropdown. Named vector: names are shown,
    # values (the real column names) are what input$group_var returns.
    grouping_choices <- setNames(grouping_cols, pretty_label(grouping_cols))
    grouping_choices <- grouping_choices[order(names(grouping_choices))]
    has_markers <- !is.null(seu()@misc$cluster_markers)
    has_contrasts <- !is.null(seu()@misc$contrasts)

    # DE tests offered by the generic (scdrake-free) Markers/Contrasts tabs.
    # All built into Seurat and fast enough to auto-compute. (poisson uses raw
    # counts; roc reports an AUC "power" instead of a p-value.)
    # NOTE: the slow per-gene tests (negbinom, MAST, DESeq2) were removed for now
    # -- they take minutes and these tabs auto-compute. Re-add here if wanted
    # (the counts-slot logic below already handles negbinom/DESeq2).
    de_test_choices <- c("Wilcoxon" = "wilcox", "Student's t" = "t", "ROC / AUC" = "roc",
                         "Logistic regression" = "LR", "Likelihood ratio (bimod)" = "bimod",
                         "Poisson (counts)" = "poisson")

    tabs <- list(
      tabPanel(
        "Clusters / Metadata",
        sidebarLayout(
          sidebarPanel(
            selectInput("group_var", "Color UMAP by:", choices = grouping_choices),
            uiOutput("category_filter_ui")
          ),
          mainPanel(
            div(style = "display: flex; justify-content: flex-start; width: 100%; overflow: hidden;",
                div(id = "umap_plot_box", style = "width: calc(min(100vh - 230px, 100vw - 800px) + 240px); height: min(100vh - 230px, 100vw - 800px);",
                    withSpinner(plotlyOutput("umap_plot", width = "100%", height = "100%"), caption = "Rendering plot...", color = "#39FF14"))
            )
          )
        )
      ),
      tabPanel(
        "Gene Expression",
        wellPanel(
          fluidRow(
            column(3, selectInput("gene_sort", "Sort genes by:", width = "100%",
                                  choices = c("Alphabetical", "Highest average expression", "Lowest average expression"))),
            column(3, selectizeInput("genes", "Gene(s):", choices = NULL, multiple = TRUE, width = "100%",
                                     options = list(placeholder = "Type one or more gene symbols..."))),
            # TWEAK: default / allowed plots-per-row for the Gene Expression grid.
            column(2, numericInput("feature_ncol", "Plots per row:", value = 2, min = 1, max = 4, width = "100%")),
            column(4, radioButtons("scale_mode", "Color scale:",
                                   choices = c("Same scale for all genes (comparable)" = "shared",
                                               "Scale adjusted per gene (all clearly visible)" = "per_gene"),
                                   selected = "shared"))
          )
        ),
        div(style = "width: 100%; display: flex; justify-content: center;",
            uiOutput("dynamic_feature_plots"))
      )
    )

    # Cluster Markers tab. If the object carries scdrake's precomputed
    # @misc$cluster_markers we show those (instant). Otherwise we compute markers
    # on demand from any Seurat object (FindMarkers, behind a Run button).
    if (has_markers) {
      tabs <- c(tabs, list(tabPanel(
        "Cluster Markers",
        wellPanel(
          fluidRow(
            column(3, selectInput("marker_resolution", "Clustering resolution:", width = "100%",
                                  choices = unique(seu()@misc$cluster_markers$source_column))),
            column(2, selectInput("marker_test", "Test type:", width = "100%",
                                  choices = unique(seu()@misc$cluster_markers$test_type))),
            column(2, uiOutput("marker_cluster_ui")),
            column(3, selectInput("marker_search_field", "Search in:", width = "100%",
                                  choices = c("All columns", "ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))),
            column(2, textInput("marker_search_text", "Search:", width = "100%", placeholder = "Type to search..."))
          )
        ),
        withSpinner(DTOutput("marker_table"), caption = "Loading table...", color = "#39FF14")
      )))
    } else {
      tabs <- c(tabs, list(tabPanel(
        "Cluster Markers",
        wellPanel(
          fluidRow(
            column(3, selectInput("gmarker_group", "Grouping column:", choices = grouping_choices, width = "100%")),
            column(3, uiOutput("gmarker_cluster_ui")),
            column(3, selectInput("gmarker_test", "Test type:", choices = de_test_choices, width = "100%")),
            column(3, textInput("gmarker_search_text", "Search:", width = "100%", placeholder = "Type to search..."))
          ),
          helpText("Markers of the chosen group vs all other cells. Computed automatically.")
        ),
        withSpinner(DTOutput("gmarker_table"), caption = "Computing markers...", color = "#39FF14")
      )))
    }

    # Contrasts tab. Precomputed scdrake @misc$contrasts if present, otherwise
    # differential expression between two chosen groups computed on demand.
    if (has_contrasts) {
      tabs <- c(tabs, list(tabPanel(
        "Contrasts",
        wellPanel(
          fluidRow(
            column(4, selectInput("contrast_choice", "Contrast:", width = "100%",
                                  choices = unique(paste0(seu()@misc$contrasts$target, " vs ", seu()@misc$contrasts$reference)))),
            column(2, selectInput("contrast_test", "Test type:", width = "100%",
                                  choices = unique(seu()@misc$contrasts$test_type))),
            column(3, selectInput("contrast_search_field", "Search in:", width = "100%",
                                  choices = c("All columns", "ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))),
            column(3, textInput("contrast_search_text", "Search:", width = "100%", placeholder = "Type to search..."))
          )
        ),
        withSpinner(DTOutput("contrast_table"), caption = "Loading table...", color = "#39FF14")
      )))
    } else {
      tabs <- c(tabs, list(tabPanel(
        "Contrasts",
        wellPanel(
          fluidRow(
            column(3, selectInput("gcontrast_group", "Grouping column:", choices = grouping_choices, width = "100%")),
            column(2, uiOutput("gcontrast_target_ui")),
            column(2, uiOutput("gcontrast_reference_ui")),
            column(2, selectInput("gcontrast_test", "Test type:", choices = de_test_choices, width = "100%")),
            column(3, textInput("gcontrast_search_text", "Search:", width = "100%", placeholder = "Type to search..."))
          ),
          helpText("Differential expression, target vs reference group. Computed automatically.")
        ),
        withSpinner(DTOutput("gcontrast_table"), caption = "Computing contrast...", color = "#39FF14")
      )))
    }

    # Aggregated Expression (pseudobulk): summarise the selected genes' expression
    # per level of a chosen metadata column (condition). Always available -- it
    # only needs meta.data + an expression layer, so it works on any object.
    # Controls sit in a compact row across the top so the heatmap + table below
    # get the FULL page width (and a taller plot), rather than a narrow sidebar.
    tabs <- c(tabs, list(tabPanel(
      "Aggregated Expression",
      wellPanel(
        fluidRow(
          column(3, selectInput("agg_group_var", "Aggregate by:",
                                choices = grouping_choices, width = "100%")),
          column(3, selectInput("agg_gene_sort", "Sort genes by:", width = "100%",
                                choices = c("Alphabetical", "Highest average expression", "Lowest average expression"))),
          column(4, selectizeInput("agg_genes", "Gene(s):", choices = NULL, multiple = TRUE, width = "100%",
                                   options = list(placeholder = "Type one or more gene symbols..."))),
          column(2, radioButtons("agg_fun", "Summary per group:",
                                 choices = c("Mean" = "mean", "Sum" = "sum", "Median" = "median"),
                                 selected = "mean"))
        )
      ),
      withSpinner(plotlyOutput("agg_heatmap", height = "600px"),
                  caption = "Rendering plot...", color = "#39FF14"),
      tags$hr(),
      withSpinner(DTOutput("agg_table"), caption = "Loading table...", color = "#39FF14")
    )))

    # ---- Load Test -------------------------------------------------------
    # Renders N copies of a plot set (Clusters UMAP + Gene-expression UMAP +
    # Aggregated heatmap) stacked on ONE page, so all copies render at once
    # (separate tabs would suspend the hidden ones and not stress rendering).
    # A JS timer (#lt_timer) reports how long N copies took to draw.
    tabs <- c(tabs, list(tabPanel(
      "Load Test",
      wellPanel(
        fluidRow(
          column(3, sliderInput("lt_copies", "Copies of the plot set:",
                                min = 1, max = 20, value = 3, step = 1, width = "100%")),
          column(3, selectInput("lt_group", "UMAP color / aggregate by:",
                                choices = grouping_choices, width = "100%")),
          column(6, selectizeInput("lt_genes", "Gene(s) (expression + heatmap):",
                                   choices = NULL, multiple = TRUE, width = "100%",
                                   options = list(placeholder = "Type one or more gene symbols...")))
        ),
        div(id = "lt_timer", style = "font-weight: 600; margin-top: 5px; color: #2ecc0f;",
            "Move the 'Copies' slider to render and time the plots.")
      ),
      uiOutput("lt_blocks")
    )))

    do.call(tabsetPanel, tabs)
  })

  # The per-category show/hide checkboxes under "Color UMAP by". Rebuilt every
  # time the colouring column changes; NA (unlabelled) cells become an "NA" box.
  output$category_filter_ui <- renderUI({
    req(input$group_var)
    vals <- as.character(seu()@meta.data[[input$group_var]])
    vals[is.na(vals)] <- "NA"   # un-annotated cells (e.g. pruned label columns)
    vals <- sort(unique(vals))
    checkboxGroupInput("categories_shown", "Show categories:", choices = vals, selected = vals)
  })

  # ---------------------------------------------------------------------------
  # Clusters / Metadata UMAP (one plotly scatter coloured by a metadata column)
  # ---------------------------------------------------------------------------

  # Base data for the UMAP, in a fixed (seeded) random point order. Depends ONLY
  # on the object + colouring column -- NOT on which categories are ticked -- so
  # that toggling categories doesn't rebuild it (the observer below just
  # recolours the existing plot in place).
  umap_base <- reactive({
    req(input$group_var)
    emb <- Embeddings(seu(), "umap")
    df <- data.frame(
      UMAP_1 = emb[, 1],
      UMAP_2 = emb[, 2],
      group  = as.character(seu()@meta.data[[input$group_var]]),
      stringsAsFactors = FALSE
    )
    df$group[is.na(df$group)] <- "NA"   # un-annotated cells -> explicit "NA" category
    set.seed(1)
    df <- df[sample(nrow(df)), ]        # fixed interleaved order (avoids category stacking)
    df$tooltip <- paste0("Group: ", df$group)

    # Colours must be hex -- plotly.js renders R names like "grey85" as black.
    cats <- sort(unique(df$group))
    pal <- scales::hue_pal()(length(cats))
    names(pal) <- cats
    if ("NA" %in% cats) pal[["NA"]] <- "#B3B3B3"   # unassigned cells (grey70)

    list(df = df, cats = cats, palette = pal)
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

  output$umap_plot <- renderPlotly({
    base <- umap_base()
    df <- base$df
    # Initial colours reflect the current selection but WITHOUT a reactive
    # dependency on it (isolate) -- the observer below updates them in place.
    cols <- umap_point_cols(base, isolate(input$categories_shown))
    df$col <- cols$points

    # A SINGLE marker trace (index 0) holds every point, so toggling a category
    # is just a colour restyle -- the plot never reloads. (plotly's `color=`
    # would split into one trace per category and stack them, hiding one.)
    p <- plot_ly(
      df, x = ~UMAP_1, y = ~UMAP_2,
      type = "scattergl", mode = "markers",
      marker = list(size = 4, color = ~col),   # TWEAK: UMAP point size
      text = ~tooltip, hoverinfo = "text", showlegend = FALSE
    )

    # Legend = one invisible proxy trace per category (traces 1..N), ordered by
    # frequency (most common first). Each carries one tiny (size 0.1) real point
    # (an all-NA proxy would be dropped); itemsizing keeps swatches full size.
    legend_order <- names(sort(table(df$group), decreasing = TRUE))
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
      p,
      title = list(text = paste("UMAP colored by", pretty_label(input$group_var)), x = 0, xanchor = "left",
                   font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
      xaxis = list(title = "UMAP_1", zeroline = FALSE, range = umap_range(),
                   constrain = "domain", constraintoward = "left"),
      yaxis = list(title = "UMAP_2", zeroline = FALSE, scaleanchor = "x", scaleratio = 1,
                   constrain = "domain", constraintoward = "top", range = umap_range()),
      showlegend = TRUE,
      legend = list(title = list(text = pretty_label(input$group_var)), itemsizing = "constant",
                    y = 0.5, yanchor = "middle"),
      margin = list(l = 50, r = 10, t = 40, b = 40)
    )
    plotly::config(p, displayModeBar = FALSE)
  })

  # Toggling a category recolours the EXISTING plot via plotly restyle -- no
  # re-render, no "Rendering plot..." spinner. Trace 0 = points; traces 1..N =
  # the legend swatches (frequency order), recoloured to match the selection.
  observeEvent(input$categories_shown, {
    req(input$group_var)
    base <- umap_base()
    cols <- umap_point_cols(base, input$categories_shown)
    legend_order <- names(sort(table(base$df$group), decreasing = TRUE))

    proxy <- plotlyProxy("umap_plot", session)
    # Points trace (index 0): replace the whole colour vector.
    plotlyProxyInvoke(proxy, "restyle", list(marker.color = list(cols$points)), list(0))
    # Legend swatches (traces 1..N): grey out / restore each category.
    leg_cols <- unname(cols$legend[legend_order])
    plotlyProxyInvoke(proxy, "restyle",
                      list(marker.color = as.list(leg_cols)),
                      as.list(seq_along(legend_order)))
  }, ignoreInit = TRUE)

  # ---- Square grid for Gene Expression ----
  # Each plot is a square whose side is the SMALLER of the width available per
  # column and the height available per row (viewport-aware). Using the same
  # value for the column and row tracks makes every cell a true square (so
  # plotly's equal-scale axes render correctly, single plot included), and the
  # height term keeps the whole grid on screen -- 4 plots fit without scrolling.
  output$dynamic_feature_plots <- renderUI({
    req(input$genes)
    found <- input$genes[input$genes %in% rownames(seu())]
    if (length(found) == 0) return(NULL)

    n <- length(found)
    ncol <- max(1, min(input$feature_ncol, n))
    nrow <- ceiling(n / ncol)
    shared <- identical(input$scale_mode, "shared")

    # TWEAK: these offsets tune the (fixed) square plot size. If plots feel too
    # small/large, adjust side_offset (horizontal) or top_offset (the height cap).
    # (Controls now sit in a row ABOVE the plots, not in a left sidebar, so the
    # plots use nearly the full page width -- small side_offset.)
    #
    # The size is driven by the width available per column and CAPPED to one
    # viewport height -- it does NOT depend on the number of rows, so each plot
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
      withSpinner(
        plotlyOutput(paste0("feature_plot_individual_", i), width = "100%", height = "100%"),
        caption = "Rendering plot...", color = "#39FF14"
      )
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
    if (shared) {
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
  })

  # Single shared colorbar, shown to the right of the grid in "shared" scale mode
  output$shared_colorbar <- renderPlotly({
    req(identical(input$scale_mode, "shared"))
    upper <- shared_upper()
    plot_ly(
      x = c(0, 0), y = c(0, 1),
      type = "scatter", mode = "markers",
      marker = list(
        color = c(0, upper),
        # TWEAK: keep this gradient identical to the per-plot one above so the
        # shared colorbar matches the points (grey85 low -> green high).
        colorscale = list(c(0, "#D9D9D9"), c(1, "#1E8B0D")),   # grey85 -> green
        cmin = 0, cmax = upper,
        opacity = 0,                                           # hide the dummy points
        showscale = TRUE,
        colorbar = list(title = list(text = "Expression"), thickness = 15,
                        len = 0.5, outlinewidth = 0)
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

  # One value-stable reactive per plot slot, returning the gene assigned to that
  # slot (or NA when unused). Because each returns a single value, adding a gene
  # only changes the NEW slot -- existing slots return the same gene, so their
  # plots are NOT invalidated and keep their cached render (only the added plot
  # renders). Depending on input$genes directly, as before, re-ran ALL of them.
  slot_gene <- lapply(1:max_plots, function(i) {
    force(i)
    reactive({
      fg <- found_genes()
      if (i <= length(fg)) fg[[i]] else NA_character_
    })
  })

  lapply(1:max_plots, function(i) {
    force(i)
    output[[paste0("feature_plot_individual_", i)]] <- renderPlotly({
      g <- slot_gene[[i]]()
      req(!is.na(g))

      umap_coords <- Embeddings(seu(), "umap")
      # Fetch ONLY this slot's gene, so the plot depends on `g` alone (not the
      # whole gene list) -- adding another gene won't invalidate this plot.
      # GetAssayData(layer="data")[g, ] is faster than FetchData for one gene.
      plot_df <- data.frame(
        UMAP_1 = umap_coords[, 1],
        UMAP_2 = umap_coords[, 2],
        expression = as.numeric(GetAssayData(seu(), layer = "data")[g, ]),
        check.names = FALSE
      )

      # Determine limits for Color Scale
      if (input$scale_mode == "shared") {
        upper <- shared_upper()
      } else {
        upper <- max(plot_df$expression, na.rm = TRUE)   # true per-gene maximum
        if (!is.finite(upper) || upper <= 0) upper <- 1  # guard: all-zero expression
      }

      # Draw low-expression points first so high-expression ones sit on top.
      plot_df <- plot_df[order(plot_df$expression), ]
      plot_df$tooltip <- paste0("Gene: ", g, "<br>Expression: ", round(plot_df$expression, 3))

      # Per-gene mode shows each plot's own colorbar; shared mode uses the single
      # shared colorbar to the right of the grid, so hide the per-plot one.
      show_bar <- input$scale_mode != "shared"

      # Native plot_ly (scattergl) instead of ggplot + ggplotly -- much faster.
      # cmin/cmax = 0..upper reproduces the old scale_color_gradient(limits =,
      # oob = squish): values above `upper` clamp to the max colour.
      p <- plot_ly(
        plot_df, x = ~UMAP_1, y = ~UMAP_2,
        type = "scattergl", mode = "markers",
        marker = list(
          size = 4,                                              # TWEAK: gene-expression point size
          color = ~expression,
          colorscale = list(c(0, "#D9D9D9"), c(1, "#1E8B0D")),   # TWEAK: grey85 -> green gradient
          cmin = 0, cmax = upper,
          showscale = show_bar,
          colorbar = list(title = list(text = "Expression"), thickness = 15, len = 0.5, outlinewidth = 0)
        ),
        text = ~tooltip, hoverinfo = "text"
      )

      # Same framing as the Clusters/Metadata UMAP: fixed square window, no tick
      # chrome, left-aligned title (matching font/size).
      p <- plotly::layout(
        p,
        title = list(text = g, x = 0, xanchor = "left",
                     font = list(size = 15, family = "Open Sans, verdana, arial, sans-serif")),
        xaxis = list(title = "UMAP_1", zeroline = FALSE, showline = FALSE, ticks = "",
                     range = umap_range(), constrain = "domain", constraintoward = "left"),
        yaxis = list(title = "UMAP_2", zeroline = FALSE, showline = FALSE, ticks = "",
                     scaleanchor = "x", scaleratio = 1, constrain = "domain",
                     constraintoward = "top", range = umap_range()),
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

  # agg_matrix(): the genes x group-levels matrix of aggregated expression.
  agg_matrix <- reactive({
    req(input$agg_group_var)
    validate(need(length(input$agg_genes) > 0, "Select one or more genes to aggregate."))
    genes <- input$agg_genes[input$agg_genes %in% rownames(seu())]
    validate(need(length(genes) > 0, "None of the selected genes are in this object."))

    fun <- switch(input$agg_fun, sum = sum, median = stats::median, mean)
    expr <- FetchData(seu(), vars = genes)                 # cells x genes (data layer)
    grp <- as.character(seu()@meta.data[[input$agg_group_var]])
    grp[is.na(grp)] <- "NA"                                # un-annotated cells -> "NA" group
    groups <- sort(unique(grp))

    mat <- matrix(NA_real_, nrow = length(genes), ncol = length(groups),
                  dimnames = list(genes, groups))
    # tapply groups once per gene (vs subsetting per gene x group in a nested loop).
    for (g in genes) {
      agg_vals <- tapply(expr[[g]], grp, fun)
      mat[g, names(agg_vals)] <- agg_vals
    }
    mat
  })

  # Short label for the summary function (colorbar / column headers).
  agg_fun_label <- reactive({
    switch(input$agg_fun, sum = "Sum", median = "Median", "Mean")
  })

  # Heatmap: genes (rows) x group levels (columns), grey85 -> green like the
  # Gene Expression plots. Colour scale auto-fits the aggregated values.
  output$agg_heatmap <- renderPlotly({
    mat <- agg_matrix()
    plot_ly(
      x = colnames(mat), y = rownames(mat), z = mat,
      type = "heatmap",
      colorscale = list(c(0, "#D9D9D9"), c(1, "#1E8B0D")),   # grey85 -> green
      colorbar = list(title = list(text = agg_fun_label()), thickness = 15, outlinewidth = 0),
      hovertemplate = paste0("Gene: %{y}<br>", pretty_label(input$agg_group_var),
                             ": %{x}<br>", agg_fun_label(), ": %{z:.3f}<extra></extra>")
    ) %>%
      plotly::layout(
        title = list(text = paste(agg_fun_label(), "expression by", pretty_label(input$agg_group_var)),
                     x = 0, xanchor = "left",
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
    datatable(df, rownames = FALSE,
              options = list(pageLength = 15, searching = FALSE, scrollX = TRUE,
                             processing = TRUE, language = list(processing = "Loading table...")))
  })

  # ===========================================================================
  # Cluster Markers & Contrasts tabs (interactive DT tables).
  # Both read data frames from obj@misc (scdrake format): each row is one
  # cluster/contrast whose `markers` cell holds a per-gene stats table.
  # ===========================================================================

  # Cluster dropdown for the Markers tab, numerically sorted (1,2,...,10 not 1,10,2).
  output$marker_cluster_ui <- renderUI({
    req(input$marker_resolution, input$marker_test)
    df <- seu()@misc$cluster_markers
    df <- df[df$source_column == input$marker_resolution & df$test_type == input$marker_test, ]
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

  # Helper: filter a table by the "Search in:" field + text box. "All columns"
  # matches any column; otherwise matches just the chosen column (case-insensitive).
  apply_field_search <- function(df, field, text) {
    if (is.null(text) || text == "") return(df)
    # Case-insensitive literal search: lowercase BOTH the query and the data.
    # (grepl's `ignore.case` is silently ignored when `fixed = TRUE`, so we can't
    # rely on it -- lowercasing both sides is the reliable way.)
    text <- tolower(text)
    if (field == "All columns") {
      matches <- apply(df, 1, function(row) any(grepl(text, tolower(row), fixed = TRUE)))
    } else if (field %in% colnames(df)) {
      matches <- grepl(text, tolower(as.character(df[[field]])), fixed = TRUE)
    } else {
      matches <- rep(TRUE, nrow(df))
    }
    df[matches, , drop = FALSE]
  }

  # Marker table: pick the row matching (resolution, test, cluster), pull its
  # per-gene stats table, sort by FDR, reorder columns, apply the search filter.
  # TWEAK: pageLength controls rows per page.
  output$marker_table <- renderDT({
    req(input$marker_resolution, input$marker_test, input$marker_cluster)
    df <- seu()@misc$cluster_markers
    row <- df[df$source_column == input$marker_resolution & df$test_type == input$marker_test & df$group_level == input$marker_cluster, ]
    validate(need(nrow(row) == 1, "No data for this combination."))
    markers_df <- as.data.frame(row$markers[[1]])
    markers_df <- markers_df[order(markers_df$FDR), ]
    markers_df <- reorder_id_cols(markers_df)
    markers_df <- apply_field_search(markers_df, input$marker_search_field, input$marker_search_text)
    datatable(markers_df, rownames = FALSE,
              options = list(pageLength = 15, searching = FALSE,
                             processing = TRUE, language = list(processing = "Loading table...")))
  })

  # Contrast table: same idea as the marker table, but the row is chosen by the
  # "target vs reference" pair (split back apart here) plus the test type.
  output$contrast_table <- renderDT({
    req(input$contrast_choice, input$contrast_test)
    parts <- strsplit(input$contrast_choice, " vs ")[[1]]
    df <- seu()@misc$contrasts
    row <- df[df$target == parts[1] & df$reference == parts[2] & df$test_type == input$contrast_test, ]
    validate(need(nrow(row) == 1, "No data for this combination."))
    contrast_df <- as.data.frame(row$markers[[1]])
    contrast_df <- contrast_df[order(contrast_df$FDR), ]
    contrast_df <- reorder_id_cols(contrast_df)
    contrast_df <- apply_field_search(contrast_df, input$contrast_search_field, input$contrast_search_text)
    datatable(contrast_df, rownames = FALSE,
              options = list(pageLength = 15, searching = FALSE,
                             processing = TRUE, language = list(processing = "Loading table...")))
  })

  # ===========================================================================
  # Generic (scdrake-free) Cluster Markers & Contrasts.
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

  gd_table_opts <- list(pageLength = 15, searching = FALSE, scrollX = TRUE,
                        processing = TRUE, language = list(processing = "Loading table..."))

  # ---- Generic Cluster Markers: chosen group vs all other cells ----
  output$gmarker_cluster_ui <- renderUI({
    req(input$gmarker_group)
    lvls <- sort(unique(as.character(seu()@meta.data[[input$gmarker_group]])))
    selectInput("gmarker_cluster", "Group / cluster:", choices = lvls, width = "100%")
  })

  # Auto-computes whenever the grouping column or group changes (cached until
  # then, so editing the search box only re-filters -- no recompute).
  gmarker_res <- reactive({
    req(seu(), input$gmarker_group, input$gmarker_cluster, input$gmarker_test)
    # negbinom/poisson/DESeq2 are count-based; the rest use the "data" layer.
    slot <- if (input$gmarker_test %in% c("negbinom", "poisson", "DESeq2")) "counts" else "data"
    tryCatch(
      # suppressWarnings: some tests (e.g. negbinom) emit one benign per-gene
      # convergence warning each -- hundreds of lines. Errors are still caught.
      suppressWarnings(
        Seurat::FindMarkers(seu(), ident.1 = input$gmarker_cluster, group.by = input$gmarker_group,
                            test.use = input$gmarker_test, slot = slot)
      ),
      error = function(e) e
    )
  })

  output$gmarker_table <- renderDT({
    res <- gmarker_res()
    if (inherits(res, "error")) validate(need(FALSE, paste("Could not compute markers:", conditionMessage(res))))
    out <- format_generic_de(res)
    out <- apply_field_search(out, "All columns", input$gmarker_search_text)
    datatable(out, rownames = FALSE, options = gd_table_opts)
  })

  # ---- Generic Contrasts: target group vs reference group ----
  output$gcontrast_target_ui <- renderUI({
    req(input$gcontrast_group)
    lvls <- sort(unique(as.character(seu()@meta.data[[input$gcontrast_group]])))
    selectInput("gcontrast_target", "Target (A):", choices = lvls, width = "100%")
  })
  output$gcontrast_reference_ui <- renderUI({
    req(input$gcontrast_group)
    lvls <- sort(unique(as.character(seu()@meta.data[[input$gcontrast_group]])))
    selectInput("gcontrast_reference", "Reference (B):", choices = lvls,
                selected = if (length(lvls) > 1) lvls[2] else lvls[1], width = "100%")
  })

  # Auto-computes whenever the groups or test change (cached until then, so the
  # search box only re-filters -- no recompute).
  gcontrast_res <- reactive({
    req(seu(), input$gcontrast_group, input$gcontrast_target, input$gcontrast_reference, input$gcontrast_test)
    validate(need(input$gcontrast_target != input$gcontrast_reference, "Pick two different groups."))
    slot <- if (input$gcontrast_test %in% c("negbinom", "poisson", "DESeq2")) "counts" else "data"
    tryCatch(
      # suppressWarnings: silence the benign per-gene convergence notes (e.g.
      # negbinom's theta.ml "iteration limit reached"). Errors still caught.
      suppressWarnings(
        Seurat::FindMarkers(seu(), ident.1 = input$gcontrast_target,
                            ident.2 = input$gcontrast_reference, group.by = input$gcontrast_group,
                            test.use = input$gcontrast_test, slot = slot)
      ),
      error = function(e) e
    )
  })

  output$gcontrast_table <- renderDT({
    res <- gcontrast_res()
    if (inherits(res, "error")) validate(need(FALSE, paste("Could not compute contrast:", conditionMessage(res))))
    out <- format_generic_de(res)
    out <- apply_field_search(out, "All columns", input$gcontrast_search_text)
    datatable(out, rownames = FALSE, options = gd_table_opts)
  })

  # ===========================================================================
  # Load Test: N stacked copies of a 3-plot set, all rendered at once, to
  # benchmark how the app scales. Shared data prep (computed once, drawn N
  # times) so the timer mostly reflects RENDERING cost.
  # ===========================================================================
  lt_max_copies <- 20

  # Populate the Load Test gene selector; default-select the first gene so the
  # gene/heatmap plots always render (keeps the render-count timer accurate).
  observeEvent(seu(), {
    req(seu())
    lookup <- gene_lookup()
    ordered <- lookup[order(lookup$value), , drop = FALSE]
    updateSelectizeInput(session, "lt_genes", choices = setNames(ordered$value, ordered$label),
                         selected = ordered$value[1], server = TRUE)
  })

  # Shared, once-computed data for all copies.
  lt_umap_df <- reactive({
    req(seu(), input$lt_group)
    emb <- Embeddings(seu(), "umap")
    grp <- as.character(seu()@meta.data[[input$lt_group]]); grp[is.na(grp)] <- "NA"
    cats <- sort(unique(grp))
    pal <- scales::hue_pal()(length(cats)); names(pal) <- cats
    data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
               col = unname(pal[grp]), stringsAsFactors = FALSE)
  })

  lt_heat_mat <- reactive({
    req(seu(), input$lt_group, input$lt_genes)
    genes <- input$lt_genes[input$lt_genes %in% rownames(seu())]
    req(length(genes) > 0)
    expr <- FetchData(seu(), vars = genes)
    grp <- as.character(seu()@meta.data[[input$lt_group]]); grp[is.na(grp)] <- "NA"
    groups <- sort(unique(grp))
    mat <- matrix(NA_real_, length(genes), length(groups), dimnames = list(genes, groups))
    for (g in genes) { agg_vals <- tapply(expr[[g]], grp, mean); mat[g, names(agg_vals)] <- agg_vals }
    mat
  })

  # One row of three plotly outputs per copy; N rows stacked -> all render.
  output$lt_blocks <- renderUI({
    req(input$lt_copies)
    rows <- lapply(seq_len(input$lt_copies), function(i) {
      div(style = "display: flex; gap: 10px; margin-bottom: 10px; height: 300px;",
          div(style = "flex: 1;", plotlyOutput(paste0("lt_umap_", i), height = "100%")),
          div(style = "flex: 1;", plotlyOutput(paste0("lt_gene_", i), height = "100%")),
          div(style = "flex: 1;", plotlyOutput(paste0("lt_heat_", i), height = "100%")))
    })
    do.call(tagList, rows)
  })

  # Pre-register up to lt_max_copies of each plot; each renders only when its
  # index is within the current copy count.
  lapply(seq_len(lt_max_copies), function(i) {
    force(i)
    output[[paste0("lt_umap_", i)]] <- renderPlotly({
      req(i <= input$lt_copies)
      df <- lt_umap_df()
      plot_ly(df, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
              marker = list(size = 4, color = ~col), hoverinfo = "skip") %>%
        plotly::layout(title = list(text = paste("UMAP", i), x = 0), showlegend = FALSE,
                       xaxis = list(title = "", zeroline = FALSE),
                       yaxis = list(title = "", zeroline = FALSE, scaleanchor = "x")) %>%
        plotly::config(displayModeBar = FALSE)
    })

    output[[paste0("lt_gene_", i)]] <- renderPlotly({
      req(i <= input$lt_copies)
      g <- input$lt_genes[1]
      # Always return a plot (empty if no gene) so the timer's render count is exact.
      if (is.null(g) || !(g %in% rownames(seu()))) {
        return(plotly::config(plot_ly(type = "scattergl"), displayModeBar = FALSE))
      }
      emb <- Embeddings(seu(), "umap")
      df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                       expr = as.numeric(GetAssayData(seu(), layer = "data")[g, ]))
      plot_ly(df, x = ~UMAP_1, y = ~UMAP_2, type = "scattergl", mode = "markers",
              marker = list(size = 4, color = ~expr,
                            colorscale = list(c(0, "#D9D9D9"), c(1, "#1E8B0D")), showscale = FALSE),
              hoverinfo = "skip") %>%
        plotly::layout(title = list(text = paste(g, i), x = 0),
                       xaxis = list(title = "", zeroline = FALSE),
                       yaxis = list(title = "", zeroline = FALSE, scaleanchor = "x")) %>%
        plotly::config(displayModeBar = FALSE)
    })

    output[[paste0("lt_heat_", i)]] <- renderPlotly({
      req(i <= input$lt_copies)
      mat <- lt_heat_mat()
      plot_ly(x = colnames(mat), y = rownames(mat), z = mat, type = "heatmap",
              colorscale = list(c(0, "#D9D9D9"), c(1, "#1E8B0D")), showscale = FALSE) %>%
        plotly::layout(title = list(text = paste("Heatmap", i), x = 0),
                       yaxis = list(autorange = "reversed")) %>%
        plotly::config(displayModeBar = FALSE)
    })
  })
}

shinyApp(ui, server)
