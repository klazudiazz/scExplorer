# scExplorer — Session Handoff for the Next Agent

Handoff summary of work on **scExplorer**, a Shiny app packaged as the R package
`scexplorer`, for exploring Seurat objects. Read this before touching anything —
several non-obvious gotchas below have caused repeated confusion (and one data loss).

---

## 1. What scExplorer is

A general-purpose, **upload-based Shiny app** to explore *any* Seurat object in the
browser (not just scdrake output). Four tabs:
- **Clusters / Metadata** — one plotly UMAP coloured by a chosen metadata column, with
  per-category show/hide checkboxes.
- **Gene Expression** — one UMAP per selected gene, coloured by expression.
- **Cluster Markers** / **Contrasts** — searchable DT tables, shown only if the object
  has `obj@misc$cluster_markers` / `obj@misc$contrasts` (scdrake format).

It's wrapped as an R package with a launcher `run_app()`. Design goal (important):
**must work for ANY Seurat object**, avoid baking in scdrake-specific assumptions.

---

## 2. Where everything lives (READ THIS — the #1 source of confusion)

There are **four separate copies** of the app. Editing one does NOT update the others:

| Copy | Location | Notes |
|---|---|---|
| **Local source** (the one you edit) | `/home/user/scdrake_projects/pbmc1k/scexplorer/inst/app/app.R` | On the laptop `pavilon`. This is the working copy. |
| **Pod copy** (what the user usually RUNS) | pod `~/scdrake_shiny_app/app.R` | On k8s pod `klaudia-scdrake-0`, ns `cmm-bioit-ns`. |
| **Installed package** | `/usr/local/lib/R/site-library/scexplorer/app/app.R` (or user lib) | Used by `scexplorer::run_app()`. A COPY made at install time. |
| **GitHub** | `github.com/klazudiazz/scexplorer` (**private**) | Source of truth for distribution. |

The laptop IS the environment the assistant's Bash tool runs in (`/home/user/...`).

---

## 3. The deploy / test workflow (CRITICAL)

Local edits are invisible on the pod until copied over. The loop is:

1. Agent edits `…/scexplorer/inst/app/app.R` (local).
2. **User** copies it to the pod (from the laptop terminal):
   ```bash
   kubectl cp /home/user/scdrake_projects/pbmc1k/scexplorer/inst/app/app.R \
     cmm-bioit-ns/klaudia-scdrake-0:scdrake_shiny_app/app.R
   ```
3. On the pod: **Session → Restart R** (clears the cached example object!), then in the
   **R Console** (not the Terminal): `shiny::runApp("~/scdrake_shiny_app")`.
4. Viewed via RStudio Server proxy URL (`…dyn.cloud.trusted.e-infra.cz/p/…`).

Quick tell for "am I running the new code?": the title says **"scExplorer"** and the
"Color UMAP by" dropdown shows **friendly labels + many columns**. If not, it's an old copy.

To test on the **laptop** instead (no copy needed): `shiny::runApp("/home/user/scdrake_projects/pbmc1k/scexplorer/inst/app")`.

---

## 4. User preferences (also in the assistant's memory files)

- **Runs everything on the pod (RStudio Server), not locally** — but during this session
  we also set up a local run on the laptop.
- **Runs all git commands themselves** — give them copy-paste commands; do NOT execute
  `git` via Bash (read-only `git status/diff/fetch` for diagnosis is OK).
- **Keeps mixing up prompts**: `>` = R Console (R code), `$` = bash Terminal (shell). They
  repeatedly paste R into bash and vice-versa. Remind gently.
- The app must stay **generic** (any Seurat object).

---

## 5. Environment

- **Pod** `klaudia-scdrake-0` (ns `cmm-bioit-ns`): R 4.4.3, Bioconductor 3.20, Ubuntu Jammy,
  full single-cell stack installed. This is where it reliably runs.
- **Laptop** `pavilon` (Ubuntu Jammy): R was **4.1.2 → upgraded to 4.6** this session (CRAN
  repo). Seurat 5 stack now installs there (SeuratObject 5 needs a newer Matrix than 4.1
  could get — that's why the upgrade was needed). Packages go to the user library
  `~/R/x86_64-pc-linux-gnu-library/<ver>` (site-library isn't writable — always answer
  "yes, personal library" when prompted).
- **GitHub** `klazudiazz/scexplorer` is **private** → `install_github` needs a PAT
  (`auth_token=` or `GITHUB_PAT`). User's GitHub username: `klazudiazz`.
  ⚠️ The user pasted a token into chat twice — advise revoking any exposed token.

---

## 6. The example data (a recurring source of trouble)

- The app's "Use example data" loads a **16 MB trimmed demo** (`seurat_example.rds`,
  ~400 cells) — NOT the original 2.2 GB object.
- The trimming: subset ~400 cells → `DietSeurat(layers=c("counts","data"), dimreducs="umap",
  misc=TRUE)` → trim each `@misc$cluster_markers`/`contrasts` marker table to top ~200 genes
  and **drop the pairwise `lfc_*.*` columns** (that pairwise block was ~99% of the 2.2 GB —
  the expression data itself was only ~7 MB).
- Copies on the pod: `~/inst/app/example_data/seurat_example.rds` (16 MB), the installed
  package copy (16 MB), and the renamed full object `~/scdrake_shiny_app/example_data/seurat_example_FULL.rds` (2.2 GB, renamed so it can't be picked up).
- Local copy of the 16 MB demo also at `~/Downloads/seurat_example.rds`.
- **Gotcha 1 — the cache**: `get_example_seu()` caches the object in `.seurat_example_obj`
  in globalenv for the whole R session. If the 2.2 GB was loaded earlier, it keeps serving
  it until R is restarted. Restart R after changing the example file.
- **Gotcha 2 — DietSeurat dropped SYMBOL/ENSEMBL** feature metadata, so gene labels showed
  "NA (NA)". Fixed in `gene_lookup()` (falls back to the rowname; the demo's rownames ARE
  symbols like `Pcmtd1`).
- **Gotcha 3 — the example path**: `.example_candidate_paths` (top of app.R) tries paths in
  order. First is a hardcoded `/home/rstudio/inst/app/example_data/...` (pod-specific, has a
  `>>> OTHER USERS: change this` comment), then `system.file(...)` (portable), then relative,
  then legacy. `run_app()` resolves via the package copy; `runApp("~/scdrake_shiny_app")`
  resolves via the hardcoded path.

---

## 7. What was changed this session (all in the local app.R)

- Renamed app to **scExplorer** (title + header).
- **Clusters/Metadata UMAP**: rebuilt as a single native `plot_ly` scattergl trace (fixes an
  overplotting bug where `ggplotly`'s per-category traces let one category hide another);
  legend built from invisible proxy traces (frequency-sorted); category toggling recolours
  **in place via `plotlyProxy`** (no full re-render); colours are hex (plotly.js renders R
  names like "grey85" as black); NA-annotation cells become an explicit grey "NA" category.
- **Gene Expression**: converted from `ggplot`+`ggplotly` to **native `plot_ly` scattergl**
  (big speedup). Per-slot rendering (`slot_gene`) so adding a gene only renders the new plot.
  Square grid sized to fit the viewport (4 plots without scrolling). Shared vs per-gene colour
  scale; single shared colorbar to the right in shared mode. Colour scale uses the **true max**
  (`cmin=0, cmax=upper`).
- **UMAP sizing**: `umap_plot_box` sized to match a single gene plot; both use `umap_range()`
  (fixed square window) + `constrain="domain"`/`constraintoward` so all UMAPs are identical.
- **Metadata dropdown**: now offers **all categorical + low-cardinality-numeric** columns
  (`grouping_cols` heuristic), with **friendly labels** (`pretty_label()`: underscores→spaces,
  capitalised, sorted A–Z). Value stays the real column name.
- Case-insensitive table search (`apply_field_search`: lowercases both sides; `grepl`'s
  `ignore.case` is silently ignored with `fixed=TRUE`).
- Loaded-data **status line** replaces the "Get started" picker once data loads
  (`output$data_ready` drives conditionalPanels). Upload option disabled once example chosen.
- **Gene-list loading**: (a) skip the slow `gene_avg_expr` rowMeans unless sorting by
  expression (fast Alphabetical default); (b) **CURRENT approach: a full-app "Loading…"
  overlay** (`#app-loading`, inside `output$app_body`) covers the tabs until the server sends
  the `genes_ready` custom message (gene list populated), then JS hides it. Earlier attempts
  (disabling the selectize, then a gene-bar overlay) did NOT work reliably — the app-level
  overlay is the latest and needs verifying on the pod.

---

## 8. The data-loss incident (so it isn't repeated)

The user **`rm -rf`'d the local `scexplorer/` folder** to do a "fresh from GitHub" test — but
the local commits (all session work) had **never been pushed** (pushes kept failing on
divergence because the user also edits files in GitHub's **web editor**, creating conflicting
commits). The folder + `.git` were gone.

Recovered because: the latest `app.R` survived on the **pod** (`~/scdrake_shiny_app/app.R`) and
the 16 MB demo survived in `~/Downloads`. The rest of the package (DESCRIPTION, NAMESPACE,
run_app.R, man/, README, .gitignore) was recreated by the agent, then force-pushed to GitHub.

**Lessons:** (1) never delete the only copy before confirming a push landed; (2) the user keeps
editing on the GitHub web UI AND locally → constant divergence. Tell them to pick ONE place.
(3) A GitHub web edit had **broken the DESCRIPTION** `Authors@R` (missing quote), which failed
every `install_github` — likely the same error a collaborator hit. It's fixed now.

---

## 9. Known open items / to verify

- **Verify the new app-level loading overlay** actually shows and clears on the pod (redeploy
  needed — user was mid-test). If it clears too early/late, adjust when `genes_ready` fires
  (the CSS `appAutoHide` fallback is 20 s).
- **Point size** in the native gene plots is `size=4` (matches the UMAP); the old ggplot used
  `0.6`. If it looks too big vs before, lower it (`# TWEAK: gene-expression point size`).
- GitHub is **private** → decide public vs private (affects `install_github` for collaborators;
  README documents the token path).
- Consider **server-side DT** (`server=TRUE`) for big marker/contrast tables (next speedup;
  the demo's small tables are already fine).
- `patchwork`/`ggplot2`/`scales` are still `library()`'d though ggplot is no longer used for
  the plots — could trim, low priority.
- Package `DESCRIPTION` author is `Klaudia Zieba <ziebaklaudia4@gmail.com>` (ASCII, `License: MIT`,
  no LICENSE file — user deleted them).

---

## 10. app.R architecture quick-map

- Top: packages, `options(shiny.maxRequestSize)`, `.example_candidate_paths` + `get_example_seu()`,
  `pretty_label()`.
- `ui`: `tags$head` (big CSS block — theme `#39FF14`, spinner keyframes, `#app-loading` overlay,
  file-upload progress; JS — upload MutationObserver, "disable upload after example", and the
  `genes_ready` custom-message handler that hides `#app-loading`). Then title, data-source
  radios (+ conditional fileInput), status line, `withSpinner(uiOutput("app_body"))`.
- `server`:
  - `seu()` — loads uploaded or example object; `output$data_ready` / `output$data_status`.
  - `gene_lookup()`, `gene_avg_expr()`, `found_genes()`, `shared_upper()`, `umap_range()`.
  - gene-list `observeEvent` (sorts, `updateSelectizeInput`, sends `genes_ready`).
  - `output$app_body` — builds the tabs; `grouping_cols` heuristic + `grouping_choices`
    (friendly labels); wraps tabs in the `#app-loading` overlay.
  - `output$category_filter_ui`, `strip_plotly_axis_chrome()` (now largely unused after the
    native rewrite), `output$umap_plot` (native + plotlyProxy toggle observer), the
    Gene-Expression grid `output$dynamic_feature_plots`, `output$shared_colorbar`, the
    `slot_gene` reactives + `feature_plot_individual_i` renderers (native plot_ly),
    marker/contrast tables + helpers.

Search the file for `# TWEAK:` for the intended customisation points.
