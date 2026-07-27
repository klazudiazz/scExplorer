#' Launch the Seurat Object Explorer
#'
#' Starts the interactive Shiny app for exploring a Seurat object. The app UI
#' and server live in `inst/app/app.R`; this function just points
#' [shiny::runApp()] at it.
#'
#' @param ... Additional arguments passed to [shiny::runApp()], e.g.
#'   `launch.browser = TRUE`, `host`, or `port`.
#'
#' @return Called for its side effect (runs the app); does not return a value.
#'
#' @examples
#' \dontrun{
#'   run_app()
#' }
#'
#' @importFrom shiny runApp
#' @export
run_app <- function(...) {
  app_dir <- system.file("app", package = "scexplorer")
  if (identical(app_dir, "")) {
    stop("Could not find the app directory. Try re-installing the 'scexplorer' package.",
         call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
