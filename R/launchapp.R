#' Launch the BRAINscape Shiny App
#'
#' This function launches the Shiny app for visualizing eQTL, DESeq2, and GWAS results.
#' @export
launchApp <- function() {
  app_dir <- system.file("shinyApp", package = "BRAINscape")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing `BRAINscape`.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
