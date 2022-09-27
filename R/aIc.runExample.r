#' \code{aIc.runExample} loads the associated shiny app
#' This will load the selex example dataset with the default group sizes,
#' the user can upload their own local dataset and adjust groups accordingly.
#'
#' @return No return value, but instead opens a shiny connection to your 
#' default web browser with the selex dataset as an example.
#'
#' @author Greg Gloor
#'
#' @examplesIf interactive()
#' library(aIc)
#' aIc.runExample()
#' @importFrom shiny runApp
#' @export
aIc.runExample <- function() {
  appDir <- system.file("amIcomp", package = "aIc")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `aIc`.", call. = FALSE)
  }
  shiny::runApp(appDir)
}

