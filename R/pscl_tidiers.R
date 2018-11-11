#' Tidying methods for zero-altered models from pscl
#'
#' These methods tidy the coefficients of zero-altered (\code{zeroinfl} and \code{hurdle} models from the \pkg{pscl} package
#'
#' @param x An object of class \code{zeroinfl} or \code{hurdle}
#'
#' @return A tibble.
#'
#' @name pscl_tidiers
#' @aliases zeroinfl_tidiers
#' @aliases hurdle_tidiers
#'
#' @examples
#'
#' if (require("pscl")) {
#'    data(bioChemists)
#'    zi1 <- zeroinfl(art ~ ., data=bioChemists)
#'    h1 <- hurdle(art ~ ., data=bioChemists)
#' }
