## from David Luke Thiessen, 
#' Tidying methods for lqmm models (EXPERIMENTAL)
#'
#' These methods, suggested by
#' David Luke Thiessen on \href{https://stackoverflow.com/questions/72514230/is-it-possible-to-use-lqmm-with-a-mira-object}{Stack Exchange}, provide support for linear quantile mixed models. They have not been carefully tested - please
#' check output carefully and report problems!
#'
#' @inheritParams lme4_tidiers
#'
#' @export
tidy.lqmm <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  as_tidy_tibble(data.frame(
    estimate = coef(x),
    std.error = sqrt(
      diag(summary(x, covariance = TRUE, 
                   R = 50)$Cov[names(coef(x)),
                               names(coef(x))]))))
}

#' @rdname tidy.lqmm
#' @export
glance.lqmm <- function(x, ...) {
  as_glance_tibble(
    logLik = as.numeric(stats::logLik(x)),
    df.residual = summary(x)$rdf,
    nobs = stats::nobs(x),
    na_types = "rii")
}
