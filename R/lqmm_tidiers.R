## from David Luke Thiessen, https://stackoverflow.com/questions/72514230/is-it-possible-to-use-lqmm-with-a-mira-object
tidy.lqmm <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  as_tidy_tibble(data.frame(
    estimate = coef(x),
    std.error = sqrt(
      diag(summary(x, covariance = TRUE, 
                   R = 50)$Cov[names(coef(x)),
                               names(coef(x))]))))
}
glance.lqmm <- function(x, ...) {
  as_glance_tibble(
    logLik = as.numeric(stats::logLik(x)),
    df.residual = summary(x)$rdf,
    nobs = stats::nobs(x),
    na_types = "rii")
}
