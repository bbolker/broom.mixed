#' Tidying methods for gamlss objects
#'
#' Tidying methods for "gamlss" objects from the gamlss package.
#'
#' @param x A "gamlss" object
#' @param quick Whether to perform a fast version, and return only the coefficients
#' @param conf.int whether to return confidence intervals
#' @param ... arguments passed to \code{confint.gamlss}
#' @inheritParams lme4_tidiers
#'
#' @importFrom stats na.omit
#' @name gamlss_tidiers
#'
#' @template boilerplate
#'
#' @return A tibble with one row for each coefficient, containing columns:
#'   \item{parameter}{type of coefficient being estimated: \code{mu}, \code{sigma}, \code{nu}, or \code{tau}}
#'   \item{term}{term in the model being estimated and tested}
#'   \item{estimate}{estimated coefficient}
#'   \item{std.error}{standard error}
#'   \item{statistic}{t-statistic}
#'   \item{p.value}{two-sided p-value}
#' 
#' @examples
#' if (requireNamespace("gamlss", quietly = TRUE) &&
#'     requireNamespace("gamlss.data", quietly = TRUE)) {
#'     data(abdom, package="gamlss.data")
#'     \dontrun{
#'          mod <- gamlss(y~pb(x), sigma.fo=~pb(x), family=BCT,
#'                        data=abdom, method=mixed(1,20))
#'     }
#'     ## load stored object
#'     mod <- readRDS(system.file("extdata", "gamlss_example.rds",
#'                    package="broom.mixed"))
#'     tidy(mod)
#' }
#'
#' @export
tidy.gamlss <- function(x, quick = FALSE, conf.int = FALSE, conf.level = 0.95, ...) {
  if (quick) {
    co <- stats::coef(x)
    return(dplyr::tibble(term = names(co), estimate = unname(co)))
  }

  # need gamlss for summary to work
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("gamlss package not installed, cannot tidy gamlss")
  }

  # use capture.output to prevent summary from being printed to screen
  utils::capture.output(s <- summary(x, type = "qr"))

  # tidy the coefficients much as would be done for a linear model
  nn <- c("estimate", "std.error", "statistic", "p.value")
  ret <- fix_data_frame(s, nn, newcol="term")

  if (conf.int) {
    cilist <- lapply(x$parameters,
        function(w) confint(x, what=w, level=conf.level, ...))
    cimat <- do.call(rbind, cilist)
    ret <- bind_cols(ret, conf.low=cimat[,1], conf.high=cimat[,2])
  }
  ## add parameter types
  coefs <- purrr::map(x[paste0(x$parameters,".coefficients")],na.omit)
  parameters <- rep(x$parameters,
                    vapply(coefs,length,1L))
  bind_cols(parameter = parameters, ret)
}

#' @export
glance.gamlss <- function (x, ...) {
    ret <- list()
    ll <- function(x) c(logLik(x))
    for (cc in c("ll","AIC","BIC","deviance","df.residual","nobs")) {
        ret[[cc]] <- do.call(cc,list(x))
    }
    ret <- (ret
        %>% as_tibble()
        %>% mutate(df = x$df.fit)
        %>% dplyr::relocate(df, .before="df.residual")
    )
    return(ret)
}
