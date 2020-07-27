##' Tidying methods for TMB models
##'
##' @param x An object of class \code{TMB} (you may need to use
##' \code{class(obj) <- "TMB"} on your results from TMB)
##' @inheritParams glmmTMB_tidiers
##' @param effect which effects should be returned?
##' @param conf.method method for computing confidence intervals
##' @param ... additional arguments passed to confint function (tmbroot, tmbprofile)
##' @importFrom TMB sdreport
##' @importFrom TMB tmbroot
##' @importFrom splines backSpline interpSpline
## FIXME: retrieving stored objects doesn't work well ...
##' @examples
##' if (require("TMB")) {
##'     runExample("simple",thisR=TRUE)
##'     class(obj) <- "TMB"
##'     tidy(obj,conf.int=TRUE,conf.method="wald")
##'     tidy(obj,conf.int=TRUE,conf.method="uniroot")
##'     tidy(obj,conf.int=TRUE,conf.method="profile")
##' }
##' @export
tidy.TMB <- function(x, effect = c("fixed", "random"),
                     conf.int = FALSE,
                     conf.level = 0.95,
                     conf.method = c("wald", "uniroot", "profile"), ...) {
  ## R CMD check/global variables
  Estimate <- estimate <- std.error <- NULL
  sdr <- sdreport(x)
  retlist <- list()
  if ("fixed" %in% effect) {
    ss <- summary(sdr, select = "fixed") %>%
      as.data.frame() %>%
      ## FIXME: disambiguate repeated variable names better
      ## i.e. beta.1, beta.2 instead of beta, beta.1 ...
      tibble::rownames_to_column("term") %>%
      rename(estimate = Estimate, std.error = "Std. Error")
    if (conf.int) {
        if (tolower(conf.method == "wald")) {
            qval <- qnorm((1 + conf.level) / 2)
            ss <- mutate(ss,
                         conf.low = estimate - qval * std.error,
                         conf.high = estimate + qval * std.error
                         )
        } else if (conf.method == "uniroot") {
            ## FIXME: allow parm specs
            ## FIXME: avoid calling sdreport again inside tmbroot?
            ## do.call(rbind,...) because bind_rows needs named list
            tt <- do.call(
                rbind,
                lapply(seq(nrow(ss)),
                       tmbroot,
                       obj = x,
                       ...
                       )
            )
            ss$conf.low <- tt[, "lwr"]
            ss$conf.high <- tt[, "upr"]
        } else if (conf.method == "profile") {
            ## FIXME: allow parm specs
            ## FIXME: repeated var names?
            all_vars <- names(x$env$last.par.best)[-x$env$random]
            prof0 <- purrr:::map_dfr(seq_along(all_vars),
                                     ~ setNames(tmbprofile(x,name=.,trace=FALSE),c("focal","value")),
                                     .id="param")
            prof1 <- (prof0
                %>% group_by(param)
                %>% mutate(zeta=sqrt(2*(value-min(value))),
                           branch=ifelse(cumsum(zeta==0)<1, "lwr", "upr"))
                %>% ungroup()
            )
            bad_prof_flag <- FALSE
            critval <- qnorm((1+conf.level)/2)
            interp_fun <- function(dd) {
                bakspl <-tryCatch(backSpline(
                    forspl <- interpSpline(dd$focal, dd$zeta, na.action=na.omit)),
                    error=function(e)e)
                if (inherits(bakspl, "error")) {
                    bad_prof_flag <<- TRUE  ## set value globally
                    res <- approx(dd$zeta, dd$focal, xout=critval)$y
                } else {
                    res <- predict(bakspl, critval)$y
                }
                return(res)
            }
            tt <- prof1 %>% group_by(param, branch) %>% unique() %>% summarise(v=interp_fun(.data))
            ss$conf.low <- filter(tt, branch=="lwr") %>% pull(v)
            ss$conf.high <- filter(tt, branch=="upr") %>% pull(v)
        } else {
            stop(sprintf("conf.method=%s not implemented", conf.method))
        }
    } ## if conf.int
  }
  retlist$fixed <- ss
  ret <- dplyr::bind_rows(retlist, .id = "type")
  ## FIXME: add statistic (mean/stderr)
  ##   and p-value (2*pnorm(abs(statistic),lower.tail=FALSE))
  ##   if requested
  ##
  return(ret)
}
