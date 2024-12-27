##' Tidying methods for TMB models
##'
##' @param x An object of class \code{TMB} (you may need to use
##' \code{class(obj) <- "TMB"} on your results from TMB)
##' @inheritParams glmmTMB_tidiers
##' @param effects which effects should be returned?
##' @param conf.method method for computing confidence intervals
##' @param ... additional arguments passed to confint function (tmbroot, tmbprofile)
##' @importFrom stats approx predict
##' @importFrom splines backSpline interpSpline
##' @importFrom furrr future_map_dfr
##' @importFrom forcats fct_inorder
## FIXME: retrieving stored objects doesn't work well ...
## but we also don't want to try compiling as part of example!
##' @examples
##' if (require("TMB")) {
##' 
##'     \dontrun{
##'        runExample("simple",thisR=TRUE)
##'        class(obj) <- "TMB"
##'        tidy(obj,conf.int=TRUE,conf.method="wald")
##'     }
##'     \dontrun{tidy(obj,conf.int=TRUE,conf.method="uniroot")}
##'     \dontrun{tidy(obj,conf.int=TRUE,conf.method="profile")}
##' }
##' @export
tidy.TMB <- function(x, effects = c("fixed", "random"),
                     conf.int = FALSE,
                     conf.level = 0.95,
                     conf.method = c("wald", "uniroot", "profile"), ...) {

  assert_dependency("TMB")

  ## FIXME: implement CIs/profiling/etc. for random effects (and error out on a stub in the meantime)
  ## R CMD check/global variables
  branch <- v <- param <- value <- zeta <- Estimate <- estimate <- std.error <- NULL

  conf.method <- match.arg(conf.method)

  sdr <- TMB::sdreport(x)
  retlist <- list()
  if ("fixed" %in% effects) {
    ss <- summary(sdr, select = "fixed") %>%
      as.data.frame() %>%
      ## FIXME: disambiguate repeated variable names better
      ## i.e. beta.1, beta.2 instead of beta, beta.1 ...
      tibble::rownames_to_column("term") %>%
        rename(estimate = Estimate, std.error = "Std. Error")
    all_vars <- names(x$env$last.par.best)
    if (!is.null(rnd <- x$env$random)) {
        all_vars <- all_vars[-rnd]
    }
    if (conf.int) {
        if (tolower(conf.method) == "wald") {
            qval <- qnorm((1 + conf.level) / 2)
            ss <- mutate(ss,
                         conf.low = estimate - qval * std.error,
                         conf.high = estimate + qval * std.error
                         )
        } else if (conf.method == "uniroot") {
            ## FIXME: allow parm specs
            ## FIXME: avoid calling sdreport again inside tmbroot?
            tt <- furrr::future_map_dfr(seq_along(all_vars),
                                  ~ TMB::tmbroot(x, ...))
            ss <- mutate(ss,
                         conf.low = tt[["lwr"]],
                         conf.high = tt[["upr"]])
        } else if (conf.method == "profile") {
            ## FIXME: allow parm specs
            ## FIXME: repeated var names?
            ## FIXME: tracing/quietly/etc?
            prof0 <- furrr::future_map_dfr(seq_along(all_vars),
                                           ~ (TMB::tmbprofile(x,name=.,trace=FALSE)
                                               %>% setNames(c("focal","value"))),
                                     .id="param")
            prof1 <- (prof0
                %>% mutate(across(param, forcats::fct_inorder)) ## keep order consistent below!
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
                    ## may fail (if not at least two non-NA values ...)
                    res <- try(approx(dd$zeta, dd$focal, xout=critval)$y,
                               silent = TRUE)
                    if (inherits(res, "try-error")) {
                        res <- NA_real_
                    }
                } else {
                    res <- predict(bakspl, critval)$y
                }
                return(res)
            }
            ttd <- (prof1
                %>% group_by(param, branch)
                %>% distinct()
                %>% summarise(v=interp_fun(.data), .groups="drop")
            )
            tt <- (ttd %>% full_join(with(ttd, tidyr::crossing(param, branch)),
                                     by = c("param", "branch"))
                %>% arrange(branch, param)
            )
            ss$conf.low <- filter(tt, branch=="lwr") %>% pull(v)
            ss$conf.high <- filter(tt, branch=="upr") %>% pull(v)
        } else {
            stop(sprintf("conf.method=%s not implemented", conf.method))
        }
    } ## if conf.int
  } ## if 'fixed'
  retlist$fixed <- ss
  ## FIXME: find a way to join rather than bind rows? (ordering fragility)
  ret <- dplyr::bind_rows(retlist, .id = "type")
  ## FIXME: add statistic (mean/stderr)
  ##   and p-value (2*pnorm(abs(statistic),lower.tail=FALSE))
  ##   if requested
  ##
  return(ret)
}

#' @export
glance.TMB <- function(x, nobs = NA, ...) {
    assert_dependency("TMB")
    pars <- x$env$last.par.best
    random <- x$env$random
    if (!is.null(random)) pars <- pars[-random]
    npar <- length(pars)
    dev <- x$fn(pars)
    loglik <- -dev/2
    AIC <-  dev + 2*npar
    BIC <- dev + npar*log(nobs)
    tibble(df = npar, loglik, AIC, BIC)
}
