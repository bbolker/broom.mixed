#' Tidying methods for MCMC (Stan, JAGS, etc.) fits
#'
#' @param x a model fit to be converted to a data frame
#' @param pars (character) specification of which parameters to include
#' @param estimate.method method for computing point estimate ("mean" or "median")
#' @param effects which effects (fixed, random, etc.) to return
#' @param robust use mean and standard deviation (if FALSE) or median and mean absolute deviation (if TRUE) to compute point estimates and uncertainty?
#' @param scales scales on which to report results
#' @param conf.int (logical) include confidence interval?
#' @param conf.level probability level for CI
#' @param conf.method method for computing confidence intervals
#' ("quantile" or "HPDinterval")
#' @param drop.pars Parameters not to include in the output (such
#' as log-probability information)
#' @param rhat,ess (logical) include Rhat and/or effective sample size estimates?
#' @param index Add index column, remove index from term. For example,
#' \code{term a[13]} becomes \code{term a} and \code{index 13}.
#' @param ... mostly unused; for \code{tidy.MCMCglmm}, these represent options
#' passed through to \code{tidy.mcmc} (e.g. \code{robust}, \code{conf.int}, \code{conf.method}, ...)
#'
#' @name mcmc_tidiers
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows bind_cols
#'
#' @examples
#'
#' # Using example from "RStan Getting Started"
#' # https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
#'
#' model_file <- system.file("extdata", "8schools.stan", package = "broom.mixed")
#' schools_dat <- list(J = 8,
#'                     y = c(28,  8, -3,  7, -1,  1, 18, 12),
#'                     sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
#' ## original model
#' \dontrun{
#'     set.seed(2015)
#'     rstan_example <- rstan::stan(file = model_file, data = schools_dat,
#'                          iter = 1000, chains = 2, save_dso = FALSE)
#' }
#' if (require(rstan)) {
#'    ## load stored object
#'    rstan_example <- readRDS(system.file("extdata", "rstan_example.rds", package = "broom.mixed"))
#'    tidy(rstan_example)
#'    tidy(rstan_example, conf.int = TRUE, pars = "theta")
#'    td_mean <- tidy(rstan_example, conf.int = TRUE)
#'    td_median <- tidy(rstan_example, conf.int = TRUE, estimate.method = "median")
#' 
#'   if (require(dplyr) && require(ggplot2)) {
#'     tds <- rbind(mutate(td_mean, method = "mean"),
#'              mutate(td_median, method = "median")) %>%
#'        mutate(type=ifelse(grepl("^theta",term),"theta",
#'             ifelse(grepl("^eta",term),"eta",
#'                   "other")))
#'
#'      ggplot(tds, aes(estimate, term)) +
#'       geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),height=0) +
#'       geom_point(aes(color = method))+
#'       facet_wrap(~type,scale="free",ncol=1)
#'  } ## require(dplyr,ggplot2)
#' } ## require(rstan)
#'
#' @importFrom stats median sd
#' @importFrom coda HPDinterval as.mcmc
#' @export
tidyMCMC <- function(x,
                     pars,
                     robust = FALSE,
                     conf.int = FALSE,
                     conf.level = 0.95,
                     conf.method = c("quantile", "HPDinterval"),
                     drop.pars = c("lp__", "deviance"),
                     rhat = FALSE,
                     ess = FALSE,
                     index = FALSE,
                     ...) {
  conf.method <- match.arg(conf.method)

  stan <- inherits(x, "stanfit")
  ss <- if (stan) as.matrix(x, pars = pars) else as.matrix(x)
  ss <- ss[, !colnames(ss) %in% drop.pars, drop = FALSE] ## drop log-probability info
  if (!missing(pars) && !stan) {
    if (length(badpars <- which(!pars %in% colnames(ss))) > 0) {
      stop("unrecognized parameters: ", pars[badpars])
    }
    ss <- ss[, pars]
  }

  m <- if (robust) colMeans(ss) else apply(ss, 2, median)

  stdfun <- if (robust) stats::mad else stats::sd
  ret <- dplyr::data_frame(
    term = names(m),
    estimate = m,
    std.error = apply(ss, 2, stdfun)
  )

  ## Extract indices and remove [] if requested
  if (index) {
    ret$index <- as.integer(stringr::str_match(names(m), "\\[(\\d+)\\]")[, 2])
    ret$term <- sub("\\[\\d+\\]", "", names(m))
  }

  if (conf.int) {
    levs <- c((1 - conf.level) / 2, (1 + conf.level) / 2)

    ci <- switch(conf.method,
      quantile = t(apply(ss, 2, stats::quantile, levs)),
      HPDinterval(as.mcmc(ss), prob = conf.level)
    ) %>%
      as.data.frame()

    names(ci) <- c("conf.low", "conf.high")
    ret <- bind_cols(ret, ci)
  }

  if (!stan) {
      if (rhat) warning("ignoring 'rhat' (only available for stanfit objects)")
      if (ess) {
          ess_vals <- coda::effectiveSize(x)
          ## FIXME: deal with this
          if (!missing(pars)) warning("pars ignored in ESS computation - might break something?")
          ret$ess <- as.integer(round(ess_vals))
      }
  } else {
      summ <- rstan::summary(x, pars = pars, probs = NULL)$summary[, c("Rhat", "n_eff"), drop = FALSE]
      summ <- summ[!dimnames(summ)[[1L]] %in% drop.pars, , drop = FALSE]
    if (rhat) ret$rhat <- summ[, "Rhat"]
    if (ess) ret$ess <- as.integer(round(summ[, "n_eff"]))
  }
  return(fix_data_frame(ret) %>% reorder_cols())
}


##' @rdname mcmc_tidiers
##' @importFrom coda as.mcmc
##' @export
tidy.rjags <- function(x,
                       estimate.method = "mean",
                       conf.int = FALSE,
                       conf.level = 0.95,
                       conf.method = "quantile",
                       ...) {
  tidyMCMC(as.mcmc(x$BUGS),
    estimate.method, conf.int, conf.level,
    conf.method,
    drop.pars = "deviance"
  )
}

##' @rdname mcmc_tidiers
##' @export
tidy.stanfit <- tidyMCMC

##' @rdname mcmc_tidiers
##' @export
tidy.mcmc <- tidyMCMC
