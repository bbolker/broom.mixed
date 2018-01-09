#' Tidying methods for MCMC (Stan, JAGS, etc.) fits
#'
#' @param x a model fit to be converted to a data frame
#' @param pars (character) specification of which parameters to include
#' @param estimate.method method for computing point estimate ("mean" or "median")
#' @param effects which effects (fixed, random, etc.) to return
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
#' @param ... unused
#' 
#' @name mcmc_tidiers
#' @importFrom broom tidy
#' 
#' @examples
#' 
#' # Using example from "RStan Getting Started"
#' # https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
#' 
#' model_file <- system.file("example_data", "8schools.stan", package = "broom.mixed")
#' schools_dat <- list(J = 8, 
#'                     y = c(28,  8, -3,  7, -1,  1, 18, 12),
#'                     sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
#' \dontrun{
#'   if (require(rstan)) {
#'       set.seed(2015)
#'       rstan_example <- rstan::stan(file = model_file, data = schools_dat, 
#'                          iter = 1000, chains = 2, save_dso = FALSE)
#'      }
#' }
#' rstan_example <- readRDS(system.file("example_data", "rstan_example.rds", package = "broom.mixed"))
#' if (require(broom)) {
#'    tidy(rstan_example)
#'    tidy(rstan_example, conf.int = TRUE, pars = "theta")
#'   
#'    td_mean <- tidy(rstan_example, conf.int = TRUE)
#'    td_median <- tidy(rstan_example, conf.int = TRUE, estimate.method = "median")
#'
#'    if (require(dplyr) && require(ggplot2)) {
#'       tds <- rbind(mutate(td_mean, method = "mean"),
#'                mutate(td_median, method = "median")) %>%
#'          mutate(type=ifelse(grepl("^theta",term),"theta",
#'               ifelse(grepl("^eta",term),"eta",
#'                     "other")))
#'   
#'        ggplot(tds, aes(estimate, term)) +
#'         geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),height=0) +
#'         geom_point(aes(color = method))+
#'         facet_wrap(~type,scale="free",ncol=1)
#'    } ## require(dplyr,ggplot2)
#' } ## if require(broom)
#' 
#' if (requireNamespace("MCMCglmm", quietly = TRUE)) {
#'     data(PlodiaPO,package="MCMCglmm")  
#'     model1 <- MCMCglmm::MCMCglmm(PO~1, random=~FSfamily, data=PlodiaPO, verbose=FALSE)
#'     tidy(model1)
#' }
#' @importFrom stats median sd
#' @importFrom coda HPDinterval as.mcmc
#' @export
tidyMCMC <- function(x,
                     pars,
                     estimate.method = c("mean","median"),
                     conf.int = FALSE,
                     conf.level = 0.95,
                     conf.method = c("quantile","HPDinterval"),
                     drop.pars = c("lp__","deviance"),
                     rhat = FALSE,
                     ess = FALSE,
                     index = FALSE,
                     ...) {
    
    estimate.method <- match.arg(estimate.method)
    conf.method <- match.arg(conf.method)
    

    stan <- inherits(x, "stanfit")
    ss <- if (stan) as.matrix(x, pars = pars) else as.matrix(x)
    ss <- ss[, !colnames(ss) %in% drop.pars, drop = FALSE]  ## drop log-probability info
    if (!missing(pars) && !stan) {
        if (length(badpars <- which(!pars %in% colnames(ss))) > 0) {
            stop("unrecognized parameters: ", pars[badpars])
        }
        ss <- ss[, pars]
    }
    
    m <- switch(estimate.method,
                mean = colMeans(ss),
                median = apply(ss, 2, median))

    # Extract indexes and remove [] if requested
    if (index){
      ret <- data.frame(term0 = sub("\\[\\d+\\]", "", names(m)),
                        index = as.integer(stringr::str_match(names(m), "\\[(\\d+)\\]")[,2]),
                        estimate = m,
                        std.error = apply(ss, 2, stats::sd))
      
    } else {
      ret <- data.frame(estimate = m,
                        std.error = apply(ss, 2, stats::sd))
    }

    if (conf.int) {
        levs <- c((1 - conf.level) / 2, (1 + conf.level) / 2)

        ci <- switch(conf.method,
                     quantile = t(apply(ss, 2, stats::quantile, levs)),
                     HPDinterval(as.mcmc(ss), prob = conf.level))

        colnames(ci) <- c("conf.low", "conf.high")
        ret <- data.frame(ret, ci)
    }
    
    if (rhat || ess) {
        if (!stan) warning("ignoring 'rhat' and 'ess' (only available for stanfit objects)")
        summ <- rstan::summary(x, pars = pars, probs = NULL)$summary[, c("Rhat", "n_eff"), drop = FALSE]
        summ <- summ[!dimnames(summ)[[1L]] %in% drop.pars,, drop = FALSE]
        if (rhat) ret$rhat <- summ[, "Rhat"]
        if (ess) ret$ess <- as.integer(round(summ[, "n_eff"]))
    }
    return(fix_data_frame(ret))
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
             conf.method, drop.pars = "deviance")
}

##' @rdname mcmc_tidiers
##' @export
tidy.stanfit <- tidyMCMC

##' @rdname mcmc_tidiers
##' @export
tidy.mcmc <- tidyMCMC

##' @rdname mcmc_tidiers
##' @importFrom plyr ldply
##' @export
tidy.MCMCglmm <- function(x,effects="fixed",scales,...) {
    if (!missing(scales)) stop("tidy.MCMCglmm doesn't yet implement scales")
    ## FIXME: allow scales= parameter to get varcov on sd/corr scale?
    clist <- c(fixed="Sol",ran_pars="VCV",ran_vals="Liab")
    comp <- clist[effects]
    ## override MCMCglmm internal component names
    ## FIXME:: have to work harder to retrieve group/term information
    ##  about random parameters
    ## individual components are mcmc objects: call tidy on them
    return(plyr::ldply(setNames(x[comp],effects),
                 tidy,...,.id="effect"))
}
