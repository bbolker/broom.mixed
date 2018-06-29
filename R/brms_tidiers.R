#' Tidying methods for a brms model
#' 
#' These methods tidy the estimates from 
#' \code{\link[brms:brmsfit-class]{brmsfit-objects}}
#' (fitted model objects from the \pkg{brms} package) into a summary.
#'
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#'
#' @seealso \code{\link[brms]{brms}}, \code{\link[brms]{brmsfit-class}}
#'
#' @name brms_tidiers
#' 
#' @param x Fitted model object from the \pkg{brms} package. See 
#'   \code{\link[brms]{brmsfit-class}}.
#' @examples
#'  ## library(brms)
#'  ## fit <- brm(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars, 
#'  ##           iter = 500, chains = 2)
#'  fit <- readRDS(system.file("extdata","brmsfit_example.rds",package="broom.mixed"))
#'  if (require("broom")) {
#'    tidy(fit)
#'    tidy(fit, parameters = "^sd_", conf.int = FALSE)
#'    tidy(fit, effects = "fixed")
#'    tidy(fit, effects = "ran_vals")
#'    tidy(fit, effects = "ran_pars", robust = TRUE)
#' 
#'   # glance method
#'    glance(fit)
#'   \dontrun{
#'      glance(fit, looic = TRUE, cores = 1)
#'   }
#'
#' }
#'  
NULL

#' @rdname brms_tidiers
#' @param parameters Names of parameters for which a summary should be 
#'   returned, as given by a character vector or regular expressions.
#'   If \code{NA} (the default) summarized parameters are specified
#'   by the \code{effects} argument.
#' @param effects One of \code{"all"}, \code{"fixed"}, 
#'   \code{"ran_vals"}, or \code{"ran_pars"} (can be abbreviated). 
#'   See the Value section for details.
#' @param robust Whether to use median and median absolute deviation of
#' the posterior distribution, rather
#'   than mean and standard deviation, to derive point estimates and uncertainty
#' @param conf.int If \code{TRUE} columns for the lower (\code{conf.low})
#' and upper bounds (\code{conf.high}) of posterior uncertainty intervals are included.
#' @param prob Defines the range of the posterior uncertainty conf.int,
#'  such that \code{100 * prob}\% of the parameter's posterior distribution 
#'  lies within the corresponding interval. 
#'  Only used if \code{conf.int = TRUE}.
#' @param looic Should the LOO Information Criterion (and related info) be
#'   included? See \code{\link[rstanarm]{loo.stanreg}} for details. Note: for
#'   models fit to very large datasets this can be a slow computation.

#' @param ... Extra arguments, not used
#' 
#' @return 
#' When \code{parameters = NA}, the \code{effects} argument is used
#' to determine which parameters to summarize.
#' 
#' Generally, \code{tidy.brmsfit} returns
#' one row for each coefficient, with at least three columns:
#' \item{term}{The name of the model parameter.}
#' \item{estimate}{A point estimate of the coefficient (mean or median).}
#' \item{std.error}{A standard error for the point estimate (sd or mad).}
#' 
#' When \code{effects = "fixed"}, only population-level
#' effects are returned.
#' 
#' When \code{effects = "ran_vals"}, only group-level effects are returned.
#' In this case, two additional columns are added:
#' \item{group}{The name of the grouping factor.}
#' \item{level}{The name of the level of the grouping factor.}
#' 
#' Specifying \code{effects = "ran_pars"} selects the
#' standard deviations and correlations of the group-level parameters.
#' 
#' If \code{conf.int = TRUE}, columns for the \code{lower} and 
#' \code{upper} bounds of the posterior conf.int computed.
#'
#' @note The names \sQuote{fixed}, \sQuote{ran_pars}, and \sQuote{ran_vals}
#' (corresponding to "non-varying", "hierarchical", and "varying" respectively
#' in previous versions of the package), while technically inappropriate in
#' a Bayesian setting where "fixed" and "random" effects are not well-defined,
#' are used for compatibility with other (frequentist) mixed model types.
#' @export
tidy.brmsfit <- function(x, parameters = NA, 
                         effects = c("all", "fixed", "ran_vals", "ran_pars"), 
                         robust = FALSE, conf.int = TRUE, prob = 0.9, ...) {
    use_effects <- anyNA(parameters) 
    if (use_effects) {
        effects <- match.arg(effects)
        if (effects == "all") {
           parameters <- NA 
        } else if (effects == "fixed") {
           parameters <- "^b_"
        } else if (effects == "ran_vals") {
           parameters <- "^r_"
        } else if (effects == "ran_pars") {
           parameters <- c("^sd_", "^cor_")
        }
    }
    samples <- brms::posterior_samples(x, parameters)
    if (is.null(samples)) {
        stop("No parameter name matches the specified pattern.",
             call. = FALSE)
    }
    out <- data.frame(term = names(samples), stringsAsFactors = FALSE)
    if (use_effects) {
        if (effects == "fixed") {
            out$term <- gsub("^b_", "", out$term)
        } else if (effects == "ran_vals") {
            out$term <- gsub("^r_", "", out$term)
            out$group <- gsub("\\[.*", "", out$term)
            out$level <- gsub(".*\\[|,.*", "", out$term)
            out$term <- gsub(".*,|\\]", "", out$term)
        }
        # no renaming if effects %in% c("all", "ran_pars")
    }
    if (robust) {
        out$estimate <- apply(samples, 2, stats::median)
        out$std.error <- apply(samples, 2, stats::mad)
    } else {
        out$estimate <- apply(samples, 2, base::mean)
        out$std.error <- apply(samples, 2, stats::sd)
    }
    if (conf.int) {
        stopifnot(length(prob) == 1L)
        probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
        out[, c("conf.low", "conf.high")] <- 
            t(apply(samples, 2, stats::quantile, probs = probs))
    }
    out
}

#' @rdname brms_tidiers
#' @export
glance.brmsfit <- function(x, looic = FALSE, ...) {
    ## defined in rstanarm_tidiers.R
    glance_stan(x, looic=looic, type="brmsfit", ...)
}

#' @rdname brms_tidiers
#' @export
augment.brmsfit <- function() {}


    
