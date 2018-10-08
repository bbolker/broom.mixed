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
#'  ## original model
#'  \dontrun{
#'     brms_crossedRE <- brm(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
#'            iter = 500, chains = 2)
#'  }
#'  if (require("brms")) {
#'    ## load stored object
#'    load(system.file("extdata", "brms_example.rda", package="broom.mixed"))
#'
#'    fit <- brms_crossedRE
#'    tidy(fit)
#'    tidy(fit, parameters = "^sd_", conf.int = FALSE)
#'    tidy(fit, effects = "fixed")
#'    tidy(fit, effects = "ran_vals")
#'    tidy(fit, effects = "ran_pars", robust = TRUE)
#'    # glance method
#'    glance(fit)
#'    ## this example will give a warning that it should be run with
#'    ## reloo=TRUE; however, doing this will fail
#'    ## because the \code{fit} object has been stripped down to save space
#'    suppressWarnings(glance(fit, looic = TRUE, cores = 1))
#'    head(augment(fit))
#' }
#'
NULL
## examples for all methods (tidy/glance/augment) included in the same
##  block so we can surround them with a single "if (require(brms))" block

#' @rdname brms_tidiers
#' @param parameters Names of parameters for which a summary should be
#'   returned, as given by a character vector or regular expressions.
#'   If \code{NA} (the default) summarized parameters are specified
#'   by the \code{effects} argument.
#' @param effects A character vector including one or more of \code{"fixed"},
#'   \code{"ran_vals"}, or \code{"ran_pars"}.
#'   See the Value section for details.
#' @param robust Whether to use median and median absolute deviation of
#' the posterior distribution, rather
#'   than mean and standard deviation, to derive point estimates and uncertainty
#' @param conf.int If \code{TRUE} columns for the lower (\code{conf.low})
#' and upper bounds (\code{conf.high}) of posterior uncertainty intervals are included.
#' @param conf.level Defines the range of the posterior uncertainty conf.int,
#'  such that \code{100 * conf.level}\% of the parameter's posterior distributio
#'  lies within the corresponding interval.
#'  Only used if \code{conf.int = TRUE}.
#' @param conf.method method for computing confidence intervals
#' ("quantile" or "HPDinterval")
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
                         effects = c("fixed", "ran_pars"),
                         robust = FALSE, conf.int = TRUE,
                         conf.level = 0.95,
                         conf.method = c("quantile", "HPDinterval"),
                         ...) {
  use_effects <- anyNA(parameters)
  conf.method <- match.arg(conf.method)
  is.multiresp <- length(x$formula$forms)>1
  mkRE <- function(x) {
    sprintf("^(%s)", paste(unlist(x), collapse = "|"))
  }
  if (use_effects) {
    prefs_LA <- list(
      fixed = "b_", ran_vals = "r_",
      ## don't want to remove these pieces, so use lookahead
      ran_pars = sprintf("(?=(%s))", c("sd_", "cor_", "sigma"))
    )
    prefs <- list(
      fixed = "b_", ran_vals = "r_",
      ## no lookahead (doesn't work with grep[l])
      ran_pars = c("sd_", "cor_", "sigma")
    )
    pref_RE <- mkRE(prefs[effects])
    parameters <- pref_RE
  }
  samples <- brms::posterior_samples(x, parameters)
  if (is.null(samples)) {
    stop("No parameter name matches the specified pattern.",
      call. = FALSE
    )
  }
  terms <- names(samples)
  if (use_effects) {
    res_list <- list()
    fixed.only <- identical(effects, "fixed")
    if ("fixed" %in% effects) {
      ## empty tibble: NA columns will be filled in as appropriate
      nfixed <- sum(grepl(prefs[["fixed"]], terms))
      res_list$fixed <- as_tibble(matrix(nrow = nfixed, ncol = 0))
    }
    grpfun <- function(x) {
        if (grepl("sigma",x[[1]])) "Residual" else x[[2]]
    }
    if ("ran_pars" %in% effects) {
      rterms <- grep(mkRE(prefs$ran_pars), terms, value = TRUE)
      ss <- strsplit(rterms, "__")
      pp <- "^(cor|sd)(?=(_))"
      nodash <- function(x) gsub("^_", "", x)
      ss2 <- lapply(
        ss,
        function(x) {
          if (!is.na(pref <- stringr::str_extract(x[1], pp))) {
            return(c(pref, nodash(stringr::str_remove(x[1], pp)), x[-1]))
          }
          return(x)
        }
      )
      sep <- getOption("broom.mixed.sep1")
      termfun <- function(x) {
        if (grepl("^sigma",x[[1]])) {
            paste("sd", "Observation", sep = sep)
        } else {
          paste(x[[1]],
            paste(x[3:length(x)], collapse = "."),
            sep = sep
          )
        }
      }
      res_list$ran_pars <-
        data_frame(
          group = sapply(ss2, grpfun),
          term = sapply(ss2, termfun)
        )
    }
    if ("ran_vals" %in% effects) {
      rterms <- grep(mkRE(prefs$ran_vals), terms, value = TRUE)
      ss <- strsplit(rterms, "(_+|\\[|\\]|,)")
      termfun <- function(x) x[[length(x)]]
      grpfun <- function(x) x[[2]]
      levelfun <- function(x) x[[length(x)-1]]
      res_list$ran_vals <-
        data_frame(
          group = sapply(ss, grpfun),
          term = sapply(ss, termfun),
          level = sapply(ss, levelfun)
        )
    }
    out <- dplyr::bind_rows(res_list, .id = "effect")
    v <- if (fixed.only) seq(nrow(out)) else is.na(out$term)
    newterms <- stringr::str_remove(terms[v], mkRE(prefs[c("fixed")]))
    if (fixed.only) {
      out$term <- newterms
    } else {
      out$term[v] <- newterms
    }

    ## prefixes already removed for ran_vals; don't remove for ran_pars
  } else {
    ## if !use_effects
    out <- data_frame(term = names(samples))
  }
  pointfun <- if (robust) stats::median else base::mean
  stdfun <- if (robust) stats::mad else stats::sd
  out$estimate <- apply(samples, 2, pointfun)
  out$std.error <- apply(samples, 2, stdfun)
  if (conf.int) {
    stopifnot(length(conf.level) == 1L)
    probs <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
    if (conf.method == "HPDinterval") {
      stop("HPDinterval not implemented yet")
    } else {
      out[, c("conf.low", "conf.high")] <-
        t(apply(samples, 2, stats::quantile, probs = probs))
    }
  }
  ## figure out component
  out$component <- dplyr::case_when(grepl("^zi",out$term) ~ "zi",
                                    ## ??? is this possible in brms models
                                    grepl("^disp",out$term) ~ "disp",
                                    TRUE ~ "cond")

  out <- reorder_cols(out)
  
  return(out)
}

#' @rdname brms_tidiers

#' @export
glance.brmsfit <- function(x, looic = FALSE, ...) {
  ## defined in rstanarm_tidiers.R
  glance_stan(x, looic = looic, type = "brmsfit", ...)
}

#' @rdname brms_tidiers
#' @param data data frame
#' @param newdata new data frame
#' @param se.fit return standard errors of fit?
#' @export
augment.brmsfit <- function(x, data = stats::model.frame(x), newdata = NULL,
                            se.fit = TRUE, ...) {
  ## can't use augment_columns because residuals.brmsfit returns
  ## a 4-column matrix (because summary=TRUE by default, no way
  ## to suppress this within augment_columns)
  ## ... add resids.arg to augment_columns?
  args <- list(x, se.fit = se.fit)
  if (!missing(newdata)) args$newdata <- newdata
  ## FIXME: influence measures??
  ## allow optional arguments to augment, e.g. pred.type,
  ## residual.type, re.form ...
  pred <- do.call(stats::predict, args)
  ret <- dplyr::data_frame(.fitted = pred[, "Estimate"])
  if (se.fit) ret[[".se.fit"]] <- pred[, "Est.Error"]
  if (is.null(newdata)) {
    ret[[".resid"]] <- stats::residuals(x)[, "Estimate"]
  }
  return(ret)
}
