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
#'  if (require("brms") && .Platform$OS.type!="windows") {
#'    ## too slow on Windows, skip (>5 seconds on r-devel-windows)
#'    ## load stored object
#'    load(system.file("extdata", "brms_example.rda", package="broom.mixed"))
#'
#'    fit <- brms_crossedRE
#'    tidy(fit)
#'    tidy(fit, parameters = "^sd_", conf.int = FALSE)
#'    tidy(fit, effects = "fixed", conf.method="HPDinterval")
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
#' @param fix.intercept rename "Intercept" parameter to "(Intercept)", to match
#' behaviour of other model types?
#' @param looic Should the LOO Information Criterion (and related info) be
#'   included? See \code{\link[rstan]{loo.stanfit}} for details. (This
#'   can be slow for models fit to large datasets.)
#' @param ... Extra arguments, not used
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
#' @note At present, the components of parameter estimates are separated by parsing the column names of \code{posterior_samples} (e.g. \code{r_patient[1,Intercept]} for the random effect on the intercept for patient 1, or \code{b_Trt1} for the fixed effect \code{Trt1}. We try to detect underscores in parameter names and warn, but detection may be imperfect.
#' @export
tidy.brmsfit <- function(x, parameters = NA,
                         effects = c("fixed", "ran_pars"),
                         robust = FALSE, conf.int = TRUE,
                         conf.level = 0.95,
                         conf.method = c("quantile", "HPDinterval"),
                         fix.intercept = TRUE,
                         ...) {

  if (!requireNamespace("brms", quietly=TRUE)) {
      stop("can't tidy brms objects without brms installed")
  }
  xr <- brms::restructure(x)
  has_ranef <- nrow(xr$ranef)>0  
  if (any(grepl("_", rownames(fixef(x)))) ||
        (has_ranef && any(grepl("_", names(ranef(x)))))) {
      warning("some parameter names contain underscores: term naming may be unreliable!")
  }
  use_effects <- anyNA(parameters)
  conf.method <- match.arg(conf.method)
  is.multiresp <- length(x$formula$forms)>1
  ## make regular expression from a list of prefixes
  mkRE <- function(x,LB=FALSE) {
      pref <- "(^|_)"
      if (LB) pref <- sprintf("(?<=%s)",pref)
      sprintf("%s(%s)", pref, paste(unlist(x), collapse = "|"))
  }
  ## NOT USED:  could use this (or something like) to
  ##  obviate need for gsub("_","",str_extract(...)) pattern ...  
  prefs_LB <- list(
      fixed = "b_", ran_vals = "r_",
      ## don't want to remove these pieces, so use look*behind*
      ran_pars =   sprintf("(?<=(%s))", c("sd_", "cor_", "sigma")),
      components = sprintf("(?<=%s)", c("zi_","disp_"))
    )
    prefs <- list(
      fixed = "b_", ran_vals = "r_",
      ## no lookahead (doesn't work with grep[l])
      ran_pars = c("sd_", "cor_", "sigma"),
      components = c("zi_", "disp_")
    )
    pref_RE <- mkRE(prefs[effects])
  if (use_effects) {
    ## prefixes distinguishing fixed, random effects

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
      if (is.multiresp) {
        if ("ran_pars" %in% effects && any(grepl("^sd",terms))) {
           warning("ran_pars response/group tidying for multi-response models is currently incorrect")
        }
        ## FIXME: unfinished attempt to fix GH #39
        ## extract response component from terms
        ## resp0 <- strsplit(terms, "_+")
        ## resp1 <- sapply(resp0,
        ##          function(x) if (length(x)==2) x[2] else x[length(x)-1])
        ## ## put the pieces back together
        ## t0 <- lapply(resp0,
        ##          function(x) if (length(x)==2) x[1] else x[-(length(x)-1)])
        ## t1 <- lapply(t0,
        ##          function(x)     
        ##              case_when(
        ##                  x[[1]]=="b"  ~ sprintf("b%s",x[[2]]),
        ##                  x[[2]]=="sd" ~ sprintf("sd_%s__%s",x[[2]],x[[3]]),
        ##                  x[[3]]=="cor" ~ sprintf("cor_%s_%s_%s_%s",
        ##                                          x[[2]],x[[3]],x[[4]],x[[5]])
        ##              ))
        ## resp0 <- stringr::str_extract_all(terms, "_[^_]+")
        ## resp1 <- lapply(resp0, gsub, pattern= "^_", replacement="")
        response <- gsub("^_","",stringr::str_extract(terms,"_[^_]+"))
        terms <- sub("_[^_]+","",terms)
    }
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
      ##  split the first term (cor/sd) into tag + group
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
            ## re-attach remaining terms
            paste(x[[1]],
                  paste(x[3:length(x)], collapse = "."),
                  sep = sep
          )
        }
      }
      res_list$ran_pars <-
        dplyr::tibble(
          group = sapply(ss2, grpfun),
          term = sapply(ss2, termfun)
        )
    }
    if ("ran_vals" %in% effects) {
      rterms <- grep(mkRE(prefs$ran_vals), terms, value = TRUE)
      
      vals <- stringr::str_match_all(rterms, "_(.+?)\\[(.+?),(.+?)\\]")

      res_list$ran_vals <-
        dplyr::tibble(
          group = plyr::laply(vals, function (v) { v[[2]] }),
          term = plyr::laply(vals, function (v) { v[[4]] }),
          level = plyr::laply(vals, function (v) { v[[3]] })
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
    if (is.multiresp) {
        out$response <- response
    }
    ## prefixes already removed for ran_vals; don't remove for ran_pars
  } else {
    ## if !use_effects
    out <- dplyr::tibble(term = names(samples))
  }
  pointfun <- if (robust) stats::median else base::mean
  stdfun <- if (robust) stats::mad else stats::sd
  out$estimate <- apply(samples, 2, pointfun)
  out$std.error <- apply(samples, 2, stdfun)
  if (conf.int) {
    stopifnot(length(conf.level) == 1L)
    probs <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
    if (conf.method == "HPDinterval") {
        cc <- coda::HPDinterval(coda::as.mcmc(samples), prob=conf.level)
    } else {
        cc <- t(apply(samples, 2, stats::quantile, probs = probs))
    }
    out$conf.low <- cc[,1]
    out$conf.high <- cc[,2]
  }
  ## figure out component
  out$component <- dplyr::case_when(grepl("(^|_)zi",out$term) ~ "zi",
                                    ## ??? is this possible in brms models
                                    grepl("^disp",out$term) ~ "disp",
                                    TRUE ~ "cond")

  out$term <- stringr::str_remove(out$term,mkRE(prefs[["components"]],
                                                LB=TRUE))
  if (fix.intercept) {
      ## use lookahead/lookbehind: replace Intercept with word boundary
      ## or underscore before/after by (Intercept) - without removing
      ## underscores!
      out$term <- stringr::str_replace(out$term,
                                        "(?<=(\\b|_))Intercept(?=(\\b|_))",
                                        "(Intercept)")
  }
  out <- reorder_cols(out)
  return(out)
}


#' @importFrom stats quantile
#' @export
sigma.brmsfit <- function (object, ...)  {
    if (!("sigma" %in% names(object$fit)))
        return(1)
    if (!requireNamespace("rstanarm")) {
        warning("need to install rstanarm to use extract sigma from brms fits")
        return(NA)
    }
    stats::quantile(as.data.frame(object$fit)[["sigma"]], probs=0.5)
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
  ret <- dplyr::tibble(.fitted = pred[, "Estimate"])
  if (se.fit) ret[[".se.fit"]] <- pred[, "Est.Error"]
  if (is.null(newdata)) {
    ret[[".resid"]] <- stats::residuals(x)[, "Estimate"]
    ret <- dplyr::bind_cols(as_tibble(data), ret)
  } else {
    ret <- dplyr::bind_cols(as_tibble(newdata), ret)
  }
  return(ret)
}
