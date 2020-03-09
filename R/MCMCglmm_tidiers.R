## n.b. this stuff is not used ATM; see tidy.MCMCglmm at the *end* of this file
## (move to MCMC_misc?)
## author: Josh Wiley
## originally from https://github.com/JWiley/postMCMCglmm

#' Extract the parameter names from an \code{MCMCglmm} object
#'
#' Simple function to extract the fixed and random effects
#' parameter names from an \code{MCMCglmm} object.
#'
#' @param object An \code{MCMCglmm} object
#' @param \dots not used
#' @return A list with two elements:
#'   \item{fixed}{A character vector of the fixed effects parameter names}
#'   \item{random}{A character vector of the random effects parameter names}
#' @export
#' @seealso \code{\link{fixef.MCMCglmm}}, \code{\link{ranef.MCMCglmm}}
#' @keywords internal
#' @importFrom stats terms
#' @examples
#' \dontrun{
#' # a simple MCMCglmm model
#' if (require(MCMCglmm)) {
#'   data(PlodiaPO)
#'   m <- MCMCglmm(PO ~ 1, random = ~FSfamily, data = PlodiaPO, verbose = FALSE, pr = TRUE)
#' }
#' # extract the parameter names
#' paramNamesMCMCglmm(m)
#' }
paramNamesMCMCglmm <- function(object, ...) {
  fNames <- as.character(attr(terms(object$Fixed$formula), "variables"))[-c(1, 2)]

  rNames <- NULL

  if (!is.null(object$Random$formula)) {
    rNames <- as.character(attr(terms(object$Random$formula), "variables"))[-c(1)]
  }

  if ("(Intercept)" %in% colnames(object[["Sol"]])) {
    fNames <- c("(Intercept)", fNames)
  }

  list(fixed = fNames, random = rNames)
}

#' Extract the levels of factors used for random effects in \code{MCMCglmm} objects
#'
#' @param object An \code{MCMCglmm} model object
#' @param data The dataset used for the model
#' @param \dots Not currently used
#' @export
#' @seealso \code{\link{paramNamesMCMCglmm}}, \code{\link{ranef.MCMCglmm}}
#' @examples
#' \dontrun{
#' # a simple MCMCglmm model
#' data(PlodiaPO)
#' m <- MCMCglmm(PO ~ 1, random = ~FSfamily, data = PlodiaPO, verbose = FALSE)
#'
#' # extract the random effects levels
#' ranefLevels(m, PlodiaPO)
#' }
ranefLevels <- function(object, data, ...) {
  n <- paramNamesMCMCglmm(object)$random
  res <- lapply(n, function(n) {
    levels(data[, n])
  })
  names(res) <- n
  return(res)
}

#' Internal function to extract the fixed or random effects from an \code{MCMCglmm} object
#'
#' Extracts the fixed or random effects portions from an MCMCglmm object.
#' Note for the random, these are the estimates themselves,
#' not the variability in the estimates. The \code{use} options let you
#' get either just the posterior mean or all the posterior samples.
#'
#' @param object An \code{MCMCglmm} object
#' @param use A character string indicating whether to return \dQuote{all} the
#'   posterior samples (the default) or only the \dQuote{mean} of them.
#' @param which A character string indicating whether to return the
#'   \dQuote{fixed} (the default) or \dQuote{random} effects.
#' @param \dots Not currently used.
#' @return A matrix of the posterior samples or means for the fixed or random effects.
#' @keywords internal
#' @seealso \code{\link{fixef.MCMCglmm}}, \code{\link{ranef.MCMCglmm}}
#' @rdname extractEffects
.extractEffects <- function(object, use = c("all", "mean"),
                            which = c("fixed", "random"), ...) {
  use <- match.arg(use)
  which <- match.arg(which)

  b <- as.matrix(object[["Sol"]])

  eff <- switch(which,
    fixed = {
      # this does not work because the intercept term contains
      # special characters (), that screw with the regular expressions
      # cannot use fixed matching because factor variables can
      # expand with arbitrarily named levels
      # unlist(lapply(paramNamesMCMCglmm(object)$fixed, function(n) {
      #  grep(paste0("^", n, ".*$"), colnames(b), value = TRUE)
      # }))
      object$X@Dimnames[[2]]
    },
    random = {
      regex <- paste(paramNamesMCMCglmm(object)$random, collapse = "|")
      regex <- paste0("^(", regex, ")\\..*$")
      grep(regex, colnames(b), value = TRUE)
    }
  )

  b <- b[, eff, drop = FALSE]

  switch(use,
    all = t(b),
    mean = as.matrix(colMeans(b))
  )
}

#' Extract fixed effects from an \code{MCMCglmm} object
#'
#' Function designed to extract the fixed effects from an
#' \code{MCMCglmm} model object. Can either extract all samples from the
#' fixed effects posteriors or return the posterior means.
#'
#' @param object An \code{MCMCglmm} model object to extract the effects from
#' @param use A character string indicating whether to extract
#'   all posterior samples or the mean of the posteriors. Defaults to
#'   "all".
#' @param \dots Arguments passed on to the worker function.
#' @return A matrix of the fixed effects
#' @importFrom nlme fixef
#' @method fixef MCMCglmm
#' @export
#' @seealso \code{\link{ranef.MCMCglmm}}
#' @examples
#' \dontrun{
#' # a simple MCMCglmm model
#' data(PlodiaPO)
#' m <- MCMCglmm(PO ~ 1, random = ~FSfamily, data = PlodiaPO, verbose = FALSE)
#'
#' # only extract average fixed effects
#' fixef(m, use = "mean")
#'
#' # histogram of posterior samples of fixed effects
#' hist(fixef(m))
#' # matches the mean
#' rowMeans(fixef(m))
#' }
fixef.MCMCglmm <- function(object, use = c("all", "mean"), ...) {
  .extractEffects(object = object, use = use, which = "fixed", ...)
}

#' Extract random effects from an \code{MCMCglmm} object
#'
#' Function designed to extract the random effects from an
#' \code{MCMCglmm} model object. Can either extract all samples from the
#' random effects posteriors or return the posterior means.
#'
#' @param object An \code{MCMCglmm} model object to extract the effects from
#' @param use A character string indicating whether to extract
#'   all posterior samples or the mean of the posteriors. Defaults to
#'   "all".
#' @param \dots Arguments passed on to the worker function.
#' @return A matrix of the fixed effects
#' @importFrom nlme ranef
#' @method ranef MCMCglmm
#' @export
#' @seealso \code{\link{fixef.MCMCglmm}}
#' @examples
#' \dontrun{
#' # a simple MCMCglmm model
#' data(PlodiaPO)
#' m <- MCMCglmm(PO ~ 1, random = ~FSfamily, data = PlodiaPO, pr = TRUE, verbose = FALSE)
#'
#' # only extract average fixed effects
#' head(ranef(m, use = "mean"))
#'
#' # histogram of posterior samples of fixed effects
#' hist(ranef(m)[1, ])
#' # matches the mean
#' rowMeans(ranef(m)[1:6, ])
#' }
ranef.MCMCglmm <- function(object, use = c("all", "mean"), ...) {
  .extractEffects(object = object, use = use, which = "random", ...)
}


#' Extract standard deviation of "random" effects from an \code{MCMCglmm} object
#'
#' Function designed to extract the standard deviation of the
#' random effects from an \code{MCMCglmm} model object.
#' Note that this is not the same as the posterior distribution of
#' (co)variance matrices. It is based on the posterior distribution
#' of the random effects. This also means it requires \code{pr=TRUE}
#' to be set in the model for the information to be saved. Can optionally
#' return standard deviation of random effects after back transforming to
#' the response metric. Currently probabilities, but only for ordinal family
#' models (\code{family="ordinal"}).
#'
#' @param object An \code{MCMCglmm} model object to extract the effects from
#' @param which A list of random effects to extract or their numeric positions
#'   If there are two numbers in a list, effects are simulataneous.
#' @param type A character string indicating whether to calculate the standard
#'   deviation on the linear predictor metric, \sQuote{lp} or
#'   response, \sQuote{response}.
#' @param \dots Not currently used.
#' @return A list of class postMCMCglmmRE with means (\code{M}) and individual estimates (\code{Data})
#' @export
#' @examples
#' \dontrun{
#' # a simple MCMCglmm model
#' data(PlodiaPO)
#' PlodiaPO <- within(PlodiaPO, {
#'   PO2 <- cut(PO, quantile(PO, c(0, .33, .66, 1)))
#'   plate <- factor(plate)
#' })
#'
#' m <- MCMCglmm(PO2 ~ 1,
#'   random = ~ FSfamily + plate,
#'   family = "ordinal", data = PlodiaPO,
#'   prior = list(
#'     R = list(V = 1, fix = 1),
#'     G = list(
#'       G1 = list(V = 1, nu = .002),
#'       G2 = list(V = 1, nu = .002)
#'     )
#'   ), verbose = FALSE, thin = 1, pr = TRUE
#' )
#'
#' # summary of the model
#' summary(m)
#'
#' # examples of extracting standard deviations of
#' # different random effects on the linear predictor metric
#' # or after transformation to probabilities (only for ordinal)
#' stdranef(m, which = list(1), type = "lp")
#' stdranef(m, which = list(2), type = "lp")
#' stdranef(m, which = list(1, 2, c(1, 2)), type = "lp")
#' stdranef(m, type = "lp")
#'
#' ## error because no 3rd random effect
#' # stdranef(m, which = list(1, 2, 3), type = "lp")
#'
#' stdranef(m, which = list("FSfamily", "plate"), type = "lp")
#'
#' # mean standard deviations on the probability metric
#' # also the full distributions, if desired in the Data slot.
#' res <- stdranef(m, type = "response")
#' res$M # means
#' hist(res$Data$FSfamily[, 1]) # histogram
#' }
stdranef <- function(object, which, type = c("lp", "response"), ...) {
  type <- match.arg(type)

  if (is.null(object$Z)) stop("Z matrix must be saved")
  ## z <- object$Z
  z <- diag(ncol(object$Z))
  colnames(z) <- colnames(object$Z)

  re <- paramNamesMCMCglmm(object)$random

  if (missing(which)) which <- c(re, list(re))

  stopifnot(is.list(which))

  if (is.numeric(unlist(which))) {
    stopifnot(all(unlist(which) %in% seq_along(re)))
    which <- lapply(which, function(i) re[i])
  } else {
    stopifnot(all(unlist(which) %in% re))
  }

  index <- lapply(which, function(n) {
    n <- paste(n, collapse = "|")
    regex <- paste0("^(", n, ")\\..*$")
    index <- grep(regex, colnames(z))
    return(index)
  })

  # tricky, note that coefficients and predicted probabilities
  # do not come out with the same dimensions, the are transposed
  # samples on columns for coefs, samples on rows for predicted probs
  res <- switch(type,
    lp = {
      yhat <- lapply(index, function(i) ranef(object, use = "all")[i, ])
      lapply(yhat, function(m) matrix(apply(m, 2, sd)))
    }, response = {
      yhat <- lapply(index, function(i) {
        tmp <- z
        if (length(i) < ncol(tmp)) {
          tmp[, -i] <- 0L # zero out all random effects we are not interested in
        }
        predict2.MCMCglmm(object, X = NULL, Z = tmp, use = "all", type = type)
      })
      lapply(yhat, function(m) sapply(m, function(n) apply(n, 1, sd)))
    }
  )

  names(res) <- which

  M <- do.call(rbind, lapply(res, colMeans))

  finalres <- list(M = M, Data = res)
  class(finalres) <- "postMCMCglmmRE"

  return(finalres)
}

##' @rdname mcmc_tidiers
##' @examples
##' if (require("MCMCglmm")) {
##'   ## original model
##'   \dontrun{
##'       mm0 <- MCMCglmm(Reaction ~ Days,
##'                  random = ~Subject, data = sleepstudy,
##'                  nitt=4000,
##'                  pr = TRUE
##'              )
##'    }
##'    ## load stored object
##'    load(system.file("extdata","MCMCglmm_example.rda",
##'                                      package="broom.mixed"))
##'    tidy(mm0)
##'    tidy(mm1)
##'    tidy(mm2)
##'    tail(tidy(mm0,effects="ran_vals"))
##' }
##' @export
tidy.MCMCglmm <- function(x, effects = c("fixed", "ran_pars"),
                          scales = NULL, ## c("sdcor","vcov",NA),
                          ...) {
  ## FIXME: allow scales= parameter to get varcov on sd/corr scale?
  clist <- c(fixed = "Sol", ran_pars = "VCV", ran_vals = "Sol")
  comp <- clist[effects]
  if (!is.null(scales)) {
    if (length(scales) != length(effects)) {
      stop(
        "if scales are specified, values (or NA) must be provided ",
        "for each effect"
      )
    }
  }

  ## FIXME: override MCMCglmm internal component names
  ## FIXME:: have to work harder to retrieve group/term information
  ##  about random parameters
  ## individual components are mcmc objects: call tidy on them
  retList <- (purrr::map(x[comp], tidy, ...)
  %>%
    setNames(effects))
  ## see .extractEffects
  fnames <- x$X@Dimnames[[2]]

  if ("fixed" %in% effects) {
    retList$fixed <- filter(retList$fixed, term %in% fnames)
  }
  if ("ran_vals" %in% effects) {
    retList$ran_vals <- filter(retList$ran_vals, !(term %in% fnames))
    if (nrow(retList$ran_vals) == 0) {
      stop("for tidying random effects values, must run MCMglmm with pr=TRUE")
    }
    ss <- strsplit(retList$ran_vals$term, "\\.")
    retList$ran_vals$level <- sapply(ss, utils::tail, 1)
    retList$ran_vals$group <- sapply(ss, function(x) x[length(x) - 1])
    retList$ran_vals$term <- sapply(
      ss,
      function(x) if (length(x) == 3) x[1] else "(Intercept)"
    )
  }

  if ("ran_pars" %in% effects) {
    ss <- strsplit(retList$ran_pars$term, "(:|\\.)")
    if (any(sapply(ss, length) > 4)) {
      warning("surprising term names: models with variable names containing colons or dots may work unreliably")
    }
    ss2 <- lapply(
      ss,
      function(x) {
        if (x[1] == "units") {
          return(c(rep("Observation", 2), "Residual"))
        }
        if (length(x) == 1) {
          return(c(rep("(Intercept)", 2), x[1]))
        }
        ## reconstruct variance term if necessary
        if (length(x) == 2) {
          return(c(rep(x[1], 2), x[2]))
        }
        return(c(sort(x[1:2]), x[3]))
      }
    )
    ok <- !duplicated(ss2)
    ss2 <- ss2[ok]
    if (is.null(scales)) {
      rscale <- "vcov"
    } else {
      rscale <- scales[effects == "ran_pars"]
    }
    if (!rscale %in% c("sdcor", "vcov")) {
      stop(sprintf("unrecognized ran_pars scale %s", sQuote(rscale)))
    }
    if (rscale != "vcov") {
      stop("only vcov scale implemented")
      ## ugh: to implement sdcor we have to go back to the samples
      ## and convert them ...
    }
    pref <- sapply(ss2, function(x) ifelse(x[1] == x[2], "var", "cov"))
    group <- sapply(ss2, function(x) x[[3]])
    term <- paste(pref,
      sapply(ss2, function(x) {
        ifelse(x[1] == x[2], x[1],
          paste(x[1:2], collapse = ".")
        )
      }),
      sep = getOption("broom.mixed.sep1")
    )
    retList$ran_pars <- retList$ran_pars[ok, ]
    retList$ran_pars$term <- term
    retList$ran_pars$group <- group
  }
  ret <- (retList
  %>%
    bind_rows(.id = "effect")
    %>%
    reorder_cols()
  )
  return(ret)
}
