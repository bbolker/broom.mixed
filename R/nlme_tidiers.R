#' Tidying methods for mixed effects models
#'
#' These methods tidy the coefficients of mixed effects models
#' of the \code{lme} class from functions of the \code{nlme} package.
#'
#' @param x An object of class \code{lme}, such as those from \code{lme}
#' or \code{nlme}
#'
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#'
#' @name nlme_tidiers
#'
#' @inheritParams tidy.merMod
#'
#' @examples
#'
#' if (require("nlme") && require("lme4")) {
#'     data("sleepstudy", package="lme4")
#'     ## original model
#'     \dontrun{
#'          lmm1 <- lme(Reaction ~ Days, random=~ Days|Subject, sleepstudy)
#'     }
#'     ## load stored object
#'     load(system.file("extdata","nlme_example.rda", package="broom.mixed"))
#'     tidy(lmm1)
#'     tidy(lmm1, effects = "fixed")
#'     tidy(lmm1, conf.int = TRUE)
#'     tidy(lmm1, effects = "ran_pars")
#'     tidy(lmm1, effects = "ran_vals")
#'     tidy(lmm1, effects = "ran_coefs")
#'     head(augment(lmm1, sleepstudy))
#'     glance(lmm1)
#'
#'     startvec <- c(Asym = 200, xmid = 725, scal = 350)
#'     nm1 <- nlme(circumference ~ SSlogis(age, Asym, xmid, scal),
#'                   data = Orange,
#'                   fixed = Asym + xmid + scal ~1,
#'                   random = Asym ~1,
#'                   start = startvec)
#'     tidy(nm1)
#'     tidy(nm1, effects = "fixed")
#'     head(augment(nm1, Orange))
#'     glance(nm1)
#'
#'     gls1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
#'                          correlation = corAR1(form = ~ 1 | Mare))
#'     tidy(gls1)
#'     glance(gls1)
#'     head(augment(gls1))
#' }
#'
#' @rdname nlme_tidiers
#'
#' @param effects One or more of "var_model", "ran_pars", "fixed", "ran_vals",
#'   and/or "ran_coefs".
#'
#' @return \code{tidy} returns one row for each estimated effect, either
#' random or fixed depending on the \code{effects} parameter. If
#' \code{effects = "ran_vals"} (or \code{"ran_pars"}), it contains the columns
#'   \item{group}{the group within which the random effect is being estimated}
#'   \item{level}{level within group}
#'   \item{term}{term being estimated}
#'   \item{estimate}{estimated coefficient}
#'   \item{estimated}{This column is only included if some parameters are fixed.
#'     TRUE if the residual error is estimated and FALSE if the residual error
#'     is fixed.}
#'
#' If \code{effects="fixed"}, \code{tidy} returns the columns
#'   \item{term}{fixed term being estimated}
#'   \item{estimate}{estimate of fixed effect}
#'   \item{std.error}{standard error}
#'   \item{statistic}{t-statistic}
#'   \item{p.value}{P-value computed from t-statistic}
#'
#' If \code{effects="var_model"} (the \code{weights} argument to the model),
#' \code{tidy} returns the columns defined in the help for \code{tidy.varFunc}.
#'
#' @importFrom nlme getVarCov intervals
#' @import dplyr
#' @importFrom tidyr replace_na
## @importFrom dplyr tibble select full_join
#'
#' @export
tidy.lme <- function(x, effects = c("var_model", "ran_pars", "fixed"),
                     scales = NULL,
                     conf.int = FALSE,
                     conf.level = 0.95,
                     ...) {

  ## R CMD global var check
  lower <- upper <- NULL

  effect_names <- c("var_model", "ran_pars", "fixed", "ran_vals", "ran_coefs")
  if (length(miss <- setdiff(effects, effect_names)) > 0) {
    stop("unknown effect type ", miss)
  }

  ## R CMD check false positives
  term <- estimate <- .id <- level <- std.error <- NULL

  if (!is.null(scales)) {
    if (length(scales) != length(effects)) {
      stop(
        "if scales are specified, values (or NA) must be provided ",
        "for each effect"
      )
    }
  }

  ret_list <- list()

  if ("fixed" %in% effects) {
    # return tidied fixed effects
    ret <- summary(x)[["tTable"]] %>%
      data.frame(check.names = FALSE) %>%
      rename_regex_match() %>%
      tibblify("term")
    if (conf.int) {
      cifix <- intervals(x, which = "fixed")[["fixed"]] %>%
        data.frame() %>%
        dplyr::select(lower, upper) %>%
        setNames(c("conf.low", "conf.high")) %>%
        tibble::rownames_to_column("term")
      ret <- dplyr::full_join(ret, cifix, by = "term")
    }

    ran_effs <- sprintf("ran_%s", c("pars", "vals", "coefs"))

    ret_list$fixed <- ret %>%
      reorder_frame()
  }


  if ("ran_pars" %in% effects) {
    if (is.null(scales)) {
      rscale <- "sdcor"
    } else {
      rscale <- scales[effects == "ran_pars"]
    }
    if (!rscale %in% c("sdcor", "vcov")) {
      stop(sprintf("unrecognized ran_pars scale %s", sQuote(rscale)))
    }
    grplen <- attr(x$modelStruct$reStruct, "plen")
    multilevel <- (length(grplen) > 1)
    nonlin <- inherits(x, "nlme")

    if (multilevel || nonlin) {
      warning(
        "ran_pars not yet implemented for ",
        if (multilevel) {
          "multiple levels of nesting"
        } else {
          "nonlinear models"
        }
      )
      ret <- dplyr::tibble()
    } else {
      vc <- nlme::getVarCov(x)
      ran_prefix <- switch(rscale,
        vcov = c("var", "cov"),
        sdcor = c("sd", "cor")
      )
      ## construct appropriate sd/cor names from dimnames
      nmvec <- outer(
        colnames(vc), rownames(vc),
        function(x, y) {
          ifelse(x == y,
            sprintf(
              "%s_%s",
              ran_prefix[1],
              x
            ),
            sprintf(
              "%s_%s.%s",
              ran_prefix[2],
              x, y
            )
          )
        }
      )
      lwrtri <- function(x) {
        x[lower.tri(x, diag = TRUE)]
      }
      nmvec <- c(
        lwrtri(nmvec),
        sprintf("%s_%s", ran_prefix[1], "Observation")
      )
      grpnames <- c(rep(names(grplen), grplen), "Residual")
      if (rscale == "vcov") {
        vals <- c(lwrtri(vc), sigma(x)^2)
      } else {
        vals <- cov2cor(vc)
        diag(vals) <- sqrt(diag(vc))
        vals <- c(lwrtri(vals), sigma(x))
      }
      ret <- dplyr::tibble(
        effect = "ran_pars",
        group = grpnames,
        term = c(nmvec),
        estimate = c(vals)
      )

      if (conf.int) {
        ii <- intervals(x, which = "var-cov")$reStruct
        trfun <- function(z) {
          nm <- rownames(z)
          nm <- sub(
            "\\(", "_",
            sub(
              ")$", "",
              sub(",", ".", nm)
            )
          )
          ## ugh, swap order of vars in cor term
          corterms <- grepl("^cor", nm)
          re <- "_([^\\.]+)\\.(.+)$"
          nm[corterms] <-
            sub(re, "_\\2.\\1", nm[corterms],
              perl = TRUE
            )
          return(dplyr::tibble(
            term = nm, conf.low = z[, "lower"],
            conf.high = z[, "upper"]
          ))
        }
        ci <- dplyr::bind_rows(lapply(ii, trfun), .id = "group")
        if (rscale != "sdcor") {
          ## FIXME: transform/recompute from scratch
          warning("confidence intervals for ran pars only available on sdcor scale")
          ci$conf.low <- ci$conf.high <- NA
        }
        ## FIXME: also do confint on residual
        ret <- dplyr::full_join(ret, ci, by = c("group", "term"))
      }
    } ## if not multi-level model
    if (attr(x$modelStruct, 'fixedSigma')) {
      mask_residual <- ret$group == "Residual"
      if (sum(mask_residual) != 1) {
        stop("More than one residual estimate found, please report this as a bug") # nocov
      }
      ret$estimated <- !mask_residual
    }
    ret_list$ran_pars <- ret
  }

  if ("ran_vals" %in% effects) {
    ret_list$ran_vals <-
      ranef(x) %>% as.data.frame()
  }
  if ("ran_coefs" %in% effects) {
    ret_list$ran_coefs <-
      stats::coef(x) %>%
      tibble::rownames_to_column("level") %>%
      tidyr::gather(key = term, value = estimate, -level)
    ## FIXME: group?
  }
  if ("var_model" %in% effects) {
    # x$modelStruct$varStruct will be NULL if the model does not use a
    # varStruct, so this will return an empty tibble.  So, we don't need to test
    # that varStruct is actually there.
    ret_list$var_model <- tidy(x$modelStruct$varStruct)
  }
  ret <- bind_rows(ret_list, .id = "effect") %>%
      dplyr::select(any_of(c("effect", "group", "level", "term",
                             "estimate", "estimated", "std.error", "df",
                             "statistic", "p.value", "conf.low", "conf.high")))
  if ("estimated" %in% names(ret)) {
    # NA in the estimated column at this point indicates that part of the mdoel
    # was estimated, but the "estimated" column was not set in that sub-part of
    # the tidier.
    ret$estimated <- tidyr::replace_na(ret$estimated, TRUE)
  }
  return(ret)
}

#' @rdname nlme_tidiers
#'
#' @param data original data this was fitted on; if not given this will
#' attempt to be reconstructed
#' @param newdata new data to be used for prediction; optional
#'
#' @template augment_NAs
#'
#' @return \code{augment} returns one row for each original observation,
#' with columns (each prepended by a .) added. Included are the columns
#'   \item{.fitted}{predicted values}
#'   \item{.resid}{residuals}
#'   \item{.fixed}{predicted values with no random effects}
#'
#'
#' @importFrom broom augment augment_columns
#' @export
augment.lme <- function(x, data = x$data, newdata, ...) {
  if (is.null(data)) {
    stop("augment.lme must be called with an explicit 'data' argument for this (n)lme fit  because of an inconsistency in nlme.")
  }
  # move rownames if necessary
  if (missing(newdata)) {
    newdata <- NULL
  }
  ret <- augment_columns(x, data, newdata, se.fit = NULL)

  # add predictions with no random effects (population means)
  predictions <- stats::predict(x, level = 0)
  if (length(predictions) == nrow(ret)) {
    ret$.fixed <- predictions
  }

  return(tibblify(ret, var = NULL))
}


#' @importFrom dplyr across mutate
#' @export
as.data.frame.ranef.lme <- function(x, row.names, optional=TRUE,
                                    stringsAsFactors = FALSE, ...) {
  where <- NULL
  ## see https://github.com/r-lib/tidyselect/issues/201
  ## this is better than utils::globalVariables() which is global ...

    group <- term <- level <- estimate <- NULL ## NSE arg checking
    if (!missing(row.names)) stop(sQuote("row.names"),
                                  "  argument not implemented")
    if (length(list(...)>0)) warning("additional arguments ignored")
    ## see ?as.data.frame: optional==FALSE corresponds to check.names==TRUE
    ## "check_unique" is the tibble default, "universal" behaves like
    ##  check.names == TRUE
    name_repair <- if (optional) "check_unique" else "universal"
    melt <- function(x) purrr::map_dfr(as.list(x),
                 ~tibble(level = rownames(x), estimate = .), .id = "term")
    grps <- attr(x, "grpNames")
    if (length(grps)==1) x <- list(x)
    mm <- purrr::map(x, melt)
    res <- purrr::map2_dfr(.x = mm, .y = grps, ~mutate(.x, group = .y)) %>%
      dplyr::select(group, term, level, estimate)
    if (stringsAsFactors) {
      res <- res %>% dplyr::mutate(across(where(is.character), as.factor))
    }
    return(res)
}

#' @rdname nlme_tidiers
#'
#' @param ... extra arguments (not used)
#'
#' @return \code{glance} returns one row with the columns
#'   \item{sigma}{the square root of the estimated residual variance}
#'   \item{logLik}{the data's log-likelihood under the model}
#'   \item{AIC}{the Akaike Information Criterion}
#'   \item{BIC}{the Bayesian Information Criterion}
#'   \item{deviance}{returned as NA. To quote Brian Ripley on R-help
#' \url{https://stat.ethz.ch/pipermail/r-help/2006-May/104744.html},
#'  "McCullagh & Nelder (1989) would be the authorative [sic] reference, but the 1982
#' first edition manages to use 'deviance' in three separate senses on one
#' page." }
#'
#' @export
glance.lme <- function(x, ...) {
  finish_glance(x = x)
}

#' @rdname nlme_tidiers
#' @export
tidy.gls <- function(x,
                     conf.int = FALSE,
                     conf.level = 0.95,
                     ...) {
  . <- Value <- Std.Error <- `t-value` <- `p-value` <- NULL ## glob var checks
  ret <- (summary(x)[["tTable"]]
      %>% as.data.frame() ## have to convert to df *first*
      %>% tibble::rownames_to_column(var = "term")
      %>% as_tibble()
      %>%  dplyr::rename(
                      estimate = Value,
                      std.error = Std.Error,
                      statistic = `t-value`,
                      p.value = `p-value`
                  )
  )
  if (conf.int) {
      cc <- (confint(x, level=conf.level)
          %>% as.data.frame()
          %>% setNames(c("conf.low","conf.high"))
      )
      ret <- dplyr::bind_cols(ret,cc)
  }
  return(ret)
}

#' @export
glance.gls <- function(x, ...) {
  ss <- summary(x)
  with(
    ss,
    tibble::tibble(sigma,
      df = dims[["p"]],
      logLik,
      AIC,
      BIC,
      df.residual = dims[["N"]] - dims[["p"]]
    )
  )
}

#' @rdname nlme_tidiers
#' @export
augment.gls <- function(x, data = nlme::getData(x), newdata, ...) {

  # move rownames if necessary
  if (missing(newdata)) {
    newdata <- NULL
  }
  ret <- augment_columns(x, data, newdata, se.fit = NULL)
  ret
}

# varFunc tidiers ####

#' Tidy variance structure for the \code{nlme} package.
#'
#' Returns a tibble with the following columns:
#' \itemize{
#' \item{group}{type of varFunc, along with the right hand side of the formula
#'   in parentheses e.g. \code{"varExp(age | Sex)"}.}
#' \item{term}{terms included in the formula of the variance model, specifically
#'   the names of the coefficients.  If the value is fixed, it will be appended
#'   with \code{" ; fixed"}.}
#' \item{estimate}{estimated coefficient}
#' \item{estimated}{This column is only included if some parameters are fixed.
#'   TRUE if the parameter is estimated and FALSE if the parameter is fixed.}
#' }
#'
#' @param x An object of class \code{varFunc}, such as those used as the
#'   \code{weights} argument from the \code{nlme} package
#' @param ... Ignored
#' @return If the \code{varFunc} is uninitialized or has no parameters, the
#'   function will return an empty tibble.  Otherwise, it will return a tibble with
#'   names described in the details section.
#' @examples
#' \dontrun{
#' if (require("nlme")) {
#' ChickWeight_arbitrary_group <- datasets::ChickWeight
#' ChickWeight_arbitrary_group$group_arb_n <-
#'   1 + (
#'     as.integer(ChickWeight_arbitrary_group$Chick) >
#'     median(as.integer(ChickWeight_arbitrary_group$Chick))
#'   )
#' ChickWeight_arbitrary_group$group_arb <- c("low", "high")[ChickWeight_arbitrary_group$group_arb_n]
#'
#' fit_with_fixed <-
#'   lme(
#'     weight ~ Diet * Time,
#'     random = ~Time | Chick,
#'     data =ChickWeight_arbitrary_group,
#'     weights=varIdent(fixed=c("low"=5), form=~1|group_arb)
#'   )
#' # Show all parameters
#' tidy(fit_with_fixed)
#' # Exclude fixed parameters
#' tidy(fit_with_fixed) %>%
#'   filter(across(any_of("estimated"), ~.x))
#' }
#' }
#' @export
#' @importFrom tibble tibble
#' @importFrom stats as.formula
tidy.varFunc <- function(x, ...) {
  aux <- coef(x, unconstrained = FALSE, allCoef = TRUE)
  if (length(aux) == 0) {
    warning(
      "Variance function structure of class", class(x)[1],
      "with no parameters, or uninitialized"
    )
    return(tibble::tibble())
  }
  x_formula <- stats::as.formula(x)
  if (!is.null(x_formula)) {
    group_text <- sprintf("%s(%s)", class(x)[1], as.character(x_formula)[2])
  } else {
    group_text <- class(x)[1]
  }
  term_text <- names(aux)
  ret <-
    tibble::tibble(
      #effect="var_model",
      group=group_text,
      term=term_text,
      estimate=aux
    )
  if (any(attr(x, "whichFix"))) {
    # Detect fixed parameters based on name.  Name order may differ from the
    # order in the 'whichFix' attribute, so matching must be done by name.
    ret$estimated <- !(ret$term %in% names(attr(x, "fixed")))
  }
  ret
}

#' @rdname tidy.varFunc
#' @export
tidy.varComb <- function(x, ...) {
  # A varComb object is a named list ("A", "B", "C", ...) of varFunc objects
  ret <-
    dplyr::bind_rows(lapply(
      X=x,
      FUN=tidy,
      ...
    ))
  ret
}
