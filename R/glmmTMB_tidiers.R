#' Tidying methods for glmmTMB models
#'
#' These methods tidy the coefficients of mixed effects models, particularly
#' responses of the \code{merMod} class
#'
#' @param x An object of class \code{merMod}, such as those from \code{lmer},
#' \code{glmer}, or \code{nlmer}
#'
#' @return All tidying methods return a \code{tibble}.
#' The structure depends on the method chosen.
#'
#' @name glmmTMB_tidiers
#'
#' @examples
#' if (require("glmmTMB") && require("lme4")) {
#'     data("sleepstudy",package="lme4")
#'     ## original model:
#'     \dontrun{
#'         lmm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
#'     }
#'     ## load stored object
#'     load(system.file("extdata","glmmTMB_example.rda",package="broom.mixed"))
#'     tidy(lmm1)
#'     tidy(lmm1, effects = "fixed")
#'     tidy(lmm1, effects = "fixed", conf.int=TRUE)
#'     tidy(lmm1, effects = "fixed", conf.int=TRUE, conf.method="uniroot")
#'     ## FIX: tidy(lmm1, effects = "ran_vals", conf.int=TRUE)
#'     head(augment(lmm1, sleepstudy))
#'     glance(lmm1)
#'
#'     ## original model:
#'     ##  glmm1 <- glmmTMB(incidence/size ~ period + (1 | herd),
#'     ##                  data = cbpp, family = binomial, weights=size)
#'     tidy(glmm1)
#'     tidy(glmm1, effects = "fixed")
#'     tidy(glmm1, effects = "fixed", exponentiate=TRUE)
#'     tidy(glmm1, effects = "fixed", conf.int=TRUE, exponentiate=TRUE)
#'     head(augment(glmm1, cbpp))
#'     head(augment(glmm1, cbpp, type.residuals="pearson"))
#'     glance(glmm1)
#' \dontrun{
#' ## profile CIs - a little bit slower but more accurate
#'     tidy(glmm1, effects = "fixed", conf.int=TRUE, exponentiate=TRUE, conf.method="profile")
#' }
#' }
NULL


#' @rdname glmmTMB_tidiers
#'
#' @inheritParams lme4_tidiers
#' @param effects A character vector including one or more of "fixed" (fixed-effect parameters), "ran_pars" (variances and covariances or standard deviations and correlations of random effect terms) or "ran_vals" (conditional modes/BLUPs/latent variable estimates)
#' @param component which component to extract (e.g. \code{cond} for conditional effects (i.e., traditional fixed effects); \code{zi} for zero-inflation model; \code{disp} for dispersion model
#' @param conf.int whether to include a confidence interval
#' @param conf.level confidence level for CI
#' @param conf.method method for computing confidence intervals (see \code{\link[lme4]{confint.merMod}})
#' @param scales scales on which to report the variables: for random effects, the choices are \sQuote{"sdcor"} (standard deviations and correlations: the default if \code{scales} is \code{NULL}) or \sQuote{"varcov"} (variances and covariances). \code{NA} means no transformation, appropriate e.g. for fixed effects; inverse-link transformations (exponentiation
#' or logistic) are not yet implemented, but may be in the future.
#' @param ran_prefix a length-2 character vector specifying the strings to use as prefixes for self- (variance/standard deviation) and cross- (covariance/correlation) random effects terms
#'
#' @return \code{tidy} returns one row for each estimated effect, either
#' with groups depending on the \code{effects} parameter.
#' It contains the columns
#'   \item{group}{the group within which the random effect is being estimated: \code{NA} for fixed effects}
#'   \item{level}{level within group (\code{NA} except for modes)}
#'   \item{term}{term being estimated}
#'   \item{estimate}{estimated coefficient}
#'   \item{std.error}{standard error}
#'   \item{statistic}{t- or Z-statistic (\code{NA} for modes)}
#'   \item{p.value}{P-value computed from t-statistic (may be missing/NA)}
#'
#' @note zero-inflation parameters (including the intercept) are reported
#' on the logit scale
#'
#' @importFrom plyr ldply rbind.fill
#' @import dplyr
#' @importFrom tidyr gather spread
#' @importFrom nlme VarCorr ranef
#' @importFrom stats qnorm confint coef na.omit setNames
## FIXME: is it OK/sensible to import these from (priority='recommended')
## nlme rather than (priority=NA) lme4?
#'
#' @export
tidy.glmmTMB <- function(x, effects = c("ran_pars", "fixed"),
                         component = c("cond", "zi"),
                         scales = NULL, ## c("sdcor",NA),
                         ran_prefix = NULL,
                         conf.int = FALSE,
                         conf.level = 0.95,
                         conf.method = "Wald",
                         exponentiate = FALSE,
                         ...) {

  ## FIXME:  cleanup
  ##   - avoid (as.)data.frame

  ## R CMD check false positives
  condsd <- condval <- grp <- grpvar <-
        term <- estimate <- .id <- level <- std.error <- . <- NULL

  drop.missing <- function(x) x[vapply(x,length,numeric(1))>0]
  ss <- stats::coef(summary(x))
  ss <- drop.missing(ss)
  ## FIXME: warn if !missing(component) and component includes
  ##  NULL terms
  component <- intersect(component, names(ss))
  if (length(component[!component %in% c("cond", "zi")]) > 0L) {
    stop("only works for conditional and (partly for) zero-inflation components")
  }
  ss <- ss[component]
  effect_names <- c("ran_pars", "fixed", "ran_vals")
  if (!is.null(scales)) {
    if (length(scales) != length(effects)) {
      stop(
        "if scales are specified, values (or NA) must be provided ",
        "for each effect"
      )
    }
  }
  if (length(miss <- setdiff(effects, effect_names)) > 0) {
    stop("unknown effect type ", miss)
  }
  ret <- list()
  ret_list <- list()
  if ("fixed" %in% effects) {
    # return tidied fixed effects rather than random
    ret <- lapply(
      ss,
      function(x) {
        x %>%
          as.data.frame(stringsAsFactors = FALSE) %>%
          setNames(c("estimate", "std.error", "statistic", "p.value")) %>%
          tibble::rownames_to_column("term")
      }
    )
    # p-values may or may not be included
    # HACK: use the columns from the conditional component, preserving previous behaviour
    if (conf.int) {
      for (comp in component) {
        cifix <- confint(x,
          method = tolower(conf.method),
          level = conf.level,               
          component = comp,
          estimate = FALSE,
          ## conditional/zi components
          ## include random-effect parameters
          ## as well, don't want those right now ...
          parm = seq(nrow(ret[[comp]])), ...
        ) %>%
          as.data.frame(stringsAsFactors = FALSE) %>%
          setNames(c("conf.low", "conf.high"))
        ret[[comp]] <- bind_cols(
          ret[[comp]],
          cifix
        )
      }
    }
    ret_list$fixed <- bind_rows(ret, .id = "component")

    if (exponentiate) {
      vv <- intersect(c("estimate", "conf.low", "conf.high"), names(ret_list$fixed))
      ret_list$fixed <- (ret_list$fixed
          %>% mutate_at(vars(vv), ~exp(.))
          %>% mutate(std.error = std.error * estimate)
      )
    }

  }
  if ("ran_pars" %in% effects &&
    !all(sapply(VarCorr(x), is.null))) {
    ## FIXME: do something sensible about standard errors, confint

    if (is.null(scales)) {
      rscale <- "sdcor"
    } else {
      rscale <- scales[effects == "ran_pars"]
    }
    if (!rscale %in% c("sdcor", "vcov")) {
      stop(sprintf("unrecognized ran_pars scale %s", sQuote(rscale)))
    }
    ## kluge for now ...
    vv <- list()
    if ("cond" %in% component) {
      vv$cond <- VarCorr(x)[["cond"]]
      class(vv$cond) <- "VarCorr.merMod"
    }
    if ("zi" %in% component) {
      if (!is.null(vv$zi <- VarCorr(x)[["zi"]])) {
        class(vv$zi) <- "VarCorr.merMod"
      }
    }

    ret <- (
      purrr::map(vv, as.data.frame, stringsAsFactors = FALSE)
      %>%
        bind_rows(.id = "component")
        %>%
        mutate_if(., is.factor, as.character)
    )
    if (is.null(ran_prefix)) {
      ran_prefix <- switch(rscale,
        vcov = c("var", "cov"),
        sdcor = c("sd", "cor")
      )
    }

    ## DRY! refactor glmmTMB/lme4 tidiers

    ## don't try to assign as rowname (non-unique anyway),
    ## make it directly into a term column
    if (nrow(ret)>0) {
        ret[["term"]] <- apply(ret[c("var1", "var2")], 1,
                               ran_pars_name,
                               ran_prefix = ran_prefix
                               )

       ## keep only desired term, rename
       ## FIXME: should use select + tidyeval + rename ... ?
       ranpar_names <- c("component", "group", "term", "estimate")
       ret <- setNames(
          ret[c("component", "grp", "term", rscale)],
          ranpar_names
       )
    } else {
        ret <- dplyr::tibble(component=character(0),
                      group=character(0),
                      term=character(0),
                      estimate=numeric(0))
    }
    ## rownames(ret) <- seq(nrow(ret))

      if (conf.int) {
        thpar <- "theta_"
        if (utils::packageVersion("glmmTMB")<="0.2.2.0") {
             thpar <- which(names(x$obj$par)=="theta")
        }
        ciran <- (confint(x,
                          ## for next glmmTMB (> 0.2.3) can be "theta_",
                          parm = thpar,
                          method = conf.method,
                          level = conf.level,                
                          estimate = FALSE,
                          ...
      )
      %>% as_tibble()
      %>% setNames(c("conf.low", "conf.high"))
      )
      ret <- bind_cols(ret, ciran)
    }
    ret_list$ran_pars <- ret
  }

  if ("ran_vals" %in% effects) {
    ## fix each group to be a tidy data frame

      ret <- (ranef(x, condVar = TRUE)
        %>% as.data.frame(stringsAsFactors = FALSE)
        %>% dplyr::rename(
                       group = grpvar, level = grp,
                       estimate = condval, std.error = condsd
                   )
        %>% dplyr::mutate_if(is.factor, as.character)
    )

    if (conf.int) {
      if (conf.method != "Wald") {
        stop("only Wald CIs available for conditional modes")
      }

      mult <- qnorm((1 + conf.level) / 2)
      ret <- transform(ret,
        conf.low = estimate - mult * std.error,
        conf.high = estimate + mult * std.error
      )
    }

    ret_list$ran_vals <- ret
  }
  ret <- (ret_list
      %>% dplyr::bind_rows(.id = "effect")
      %>% as_tibble()
      %>% reorder_cols()
  )
  return(ret)
}

#' @rdname glmmTMB_tidiers
#'
#' @template augment_NAs
#'
#' @param data original data this was fitted on; if not given this will
#' attempt to be reconstructed
#' @param newdata new data to be used for prediction; optional

#' @return \code{augment} returns one row for each original observation,
#' with columns (each prepended by a .) added. Included are the columns
#'   \item{.fitted}{predicted values}
#'   \item{.resid}{residuals}
#'   \item{.fixed}{predicted values with no random effects}
#'
#' @export
augment.glmmTMB <- function(x, data = stats::model.frame(x), newdata=NULL,
                            ...) {
  broom::augment_columns(x, data, newdata, ...)
}

#' @rdname glmmTMB_tidiers
#'
#' @param ... extra arguments (not used)
#'
#' @return \code{glance} returns one row with the columns
#'   \item{sigma}{the square root of the estimated residual variance}
#'   \item{logLik}{the data's log-likelihood under the model}
#'   \item{AIC}{the Akaike Information Criterion}
#'   \item{BIC}{the Bayesian Information Criterion}
#'   \item{deviance}{deviance}
#'
#' @rawNamespace if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigma)
#' @export
glance.glmmTMB <- function(x, ...) {
  finish_glance(x = x)
}
