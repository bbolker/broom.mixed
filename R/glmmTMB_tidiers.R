#' Tidying methods for glmmTMB models
#' 
#' These methods tidy the coefficients of mixed effects models, particularly
#' responses of the \code{merMod} class
#' 
#' @param x An object of class \code{merMod}, such as those from \code{lmer},
#' \code{glmer}, or \code{nlmer}
#' 
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#' 
#' @name glmmTMB_tidiers
#'
#' @examples
#' if (require("broom") && require("glmmTMB") && require("lme4")) {
#'     
#'     # example regressions are from lme4 documentation
#'     ## lmm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
#'     load(system.file("extdata","glmmTMB_example.rda",package="broom.mixed"))
#'     tidy(lmm1)
#'     tidy(lmm1, effects = "fixed")
#'     tidy(lmm1, effects = "fixed", conf.int=TRUE)
#'     \dontrun{
#'        ## FIXME: fails on stored model fit?
#'        tidy(lmm1, effects = "fixed", conf.int=TRUE, conf.method="uniroot")
#'     }
#'     ## tidy(lmm1, effects = "ran_modes", conf.int=TRUE)
#'     head(augment(lmm1, sleepstudy))
#'     glance(lmm1)
#'     
#'     ## glmm1 <- glmmTMB(incidence/size ~ period + (1 | herd),
#'     ##                  data = cbpp, family = binomial, weights=size)
#'     tidy(glmm1)
#'     tidy(glmm1, effects = "fixed")
#'     head(augment(glmm1, cbpp))
#'     head(augment(glmm1, cbpp, type.residuals="pearson"))
#'     glance(glmm1)
#' }
NULL


#' @rdname glmmTMB_tidiers
#'
#' @param effects A character vector including one or more of "fixed" (fixed-effect parameters), "ran_pars" (variances and covariances or standard deviations and correlations of random effect terms) or "ran_modes" (conditional modes/BLUPs/latent variable estimates)
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
#'   \item{group}{the group within which the random effect is being estimated: \code{"fixed"} for fixed effects}
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
tidy.glmmTMB <- function(x, effects = c("ran_pars","fixed"),
                         component = c("cond","zi"),
                         scales = NULL, ## c("sdcor",NA),
                         ran_prefix=NULL,
                         conf.int = FALSE,
                         conf.level = 0.95,
                         conf.method = "Wald",
                        ...) {

    ## R CMD check false positives
    term <- estimate <- .id <- level <- std.error <- NULL

    ss <- stats::coef(summary(x))
    ss <- ss[!sapply(ss,is.null)]
    ## FIXME: warn if !missing(component) and component includes
    ##  NULL terms
    component <- intersect(component,names(ss))
    if (length(component[!component %in% c("cond", "zi")]) > 0L) {
        stop("only works for conditional and (partly for) zero-inflation components")
    }
    effect_names <- c("ran_pars", "fixed", "ran_modes")
    if (!is.null(scales)) {
        if (length(scales) != length(effects)) {
            stop("if scales are specified, values (or NA) must be provided ",
                 "for each effect")
        }
    }
    if (length(miss <- setdiff(effects,effect_names))>0)
        stop("unknown effect type ",miss)
    base_nn <- c("estimate", "std.error", "statistic", "p.value")
    ret <- list()
    ret_list <- list()
    if ("fixed" %in% effects) {
        # return tidied fixed effects rather than random
        ret <- lapply(ss,
                      function(x) {
            x %>% as.data.frame %>%
                setNames(c("estimate", "std.error", "statistic", "p.value")) %>%
                tibble::rownames_to_column("term")
        })
        # p-values may or may not be included
        # HACK: use the columns from the conditional component, preserving previous behaviour
        nn <- base_nn[1:ncol(ret$cond)]
        if (conf.int) {
            for (comp in component) {
                cifix <- confint(x,method=tolower(conf.method),
                                 component=comp,
                                 estimate=FALSE,
                                 ## conditional/zi components
                                 ## include random-effect parameters
                                 ## as well, don't want those right now ...
                                 parm=seq(nrow(ret[[comp]])),...) %>%
                    as.data.frame %>%
                    setNames(c("conf.low","conf.high"))
                ret[[comp]] <- bind_cols(
                    tibble::rownames_to_column(ret[[comp]],var="term"),
                    cifix)
            }
            nn <- c(nn,"conf.low","conf.high")
        }
        if ("ran_pars" %in% effects || "ran_modes" %in% effects) {
            ret <- lapply(ret, function(x) data.frame(x,group="fixed"))
            nn <- c(nn,"group")
        }
        ret_list$fixed <- plyr::ldply(ret, .fun = reorder_frame,
                                      .id = "component")
    }
    if ("ran_pars" %in% effects &&
        !all(sapply(VarCorr(x),is.null))) {
        if (is.null(scales)) {
            rscale <- "sdcor"
        } else rscale <- scales[effects=="ran_pars"]
        if (!rscale %in% c("sdcor","vcov"))
            stop(sprintf("unrecognized ran_pars scale %s",sQuote(rscale)))
        ## kluge for now ...
        vv <- list()
        if ("cond" %in% component){
          vv$cond <- VarCorr(x)[["cond"]]
          class(vv$cond) <- "VarCorr.merMod"
        }
        if ("zi" %in% component) {
          if (!is.null(vv$zi <- VarCorr(x)[["zi"]])) {
              class(vv$zi) <- "VarCorr.merMod"
          }
        }
    
        ret <- plyr::ldply(vv, as.data.frame, .id = "component")
        ret[] <- lapply(ret, function(x) if (is.factor(x))
                                                 as.character(x) else x)
        if (is.null(ran_prefix)) {
            ran_prefix <- switch(rscale,
                                 vcov=c("var","cov"),
                                 sdcor=c("sd","cor"))
        }
        pfun <- function(x) {
            v <- na.omit(unlist(x))
            if (length(v)==0) v <- "Observation"
            p <- paste(v,collapse=".")
            if (!identical(ran_prefix,NA)) {
                p <- paste(ran_prefix[length(v)],p,sep="_")
            }
            return(p)
        }

        ## DRY! try to refactor glmmTMB/lme4 tidiers
        
        ## don't try to assign as rowname (non-unique anyway),
        ## make it directly into a term column
        ret[["term"]] <- apply(ret[c("var1","var2")],1,pfun)

        ## keep only desired term, rename
        ret <- setNames(ret[c("grp","term",rscale)],
                        c("group","term","estimate"))

        ## rownames(ret) <- seq(nrow(ret))

        if (conf.int) {
            ciran <- confint(x,parm="theta_",method=conf.method,...)
            ret <- data.frame(ret,ciran)
            nn <- c(nn,"conf.low","conf.high")
        }
        ret_list$ran_pars <- ret
    }
    
    if ("ran_modes" %in% effects) {
        ## fix each group to be a tidy data frame

        nn <- c("estimate", "std.error")
        re <- ranef(x,condVar=TRUE)
        getSE <- function(x) {
            v <- attr(x,"postVar")
            setNames(as.data.frame(sqrt(t(apply(v,3,diag)))),
                     colnames(x))
        }
        fix <- function(g,re,.id) {
             newg <- broom::fix_data_frame(g, newnames = colnames(g), newcol = "level")
             # fix_data_frame doesn't create a new column if rownames are numeric,
             # which doesn't suit our purposes
             newg$level <- rownames(g)
             newg$type <- "estimate"

             newg.se <- getSE(re)
             newg.se$level <- rownames(re)
             newg.se$type <- "std.error"

             data.frame(rbind(newg,newg.se),.id=.id,
                        check.names=FALSE)
                        ## prevent coercion of variable names
        }

        mm <- do.call(rbind,Map(fix,coef(x),re,names(re)))

        ## block false-positive warnings due to NSE
        type <- spread <- est <- NULL
        mm %>% gather(term, estimate, -.id, -level, -type) %>%
            spread(type,estimate) -> ret

        ## FIXME: doesn't include uncertainty of population-level estimate

        if (conf.int) {
            if (conf.method != "Wald")
                stop("only Wald CIs available for conditional modes")

            mult <- qnorm((1+conf.level)/2)
            ret <- transform(ret,
                             conf.low=estimate-mult*std.error,
                             conf.high=estimate+mult*std.error)
        }

        ret <- dplyr::rename(ret,grp=.id)
        ret_list$ran_modes <- ret
    }
    ## use ldply to get 'effect' added as a column
    return(plyr::ldply(ret_list,identity,.id="effect"))

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
augment.glmmTMB <- function(x, data = stats::model.frame(x), newdata,
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
    finish_glance(x=x)
}
