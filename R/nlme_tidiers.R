#' Tidying methods for mixed effects models
#' 
#' These methods tidy the coefficients of mixed effects models
#' of the \code{lme} class from functions  of the \code{nlme} package.
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
#' if (require("broom") && require("nlme") & require("lme4")) {
#'     # example regressions are from lme4 documentation
#'     lmm1 <- lme(Reaction ~ Days, random=~ Days|Subject, sleepstudy)
#'     tidy(lmm1)
#'     tidy(lmm1, effects = "fixed")
#'     tidy(lmm1, conf.int = TRUE)
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
#' }
#' 
#' @rdname nlme_tidiers
#' 
#' @param effects Either "random" (default) or "fixed"
#' 
#' @return \code{tidy} returns one row for each estimated effect, either
#' random or fixed depending on the \code{effects} parameter. If
#' \code{effects = "random"}, it contains the columns
#'   \item{group}{the group within which the random effect is being estimated}
#'   \item{level}{level within group}
#'   \item{term}{term being estimated}
#'   \item{estimate}{estimated coefficient}
#' 
#' If \code{effects="fixed"}, \code{tidy} returns the columns
#'   \item{term}{fixed term being estimated}
#'   \item{estimate}{estimate of fixed effect}
#'   \item{std.error}{standard error}
#'   \item{statistic}{t-statistic}
#'   \item{p.value}{P-value computed from t-statistic}
#' 
#' @importFrom plyr ldply
#' @importFrom nlme getVarCov intervals
#' @import dplyr
## @importFrom dplyr data_frame select full_join
#' 
#' @export
tidy.lme <- function(x, effects = c("ran_pars", "fixed"),
                     scales = NULL,
                     conf.int = FALSE,
                     conf.level = 0.95,
                     ...) {

    ## R CMD global var check
    lower <- upper <- NULL
    
    effect_names <- c("ran_pars", "fixed", "ran_modes", "ran_coefs")
    if (length(miss <- setdiff(effects,effect_names))>0)
        stop("unknown effect type ",miss)

    ## R CMD check false positives
    term <- estimate <- .id <- level <- std.error <- NULL

    if (!is.null(scales)) {
        if (length(scales) != length(effects)) {
            stop("if scales are specified, values (or NA) must be provided ",
                 "for each effect")
        }
    }

    ret_list <- list()
    
    if ("fixed" %in% effects) {
        # return tidied fixed effects
        ret <- summary(x)[["tTable"]] %>% data.frame(check.names=FALSE) %>%
                        rename_regex_match() %>%
                            tibble::rownames_to_column("term")
        if (conf.int) {
             cifix <- intervals(x,which="fixed")[["fixed"]] %>%
                data.frame() %>%
                dplyr::select(lower,upper) %>%
                setNames(c("conf.low","conf.high")) %>%
                tibble::rownames_to_column("term")
             ret <- dplyr::full_join(ret,cifix,by="term")
        }

        ran_effs <- sprintf("ran_%s",c("pars","modes","coefs"))
        if (any(purrr::map_lgl(ran_effs, ~. %in% effects))) {
            ## add group="fixed" to tidy table for fixed effects
            ret <- mutate(ret,effect="fixed",group="fixed")
        }

        ret_list$fixed <- ret %>%
            reorder_frame()
    }

    
    if ("ran_pars" %in% effects) {
        if (is.null(scales)) {
            rscale <- "sdcor"
        } else rscale <- scales[effects=="ran_pars"]
        if (!rscale %in% c("sdcor","vcov"))
            stop(sprintf("unrecognized ran_pars scale %s",sQuote(rscale)))
        grplen <- attr(x$modelStruct$reStruct,"plen")
        multilevel <- (length(grplen)>1)
        nonlin <- inherits(x,"nlme")

        if (multilevel || nonlin) {
            warning("ran_pars not yet implemented for ",
                    if (multilevel) {
                        "multiple levels of nesting"
                    } else {
                        "nonlinear models"
                    })
            ret <- dplyr::data_frame()
        } else {
            vc <- nlme::getVarCov(x)
            ran_prefix <- switch(rscale,
                                 vcov=c("var","cov"),
                                 sdcor=c("sd","cor"))
            ## construct appropriate sd/cor names from dimnames
            nmvec <- outer(colnames(vc),rownames(vc),
                           function(x,y) { ifelse(x==y,
                                      sprintf("%s_%s",
                                        ran_prefix[1],
                                        x),
                                sprintf("%s_%s.%s",
                                        ran_prefix[2],
                                        x,y))
                                })
            lwrtri <- function(x) { x[lower.tri(x,diag=TRUE)] }
            nmvec <- c(lwrtri(nmvec),
                       sprintf("%s_%s",ran_prefix[1],"Observation"))
            grpnames <- c(rep(names(grplen),grplen),"Residual")
            if (rscale=="vcov") {
                vals <- c(lwrtri(vc),sigma(x)^2)
            } else {
                vals <- cov2cor(vc)
                diag(vals) <- sqrt(diag(vc))
                vals <- c(lwrtri(vals),sigma(x))
            }
            ret <- dplyr::data_frame(effect="ran_pars",
                                     group=grpnames,
                                     term=c(nmvec),
                                     estimate=c(vals))
            
            if (conf.int) {
                ii <- intervals(x,which="var-cov")$reStruct
                trfun <- function(z) {
                    nm <- rownames(z)
                    nm <- sub("\\(","_",
                              sub(")$","",
                                  sub(",",".",nm)))
                    ## ugh, swap order of vars in cor term
                    corterms <- grepl("^cor",nm)
                    re <- "_([^\\.]+)\\.(.+)$"
                    nm[corterms] <-
                        sub(re,"_\\2.\\1",nm[corterms],
                            perl=TRUE)
                    return(dplyr::data_frame(term=nm,conf.low=z[,"lower"],
                                             conf.high=z[,"upper"]))
                }
                ci <- dplyr::bind_rows(lapply(ii,trfun),.id="group")
                if (rscale!="sdcor") {
                    ## FIXME: transform/recompute from scratch
                    warning("confidence intervals for ran pars only available on sdcor scale")
                    ci$conf.low <- ci$conf.high <- NA
                }
                ## FIXME: also do confint on residual
                ret <- dplyr::full_join(ret,ci,by=c("group","term"))
            }
        } ## if not multi-level model
        ret_list$ran_pars <- ret
    }
    
    if ("ran_modes" %in% effects) {

        ret_list$ran_modes <-
            ranef(x) %>%
            tibble::rownames_to_column("level") %>%
            tidyr::gather(key=term,value=estimate,-level)
        ## FIXME: group?
    }
    if ("ran_coefs" %in% effects) {
        ret_list$ran_coefs <-
            stats::coef(x) %>%
            tibble::rownames_to_column("level") %>%
            tidyr::gather(key=term,value=estimate,-level)
        ## FIXME: group?
    }

    ret <-  bind_rows(ret_list)
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
    if (is.null(data)){
        stop("augment.lme must be called with an explicit 'data' argument for this (n)lme fit  because of an inconsistency in nlme.")
    }
    # move rownames if necessary
    if (missing(newdata)) {
        newdata <- NULL
    }
    ret <- augment_columns(x, data, newdata, se.fit = NULL)
    
    # add predictions with no random effects (population means)
    predictions <- stats::predict(x, level=0)
    if (length(predictions) == nrow(ret)) {
        ret$.fixed <- predictions
    }
    
    unrowname(ret)
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
#' \url{https://stat.ethz.ch/pipermail/r-help/2006-May/104744.html}
#'  McCullagh & Nelder (1989) would be the authorative [sic] reference, but the 1982
#' first edition manages to use 'deviance' in three separate senses on one
#' page. }
#'
#' @export
glance.lme <- function(x) {
    finish_glance(x=x)
}

#' @rdname nlme_tidiers
#' @export
tidy.gls <- function(x, 
                     conf.int = FALSE,
                     conf.level = 0.95,
                     ...) {
    . <- Value <- Std.Error <- `t-value` <- `p-value` <- NULL ## glob var checks
    summary(x)[["tTable"]] %>%
        as.data.frame() %>% ## have to convert to df *first*
        tibble::rownames_to_column(var="term") %>%
        as_tibble() %>%
        dplyr::rename(estimate=Value,
                      std.error=Std.Error,
                      statistic=`t-value`,
                      p.value=`p-value`)
        
}

#' @export
glance.gls <- function(x) {
    ss <- summary(x)
    with(ss,
         tibble::tibble(sigma,
                        df=dims[["p"]],
                        logLik,
                        AIC,
                        BIC,
                        df.residual=dims[["N"]]-dims[["p"]])
         )
}
