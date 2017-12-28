#' Tidying methods for mixed effects models
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
#' @name lme4_tidiers
#'
#' @examples
#' 
#' if (require("lme4")) {
#'     # example regressions are from lme4 documentation
#'     lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'     tidy(lmm1)
#'     tidy(lmm1, effects = "fixed")
#'     tidy(lmm1, effects = "fixed", conf.int=TRUE)
#'     tidy(lmm1, effects = "fixed", conf.int=TRUE, conf.method="profile")
#'     ## pp <- profile(lmm1)
#'     lmm1_prof <- readRDS(system.file("example_data","lmm1_prof.rds",
#'                                      package="broom.mixed"))
#'     tidy(lmm1, conf.int=TRUE, conf.method="profile", profile=lmm1_prof)
#'     ## conditional modes (group-level deviations from population-level estimate)
#'     tidy(lmm1, effects = "ran_modes", conf.int=TRUE)
#'     ## coefficients (group-level estimates)
#'     tidy(lmm1, effects = "ran_coefs")
#'     head(augment(lmm1, sleepstudy))
#'     glance(lmm1)
#'     
#'     glmm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#'                   data = cbpp, family = binomial)
#'     tidy(glmm1)
#'     tidy(glmm1, effects = "fixed")
#'     head(augment(glmm1, cbpp))
#'     glance(glmm1)
#'     
#'     startvec <- c(Asym = 200, xmid = 725, scal = 350)
#'     nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
#'                   Orange, start = startvec)
#'     tidy(nm1)
#'     tidy(nm1, effects = "fixed")
#'     head(augment(nm1, Orange))
#'     glance(nm1)
#' }
NULL

confint.rlmerMod <- function(x, parm,
                             method="Wald", ...) {
    cc <- class(x)
    class(x) <- "merMod"
    if (method!="Wald") {
        warning("only Wald method implemented for rlmerMod objects")
    }
    return(confint(x,parm=parm, method="Wald", ...))
}
    
#' @rdname lme4_tidiers
#'
#' @param debug print debugging output?
#' @param effects A character vector including one or more of "fixed" (fixed-effect parameters); "ran_pars" (variances and covariances or standard deviations and correlations of random effect terms); "ran_modes" (conditional modes/BLUPs/latent variable estimates); or "ran_coefs" (predicted parameter values for each group, as returned by \code{\link[lme4]{coef.merMod}})
#' @param conf.int whether to include a confidence interval
#' @param conf.level confidence level for CI
#' @param conf.method method for computing confidence intervals (see \code{lme4::confint.merMod})
#' @param scales scales on which to report the variables: for random effects, the choices are \sQuote{"sdcor"} (standard deviations and correlations: the default if \code{scales} is \code{NULL}) or \sQuote{"vcov"} (variances and covariances). \code{NA} means no transformation, appropriate e.g. for fixed effects; inverse-link transformations (exponentiation or logistic) are not yet implemented, but may be in the future.
#' @param ran_prefix a length-2 character vector specifying the strings to use as prefixes for self- (variance/standard deviation) and cross- (covariance/correlation) random effects terms
#' @param profile pre-computed profile object, for speed when using \code{conf.method="profile"}
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
#' @importFrom plyr ldply 
#' @importFrom dplyr mutate bind_rows data_frame
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather spread
#' @importFrom purrr map
#' @importFrom nlme VarCorr ranef
#' @importFrom methods is
#' @importFrom broom fix_data_frame
#' 
#' @export
tidy.merMod <- function(x, effects = c("ran_pars","fixed"),
                        scales = NULL, ## c("sdcor","vcov",NA),
                        ran_prefix=NULL,
                        conf.int = FALSE,
                        conf.level = 0.95,
                        conf.method = "Wald",
                        profile = NULL,
                        debug=FALSE,
                        ...) {

    ## R CMD check false positives
    term <- estimate <- .id <- level <- std.error <- NULL

    effect_names <- c("ran_pars", "fixed", "ran_modes", "ran_coefs")
    if (!is.null(scales)) {
        if (length(scales) != length(effects)) {
            stop("if scales are specified, values (or NA) must be provided ",
                 "for each effect")
        }
    }
    if (length(miss <- setdiff(effects,effect_names))>0)
        stop("unknown effect type ",miss)
    base_nn <- c("estimate", "std.error", "statistic", "p.value")
    ret_list <- list()
    if ("fixed" %in% effects) {
        # return tidied fixed effects rather than random
        ss <- summary(x)
        ret <- stats::coef(ss)
        if (debug) {
            cat("output from coef(summary(x)):\n")
            print(coef(ss))
        }
        if (is(x,"merModLmerTest")) {
            ret <- ret[,!colnames(ret) %in% "df"]
        }            

        # p-values may or may not be included
        nn <- base_nn[1:ncol(ret)]

        if (conf.int && conf.method=="profile" && !is.null(profile)) {
            p <- profile
        } else p <- x

        if (conf.int) {
            if (is(x,"merMod") || is(x,"rlmerMod")) {
                cifix <- confint(p,parm="beta_",method=conf.method,...)
            } else {
                ## ?? for glmmTMB?  check ...
                cifix <- confint(p,...)
            }
            ret <- data.frame(ret,cifix,stringsAsFactors=FALSE)
            nn <- c(nn,"conf.low","conf.high")
        }
        if ("ran_pars" %in% effects || "ran_modes" %in% effects) {
            ret <- data.frame(ret,group="fixed",stringsAsFactors=FALSE)
            nn <- c(nn,"group")
        }
        ret_list$fixed <-
            fix_data_frame(ret, newnames = nn)
    }
    if ("ran_pars" %in% effects) {
        if (is.null(scales)) {
            rscale <- "sdcor"
        } else rscale <- scales[effects=="ran_pars"]
        if (!rscale %in% c("sdcor","vcov"))
            stop(sprintf("unrecognized ran_pars scale %s",sQuote(rscale)))
        vc <- VarCorr(x)
        if (!is(x,"merMod") && grepl("^VarCorr",class(vc)[1])) {
            if (!is(x,"rlmerMod")) {
                ## hack: attempt to augment glmmADMB (or other)
                ##   values so we can use as.data.frame.VarCorr.merMod
                vc <- lapply(vc,
                             function(x) {
                    attr(x,"stddev") <- sqrt(diag(x))
                    attr(x,"correlation") <- stats::cov2cor(x)
                    x
                })
                attr(vc,"useScale") <- (stats::family(x)$family=="gaussian")
            }
            class(vc) <- "VarCorr.merMod"
        }
        ret <- as.data.frame(vc)
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

        rownames(ret) <- paste(apply(ret[c("var1","var2")],1,pfun),
                                ret[,"grp"], sep = ".")

        ## keep only desired term, rename
        ret <- setNames(ret[c("grp",rscale)],
                        c("group","estimate"))

        if (conf.int) {
            ciran <- confint(p,parm="theta_",method=conf.method,...) %>%
                as.data.frame %>%
                setNames(c("conf.low","conf.high"))
            ret <- data.frame(ret,ciran)
        }
        ret_list$ran_pars <- ret %>% tibble::rownames_to_column("term")
    }

    
    getSE <- function(x) {
        v <- attr(x,"postVar")
        re_sd <- sqrt(t(apply(v,3,diag)))
        ## for single random effect term, need to re-transpose
        if (nrow(re_sd)==1) re_sd <- t(re_sd)
        r <- setNames(data.frame(re_sd,row.names=rownames(x)),
                      colnames(x))
        return(r)
    }
    fix_ran_modes <- function(g, do_SE=FALSE) {
        r <- g %>%
            tibble::rownames_to_column("level") %>%
            gather(term,estimate,-level)
g
        if (do_SE) {
            newg.se <- getSE(g) %>%
                tibble::rownames_to_column("level") %>%
                gather(term,std.error,-level)

            r <- full_join(r,newg.se,by=c("level","term"))
        }
        return(r)
    }

    if ("ran_modes" %in% effects) {
        ## fix each group to be a tidy data frame

        ret <- ranef(x,condVar=TRUE)  %>%
            purrr::map(fix_ran_modes, do_SE=TRUE) %>%
            bind_rows(.id="group")

        if (conf.int) {
            if (conf.method != "Wald")
                stop("only Wald CIs available for conditional modes")

            mult <- qnorm((1+conf.level)/2)
            ret <- ret %>% mutate(
                               conf.low=estimate-mult*std.error,
                               conf.high=estimate+mult*std.error)
        }

        ret_list$ran_modes <- ret
    }
    if ("ran_coefs" %in% effects) {
        ret <- coef(x) %>%
            purrr::map(fix_ran_modes, do_SE=FALSE) %>%
            bind_rows(.id="group")
        
        if (conf.int) {
            warning("CIs not available for random-effects coefficients: returning NA")
            ret <- ret %>% mutate(conf.low=NA,conf.high=NA)
        }


        ret_list$ran_coefs <- ret
    }

    return(dplyr::bind_rows(ret_list,.id="effect"))
}

#' @rdname lme4_tidiers
#' @export
tidy.rlmerMod <- tidy.merMod

#' @rdname lme4_tidiers
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
#' Also added for "merMod" objects, but not for "mer" objects,
#' are values from the response object within the model (of type
#' \code{lmResp}, \code{glmResp}, \code{nlsResp}, etc). These include \code{".mu",
#' ".offset", ".sqrtXwt", ".sqrtrwt", ".eta"}.
#'
#' @importFrom broom augment augment_columns
#' @export
augment.merMod <- function(x, data = stats::model.frame(x), newdata, ...) {    
    # move rownames if necessary
    if (missing(newdata)) {
        newdata <- NULL
    }
    ret <- augment_columns(x, data, newdata, se.fit = NULL, ...)
    
    # add predictions with no random effects (population means)
    predictions <- stats::predict(x, re.form = NA)
    # some cases, such as values returned from nlmer, return more than one
    # prediction per observation. Not clear how those cases would be tidied
    if (length(predictions) == nrow(ret)) {
        ret$.fixed <- predictions
    }

    # columns to extract from resp reference object
    # these include relevant ones that could be present in lmResp, glmResp,
    # or nlsResp objects

    respCols <- c("mu", "offset", "sqrtXwt", "sqrtrwt", "weights",
                  "wtres", "gam", "eta")
    cols <- lapply(respCols, function(n) x@resp[[n]])
    names(cols) <- paste0(".", respCols)
    cols <- as.data.frame(compact(cols))  # remove missing fields
    
    cols <- insert_NAs(cols, ret)
    if (length(cols) > 0) {
        ret <- cbind(ret, cols)
    }

    return(unrowname(ret))
}


#' @rdname lme4_tidiers
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
#' @importFrom broom glance finish_glance
#' @export
glance.merMod <- function(x, ...) {
    # We cannot use stats::sigma or lme4::sigma here, even in an
    # if statement, since that leads to R CMD CHECK warnings on 3.2
    # or dev R, respectively
    ret <- unrowname(data.frame(sigma = sigma(x)))
    broom::finish_glance(ret, x)
}

##' Augmentation for random effects (for caterpillar plots etc.)
##' 
##' @param x ranef (conditional mode) information from an lme4 fit, using \code{ranef(.,condVar=TRUE)}
##' @param ci.level level for confidence intervals
##' @param reorder reorder levels by conditional mode values?
##' @param order.var numeric or character: which variable to use for ordering levels?
##' @param \dots additional arguments (unused: for generic consistency)
##' @importFrom reshape2 melt
##' @importFrom plyr ldply
##' @examples
##' if (require("lme4")) {
##'    fit <- lmer(Reaction~Days+(Days|Subject),sleepstudy)
##'    rr <- ranef(fit,condVar=TRUE)
##'    aa <- broom::augment(rr)
##'    ## Q-Q plot:
##'    if (require(ggplot2) && require(dplyr)) {
##'       g0 <- ggplot(aa,aes(estimate,qq,xmin=lb,xmax=ub))+
##'           geom_errorbarh(height=0)+
##'           geom_point()+facet_wrap(~variable,scale="free_x")
##'       ## regular caterpillar plot:
##'       g1 <- ggplot(aa,aes(estimate,level,xmin=lb,xmax=ub))+
##'          geom_errorbarh(height=0)+
##'          geom_vline(xintercept=0,lty=2)+
##'          geom_point()+facet_wrap(~variable,scale="free_x")
##'       ## emphasize extreme values
##'       aa2 <- aa  %>% group_by(grp,level) %>%
##'             mutate(keep=any(estimate/std.error>2))
##'        ## Update caterpillar plot with extreme levels highlighted
##'        ##  (highlight all groups with *either* extreme intercept *or*
##'        ##   extreme slope)
##'       ggplot(aa2,aes(estimate,level,xmin=lb,xmax=ub,colour=factor(keep)))+
##'          geom_errorbarh(height=0)+
##'          geom_vline(xintercept=0,lty=2)+
##'          geom_point()+facet_wrap(~variable,scale="free_x")+
##'          scale_colour_manual(values=c("black","red"), guide=FALSE)
##'    }
##' }
##' @importFrom stats ppoints
##' @export 
augment.ranef.mer <- function(x,
                                 ci.level=0.9,
                                 reorder=TRUE,
                                 order.var=1, ...) {

    variable <- level <- estimate <- std.error <- NULL
    tmpf <- function(z) {
        if (is.character(order.var) && !order.var %in% names(z)) {
            order.var <- 1
            warning("order.var not found, resetting to 1")
        }
        ## would use plyr::name_rows, but want levels first
        zz <- data.frame(level=rownames(z),z,check.names=FALSE)
        if (reorder) {
            ## if numeric order var, add 1 to account for level column
            ov <- if (is.numeric(order.var)) order.var+1 else order.var
            zz$level <- reorder(zz$level, zz[,order.var+1], FUN=identity)
        }
        ## Q-Q values, for each column separately
        qq <- c(apply(z,2,function(y) {
                  qnorm(stats::ppoints(nrow(z)))[order(order(y))]
              }))
        rownames(zz) <- NULL
        pv   <- attr(z, "postVar")
        cols <- 1:(dim(pv)[1])
        se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
        ## n.b.: depends on explicit column-major ordering of se/melt
        zzz <- cbind(melt(zz,id.vars="level",value.name="estimate"),
                     qq=qq,std.error=se)
        ## reorder columns:
        subset(zzz,select=c(variable, level, estimate, qq, std.error))
    }
    dd <- ldply(x,tmpf,.id="grp")
    ci.val <- -qnorm((1-ci.level)/2)
    transform(dd,
              ## p=2*pnorm(-abs(estimate/std.error)), ## 2-tailed p-val
              lb=estimate-ci.val*std.error,
              ub=estimate+ci.val*std.error)
}