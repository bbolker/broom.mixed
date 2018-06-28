##' Tidying methods for TMB models
##'
##' @param x An object of class \code{TMB} (you may need to use
##' \code{class(obj) <- "TMB"} on your results from TMB)
##' @inheritParams glmmTMB_tidiers
##' @param effect which effects should be returned?
##' @param conf.method method for computing confidence intervals
##' @importFrom TMB sdreport
##' @importFrom glmmTMB tmbroot
##' @examples
##' \dontrun{
##' runExample("simple",thisR=TRUE)
##' class(obj) <- "TMB"
##' tidy(obj,conf.int=TRUE,conf.method="wald")
##' tidy(obj,conf.int=TRUE,conf.method="uniroot")
##' }
##' @export
tidy.TMB <- function(x,effect=c("fixed","random"),
                     conf.int=FALSE,
                     conf.level = 0.95,
                     conf.method=c("wald","uniroot","profile")) {
    ## R CMD check/global variables
    Estimate <- estimate <- std.error <- NULL
    sdr <- sdreport(x)
    retlist <- list()
    if ("fixed" %in% effect) {
        ss <- summary(sdr,select="fixed") %>%
            as.data.frame %>%
            ## FIXME: disambiguate repeated variable names better
            ## i.e. beta.1, beta.2 instead of beta, beta.1 ...
            tibble::rownames_to_column("term") %>%
            rename(estimate=Estimate,std.error="Std. Error")
        if (conf.int) {
            if (tolower(conf.method=="wald")) {
                qval <- qnorm((1+conf.level)/2)
                ss <- mutate(ss,
                             conf.low=estimate-qval*std.error,
                             conf.high=estimate+qval*std.error)
            } else if (conf.method=="uniroot") {
                ## FIXME: allow parm specs
                ## FIXME: avoid calling sdreport again inside tmbroot?
                ## do.call(rbind,...) because bind_rows needs named list
                tt <- do.call(rbind,
                              lapply(seq(nrow(ss)),
                                     glmmTMB::tmbroot,obj=x))
                ss$conf.low <- tt[,"lwr"]
                ss$conf.high <- tt[,"upr"]
            } else {
                ##  profile:
                ##  call TMB:::confint.TMB(TMB::tmbprofile(x,params))
                stop(sprintf("conf.method=%s not implemented",conf.method))
            }
        }
    }
    retlist$fixed <- ss
    ret <- dplyr::bind_rows(retlist,.id="type")
    ## FIXME: add statistic (mean/stderr)
    ##   and p-value (2*pnorm(abs(statistic),lower.tail=FALSE))
    ##   if requested
    ##

    return(ret)
}

