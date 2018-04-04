library(TMB)
## devtools::install_github("bbolker/broom.mixed")
runExample("simple",thisR=TRUE)
sdreport(obj)
tidy.TMB <- function(object,effect=c("fixed","random"),se=TRUE,
                     conf.int=FALSE,
                     conf.method=c("wald","uniroot","profile")) {
    sdr <- sdreport(object)
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
                ## FIXME: allow alpha spec
                qval <- qnorm(0.975)
                ss <- mutate(ss,conf.low
            }
    }

    }
    retlist$fixed <- ss
    ret <- dplyr::bind_rows(retlist,.id="type")
    ## FIXME: add statistic (mean/stderr)
    ##   and p-value (2*pnorm(abs(statistic),lower.tail=FALSE))
    ##   if requested
    ##
    ## ## FIXME: get confidence intervals if requested
    ##  Wald: find qnorm(0.975)*c(-1,1)*std.error+estimate
    ##  uniroot: call glmmTMB:::tmbroot
    ##  profile: call TMB:::confint.TMB(TMB::tmbprofile(object,params))
    return(ret)
}

class(object) <- "TMB"
tidy(object)
