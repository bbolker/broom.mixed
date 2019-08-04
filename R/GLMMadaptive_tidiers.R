tidy.MixMod <- function(x, effects = c("fixed"), conf.int=TRUE,
                        conf.level=0.95) {
    estimate <- std.err <- NULL ## R CMD check/NSE
    effects <- match.arg(effects)
    ret <- (coef(summary(x))
        %>% as.data.frame()
        %>% tibble::rownames_to_column("term")
        %>% rename(estimate="Estimate",
                   std.error="Std.Err",
                   statistic="z-value",
                   p.value="p-value")
    )
    if (conf.int) {
        qq <- qnorm((1+c(-1,1)*conf.level)/2)
        ret <- (ret
            %>% mutate(conf.low=estimate+qq[1]*std.error,
                       conf.high=estimate+qq[2]*std.error)
        )
    }
    return(ret)
}
