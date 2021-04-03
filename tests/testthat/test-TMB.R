library(broom.mixed)
library(testthat)
if (require("TMB")) {
    if (identical(Sys.getenv("NOT_CRAN"), "true")) {
        cc <- capture.output(runExample("simple",thisR=TRUE))
        class(obj) <- "TMB"
        test_that("basic tidying", {
            t0 <- tidy(obj)
            expect_equal(dim(t0),c(4,4))
            expect_equal(names(t0), c("type", "term", "estimate", "std.error"))
        })
        test_that("confint", {          
            t1 <- tidy(obj,conf.int=TRUE,conf.method="wald")
            t2 <- tidy(obj,conf.int=TRUE,conf.method="uniroot")
            ## CIs *not* identical, but should be close ...
            expect_false(identical(t1$conf.low, t2$conf.low))
            expect_equal(t1$conf.low, t2$conf.low, tolerance=1e-2)
            ## t3 <- tidy(obj,conf.int=TRUE,conf.method="profile")
        })
    } ## not on CRAN
} ## if require(TMB)
        
