stopifnot(require("testthat"), require("broom.mixed"))


if (require(glmmTMB, quietly = TRUE)) {
    load(system.file("extdata","glmmTMB_example.rda",package="broom.mixed",
                     mustWork=TRUE))

    context("glmmTMB models")
    
    test_that("components included for zi models", {
        td <- tidy(zipm3)
        check_tidy(td, 29, 8,
                   c("effect", "component", "group", "term",
                     "estimate", "std.error", "statistic", "p.value"))
    })
}
