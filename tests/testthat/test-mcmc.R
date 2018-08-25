## test tidy and glance methods from mcmc_tidiers
stopifnot(require("testthat"), require("broom.mixed"))

context("mcmc tidiers")

## HACK: need to find the right generic
tidy <- broom.mixed:::tidy.mcmc

if (suppressPackageStartupMessages(require(rstan, quietly = TRUE))) {
  test_that("tidy returns indexes if requested on rstanarm fits", {

    # Make sure that (inst/)extdata/run_examples.R was run to generate rds
    rstan_example <- readRDS(system.file("extdata", "rstan_example.rds", package = "broom.mixed"))
    # check_tidy from helper-checkers
    td <- tidy(rstan_example)
    check_tidy(td, 18, 3, c("term", "estimate", "std.error"))

    td <- tidy(rstan_example, index = TRUE)
    check_tidy(td, 18, 4, c("term", "index", "estimate", "std.error"))

    td <- tidy(rstan_example, drop.pars = NULL)
    expect_equal(td[19, ][["term"]], "lp__")

    td <- tidy(rstan_example, conf.int = TRUE)
    check_tidy(td, 18, 5, c("term", "estimate", "std.error", "conf.low", "conf.high"))

    td <- tidy(rstan_example, rhat = TRUE)
    check_tidy(td, 18, 4, c("term", "estimate", "std.error", "rhat"))

    td <- tidy(rstan_example, ess = TRUE)
    check_tidy(td, 18, 4, c("term", "estimate", "std.error", "ess"))
  })
}
