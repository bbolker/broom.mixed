## test tidy and glance methods from mcmc_tidiers
stopifnot(require("testthat"), require("broom.mixed"))

context("mcmc tidiers")

if (suppressPackageStartupMessages(require(rstan, quietly = TRUE)))  {

  test_that("tidy returns indexes if requested on rstanarm fits", {
  
      # Make sure that (inst/)extdata/run_examples.R was run to generate rds
    rstan_example <- readRDS(system.file("extdata", "rstan_example.rds", package = "broom.mixed"))
    td <- broom::tidy(rstan_example) # For S3 only, avoid importing broom full namespace
    # check_tidy from helper-checkers
    check_tidy(td, 18, 3, c("term", "estimate", "std.error"))
    td_mcmc <- broom::tidyMCMC(rstan_example)
    expect_equal(td, td_mcmc)
    
    td <- tidyMCMC(rstan_example, index = TRUE)
    check_tidy(td, 18, 5, c("term", "term0", "index", "estimate", "std.error"))

    td <- tidyMCMC(rstan_example, drop.pars = NULL)
    expect_equal(td[19, "term"], "lp__")
    
    td <- tidyMCMC(rstan_example, conf.int = TRUE)
    check_tidy(td, 18, 5, c("term", "estimate", "std.error", "conf.low", "conf.high"))

    td <- tidyMCMC(rstan_example, rhat = TRUE)
    check_tidy(td, 18, 4, c("term", "estimate", "std.error", "rhat"))

    td <- tidyMCMC(rstan_example, ess = TRUE)
    check_tidy(td, 18, 4, c("term", "estimate", "std.error", "ess"))
  })
}
