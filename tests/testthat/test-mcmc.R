## test tidy and glance methods from mcmc_tidiers
stopifnot(require("testthat"), require("broom.mixed"))

context("mcmc tidiers")

if (suppressPackageStartupMessages(require(rstan, quietly = TRUE)))  {

  test_that("tidy returns indexes if requested on rstanarm fits", {
  
      # Make sure that (inst/)extdata/run_examples.R was run to generate rds
    rstan_example <- readRDS(system.file("extdata", "rstan_example.rds", package = "broom.mixed"))
    td <- broom::tidy(rstan_example) # For S3 only, avoid importing broom full namespace
    expect_equal(names(td), c("term", "estimate", "std.error"))
    expect_equal(nrow(td), 18)
    td_mcmc <- broom::tidyMCMC(rstan_example)
    expect_equal(td, td_mcmc)
    
    td <- tidyMCMC(rstan_example, index = TRUE)
    expect_equal(names(td),
                 c("term", "term0", "index", "estimate", "std.error"))
    expect_equal(nrow(td), 18)
    td <- tidyMCMC(rstan_example, drop.pars = NULL)
    expect_equal(td[19, "term"], "lp__")
    td <- tidyMCMC(rstan_example, conf.int = TRUE)
    expect_equal(names(td),
                 c("term", "estimate", "std.error", "conf.low", "conf.high"))
    expect_equal(nrow(td), 18)
    td <- tidyMCMC(rstan_example, rhat = TRUE)
    expect_equal(names(td), c("term", "estimate", "std.error", "rhat"))
    expect_equal(nrow(td), 18)
    td <- tidyMCMC(rstan_example, ess = TRUE)
    expect_equal(names(td), c("term", "estimate", "std.error", "ess"))
    expect_equal(nrow(td), 18)
  })
}
