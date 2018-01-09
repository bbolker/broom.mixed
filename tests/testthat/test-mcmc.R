## test tidy and glance methods from mcmc_tidiers
stopifnot(require("testthat"), require("broom.mixed"))

context("mcmc tidiers")

if (suppressPackageStartupMessages(require(rstan, quietly = TRUE)))  {
  rstan_tests <- function(rstan_example) {
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
  }
  
  test_that("tidy returns indexes if requested on rstanarm fits (not CRAN)", {
    skip_on_cran()
    set.seed(2016)
    model_file <-
      system.file("example_data", "8schools.stan", package = "broom.mixed")
    schools_dat <- list(
      J = 8,
      y = c(28,  8, -3,  7, -1,  1, 18, 12),
      sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
    )
    rstan_example_notcran <- stan(
      file = model_file,
      data = schools_dat,
      iter = 100,
      chains = 2
    )
    # Sanity check for rds versions
    # TODO: does not work currently, anyone knows a more portable way to check
    # that the compiled model result has been updated?
##    rstan_example <- 
##      readRDS(system.file("example_data", "rstan_example.rds", package = "broom.mixed"))
##    expect_equal(rstan_example, rstan_example_notcran)
    rstan_tests(rstan_example_notcran)
  })
  
  test_that("tidy returns indexes if requested on rstanarm fits (CRAN)", {
    rstan_example <- readRDS(system.file("example_data", "rstan_example.rds", package = "broom.mixed"))
    rstan_tests(rstan_example)
  })
}
