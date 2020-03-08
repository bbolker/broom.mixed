
stopifnot(require("testthat"), require("broom.mixed"))

if (require(brms, quietly = TRUE) && require(rstanarm, quietly=TRUE)) {
  load(system.file("extdata", "brms_example.rda",
    package = "broom.mixed",
    mustWork = TRUE
    ))

  context("brms models")

  ## GH #87
  expect_true(all(c("conf.low","conf.high") %in%
                  names(tidy(brms_multi,conf.int=TRUE))))
}
