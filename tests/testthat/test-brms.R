stopifnot(require("testthat"), require("broom.mixed"))

if (require(brms, quietly = TRUE) && require(rstanarm, quietly=TRUE)) {
  load(system.file("extdata", "brms_example.rda",
    package = "broom.mixed",
    mustWork = TRUE
    ))

  ## GH #87
  tt <- suppressWarnings(tidy(brms_multi,conf.int=TRUE))
  expect_true(all(c("conf.low","conf.high") %in%
                  names(tt)))

  ## GH #101
  gg <- glance(brms_noran)
  expect_equal(names(gg),c("algorithm","pss","nobs","sigma"))

  ## Check the descriptive columns of tidy summaries
  ### brms_RE
  expected <- tibble::tribble(
     ~effect, ~component,     ~group,                         ~term,
     "fixed",     "cond",         NA,                 "(Intercept)",
     "fixed",     "cond",         NA,                  "Days_extra",
  "ran_pars",     "cond",  "Subject",             "sd__(Intercept)",
  "ran_pars",     "cond",  "Subject",              "sd__Days_extra",
  "ran_pars",     "cond",  "Subject", "cor__(Intercept).Days_extra",
  "ran_pars",     "cond", "Residual",             "sd__Observation"
  )
  observed <- suppressWarnings(tidy(brms_RE))
  expect_equal(observed[, 1:4], expected)
  ### brms_noran
  expected <- tibble::tribble(
     ~effect, ~component,     ~group,             ~term,
     "fixed",     "cond",         NA,     "(Intercept)",
     "fixed",     "cond",         NA,              "wt",
  "ran_pars",     "cond", "Residual", "sd__Observation"
  )
  observed <- suppressWarnings(tidy(brms_noran))
  expect_equal(observed[, 1:4], expected)
  ### brms_brm_fit4
  expected <- tibble::tribble(
  ~effect, ~component, ~group,         ~term,
  "fixed",     "cond",     NA, "(Intercept)",
  "fixed",     "cond",     NA,           "x"
  )
  observed <- suppressWarnings(tidy(brms_brm_fit4))
  expect_equal(observed[, 1:4], expected)

}
