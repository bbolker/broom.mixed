## test tidy and glance methods from rstanarm_tidiers.R
stopifnot(require("testthat"), require("broom.mixed"), require("broom"))

if (suppressPackageStartupMessages(require(rstanarm, quietly = TRUE))) {
  load(system.file("extdata", "rstanarm_example.rda", package = "broom.mixed"))
  ## fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
  ##   iter = 200, chains = 2)

  context("rstanarm tidiers")
  test_that("tidy works on rstanarm fits", {
    td1 <- tidy(fit)
    td2 <- tidy(fit, effects = "ran_vals")
    td3 <- tidy(fit, effects = "ran_pars")
    td4 <- tidy(fit, effects = "auxiliary")
    expect_equal(colnames(td1), c("term", "estimate", "std.error", "group"))
  })

  test_that("tidy with multiple 'effects' selections works on rstanarm fits", {
    td1 <- tidy(fit, effects = c("ran_vals", "auxiliary"))
    expect_true(all(c("sigma", "mean_PPD") %in% td1$term))
    expect_equal(colnames(td1), c("term", "estimate", "std.error", "level", "group"))
  })

  test_that("conf.int works on rstanarm fits", {
    td1 <- tidy(fit, conf.int = TRUE, conf.level = 0.8)
    td2 <- tidy(fit, effects = "ran_vals", conf.int = TRUE, conf.level = 0.5)
    td3 <- tidy(fit, conf.int = TRUE, conf.level = 0.95)
    nms <- c("level", "group", "term", "estimate", "std.error", "conf.low", "conf.high")
    expect_equal(colnames(td2), nms)
    ## FIXME: why NA values in std.error/conf.low/conf.high here?
    expect_true(all(is.na(td3$conf.low) | td3$conf.low < td1$conf.low))
    expect_true(all(is.na(td3$conf.high) | td3$conf.high > td1$conf.high))
  })

  test_that("glance works on rstanarm fits", {
    skip_on_cran()
    g1 <- glance(fit)
    g2 <- glance(fit, looic = TRUE, cores = 1, k_threshold = 0.7)
    expect_equal(colnames(g1), c("algorithm", "pss", "nobs", "sigma"))
    expect_equal(colnames(g2), c(colnames(g1), "looic", "elpd_loo", "p_loo"))
  })

  test_that("exponentiation", {
      td1 <- tidy(fit2, conf.int = TRUE, effects = "fixed")
      td1e <- tidy(fit2, conf.int = TRUE, exponentiate = TRUE, effects = "fixed")
      expect_equal(td1e$estimate, exp(td1$estimate))
      expect_equal(td1e$conf.low, exp(td1$conf.low))
      expect_equal(td1e$conf.high, exp(td1$conf.high))
      expect_equal(td1e$std.error, exp(td1$estimate)*td1$std.error)
  })

  ## GH 153
  if (requireNamespace("mice", quietly = TRUE)) {
      test_that("mice imputed data", {
          data(nhanes, package = "mice")
          imp <- mice::mice(nhanes, m = 3, print = FALSE)
          suppressWarnings(
              capture.output(ms <- with(imp,
                                        stan_glm(age ~ bmi + chl)))
          )
          expect_is(summary(mice::pool(ms)), "mipo.summary")
      })
  }

} ## rstanarm available

