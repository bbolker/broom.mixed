stopifnot(require("testthat"), require("broom.mixed"))
context("glmmTMB models")
## source("helper-checkers.R") ## done automatically in test workflow

if (require(glmmTMB, quietly = TRUE)
    ## do we need all this?
    ## &&
    ## checkDepPackageVersion(dep_pkg = "TMB",
    ##                        this_pkg = "glmmTMB",
    ##                        warn = FALSE) &&
    ## checkDepPackageVersion(dep_pkg = "Matrix",
    ##                        this_pkg = "TMB",
    ##                        warn = FALSE)
    ) {
  L <- load(system.file("extdata", "glmmTMB_example.rda",
    package = "broom.mixed",
    mustWork = TRUE
    ))
  for (obj in L) {
    assign(obj, glmmTMB::up2date(get(obj)))
  }

  test_that("components included for zi models", {
    td <- tidy(zipm3)
    check_tidy(
      td, 29, 8,
      c(
        "effect", "component", "group", "term",
        "estimate", "std.error", "statistic", "p.value"
      )
    )
  })

    test_that("tidy respects components argument", {
        tdc <- tidy(zipm3, component="cond", effects="fixed")
        check_tidy(
            tdc, 14, 7,
            c(
                "effect", "component", "term",
                "estimate", "std.error", "statistic", "p.value"
            )
        )
        tdz <- tidy(zipm3, component="zi", effects="fixed")
        check_tidy(
            tdc, 14, 7,
            c(
                "effect", "component", "term",
                "estimate", "std.error", "statistic", "p.value"
            )
        )
  })

  test_that("tidy respects conf.level", {
      tmpf <- function(cl=0.95) {
          return(tidy(zipm3,conf.int=TRUE,conf.level=cl)[1,][["conf.low"]])
      }
      expect_equal(tmpf(),-2.088147,tolerance=1e-4)
      expect_equal(tmpf(0.5),-0.7871105,tolerance=1e-4)
  })

  test_that("empty components are OK", {
      expect_equal(dim(tidy(zipm3, effects = "ran_pars", component = "zi")),
                   c(0,5))
  })

  test_that("ran_pars works", {
      expect_equal(dim(tidy(zipm3, effects = "ran_pars", component = "zi")),
                   c(0,5))
  })

  test_that("ran_vals works", {
      ## GH103
      expect_equal(dim(tidy(glmm1, effects="ran_vals")), c(15,7))
  })

  test_that("profile tidying works", {
      td <- tidy(glmm1, effects="fixed", conf.int=TRUE, conf.method="profile")
      check_tidy(
          td, 4, 9,
      c(
        "effect", "component", "term",
        "estimate", "std.error", "statistic", "p.value",
        "conf.low","conf.high"))
      expect_equal(td$conf.low,
                   c(-1.9012409337, -1.61676924, -1.8010155, -2.50085),
                   tolerance=1e-4)
  })

  test_that("confint with non-pos-def results", {

      d <- data.frame(x = rep(1,100))
      suppressWarnings(m1 <- glmmTMB(x~1, family=nbinom2, data=d))
      expect_is(tidy(m1), "tbl_df")
      res <- unlist(suppressWarnings(
          tidy(m1, conf.int=TRUE)[,c("estimate", "conf.low","conf.high")]))
      ## this test is not reliable ... not sure what else to test?
      ## expect_equal(unname(res), c(-0.196,0.196), tolerance=1e-4)
      expect_equal(unname(glance(m1)[,c("nobs","df.residual")]),
                   unname(tibble::tibble(100L, 98L)))
  })

  test_that("confint with multiple REs", {
      if (packageVersion("glmmTMB") > "1.1.3" && requireNamespace("lme4")) {
          dd <- expand.grid(r = 1:10, a = factor(1:2), b = factor(1:3),
                            f = factor(1:5), g = factor(1:6))
        dd$y <- simulate(
            seed = 101,
            ~ 1 + (a|f) + (b|g),
            newdata = dd,
            newparams = list(beta = 1,
                             theta = rep(1,9),
                             sigma = 1),
            family = gaussian)[[1]]
          res <- glmmTMB(y~ 1 + (a+0|f) + (b+0|g), data = dd)
          td <- tidy(res, conf.int = TRUE)
          check_tidy(
              td, 11, 10,
              c("effect", "component", "group", "term", "estimate", "std.error", 
                "statistic", "p.value", "conf.low", "conf.high"))
          } ## require lme4 (for simulate)
  }
  )
  if (requireNamespace("lme4")) {
      ## GH #136
      data("sleepstudy", package = "lme4")
      ## FIXME: speed up by storing this?
      ## warning because using bogus example
      suppressWarnings(fm3ZIP <-
                           glmmTMB(round(Reaction) ~ Days + (1|Subject), family=poisson,
                                   ziformula=~(1|Subject),
                                   data = sleepstudy))
      t1 <- tidy(fm3ZIP, conf.int = TRUE, component = "cond", effect = "ran_pars")
      t2 <- tidy(fm3ZIP, conf.int = TRUE, effect = "ran_pars")
      expect_identical(nrow(t1), 1L)
      expect_identical(nrow(t2), 2L)
  } ## if requireNamespace("lme4")


  test_that("vv without RE in conditional model", {
      m <- glmmTMB(mpg ~ hp, data = mtcars,
                   dispformula = ~ 1 + (1|cyl), family = gaussian)
      expect_equal(nrow(tidy(m)), 2L)
   })

  test_that("tidy reports dispersion for glmmTMB fits with an estimated scale parameter", {
      ## Families like Gamma (or nbinom2) have a genuinely estimated
      ## residual/dispersion parameter that VarCorr()'s per-component
      ## useSc attribute does not reflect (it describes glmmTMB's internal
      ## Cholesky parameterization, which is unscaled regardless of
      ## family), unlike the overall VarCorr object's own "sc" attribute.
      ## Previously this meant the dispersion estimate was silently
      ## dropped from tidy() output entirely for such models.
      skip_if_not_installed("glmmTMB")
      epil2 <- get(data("epil2", package = "glmmTMB", envir = environment()))
      epil2_red <- subset(epil2, y > 0)
      mod_tmb <- glmmTMB(y ~ trt + (1 | subject),
        data = epil2_red, family = Gamma(link = "log")
      )

      td <- tidy(mod_tmb)
      td_ran <- td[td$effect == "ran_pars", ]
      expect_equal(nrow(td_ran), 2)
      expect_equal(td_ran$group, c("subject", "Residual"))
      expect_equal(td_ran$term, c("sd__(Intercept)", "sd__Observation"))
      expect_equal(td_ran$estimate[2], sigma(mod_tmb))

      ## models without a real scale parameter (e.g. binomial) must be
      ## unaffected: no spurious Residual row
      skip_if_not_installed("lme4")
      cbpp <- get(data("cbpp", package = "lme4", envir = environment()))
      gm1 <- glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd),
        data = cbpp, family = binomial
      )
      td_gm1 <- tidy(gm1, effects = "ran_pars")
      expect_false("Residual" %in% td_gm1$group)

      ## gaussian fits already got a Residual row via the lme4-borrowed
      ## as.data.frame.VarCorr.merMod dispatch (useSc == TRUE there); must
      ## not get a second, duplicate one
      sleepstudy <- get(data("sleepstudy", package = "lme4", envir = environment()))
      fm <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
      td_fm <- tidy(fm, effects = "ran_pars")
      expect_equal(sum(td_fm$group == "Residual"), 1)

      ## NB: conf.int=TRUE for ran_pars is not checked here -- it is
      ## affected by a separate, pre-existing bug in the s_component
      ## filtering used for theta_ confidence intervals (CIs are
      ## currently misaligned/dropped for *all* glmmTMB models, not
      ## specific to dispersion handling).
   })
} ## if require(glmmTMB)
