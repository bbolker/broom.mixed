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

      ## conf.int=TRUE: the Residual row's CI must come from confint(x,
      ## parm="sigma"), correctly aligned with (not swapped with, or
      ## dropped in favour of) the random-effect term's own CI (see the
      ## "conf.int alignment" tests below for the underlying bug this
      ## depends on being fixed)
      td_ci <- tidy(mod_tmb, conf.int = TRUE)
      td_ci_ran <- td_ci[td_ci$effect == "ran_pars", ]
      expect_equal(
        as.numeric(td_ci_ran[1, c("conf.low", "conf.high")]),
        as.numeric(confint(mod_tmb, parm = "theta_", method = "Wald")[1, 1:2])
      )
      expect_equal(
        as.numeric(td_ci_ran[2, c("conf.low", "conf.high")]),
        as.numeric(confint(mod_tmb, parm = "sigma", method = "Wald")[1, 1:2])
      )
   })

  test_that("ran_pars conf.int alignment with theta_/sigma confint() (single component)", {
      ## Regression test for a pre-existing bug: safe_confint()'s
      ## s_component regex filter assumed theta_ CI row names are always
      ## prefixed with the component name (e.g. "cond."), but glmmTMB only
      ## adds that prefix when random effects span more than one
      ## component; with a single component (the common case) row names
      ## are unprefixed, so the filter matched nothing and CIs came out
      ## entirely dropped, or (worse) misaligned with the wrong ran_pars
      ## row once a Residual row was also present.
      skip_if_not_installed("lme4")
      sleepstudy <- get(data("sleepstudy", package = "lme4", envir = environment()))
      fm <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
      td <- tidy(fm, conf.int = TRUE, effects = "ran_pars")
      raw_theta <- confint(fm, parm = "theta_", method = "Wald")[, 1:2, drop = FALSE]
      raw_sigma <- confint(fm, parm = "sigma", method = "Wald")
      expect_equal(nrow(td), nrow(raw_theta) + 1)
      expect_equal(
        unname(as.matrix(td[seq_len(nrow(raw_theta)), c("conf.low", "conf.high")])),
        unname(raw_theta),
        ignore_attr = TRUE
      )
      expect_equal(
        as.numeric(td[nrow(td), c("conf.low", "conf.high")]),
        as.numeric(raw_sigma[1, 1:2])
      )

      cbpp <- get(data("cbpp", package = "lme4", envir = environment()))
      gm1 <- glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd),
        data = cbpp, family = binomial
      )
      td_gm1 <- tidy(gm1, conf.int = TRUE, effects = "ran_pars")
      raw_gm1 <- confint(gm1, parm = "theta_", method = "Wald")
      expect_equal(
        as.numeric(td_gm1[1, c("conf.low", "conf.high")]),
        as.numeric(raw_gm1[1, 1:2])
      )
   })

  test_that("ran_pars conf.int alignment with theta_/sigma confint() (multiple components)", {
      ## when random effects DO span multiple components, glmmTMB prefixes
      ## theta_ CI row names (e.g. "cond.Std.Dev...", "zi.Std.Dev..."), and
      ## the s_component filter must still correctly select only the
      ## requested component(s)
      skip_if_not_installed("lme4")
      sleepstudy <- get(data("sleepstudy", package = "lme4", envir = environment()))
      suppressWarnings(
        fm3ZIP <- glmmTMB(round(Reaction) ~ Days + (1 | Subject),
          family = poisson, ziformula = ~ (1 | Subject), data = sleepstudy
        )
      )
      td_both <- suppressWarnings(tidy(fm3ZIP, conf.int = TRUE, effects = "ran_pars"))
      expect_equal(nrow(td_both), 2)
      expect_equal(td_both$component, c("cond", "zi"))

      td_cond <- suppressWarnings(
        tidy(fm3ZIP, conf.int = TRUE, effects = "ran_pars", component = "cond")
      )
      expect_equal(nrow(td_cond), 1)
      expect_equal(td_cond$component, "cond")
      expect_equal(td_cond$estimate, td_both$estimate[1])
   })

  test_that("no Residual row/CI added when dispformula is non-trivial", {
      ## sigma(x) and confint(x, parm="sigma") are not comparable when
      ## dispformula has covariates/random effects of its own: glmmTMB's
      ## own confint() falls back to reporting the raw (untransformed,
      ## link-scale) dispformula coefficient CI instead of a natural-scale
      ## CI for sigma() in that case, so no Residual row should be added
      ## (and any pre-existing one from the lme4-borrowed VarCorr.merMod
      ## dispatch must be stripped, since it would use the same
      ## not-really-comparable sigma(x) value)
      skip_if_not_installed("lme4")
      sleepstudy <- get(data("sleepstudy", package = "lme4", envir = environment()))
      sleepstudy$grp2 <- factor(rep(1:3, length.out = nrow(sleepstudy)))
      suppressWarnings(
        fm_nt <- glmmTMB(Reaction ~ Days + (1 | Subject), data = sleepstudy,
          dispformula = ~ 1 + (1 | grp2)
        )
      )
      td <- suppressWarnings(tidy(fm_nt, conf.int = TRUE, effects = "ran_pars"))
      expect_false("Residual" %in% td$group)
      expect_equal(nrow(td), 1)
   })
} ## if require(glmmTMB)
