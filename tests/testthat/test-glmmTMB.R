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
} ## if require(glmmTMB)
