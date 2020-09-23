stopifnot(require("testthat"), require("broom.mixed"))
context("glmmTMB models")
## source("helper-checkers.R") ## done automatically in test workflow

if (require(glmmTMB, quietly = TRUE)) {
  load(system.file("extdata", "glmmTMB_example.rda",
    package = "broom.mixed",
    mustWork = TRUE
  ))

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
      d <- data.frame(x=rep(1,100))
      suppressWarnings(m1 <- glmmTMB(x~1,family=nbinom2,data=d))
      expect_is(tidy(m1),"tbl_df")
      res <- unlist(suppressWarnings(tidy(m1,conf.int=TRUE)[,c("conf.low","conf.high")]))
      expect_equal(unname(res),c(-0.196,0.196),tolerance=1e-4)
  })
  
} ## if require(glmmTMB)

