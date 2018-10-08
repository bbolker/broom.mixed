context("glmmTMB models")
stopifnot(require("testthat"), require("broom.mixed"))

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

  if (getRversion()>="3.6.0") {
      ## don't know why this fails on r-release (3.5.x);
      ## maybe getting glmmTMB from CRAN rather than GitHub??
      test_that("tidy respects conf.level", {
          tmpf <- function(cl=0.95) {
              return(tidy(zipm3,conf.int=TRUE,conf.level=cl)[1,][["conf.low"]])
          }
          expect_equal(tmpf(),-2.088147,tolerance=1e-4)
          expect_equal(tmpf(0.5),-0.7871105,tolerance=1e-4)
      })
  }

  test_that("empty components are OK", {
      expect_equal(dim(tidy(zipm3, effects = "ran_pars", component = "zi")),
                   c(0,5))
  })
} ## if require(glmmTMB)
