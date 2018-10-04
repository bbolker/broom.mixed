stopifnot(require("testthat"), require("broom.mixed"))


if (require(glmmTMB, quietly = TRUE)) {
  load(system.file("extdata", "glmmTMB_example.rda",
    package = "broom.mixed",
    mustWork = TRUE
  ))

  context("glmmTMB models")

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
}

test_that("tidy respects conf.level", {
     tmpf <- function(cl=0.95) {
         return(tidy(zipm3,conf.int=TRUE,conf.level=cl)[1,][["conf.low"]])
     }
     expect_equal(tmpf(),-2.088147,tolerance=1e-4)
     expect_equal(tmpf(0.5),-0.7871105,tolerance=1e-4)
})
