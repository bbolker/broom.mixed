stopifnot(require("testthat"), require("broom.mixed"))
## test lmerTest

if (require(lmerTest, quietly = TRUE)) {
  test_that("testing lmerTest p-values behind Douglas Bates' back", {
    lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    td <- tidy(lmm1, "fixed")
    expect_equal(td$df, c(17, 17), tolerance = 1e-3)
    check_tidy(td, 2, 7, c(
      "effect", "term", "estimate",
      "std.error", "df", "statistic", "p.value"
    ))
    td_ran <- tidy(lmm1, "ran_pars")
    check_tidy(td_ran, 4, 4, c("effect", "group", "term", "estimate"))
    expect_false(all(is.na(td_ran$estimate)))

    if (requireNamespace("pbkrtest")) {
      td_kr <- tidy(lmm1, "fixed", ddf.method = "Kenward-Roger")
      expect_equal(td_kr$df, c(17, 17), tol = 1e-4)
    }

    td_nodf <- tidy(lmm1, "fixed", ddf.method = "lme4")
    check_tidy(td_nodf, 2, 5, c("effect", "term", "std.error", "statistic"))
  })
}
