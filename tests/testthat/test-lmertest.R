stopifnot(require("testthat"), require("broom.mixed"))
## test lmerTest

if (require(lmerTest, quietly = TRUE)) {
  test_that("testing lmerTest p-values behind Douglas Bates' back", {
    lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    td <- tidy(lmm1, "fixed")
    expect_equal(td$df, c(17, 17), tolerance=1e-3)
    check_tidy(td, 2, 7, c(
      "effect", "term", "estimate",
      "std.error", "df", "statistic", "p.value"
    ))
    td_ran <- tidy(lmm1, "ran_pars")
    check_tidy(td_ran, 4, 4, c("effect", "group", "term", "estimate"))
    expect_false(all(is.na(td_ran$estimate)))

    if (requireNamespace("pbkrtest")) {
        td_kr <- tidy(lmm1, "fixed", ddf.method="Kenward-Roger")
        expect_equal(td_kr$df, c(17,17), tol=1e-4)
    }
    
    td_nodf <- tidy(lmm1, "fixed", ddf.method="lme4")
    check_tidy(td_nodf, 2, 5, c("effect", "term", "std.error", "statistic"))
  })
  test_that("Wald t intervals", {
      set.seed(101)
      ## unbalance to make K-R slightly different from Satterthwaite
      ss <- sleepstudy[sample(seq(nrow(sleepstudy)),size=round(0.9*nrow(sleepstudy))),]
      m1 <- lmer(Reaction~Days+(1|Subject),REML=TRUE,ss)
      tmpf <- function(ddfm="Satterthwaite") {
          tt <- tidy(m1,conf.int=TRUE,conf.method="Wald",ddf.method=ddfm,effect="fixed")
          unname(unlist(tt[,c("conf.low","conf.high")][1,]))
      }
      expect_equal(tmpf("Satterthwaite"), c(231.320558648089, 271.015491535434), tolerance=1e-6)
      expect_equal(tmpf("lme4"), c(232.327411469392, 270.008638714131), tolerance=1e-6)
      expect_equal(tmpf("Kenward-Roger"), c(231.331769141079, 271.004281042444), tolerance=1e-6)
  })
} ## if require(lmerTest)
