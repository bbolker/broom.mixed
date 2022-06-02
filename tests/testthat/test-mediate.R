
stopifnot(require("testthat"), require("broom.mixed"))
## test tidy method from mediation_tidiers.R

if (require(lme4, quietly = TRUE) && require(mediation, quietly = TRUE)) {
  load(system.file("extdata", "mediation_example.rda",
                   package = "broom.mixed",
                   mustWork = TRUE
  ))
  
  context("mediation models")
  
  d <- CO2
  colnames(d) <- c("id", "loc", "tx", "x", "y")
  d$tx <- as.integer(d$tx) - 1L
  fit <- lmer(y ~ tx + x + loc + (1 | id), data = d)
  med <- lmer(x ~ tx + loc + (1 | id), data = d)
  mod <- mediate(med, fit, treat = "tx", mediator = "x", sims = 20L)
  
  test_that("tidy works on multilevel mediation fits", {
    td <- tidy(mod)
    expect_equal(dim(td), c(4L, 4L))
    expect_equal(
      names(td),
      c("term", "estimate", "std.error", "p.value")
    )
    expect_equal(td$term, c("acme_0", "acme_1", "ade_0", "ade_1"))
  })
  
  test_that("conf.int adds columns and preserves term names", {
    tdci <- tidy(mod, conf.int = TRUE)
    expect_equal(
      names(tdci),
      c("term", "estimate", "std.error", "p.value", "conf.low", "conf.high")
    )
    expect_equal(tdci$term, c("acme_0", "acme_1", "ade_0", "ade_1"))
  })
}
