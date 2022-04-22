stopifnot(require("testthat"), require("broom.mixed"))

if (require(nlme, quietly = TRUE)) {
  test_that("as.data.frame.lme argument handling", {
    fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
                data = Loblolly,
                fixed = Asym + R0 + lrc ~ 1,
                random = Asym ~ 1,
                start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
    rr <- ranef(fm1)
    expect_equal(head(data.frame(rr), 1),
                 structure(list(group = "Seed", term = "Asym", level = "329",
                              estimate = -5.56546756677834),
                         row.names = 1L, class = "data.frame"))
  expect_equal(head(data.frame(rr, stringsAsFactors= TRUE), 1),
               structure(list(group = structure(1L, levels = "Seed", class = "factor"),
                              term = structure(1L, levels = "Asym", class = "factor"),
                              level = structure(13L, levels = c("301", "303", "305", "307",
                                                                "309", "311", "315", "319", "321", "323", "325", "327", "329",
                                                                "331"), class = "factor"), estimate = -5.56546756677834), row.names = 1L, class = "data.frame")
               )
    })
}

