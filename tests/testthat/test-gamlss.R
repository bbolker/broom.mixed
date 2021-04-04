library(broom.mixed)
if (require("testthat") && require("gamlss") && require("gamlss.data")) {
    data(abdom, package="gamlss.data")
    gamlss1 <- gamlss(
        y ~ pb(x),
        sigma.fo = ~ pb(x),
        family = BCT,
        data = abdom,
        method = mixed(1, 20)
    )
    aa <- augment(gamlss1, newdata = abdom)
    test_that("gamlss augment with new data", {
    expect_equal(head(aa),
                 structure(list(y = c(59L, 64L, 56L, 61L, 74L, 60L),
                                x = c(12.29, 12.29, 12.29, 12.43, 12.71, 12.71),
          .fitted = c(60.3907912240965, 
                      60.3907912240965, 60.3907912240965,
                      62.1788058184065, 65.7518835998086, 65.7518835998086)),
          row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame")))
    })
}


