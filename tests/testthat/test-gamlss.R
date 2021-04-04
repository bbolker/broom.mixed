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

    df <- structure(list(ID = structure(c(1L, 1L, 3L, 3L, 1L, 1L, 3L, 3L, 
                                          2L, 2L, 5L, 5L, 2L, 2L, 4L, 4L, 5L, 5L, 2L, 2L, 4L, 4L, 5L, 5L, 
                                          1L, 1L, 3L, 3L), .Label = c("1", "2", "3", "4", "5"), class = "factor"), 
                         score = c(0.041929622875205, 0.0825927706160165, 0.0885016980957589, 
                                   0.10072436783134, 3.40303792499053e-05, 6.34114774779728e-05, 
                                   0.199542264154635, 0.485745295078827, 0.000192903619829098, 
                                   3.53678015687278e-06, 1.73903243155792e-05, 9.49768211665946e-06, 
                                   0.160452490810036, 0.118529390854194, 0.00577127495391104, 
                                   0.00633025242790636, 0.387891969035188, 0.296403006174033, 
                                   0.0239720610077002, 0.0201663110699289, 1.69247724516557e-05, 
                                   4.27384007240083e-05, 0.0219262718916132, 0.0211082684976785, 
                                   1.30094341980175e-06, 3.20199933225434e-06, 0.106151246774257, 
                                   0.443433279410949), time = structure(c(1L, 2L, 1L, 2L, 1L, 
                                                                          2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 
                                                                          1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), .Label = c("early", "late"
                                                                                                                      ), class = "factor"), region = structure(c(2L, 2L, 2L, 2L, 
                                                                                                                                                                 3L, 3L, 3L, 3L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                                                                                                                                                 3L, 3L, 3L, 3L, 3L, 1L, 1L, 1L, 1L), .Label = c("1", "2", 
                                                                                                                                                                 "3"), class = "factor")), row.names = c(NA, -28L), class = "data.frame")

    m1 <- gamlss(score ~ time + region
               , sigma.formula = ~ region
               , data=df
               , family = BE()
               , control=gamlss.control(trace=FALSE))
    t1 <- tidy(m1)  
    test_that("basic gamlss tidying", {
        expect_equal(dim(t1),c(7,6))
        expect_equal(t1[c("parameter","term")],
                     structure(list(parameter =
                                        rep(c("mu", "sigma"), c(4,3)),
          term = c("(Intercept)", "timelate", "region2", "region3",
                   "(Intercept)", "region2", "region3")),
             row.names = c(NA, -7L),
              class = c("tbl_df", "tbl", "data.frame")))
    })
    m2 <- gamlss(score ~ time + region + random(ID)
               , sigma.formula = ~ 1
               , data=df, family = BE()
               , control = gamlss.control(n.cyc = 1e4, trace=FALSE))
    t2 <- tidy(m2)
    test_that("basic gamlss tidying 2", {
        expect_equal(dim(t2),c(5,6))
        expect_equal(t2[c("parameter","term")],
                     structure(list(parameter =
                                        rep(c("mu", "sigma"), c(4,1)),
              term = c("(Intercept)", "timelate", "region2", 
                           "region3", "(Intercept)")), row.names = c(NA, -5L),
              class = c("tbl_df", "tbl", "data.frame")))
    })

}


