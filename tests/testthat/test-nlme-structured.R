test_that("tidy handles the stored diagonal lme model", {
    models <- new.env()
    load(system.file("extdata", "nlme_example.rda", package = "broom.mixed"),
         envir = models)
    td <- tidy(models$lmm2)
    expect_equal(td$estimate[td$effect == "fixed"],
                 unname(nlme::fixef(models$lmm2)))
    expect_equal(td$group[td$effect == "ran_pars"],
                 c(rep("Subject", 3), "Residual"))
})

test_that("structured lme covariance rows have matching groups and values", {
    structures <- list(nlme::pdSymm(~ age), nlme::pdDiag(~ age),
                       nlme::pdIdent(~ age), nlme::pdCompSymm(~ age))
    for (pd_structure in structures) {
        fit <- nlme::lme(distance ~ age + Sex, data = nlme::Orthodont,
                         random = list(Subject = pd_structure))
        vc <- nlme::getVarCov(fit)
        for (scale in c("sdcor", "vcov")) {
            td <- tidy(fit, effects = "ran_pars", scales = scale)
            expect_equal(td$effect, rep("ran_pars", 4))
            expect_equal(td$group, c(rep("Subject", 3), "Residual"))
            if (scale == "sdcor") {
                expect_equal(td$term, c("sd_(Intercept)", "cor_age.(Intercept)",
                                        "sd_age", "sd_Observation"))
                expected <- c(sqrt(vc[1, 1]), vc[2, 1] / sqrt(prod(diag(vc))),
                              sqrt(vc[2, 2]), sigma(fit))
            } else {
                expect_equal(td$term, c("var_(Intercept)", "cov_age.(Intercept)",
                                        "var_age", "var_Observation"))
                expected <- c(vc[1, 1], vc[2, 1], vc[2, 2], sigma(fit)^2)
            }
            expect_equal(td$estimate, unname(expected))
        }
    }
})

test_that("diagonal lme confidence intervals match nlme", {
    fit <- nlme::lme(distance ~ age + Sex, data = nlme::Orthodont,
                     random = list(Subject = nlme::pdDiag(~ age)))
    for (level in c(0.8, 0.95)) {
        td <- tidy(fit, effects = "ran_pars", conf.int = TRUE, conf.level = level)
        ci <- nlme::intervals(fit, which = "var-cov", level = level)$reStruct$Subject
        sd_rows <- match(c("sd_(Intercept)", "sd_age"), td$term)
        expect_equal(nrow(td), 4)
        expect_equal(td$estimate[sd_rows], unname(ci[, "est."]))
        expect_equal(td$conf.low[sd_rows], unname(ci[, "lower"]))
        expect_equal(td$conf.high[sd_rows], unname(ci[, "upper"]))
        expect_true(all(is.na(td$conf.low[-sd_rows])))
        expect_true(all(is.na(td$conf.high[-sd_rows])))
    }
})
