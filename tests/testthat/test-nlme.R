
## test tidy, augment, glance methods from nlme-tidiers.R
stopifnot(require("testthat"), require("broom.mixed"))

## HACK: need to find the right generic
tidy <- broom.mixed:::tidy.lme


if (suppressPackageStartupMessages(require(nlme, quietly = TRUE))) {
    context("nlme models")
    
    d <- as.data.frame(ChickWeight)
    colnames(d) <- c("y", "x", "subj", "tx")
    fit <- lme(y ~ tx * x, random = ~x | subj, data = d)
    
    test_that("tidy works on nlme/lme fits", {
        td <- tidy(fit)
        expect_equal(names(td),
           c("effect", "group", "term", "estimate",
             "std.error", "df", "statistic", "p.value"))
    })
    
    test_that("augment works on lme4 fits with or without data", {
        au1 <- augment(fit)
        au2 <- augment(fit, d)
        expect_equal(au1,au2)
        expect_equal(dim(au1),c(578,7))
    })
    dNAs <- d
    dNAs$y[c(1, 3, 5)] <- NA
    
    test_that("augment works on lme fits with NAs and na.omit", {
        fitNAs <- lme(y ~ tx*x, random = ~x | subj, data = dNAs, 
                      na.action = "na.omit")
        au <- augment(fitNAs)
        expect_equal(nrow(au), sum(complete.cases(dNAs)))
    })
    
    
    test_that("augment works on lme fits with na.omit", {
        fitNAs <- lme(y ~ tx*x, random = ~x | subj, data = dNAs, 
                      na.action = "na.exclude")
        
        au <- augment(fitNAs, dNAs)
        
        # with na.exclude, should have NAs in the output where there were NAs in input
        expect_equal(nrow(au), nrow(dNAs))
        expect_equal(complete.cases(au), complete.cases(dNAs))
    })

    ## FIXME: weak tests - only shows that no error occurs and
    ##  the right type is returned!
    test_that("glance works on nlme fits", {
        expect_is(glance(fit),"data.frame")
    })
    
    
    testFit <- function(fit, data = NULL){
        test_that("Pinheiro/Bates fit works", {
            expect_is(tidy(fit, "fixed"),"data.frame")
            expect_is(tidy(fit),"data.frame")
            expect_is(glance(fit),"data.frame")
            if (is.null(data)) {
                expect_is(augment(fit),"data.frame")
            } else {
                expect_is(augment(fit, data),"data.frame")
            }
        })
    }
    
    testFit(lme(score ~ Machine, data = Machines, random = ~1 | Worker))
    testFit(lme(score ~ Machine, data = Machines, random = ~1 | Worker, method="ML"))
    testFit(lme(score ~ Machine, data = Machines, random = ~1 | Worker / Machine))
    testFit(lme(pixel ~ day + day ^ 2, data = Pixel, random = list(Dog = ~day, Side = ~1)))
    testFit(lme(pixel ~ day + day ^ 2 + Side, data = Pixel, 
                random = list(Dog = ~day, Side = ~1)))
    
    testFit(lme(yield ~ ordered(nitro)*Variety, data = Oats, 
                random = ~1/Block/Variety))
    # There are cases where no data set is returned in the result 
    # We can do nothing about this inconsistency but give a useful error message in augment
    fit  = nlme(conc ~ SSfol(Dose, Time, lKe, lKa, lCl), data = Theoph,
                random = pdDiag(lKe + lKa + lCl ~ 1))
    test_that(
        "Fit without data in returned structure works when data are given", {
        expect_true(testFit(fit, Theoph))
        })
    # When no data are passed, a meaningful message is issued
    expect_error(augment(fit), "explicit")
}
