stopifnot(require("testthat"), require("broom.mixed"))
## test tidy, augment, glance methods from lme4-tidiers.R

## HACK: need to make sure we find the right generic!
tidy <- broom.mixed:::tidy.merMod

if (require(lme4, quietly = TRUE)) {
    context("lme4 models")
    
    d <- as.data.frame(ChickWeight)
    colnames(d) <- c("y", "x", "subj", "tx")
    fit <<- lmer(y ~ tx*x + (x | subj), data=d)
    
    test_that("tidy works on lme4 fits", {
        td <- tidy(fit)
        expect_equal(dim(td),c(12,6))
        expect_equal(names(td),
             c("effect", "group", "term", "estimate",
               "std.error", "statistic"))
        expect_equal(td$term,
                     c("(Intercept)", "tx2", "tx3", "tx4", "x",
                       "tx2:x", "tx3:x", "tx4:x",
                       "sd_(Intercept)", "sd_x",
                       "cor_(Intercept).x", "sd_Observation"))
    })
    
    test_that("tidy/glance works on glmer fits",{
      gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                  cbpp, binomial, nAGQ = 0)
      ggm <- glance(gm)
      expect_equal(names(ggm), c("sigma", "logLik", "AIC", "BIC", "deviance", "df.residual"))
      td <- tidy(gm)
      expect_equal(names(td),
                   c("effect", "group", "term", "estimate",
                     "std.error", "statistic", "p.value"))
      td_ran <- tidy(gm, "ran_pars")
      expect_equal(names(td_ran), c("effect","group", "term", "estimate"))
    })

    test_that("tidy works on non-linear fits",{
      startvec <- c(Asym = 200, xmid = 725, scal = 350)
      # use nAGQ = 0 to avoid warnings
      nm <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
                  Orange, start = startvec, nAGQ = 0L)
      gnm <- glance(nm)
      expect_equal(names(gnm), c("sigma", "logLik", "AIC", "BIC", "deviance", "df.residual"))
      td <- tidy(nm)
      expect_equal(names(td),
                   c("effect", "group", "term", "estimate",
                     "std.error", "statistic"))
      td_ran <- tidy(nm, "ran_pars")
      expect_equal(names(td_ran), c("effect","group", "term", "estimate"))
    })
    
    test_that("scales works", {
        t1 <- tidy(fit,effects="ran_pars")
        t2 <- tidy(fit,effects="ran_pars",scales="sdcor")
        expect_equal(t1$estimate,t2$estimate)
        expect_error(tidy(fit,effects="ran_pars",scales="varcov"),
                     "unrecognized ran_pars scale")
        t3  <- tidy(fit,effects="ran_pars",scales="vcov")
        expect_equal(t3$estimate[c(1,2,4)],
                     t2$estimate[c(1,2,4)]^2)
        expect_error(tidy(fit,scales="vcov"),
                     "must be provided for each effect")
              })
    
    test_that("tidy works with more than one RE grouping variable", {
       dd <- expand.grid(f=factor(1:10),g=factor(1:5),rep=1:3)
       dd$y <- suppressMessages(simulate(~(1|f)+(1|g),newdata=dd,
                        newparams=list(beta=1,theta=c(1,1)),
                        family=poisson, seed=101))[[1]]
       gfit <- glmer(y~(1|f)+(1|g),data=dd,family=poisson)
       tnames <- as.character(tidy(gfit,effects="ran_pars")$term)
       expect_equal(tnames,rep("sd_(Intercept)",2))
   })
              
    test_that("augment works on lme4 fits with or without data", {
        au1 <- augment(fit)
        au2 <- augment(fit, d)
        ## FIXME: columns not ordered the same??
        expect_equal(au1,au2[names(au1)])
    })

    dNAs <- d
    dNAs$y[c(1, 3, 5)] <- NA
    
    test_that("augment works on lme4 fits with NAs", {
        fitNAs <- lmer(y ~ tx*x + (x | subj), data = dNAs)
        au <- augment(fitNAs)
        expect_equal(nrow(au), sum(complete.cases(dNAs)))
    })
    
    test_that("augment works on lme4 fits with na.exclude", {
        fitNAs <- lmer(y ~ tx*x + (x | subj), data = dNAs, na.action = "na.exclude")
        
        #expect_error(suppressWarnings(augment(fitNAs)))
        au <- augment(fitNAs, dNAs)
        
        # with na.exclude, should have NAs in the output where there were NAs in input
        expect_equal(nrow(au), nrow(dNAs))
        expect_equal(complete.cases(au), complete.cases(dNAs))
    })

    test_that("glance works on lme4 fits", {
        g <- glance(fit)
        expect_equal(dim(g),c(1,6))
    })

    test_that("ran_modes works", {
        fm1 <- lmer(Reaction~Days+(1|Subject),sleepstudy)
        fm2 <- lmer(Reaction~Days+(Days|Subject),sleepstudy)
        td1 <- tidy(fm1,"ran_modes")
        td2 <- tidy(fm2,"ran_modes")
        expect_equal(dim(td1),c(18,6))
        expect_equal(dim(td2),c(36,6))
    })

}
