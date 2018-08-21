## run in pkg root
if (sub(".*/","",getwd())!="broom.mixed") {
    stop("should run from package root directory")
}

## FIXME: should disaggregate/implement a Makefile!
run_brms <- FALSE ## slow, only do if necessary


save_file <- function(..., pkg, type="rda") {
  f <- file.path("inst","extdata",sprintf("%s_example.%s",pkg,type))
  if (type=="rda") {
    save(...,file=f)
  } else {
    saveRDS(...,file=f)
  }
  invisible(NULL)
}

run_pkg <- function(pkg,e) {
  if (require(pkg,character.only=TRUE)) {
    eval(e)
    return(TRUE)
  } else {
    cat(sprintf("%s examples not run\n",pkg))
    return(FALSE)
  }
}

run_pkg("nlme",
{
    data(sleepstudy,package="lme4")
    lmm0 <- lme(Reaction ~ Days, random = ~ 1| Subject, sleepstudy)
    lmm0ML <- lme(Reaction ~ Days, random = ~ 1| Subject, sleepstudy,
                  method="ML")
    lmm1 <- lme(Reaction ~ Days, random = ~ Days | Subject, sleepstudy)
    ## doesn't work yet
    lmm2 <- lme(Reaction ~ Days, random = list(Subject=pdDiag(~Days)), sleepstudy)
    save_file(lmm0, lmm0ML, lmm1, lmm2, pkg="nlme", type = "rda")
})

run_pkg("lme4",
{
    lmm0 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
    lmm0ML <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy,REML=FALSE)
    lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    lmm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days | Subject), sleepstudy)
    lmm1_prof <- profile(lmm1)
    save_file(lmm1_prof, lmm0, lmm0ML, lmm1, lmm2, pkg="lme4", type = "rda")
})


run_pkg("rstan",
{
          rstan_options(auto_write=TRUE)
          model_file <- system.file("extdata", "8schools.stan", package = "broom.mixed")
          schools_dat <- list(J = 8, 
                              y = c(28,  8, -3,  7, -1,  1, 18, 12),
                              sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
          set.seed(2015)
          rstan_example <- stan(file = model_file, data = schools_dat, 
                                iter = 1000, chains = 2, save_dso = FALSE)
          save_file(rstan_example,pkg="rstan",type="rds")
        })

run_pkg("glmmTMB",
        {
          data(sleepstudy,package="lme4")
          data(cbpp,package="lme4")
          lmm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
          glmm1 <- glmmTMB(incidence/size ~ period + (1 | herd),
                           data = cbpp, family = binomial, weights=size)
          save_file(lmm1,glmm1,pkg="glmmTMB",type="rda")
        })

run_pkg("glmmADMB",
        {
          data(sleepstudy,package="lme4")
          lmm1 <- glmmadmb(Reaction ~ Days + (Days | Subject), sleepstudy,
                           family="gaussian")
          save_file(lmm1,pkg="glmmADMB",type="rda")
        })

run_pkg("MCMCglmm",
{
    data("sleepstudy",package="lme4")
    mm0 <- MCMCglmm(Reaction ~ Days, random = ~ Subject, data=sleepstudy,
                    pr=TRUE)
    mm1 <- MCMCglmm(Reaction ~ Days, random = ~us(1+Days):Subject,
                    ## parameter-expanded priors
                    ## V is 2x2 identity wlog
                    ## t(2) with standard dev of 2
                    prior=list(G=list(list(nu=2,V=diag(2),
                                           alpha.mu=rep(0,2),
                                           alpha.V=diag(rep(4,2))))),
                    data=sleepstudy, pr=TRUE)
    mm2 <- MCMCglmm(Reaction ~ Days, random = ~idh(1+Days):Subject, data=sleepstudy, pr=TRUE)
    save_file(mm0, mm1, mm2, pkg="MCMCglmm", type = "rda")
})

## slowest stuff last

run_pkg("rstanarm",
        {
          fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars, 
                            iter = 300, chains = 2)
          save_file(fit, pkg="rstanarm", type="rds")
        })


if (run_brms) {
    run_pkg("brms",
        {
            fit <- brm(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars, 
                       iter = 500, chains = 2,
                       ## ugh: compiled object gets stored as part of fitted
                       ##  object - not cached. So it just blows up the
                       ##  size of the fitted object, doesn't speed things
                       ##  up when re-running this code ...
                       save_dso= FALSE)
            save_file(fit, pkg="brmsfit", type="rds")
        })
}

