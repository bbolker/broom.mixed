## run in pkg root

save_file <- function(...,pkg,type="rda") {
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

run_pkg("lme4",
        {
          lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
          pp <- profile(lmm1)
          save_file(pp, pkg="lme4", type = "rds")
        })


run_pkg("rstan",
        {
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

run_pkg("rstanarm",
        {
          fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars, 
                            iter = 300, chains = 2)
          save_file(fit,pkg="rstanarm",type="rds")
        })


run_pkg("brms",
        {
            fit <- brm(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars, 
                       iter = 500, chains = 2, save_dso=FALSE)
            save_file(fit,pkg="brmsfit",type="rds")
        })
