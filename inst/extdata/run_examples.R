## run in pkg root
if (sub(".*/", "", getwd()) != "broom.mixed") {
  stop("should run from package root directory")
}

source(system.file("extdata","example_helpers.R",package="broom.mixed"))

## FIXME: should disaggregate/implement a Makefile!
## environment variable, command-line arguments to R script?
## slow stuff, disable for speed
run_brms <- FALSE
run_stan <- FALSE

run_pkg("nlme", {
  data(sleepstudy, package = "lme4")
  lmm0 <- lme(Reaction ~ Days, random = ~1 | Subject, sleepstudy)
  lmm0ML <- lme(Reaction ~ Days,
    random = ~1 | Subject, sleepstudy,
    method = "ML"
  )
  lmm1 <- lme(Reaction ~ Days, random = ~Days | Subject, sleepstudy)
  ## this model doesn't work yet ...erm
  lmm2 <- lme(Reaction ~ Days, random = list(Subject = pdDiag(~Days)), sleepstudy)
  startvec <- c(Asym = 200, xmid = 725, scal = 350)
  nm1 <- nlme(circumference ~ SSlogis(age, Asym, xmid, scal),
                  data = Orange,
                  fixed = Asym + xmid + scal ~1,
                  random = Asym ~1,
              start = startvec)

  gls1 <-  gls(
          model = follicles ~ sin(2 * pi * Time) + cos(2 * pi * Time),
          data = Ovary,
          correlation = corAR1(form = ~ 1 | Mare)
          )

  save_file(lmm0, lmm0ML, lmm1, lmm2, nm1, gls1, pkg = pkg, type = "rda")
})

run_pkg("lme4", {
  lmm0 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  lmm0ML <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE)
  lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  lmm2 <- lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), sleepstudy)
  lmm1_prof <- profile(lmm1)
  save_file(lmm1_prof, lmm0, lmm0ML, lmm1, lmm2, pkg = pkg, type = "rda")
})


run_pkg("glmmTMB", {
  data(sleepstudy, package = "lme4")
  data(cbpp, package = "lme4")
  lmm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
  glmm1 <- glmmTMB(incidence / size ~ period + (1 | herd),
    data = cbpp, family = binomial, weights = size
  )

  data(Salamanders, package = "glmmTMB")
  zipm3 <- glmmTMB(count ~ spp * mined + (1 | site),
    zi = ~spp * mined, Salamanders, family = "poisson"
  )

  save_file(lmm1, glmm1, zipm3, pkg = pkg, type = "rda")
})

run_pkg("gamlss", {
  data(abdom, package="gamlss")
    gamlss1 <- gamlss(
        y ~ pb(x),
        sigma.fo = ~ pb(x),
        family = BCT,
        data = abdom,
        method = mixed(1, 20)
    )
  save_file(gamlss1, pkg=pkg, type="rds")
})

run_pkg("glmmADMB", {
  data(sleepstudy, package = "lme4")
  lmm1 <- glmmadmb(Reaction ~ Days + (Days | Subject), sleepstudy,
    family = "gaussian"
  )
  save_file(lmm1, pkg = pkg, type = "rda")
})

run_pkg("R2jags", {
    set.seed(123)
    library(R2jags)
    model.file <- system.file(package = "R2jags", "model", "schools.txt")
    ## data
    J <- 8.0
    y <- c(28.4, 7.9, -2.8, 6.8, -0.6, 0.6, 18.0, 12.2)
    sd <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)
    ## setting up model
    jags.data <- list("y", "sd", "J")
    jags.params <- c("mu", "sigma", "theta")
    jags.inits <- function() {
        list("mu" = rnorm(1), "sigma" = runif(1), "theta" = rnorm(J))
    }
    ## model fitting
    jagsfit <- R2jags::jags(
                           data = list("y", "sd", "J"),
                           inits = jags.inits,
                           jags.params,
                           n.iter = 10,
                           model.file = model.file
                       )
    save_file(jagsfit, pkg = pkg, type = "rds")
})

run_pkg("MCMCglmm", {
  data("sleepstudy", package = "lme4")
  mm0 <- MCMCglmm(Reaction ~ Days,
                  random = ~Subject, data = sleepstudy,
                  nitt=4000,
                  pr = TRUE
  )
  mm1 <- MCMCglmm(Reaction ~ Days,
    random = ~us(1 + Days):Subject,
    ## parameter-expanded priors
    ## V is 2x2 identity wlog
    ## t(2) with standard dev of 2
    prior = list(G = list(list(
      nu = 2, V = diag(2),
      alpha.mu = rep(0, 2),
      alpha.V = diag(rep(4, 2))
      ))),
    nitt=4000,
    data = sleepstudy, pr = TRUE
  )
  mm2 <- MCMCglmm(Reaction ~ Days, random = ~idh(1 + Days):Subject, data = sleepstudy, pr = TRUE)
  save_file(mm0, mm1, mm2, pkg = pkg, type = "rda")
})

if (FALSE) {
    ## skip for now.
    ## FIXME:
    ## (1) retrieving stored TMB objects doesn't seem to work
    ##  > Error in .Call("getParameterOrder", data, parameters, new.env(), PACKAGE = DLL) : 
    ## > "getParameterOrder" not available for .Call() for package "simple"
    ## (2) if there's a stored object we expect to be able to augment/glance it
    ##  (in test-alltibbles.R)
    run_pkg("TMB", {
        runExample("simple",thisR=TRUE)
        class(obj) <- "TMB"
        save_file(obj, pkg=pkg, type="rds")
    })
}


## slowest stuff last

run_pkg("rstanarm", {
  fit <- stan_glmer(mpg ~ wt + (1 | cyl) + (1 + wt | gear),
    data = mtcars,
    iter = 100, chains = 2
  )
  save_file(hack_size(fit), pkg = pkg, type = "rds")
})


## we really do need small iter/number of chains to keep
##  stored object size limited, *in addition to* hack_size()
##  (which gets rid of large environments associated with Stan
##   models)
if (run_brms) {
  ## ggeffects::efc, "Sample dataset from the EUROFAMCARE project"
  efc <- readRDS(system.file("extdata","efc.rds",package="broom.mixed"))
  run_pkg("brms", {
    brms_crossedRE <- brm(mpg ~ wt + (1 | cyl) + (1 + wt | gear),
      data = mtcars,
      iter = 100, chains = 2,
      family=gaussian
      )

    brms_crossedRE <- hack_size(brms_crossedRE)
    ## save_file(brms_crossedRE,  pkg = "brms", type = "rds")

    
    data(Salamanders, package="glmmTMB")
    brms_zip <- brm(bf(count ~ spp * mined + (1 | site),
                       zi = ~ mined + (1|site)),
                    data= Salamanders,
                    family = zero_inflated_poisson,
                    iter=100, chains=2,
                    save_dso=FALSE)

    brms_zip <- hack_size(brms_zip)

    f1 <- bf(neg_c_7 ~ e42dep + c12hour + c172code)
    f2 <- bf(c12hour ~ c172code)
    brms_multi <- brm(f1 + f2 + set_rescor(FALSE), data = efc,
             chains = 1, iter = 200,
             save_dso=FALSE)
    brms_multi <- hack_size(brms_multi)

    
    f3 <- bf(neg_c_7 ~ e42dep + c12hour + c172code + (1 |ID| e15relat))
    f4 <- bf(c12hour ~ c172code + (1 |ID| e15relat))
    b14 <- brm(f3 + f4 + set_rescor(FALSE), data = efc, iter = 500, chains = 1)
    brms_multi_RE <- hack_size(b14)
    
    data(sleepstudy, package="lme4")
    ss.tmp <- sleepstudy
    ss.tmp$Days_extra <- ss.tmp$Days # To check underscores in column names are handled properly
    brms_RE <- brm(Reaction ~ Days_extra + (Days_extra|Subject), 
                   data = ss.tmp, 
                   chains = 1,
                   iter=200)
    brms_RE <- hack_size(brms_RE)
    
    save_file(brms_crossedRE, brms_zip, brms_multi,
              brms_multi_RE, brms_RE, pkg = pkg, type = "rda")
  })
}

## put rstan **after** brms:
## https://github.com/paul-buerkner/brms/issues/558
if (run_stan) {
  run_pkg("rstan", {
    rstan_options(auto_write = TRUE)
    model_file <- system.file("extdata", "8schools.stan", package = "broom.mixed")
    schools_dat <- list(
      J = 8,
      y = c(28, 8, -3, 7, -1, 1, 18, 12),
      sigma = c(15, 10, 16, 11, 9, 11, 10, 18)
    )
    set.seed(2015)
    rstan_example <- stan(
      file = model_file, data = schools_dat,
      iter = 1000, chains = 2, save_dso = FALSE
    )
    rstan_example <- hack_size(rstan_example)
    save_file(rstan_example, pkg = pkg, type = "rds")
  })
}

