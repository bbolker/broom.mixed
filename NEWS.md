# CHANGES IN broom.mixed VERSION 0.2.9.5

## NEW FEATURES

- `stanreg` tidier gains `exponentiate` argument (wish of GH #122)
- `tidy.brmsfit` gains optional `rhat` and `ess` columns (Alexey Stukalov)
- experimental support for `lqmm` models (David Luke Thiessen)

## BUG FIXES

- bug fixed for `glmmTMB` tidying with `conf.int=TRUE`, random effects in multiple model components, subset of components requested in tidy output (GH #136, Daniel Sjoberg)
- `tidy.brmsfit` works better for models with no random/group-level effects (Matthieu Bruneaux)

## USER-VISIBLE CHANGES

- `as.data.frame.ranef.lme` now processes the optional argument (see `?as.data.frame)`, so that `data.frame(ranef_object)` works

- `stanreg` tidier now checks for spurious values in `...`

# CHANGES IN broom.mixed VERSION 0.2.9.4 (2022-03-28)

- minor changes only; test tweaks for CRAN compatibility

# CHANGES IN broom.mixed VERSION 0.2.9.3 (2021-07-07)

## BUG FIXES

- improved profile robustness in `TMB` tidiers

## NEW FEATURES

- `lme` tidier gets functionality for information about variance models (use `effects = "var_model"`) (Bill Denney)

- support for models with fixed sigma values in `lme` tidier (Bill Denney)

- added `tidy` and `glance` methods for `allFit` objects from the `lme4` package

- `get_methods()` function returns a table of all available `tidy`/`glance`/`augment` methods

## USER-VISIBLE CHANGES

- improved lme tidying for random effects values

- brms tidiers no longer use deprecated `posterior_samples`

- `glance.lme4` now returns nobs (Cory Brunson)

- some tidiers are less permissive about unused arguments passed via `...`

# CHANGES IN broom.mixed VERSION 0.2.7 (2021-07-07)

## NEW FEATURES

- experimental `TMB` tidiers (the TMB package does not return an object of class TMB, so users should run `class(fit) <- "TMB"` before tidying)

## USER-VISIBLE CHANGES

- term names are no longer "sanitized" in `gamlss` tidiers (e.g. "(Intercept)" is not converted to "X.Intercept.")

- gamlss `glance` method returns `nobs` (GH #113)

## BUG FIXES

- Wald confidence intervals for `lmerTest` models now respect `ddf.method`

- `tidy.glmmTMB(.,effects="ran_vals")` fixed for `stringsAsFactors` changes in glmmTMB (GH #103)

- `tidy.gamlss` now works in a wider range of cases (GH #74)

- `tidy.brmsfit` works for models without group effects (GH #100)

# CHANGES IN broom.mixed VERSION 0.2.6 (2020-05-17)

- No improvements; compatibility with `dplyr` 1.0.0; skip examples

# CHANGES IN broom.mixed VERSION 0.2.5 (2020-04-19)

## NEW FEATURES

- `lmer` tidier gets `ddf.method` (applies only to `lmerTest` fits)

- `glmmTMB` gets `exponentiate` options

- experimental `GLMMadaptive` tidiers

## OTHER CHANGES

- fixes for updates in `tibble` package

# CHANGES IN broom.mixed VERSION 0.2.4

## NEW FEATURES

- `gls` tidier gets `confint` (GH #49)

## USER-VISIBLE CHANGES

- redundant `estimate.method` in MCMC tidiers goes away; use `robust` to compute point estimates/uncertainty via median and MAD rather than mean and SE

## BUG FIXES

- misc fixes: lme4 tidiers (confint for `ran_vals`, profile conf intervals fixed), R2jags, gamlss ...

- `ran_vals` works for glmmTMB

# CHANGES IN broom.mixed VERSION 0.2.3

## BUG FIXES

- don't ignore `conf.level` in `tidy.(merMod|glmmTMB)` (GH #30,31: @strengejacke)

- levels correct in `tidy.brmsfit` (GH #36: @strengejacke)

- `component` argument works for random effects in `glmmTMB` (GH #33: @strengejacke)

## NEW FEATURES

- `brmsfit` and `rstanarm` methods allow `conf.method="HPDinterval"`

## USER-VISIBLE CHANGES

- `tidy.brmsfit` gets component column (GH #35: @strengejacke), response column for multi-response models (GH #34: @strengejacke)

- component tags are stripped from tidied `brmsfit` objects

- "Intercept" terms in `brms` fits are re-coded as "(Intercept)" by default, for dotwhisker/cross-model compatibility; for previous behaviour, specify `fix.intercept=FALSE`

# CHANGES IN broom.mixed VERSION 0.2.2

- modify examples, for CRAN compliance

# CHANGES IN broom.mixed VERSION 0.2.1

- reduced size of stored fits for examples, for CRAN compliance

# CHANGES IN broom.mixed VERSION 0.2.0

## NEW FEATURES

- more consistent term names in `brmsfit`, `rstanreg` tidiers

- improved `tidy.MCMCglmm`

## USER-VISIBLE CHANGES

- all methods return tibbles (`tbl_df`) rather than data frames

- the value of the group variable for fixed-effect parameters has changed from `"fixed"` to `NA`

- `brmsfit` and `rstanarm` tidiers are more consistent with other tidiers (e.g. the argument for setting confidence level is `conf.level` rather than `prob`)

# CHANGES IN broom.mixed VERSION 0.0.1

## BUG FIXES

- Sorted out some of the confusion over random effect naming: `"ran_vals"` extracts conditional modes/BLUPs/varying parameters (deviations from population-level estimates), while `"ran_coefs"` extracts group-level estimates

## NEW FEATURES

- improved `nlme` tidiers

- improved `glmmTMB` tidiers (can handle some zero-inflation parameters)

- `lme4` tidiers now optionally take a pre-computed profile argument when using `conf.method="profile"`

## USER-VISIBLE CHANGES

- The default behaviour of most mixed-model tidiers has changed/will gradually be changed to the following (description modified from TJ Mahr at <https://github.com/tidymodels/broom/issues/96>):
   -  Random effect variances and covariances can now be extracted.  effects = "ran_pars" returns the standard deviations/correlations of random effects (if `scales="sdcor"` [default]) or their variances and covariances (if `scales = "varcov"`)
   - Random effects estimates are now extracted with `effects = "ran_coefs"` for the group-level estimates (previously these effects were extracted with `tidy(model, "random")`) or `effects = "ran_vals"` for the conditional modes (deviations of the group-level parameters from the population-level estimates)
   - `effects` can take a vector of values (those listed above, plus "fixed" for fixed effects). The default value is effects = c("ran_pars", "fixed") which extracts random effect variances/covariances and fixed effect estimates.

- term names for random-effect parameters no longer contain a (redundant) `group` specifier (at least for `lme4 models`); use something like `tidyr::unite(term,term,group,sep=".")` to collapse the two columns

