context("test all tidy methods return tibbles")
verbose <- FALSE

pkgs_avail <- vapply(c("testthat", "broom.mixed", "stringr", "lme4", "glmmTMB"), require,
                     quietly = TRUE,
                     character.only = TRUE,
                     FUN.VALUE = logical(1))

## FIXME: _could_ copy data objects over if we wanted to check more robustly (without requiring lme4, glmmTMB)

## harmless if rstan not available; if it is available, should (?)
##  resolve issues with brms use of rstan methods ...
## see thread at https://stat.ethz.ch/pipermail/r-package-devel/2024q3/011097.html
## for gory details
require("rstan")

if (all(pkgs_avail)) {

## objects to SKIP (from RDA files)
skip_objects <- list(nlme_example.rda = c("lmm2"),
                  ## glance method not working for this case
                  brms_example.rda = c("brms_multi","brms_multi_RE")
                  )

skip_files <- c("efc.rds")

  ## need data available for augment() methods ...

  data("sleepstudy", package = "lme4")
  data("Salamanders", package = "glmmTMB")
  efc <- readRDS(system.file("extdata", "efc.rds", package="broom.mixed"))
  ex <- list.files(system.file("extdata", package = "broom.mixed"),
                   pattern = "\\.rd"
                   )
  ex <- setdiff(ex,skip_files)


## test tidy, augment, glance methods from lme4-tidiers.R

for (e in ex) {
    p <- stringr::str_extract(e, "^[^_]+")
    ## cat(": ",p,"\n")
    if (require(p, quietly = TRUE, character.only=TRUE)) {
        f <- system.file("extdata", e, package = "broom.mixed")
        fn <- stringr::str_extract(f, "[^/]+$")
        if (verbose) cat(fn, "\n")
        if (grepl("\\.rds", e)) {
            x <- list()
            x[[1]] <- readRDS(f)
        } else {
            ## rda file
            L <- load(f)
            x <- mget(setdiff(L, skip_objects[[fn]]))
            if (p == "glmmTMB") {
               x <- lapply(x, glmmTMB::up2date)
            }
        }
        for (z in x) {
            testf <- function(fn_str, obj) {
                cc <- class(obj)[1]
                if (sum(grepl(cc, methods(fn_str))) > 0) {
                    if (verbose) cat(sprintf("found method %s for %s\n", fn_str, cc))
                    return(expect_is(get(fn_str)(obj), "tbl_df"))
                } else {
                    return(TRUE)
                }
            }
            testf("glance", z)
            testf("augment", z)
            testf("tidy", z)
        } ## loop over objects
    } ## if package available
}
} ## if all pkgs available
