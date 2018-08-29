verbose <- FALSE
stopifnot(
  require("testthat"), require("broom.mixed"), require("stringr"),
  require("lme4"), require("glmmTMB")
)

## examples to SKIP
skip_list <- list(nlme_example.rda = c("lmm2"))

## need data available for augment() methods ...
data("sleepstudy", package = "lme4")
data("Salamanders", package = "glmmTMB")
ex <- list.files(system.file("extdata", package = "broom.mixed"),
  pattern = "\\.rd"
)


## test tidy, augment, glance methods from lme4-tidiers.R

for (e in ex) {
    p <- stringr::str_extract(e, "^[^_]+")
    if (require(p,character.only=TRUE)) {
        f <- system.file("extdata", e, package = "broom.mixed")
        fn <- stringr::str_extract(f, "[^/]+$")
        if (verbose) cat(fn, "\n")
        if (grepl("\\.rds", e)) {
            x <- list()
            x[[1]] <- readRDS(f)
        } else {
            ## rda file
            L <- load(f)
            x <- mget(setdiff(L, skip_list[[fn]]))
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
