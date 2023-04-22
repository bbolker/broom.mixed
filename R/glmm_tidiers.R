##' @importFrom dplyr bind_rows tibble
##' @export
tidy.glmm <- function(x, effects = "fixed") {
    fix_nm <- names(coef(x))
    ran_nm <- x$varcomps.names
    res <- list()
    if ("fixed" %in% effects) {
        res[["fixed"]] <-
            tibble(
                term = fix_nm,
                estimate = coef(x),
                std.error = sqrt(diag(vcov(x)))[fix_nm]) |>
            mutate(statistic = estimate/std.error,
                   p.value = 2*pnorm(-abs(statistic)))
    }
    if ("ran_pars" %in% effects) {
        res[["fixed"]] <-
            tibble(
                term = ran_nm,
                estimate = x$nu,
                std.error = sqrt(diag(vcov(x)))[ran_nm]) |>
            mutate(statistic = NA_real_,
                   p.value = NA_real_)
    }
    bind_rows(res, .id = "effect")
}
