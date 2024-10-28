## https://github.com/bbolker/broom.mixed/issues/152, @phargarten2

#' Tidying VarCorr of Mixed Effect Models
#' @param x a VarCorr object
#' @param scales see \code{\link{tidy.merMod}}
#' @inheritParams tidy.merMod
#' @importFrom dplyr as_tibble
#' @export
tidy.VarCorr.lme <- function(
                             x,    
                             ## effects = c("ran_pars", "fixed"),  #effects are always "ran_pars"
                             scales = c("sdcor", "vcov"),
                             conf.int = FALSE, 
                             conf.level = 0.95, 
                             ...) {

    grp <- vcov <- var1 <- var2 <- sdcor <- estimate <- s <- NULL  ## NSE/R CMD check
    check_dots(...)
    scales <- match.arg(scales)

    if(inherits(x, "VarCorr.lme")) {
        ## Need to convert nlme object to a proper tibble
        x <- convert_VarCorr.lme(x)
    } 
    
    if (conf.int) {
        message("Can't compute confidence intervals from VarCorr")
    }

    if (scales == "vcov") {
        ests_random <- as_tibble(x) |>
            rename(group = grp, term = var1, estimate = vcov) |>
            mutate(
                term = case_when(
                    !is.na(term) & !is.na(var2) ~ paste0("cov_", term, "_", var2),
                    group == "Residual" & is.na(term) ~ paste0("var_", "Observation"),
                    TRUE ~ paste0("var_", term)
                )
            ) |>
            mutate(effect = "ran_pars", .before = "group")
        
        ests_random <- ests_random |>
            select(-var2, -sdcor)
        
    } else if(scales == "sdcor") {
        ests_random <- as_tibble(x) |>
            rename(group = grp, term = var1, estimate = sdcor) |>
            mutate(
                term = case_when(
                    !is.na(term) & !is.na(var2) ~ paste0("cor_", term, "_", var2),
                    group == "Residual" & is.na(term) ~ paste0("sd_", "Observation"),
                    TRUE ~ paste0("sd_", term)
                )
            ) |>
            mutate(effect = "ran_pars", .before = "group")
        
        ##Can estimate covariance if requested with a bit more work...
        
        total_var <- ests_random |> filter(is.na(var2)) |> summarize(s = sum(estimate^2)) |> pull(s)
        
        ests_random <- ests_random |>
            ## mutate(prop_var = case_when(
            ##            !stringr::str_detect(term, "cor_") ~ estimate^2/total_var,
            ##            TRUE ~ NA_real_
            ##        )
            ##        ) |>
            select(-var2, -vcov)
        
    }
    
    return(ests_random)
}

#' This function converts VarCorr on a nlme x to a tibble
#' i.e. VarCorr.merMod -> tibble of VarCorr.lme
#' @noRd
#' @param x  A variance-correlation component of nlme::lme object. 
#' @return A useful (?) tibble
convert_VarCorr.lme <- function(x) {
    
    sdcor <- NULL
    
    A <- as.matrix(x)
    row.residual <- stringr::str_which("Residual", rownames(A))
    
    ##Unsure how to generalize this
    corr.A <- A[-row.residual, "Corr", drop = TRUE]
    t.corr <- tibble(var = names(corr.A), corr = as.vector(corr.A))
    corr <- tibble(
        var1 = t.corr[1, "var", drop = TRUE],
        var2 = t.corr[2, "var", drop = TRUE],
        sdcor = t.corr[2, "corr", drop = TRUE]
    )
    
    tib <-  tibble(grp = NA_character_, var1 = rownames(A), var2 = NA,
                   vcov= A[ ,"Variance"], sdcor = A[ , "StdDev"]) |>
        bind_rows(corr) |>
        mutate(
            grp = case_when(
                var1 != "Residual" ~ "Subject",
                var1 == "Residual" ~ "Residual",
                TRUE ~ NA_character_
            ),
            vcov = as.numeric(vcov),
            sdcor = as.numeric(sdcor)
        )
    
    return(tib)
}

