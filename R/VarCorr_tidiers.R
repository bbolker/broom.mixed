## https://github.com/bbolker/broom.mixed/issues/152, @phargarten2

#' Tidying VarCorr of Mixed Effect Models
#' @param object a VarCorr object
#' @param scales see \code{\link{tidy.lme4}}
#' @importFrom dplyr as_tibble
#' @export
tidy.VarCorr.lme <- function(
                             object,    
                             ## effects = c("ran_pars", "fixed"),  #effects are always "ran_pars"
                             scales = c("sdcor", "vcov"),
                             conf.int = FALSE, 
                             conf.level = 0.95, 
                             ...) {

    vcov <- var1 <- NULL  ## NSE/R CMD check
    check_dots(...)
    scales <- match.arg(scales)
    
    if(inherits(object, "VarCorr.lme")) {
        ## Need to convert nlme object to a proper tibble
        object <- convert_VarCorr.MerMod_to_VarCorr.lme(object)
    } 
    
    if (conf.int) {
        cli::cli_alert_info("Can't compute confidence intervals for random effect parameters in a crossed model.")
    }

    if (scales == "vcov") {
        ests_random <- as_tibble(object) |>
            rename(group = grp, term = var1, estimate_var = vcov) |>
            mutate(
                term = case_when(
                    !is.na(term) & !is.na(var2) ~ paste0("cov_", term, "_", var2),
                    group == "Residual" & is.na(term) ~ paste0("var_", "Observation"),
                    TRUE ~ paste0("var_", term)
                )
            ) |>
            mutate(effect = "ran_pars", .before = "group")
        
        total_var <- (ests_random
            |> filter(is.na(var2))
            |> summarize(s = sum(estimate_var))
            |> pull(s)
        )
        
        ests_random <- ests_random |>
            mutate(prop_var = case_when(
                       !stringr::str_detect(term, "cov_") ~ estimate_var/total_var,
                       TRUE ~ NA_real_
                   )
                   ) |>
            select(-var2, -sdcor)
        
        ##Check
        ## print(as_tibble(object))
        
    } else if(scales == "sdcor") {
        ests_random <- as_tibble(object) |>
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
            mutate(prop_var = case_when(
                       !stringr::str_detect(term, "cor_") ~ estimate^2/total_var,
                       TRUE ~ NA_real_
                   )
                   ) |>
            select(-var2, -vcov)
        
    }
    
    return(ests_random)
}

                                        #  lmm.lme4 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
                                        #as.tibble(lme4::VarCorr(lmm.lme4))
                                        # grp      var1        var2    vcov   sdcor
                                        # <chr>    <chr>       <chr>  <dbl>   <dbl>
                                        # 1 Subject  (Intercept) NA    612.   24.7   
                                        # 2 Subject  Days        NA     35.1   5.92  
                                        # 3 Subject  (Intercept) Days    9.60  0.0656
                                        # 4 Residual NA          NA    655.   25.6  


#' This function converts VarCorr on a nlme object to a tibble form of an lme4 object
#' i.e. VarCorr.merMod -> tibble of VarCorr.lme
#' @noRd
#' @param object  A variance-correlation component of nlme::lme object. 
#' @return The tibble version of VarCorr.lme object from lme4 package.
convert_VarCorr.MerMod_to_VarCorr.lme <- function(object){
    A <- as.matrix(object)
    row.residual <- stringr::str_which("Residual", rownames(A))
    
    ##Unsure how to generalize this
    corr.A <- A[-row.residual, -(1:2), drop = TRUE]
    t.corr <- tibble(var = names(corr.A), corr = as.vector(corr.A))
    corr <- tibble(
        var1 = t.corr[1, "var", drop = TRUE],
        var2 = t.corr[2, "var", drop = TRUE],
        sdcor = t.corr[2, "corr", drop = TRUE]
    )
    
    tib <-  tibble(grp = NA_character_, var1 = rownames(A), var2 = NA, vcov= A[ ,"Variance"], sdcor = A[ , "StdDev"]) |>
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

