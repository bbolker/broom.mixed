if (requireNamespace("lme4", quietly = TRUE) &&
    requireNamespace("nlme", quietly=TRUE)) {
    library(lme4)
    library(nlme)
    devtools::load_all()
    data("sleepstudy", package="lme4")
    lmm.nlme <- lme(Reaction ~ Days, random=~ Days|Subject, sleepstudy)
                                        # > VarCorr(lmm.nlme)
                                        # Subject = pdLogChol(Days) 
                                        # Variance StdDev    Corr  
                                        # (Intercept) 612.0795 24.740241 (Intr)
                                        # Days         35.0713  5.922103 0.066 
                                        # Residual    654.9424 25.591843       

    t1 <- tidy(VarCorr(lmm.nlme), scales = "vcov")
    t2 <- tidy(VarCorr(lmm.nlme), scales = "sdcor")
    tidy(lmm.nlme)

     sleepstudy2 <- transform(sleepstudy, Group = rep(LETTERS[1:10], each = 18))
    
    lmm2 <- nlme::lme(Reaction ~ Days,
                      random = list(Subject = ~ Days, Group = ~ 1),
                      data = sleepstudy2
                      )
    tidy(lmm2)
    
}
