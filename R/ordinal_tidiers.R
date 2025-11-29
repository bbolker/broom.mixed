#' @export
predict.clmm <- function(object, ...) {
    ## hack clmm object so it looks sufficiently like a clm[m]2 object
    ##  for the predict.clm2 method to work ...
    object$location <- object$formula
    if (object$link == "logit") object$link <- "logistic"
    attr(object$location, "terms") <- object$terms
    class(object) <- c("clm2")
    predict(object, ...)
}

## predict values for every level in an ordinal response
## copied/modified from 
predict_all_clmm <- function(object, newdata, ...) {
    respvar <- attr(object$terms, "response")
    mf <- model.frame(object)
    nlev <- length(levels(mf[[respvar]]))
    if (!missing(newdata)) mf <- model.frame(object$formula, data = newdata)
    ndat <- do.call(rbind,
                    replicate(nlev, mf, simplify = FALSE))
    ndat[[respvar]] <- ordered(rep(seq(nlev), each = nrow(mf)))
    res <- matrix(predict(object, newdata = ndat), ncol=nlev)
}

#' name ordinal_tidiers
#' 
#' the \code{tidy} method for \code{clmm} objects (from the
#' \code{ordinal} package) lives in the \code{broom} package.
#'
#' @inheritParams lme4_tidiers
#' @importFrom tibble tibble
#' @importFrom stats model.frame
#' @export
augment.clmm <- function( x,
                         data = model.frame(x),
                         newdata,
                         ...) {
    
    if (!missing(newdata)) data <- newdata

    ## STUB
    ## call predict_all_clmm
    ## generate mean predictions via
    ##  sweep([pred matrix], MARGIN = 2, FUN = "*", STATS = seq(ncol[pred matrix])) %>% rowMeans
    ## get std dev similarly (sqrt(rowMeans(sweep with (1:n)^2)))
    ## residuals, pearson residuals?
}

## simulate method will look like this:
## pred_matrix <- predict_all_clmm()
## apply(pred_matrix, 1, \(p) rmultinom(1, size = 1, prob = p))

if (FALSE) {
    library(ordinal)

    fmm1 <- clmm(rating ~ temp + contact + (1|judge), data = wine)
    fmm2 <- clmm2(rating ~ temp + contact, random = judge, data = wine)

    
    mm <- predict_all_clmm(fmm1)
    stopifnot(all.equal(predict(fmm1), predict(fmm2),
                        tolerance  = 1e-6))

}
