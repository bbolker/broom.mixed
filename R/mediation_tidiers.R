#' Tidying methods for mediation analyses involving mixed effects models
#'
#' These methods tidy the coefficients of \code{mediation::mediate} output
#' (i.e., \code{mediate.mer} objects) when \code{lme4::lmer} and
#' \code{lme4::glmer} models (i.e., \code{merMod} objects) provide the input.
#'
#' @aliases tidy.mediate.mer
#' @param x an object of class \code{mediate.mer}, as from \code{mediate} using
#'   \code{lmer}, \code{glmer}, or \code{nlmer} models
#' @param \dots additional arguments (unused: for generic consistency)
#'
#' @return All tidying methods return a \code{data.frame} without rownames. The
#'   structure depends on the method chosen.
#'
#' @name mediate_tidiers
#'
#' @examples
#' if (require("lme4") && require("mediation")) {
#'     ## Borrowed from \code{help(mediation::mediate)}:
#'     ## Varying intercept for mediator 
#'     mod_m <- glmer(job_dich ~ treat + econ_hard + (1 | educ), 
#'                    family = binomial(link = "probit"), data = jobs)
#'     ## Varying intercept and slope for outcome
#'     mod_y <- glmer(work1 ~ treat + job_dich + econ_hard + (1 + treat | occp),
#'                    family = binomial(link = "probit"), data = jobs)
#'     ## Output based on mediator group ("educ")
#'     mod_med <- mediate(mod_m, mod_y, treat = "treat", 
#'                        mediator = "job_dich", sims=50, group.out="educ")
#'     ## Tidy outputs
#'     tidy(mod_m)
#'     tidy(mod_y)
#'     tidy(mod_med)
#' }
NULL

#' @rdname mediate_tidiers
#'
#' @param conf.int whether to include a confidence interval
#' @param conf.level confidence level for CI
#'
#' @return \code{tidy} returns one row for each estimated effect:
#' first the mediated effect in the control and treatment groups, respectively,
#' then the direct effect in each group.
#' It contains the columns
#'   \item{term}{term being estimated}
#'   \item{estimate}{estimated coefficient}
#'   \item{std.error}{standard error}
#'   \item{p.value}{P-value computed from t-statistic (may be missing/NA)}
#'   

#' @export
#' @seealso \code{\link[mediation]{mediate}}, \code{\link[broom]{tidy.mediate}}
tidy.mediate.mer <- function(x, conf.int = FALSE, conf.level = .95, ...) {
  
  # extract model elements as is `broom:::tidy.mediate()`
  d0 <- d1 <- z0 <- z1 <- d0.sims <- d1.sims <- z0.sims <- NULL
  z1.sims <- d0.p <- d1.p <- z0.p <- z1.p <- NULL
  s <- base::summary(x)
  nn <- c("term", "estimate", "std.error", "p.value", "conf.low", "conf.high")
  sims <- s$sims
  ci <- NULL
  co <- with(
    s,
    tibble(
      c("acme_0", "acme_1", "ade_0", "ade_1"),
      c(d0, d1, z0, z1),
      c(sd(d0.sims), sd(d1.sims), sd(z0.sims), sd(z1.sims)),
      c(d0.p, d1.p, z0.p, z1.p)
    )
  )
  
  if (conf.int) {
    low <- (1 - conf.level) / 2
    high <- 1 - low
    BC.CI <- function(theta) {
      z.inv <- length(theta[theta < mean(theta)]) / sims
      z <- qnorm(z.inv)
      U <- (sims - 1) * (mean(theta) - theta)
      top <- sum(U^3)
      under <- 6 * (sum(U^2))^{
        3 / 2
      }
      a <- top / under
      lower.inv <- pnorm(z + (z + qnorm(low)) / (1 - a * (z + qnorm(low))))
      lower2 <- lower <- quantile(theta, lower.inv)
      upper.inv <- pnorm(z + (z + qnorm(high)) / (1 - a * (z + qnorm(high))))
      upper2 <- upper <- quantile(theta, upper.inv)
      return(c(lower, upper))
    }
    ci <- with(
      x,
      sapply(list(d0.sims, d1.sims, z0.sims, z1.sims), function(x) apply(x, 1, BC.CI))
    )
    if (s$boot.ci.type != "bca") {
      CI <- function(theta) {
        return(quantile(theta, c(low, high), na.rm = TRUE))
      }
      ci <- with(
        x,
        sapply(list(d0.sims, d1.sims, z0.sims, z1.sims), function(x) apply(x, 1, CI))
      )
    }
    
    co <- cbind(co, t(ci))
  }
  
  # format tibble
  names(co) <- nn[1:ncol(co)]
  co <- reorder_cols(co)
  return(co)
}
