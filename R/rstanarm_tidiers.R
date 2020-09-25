#' Tidying methods for an rstanarm model
#'
#' These methods tidy the estimates from \code{\link[rstanarm]{stanreg-objects}}
#' (fitted model objects from the \pkg{rstanarm} package) into a summary.
#'
#'
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#'
#' @seealso \code{\link[rstan]{summary,stanfit-method}}
#'
#' @name rstanarm_tidiers
#'
#' @param x Fitted model object from the \pkg{rstanarm} package. See
#'   \code{\link[rstanarm]{stanreg-objects}}.
#' @examples
#'
#' if (require("rstanarm")) {
#' \dontrun{
#' #'     ## original model
#'     fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
#'                       iter = 300, chains = 2)
#'   }
#' ## load example data
#' fit <- readRDS(system.file("extdata", "rstanarm_example.rds", package="broom.mixed"))
#'
#'   # non-varying ("population") parameters
#'   tidy(fit, conf.int = TRUE, prob = 0.5)
#'   tidy(fit, conf.int = TRUE, conf.method = "HPDinterval", prob = 0.5)
#'
#'   # hierarchical sd & correlation parameters
#'   tidy(fit, effects = "ran_pars")
#'
#'   # group-specific deviations from "population" parameters
#'   tidy(fit, effects = "ran_vals")
#'
#'   # glance method
#'    glance(fit)
#'   \dontrun{
#'      glance(fit, looic = TRUE, cores = 1)
#'   }
#' } ## if require("rstanarm")
NULL

#' @rdname rstanarm_tidiers
#' @inheritParams brms_tidiers
#' @param conf.level See \code{\link[rstanarm]{posterior_interval.stanreg}}.
#' @param conf.int If \code{TRUE} columns for the lower (\code{conf.low}) and upper (\code{conf.high}) bounds of the
#'   \code{100*prob}\% posterior uncertainty intervals are included. See
#'   \code{\link[rstanarm]{posterior_interval.stanreg}} for details.
#'
#' @return
#' When \code{effects="fixed"} (the default), \code{tidy.stanreg} returns
#' one row for each coefficient, with three columns:
#' \item{term}{The name of the corresponding term in the model.}
#' \item{estimate}{A point estimate of the coefficient (posterior median).}
#' \item{std.error}{A standard error for the point estimate based on
#' \code{\link[stats]{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link[rstanarm]{print.stanreg}} for more details.}
#'
#' For models with group-specific parameters (e.g., models fit with
#' \code{\link[rstanarm]{stan_glmer}}), setting \code{effects="ran_vals"}
#' selects the group-level parameters instead of the non-varying regression
#' coefficients. Addtional columns are added indicating the \code{level} and
#' \code{group}. Specifying \code{effects="ran_pars"} selects the
#' standard deviations and (for certain models) correlations of the group-level
#' parameters.
#'
#' Setting \code{effects="auxiliary"} will select parameters other than those
#' included by the other options. The particular parameters depend on which
#' \pkg{rstanarm} modeling function was used to fit the model. For example, for
#' models fit using \code{\link[rstanarm]{stan_glm}} the overdispersion
#' parameter is included if \code{effects="aux"}, for
#' \code{\link[rstanarm]{stan_lm}} the auxiliary parameters include the residual
#' SD, R^2, and log(fit_ratio), etc.
#'
#' @export
tidy.stanreg <- function(x,
                         effects = "fixed",
                         conf.int = FALSE,
                         conf.level = 0.9,
                         conf.method=c("quantile","HPDinterval"),
                         ...) {
    conf.method <- match.arg(conf.method)
    effects <-
        match.arg(effects,
                  several.ok = TRUE,
                  choices = c(
                      "fixed", "ran_vals",
                      "ran_pars", "auxiliary"
                  )
                  )
    if (any(effects %in% c("ran_vals", "ran_pars"))) {
        if (!inherits(x, "lmerMod")) {
            stop("Model does not have varying ('ran_vals') or hierarchical ('ran_pars') effects.")
        }
    }

    nn <- c("estimate", "std.error")
    ret_list <- list()
    if ("fixed" %in% effects) {
        nv_pars <- names(rstanarm::fixef(x))
        ret <- cbind(
            rstanarm::fixef(x),
            rstanarm::se(x)[nv_pars]
        )

        if (inherits(x, "polr")) {
            ## also include cutpoints
            cp <- x$zeta
            se_cp <- apply(as.matrix(x, pars = names(cp)), 2, stats::mad)
            ret <- rbind(ret, cbind(cp, se_cp))
            nv_pars <- c(nv_pars, names(cp))
        }

        if (conf.int) {

            cifix <- switch(conf.method,
                            HPDinterval= {
                                m <- as.matrix(x$stanfit)
                                m <- m[,colnames(m) %in% nv_pars]
                                coda::HPDinterval(coda::as.mcmc(m),
                                                  prob=conf.level)
                            },
                            quantile=rstanarm::posterior_interval(
                                                   object = x,
                                                   pars = nv_pars,
                                                   prob = conf.level
                                               )
                            ) ## cifix
            ret <- data.frame(ret, cifix)
            nn <- c(nn, "conf.low", "conf.high")
        }
        ret_list$non_ran_vals <- fix_data_frame(ret, newnames = nn, newcol="term")
    }
    if ("auxiliary" %in% effects) {
        nn <- c("estimate", "std.error")
        parnames <- rownames(x$stan_summary)
        auxpars <- c(
            "sigma", "shape", "overdispersion", "R2", "log-fit_ratio",
            grep("mean_PPD", parnames, value = TRUE)
        )
        auxpars <- auxpars[which(auxpars %in% parnames)]
        ret <- summary(x, pars = auxpars)[, c("50%", "sd"), drop = FALSE]
        if (conf.int) {
            ints <- rstanarm::posterior_interval(x, pars = auxpars, prob = conf.level)
            ret <- data.frame(ret, ints)
            nn <- c(nn, "conf.low", "conf.high")
        }
        ret_list$auxiliary <-
            fix_data_frame(ret, newnames = nn, newcol="term")
    }
    if ("ran_pars" %in% effects) {
        ret <- (rstanarm::VarCorr(x)
            %>% as.data.frame()
            %>% mutate_if(is.factor,as.character)
        )
        rscale <- "sdcor" # FIXME
        ran_prefix <- c("sd", "cor") # FIXME
        pfun <- function(x) {
            v <- na.omit(unlist(x))
            if (length(v) == 0) v <- "Observation"
            p <- paste(v, collapse = ".")
            if (!identical(ran_prefix, NA)) {
                p <- paste(ran_prefix[length(v)], p, sep = "_")
            }
            return(p)
        }

        rownames(ret) <- paste(apply(ret[c("var1", "var2")], 1, pfun),
                               ret[, "grp"],
                               sep = "."
                               )
        ret_list$hierarchical <- fix_data_frame(ret[c("grp", rscale)],
                                                 newcol="term",
                                                newnames = c("group", "estimate"))
    }

    if ("ran_vals" %in% effects) {
        nn <- c("estimate", "std.error")
        s <- summary(x, pars = "varying") ## goes through to rstanarm
        ret <- cbind(s[, "50%"], rstanarm::se(x)[rownames(s)])

        if (conf.int) {
            ciran <- rstanarm::posterior_interval(x,
                                                  regex_pars = "^b\\[",
                                                  prob = conf.level
                                                  )
            ret <- data.frame(ret, ciran)
            nn <- c(nn, "conf.low", "conf.high")
        }

        double_splitter <- function(x, split1, sel1, split2, sel2) {
            y <- unlist(lapply(strsplit(x, split = split1, fixed = TRUE), "[[", sel1))
            unlist(lapply(strsplit(y, split = split2, fixed = TRUE), "[[", sel2))
        }
        vv <- fix_data_frame(ret, newnames = nn, newcol="term")
        nn <- c("level", "group", "term", nn)
        nms <- vv$term
        vv$term <- NULL
        lev <- double_splitter(nms, ":", 2, "]", 1)
        grp <- double_splitter(nms, " ", 2, ":", 1)
        trm <- double_splitter(nms, " ", 1, "[", 2)
        vv <- data.frame(lev, grp, trm, vv)
        ret_list$ran_vals <- fix_data_frame(vv, newnames = nn, newcol="term")
    }

    return(dplyr::bind_rows(ret_list))
}


#' @rdname rstanarm_tidiers
#'
#' @param ... For \code{glance}, if \code{looic=TRUE}, optional arguments to
#'   \code{\link[rstanarm]{loo.stanreg}}.
#' @return \code{glance} returns one row with the columns
#'   \item{algorithm}{The algorithm used to fit the model.}
#'   \item{pss}{The posterior sample size (except for models fit using
#'   optimization).}
#'   \item{nobs}{The number of observations used to fit the model.}
#'   \item{sigma}{The square root of the estimated residual variance, if
#'   applicable. If not applicable (e.g., for binomial GLMs), \code{sigma} will
#'   be given the value \code{1} in the returned object.}
#'
#'   If \code{looic=TRUE}, then the following additional columns are also
#'   included:
#'   \item{looic}{The LOO Information Criterion.}
#'   \item{elpd_loo}{The expected log predictive density (\code{elpd_loo = -2 *
#'   looic}).}
#'   \item{p_loo}{The effective number of parameters.}
#'
#' @export
glance.stanreg <- function(x, looic = FALSE, ...) {
    glance_stan(x, looic = looic, type = "stanreg", ...)
}


glance_stan <- function(x, looic = FALSE, ..., type) {
    sigma <- if (getRversion() >= "3.3.0") {
                 get("sigma", asNamespace("stats"))
             } else {
                 ## FIXME: could fail if old R & called from brms
                 ## & rstanarm not installed ...
                 get("sigma", asNamespace("rstanarm"))
             }
    if (type == "stanreg") {
        algo <- x$algorithm
        sim <- x$stanfit@sim
    } else {
        ## method is recorded for every chain; pick the first
        algo <- x$fit@stan_args[[1]][["method"]]
        sim <- x$fit@sim
    }

    ret <- dplyr::tibble(algorithm = algo)

    if (algo != "optimizing") {
        pss <- sim$n_save
        if (algo == "sampling") {
            pss <- sum(pss - sim$warmup2)
        }
        ret <- dplyr::mutate(ret, pss = pss)
    }

    ret <- mutate(ret, nobs = stats::nobs(x))
    if (length(sx <- sigma(x)) > 0) {
        ret <- dplyr::mutate(ret, sigma = sx)
    }
    if (looic) {
        if (algo == "sampling") {
            if (type == "stanreg") {
                loo1 <- rstanarm::loo(x, ...)
            } else {
                loo1 <- brms::loo(x, ...)
            }
            loo1_est <- loo1[["estimates"]]
            ret <- data.frame(
                ret,
                rbind(loo1_est[
                    c("looic", "elpd_loo", "p_loo"),
                    "Estimate"
                ])
            )
        } else {
            message("looic only available for models fit using MCMC")
        }
    }
    dplyr::as_tibble(unrowname(ret))
}
