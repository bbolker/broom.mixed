## most of these are unexported (small) functions from broom;
## could be removed if these were exported

#' check if a package is available and return informative message otherwise
#'
#' @keywords internal
assert_dependency <- function(library_name) {
  if (!requireNamespace(library_name, quietly = TRUE)) {
    stop(sprintf("Please install the %s package.", library_name))
  }
}


## https://github.com/klutometis/roxygen/issues/409
#' @importFrom broom tidy glance augment
#' @export
broom::tidy
#' @export
broom::glance
#' @export
broom::augment
#'
#' strip rownames from an object
#'
#' @param x a data frame
unrowname <- function(x) {
  rownames(x) <- NULL
  return(x)
}

## first convert to data frame, then add rownames, then tibble
tibblify <- function(x, var = "term") {
  if (is.null(var)) {
    return(dplyr::as_tibble(unrowname(x)))
  }
  ret <- (x
  %>%
    as.data.frame()
    %>%
    tibble::rownames_to_column(var)
    %>%
    dplyr::as_tibble())
  return(ret)
}

#' Remove NULL items in a vector or list
#'
#' @param x a vector or list
compact <- function(x) Filter(Negate(is.null), x)

#' insert a row of NAs into a data frame wherever another data frame has NAs
#'
#' @param x data frame that has one row for each non-NA row in original
#' @param original data frame with NAs
insert_NAs <- function(x, original) {
  indices <- rep(NA, nrow(original))
  indices[which(stats::complete.cases(original))] <- seq_len(nrow(x))
  x[indices, ]
}

## list of regex matches for mixed-effect columns -> broom names
col_matches <- list(
  estimate = "^(Estimate|Value)$",
  std.error = "Std\\. ?Error",
  df = "df",
  statistic = "(t|Z)[ -]value",
  p.value = "(Pr\\(>|[tZ]\\)|p[ -]value)"
)

## like match(), but with a table of regexes
regex_match <- function(x, table) {
  r <- sapply(
    x,
    function(z) {
      m <- vapply(col_matches, grepl, x = z, ignore.case = TRUE, logical(1))
      if (any(m)) return(which(m)) else return(NA)
    }
  )
  return(unname(r))
}

## rename columns according to regex matches
## names that are not matched are left unchanged
rename_regex_match <- function(x, table = col_matches) {
    rr <- regex_match(names(x), table)
    names(x)[!is.na(rr)] <- names(table)[na.omit(rr)]
    return(x)
}

## convert confint output to a data frame and relabel columns
cifun <- function(x, method="Wald", ddf.method=NULL, level=0.95, ...) {
  Estimate <- `Std. Error` <- NULL ## global var check
  ## compute Wald-t estimates if necessary (not handled by confint for lmerTest)
  if (!is.null(ddf.method)) {
     if (method != "Wald") warning("ddf.method ignored when conf.method != \"Wald\"")
     cc <-  as.data.frame(coef(summary(x, ddf=ddf.method)))
     if (!"df" %in% colnames(cc)) {
         mult <- qnorm((1 + level) / 2)
     } else {
         mult <- stats::qt((1+level)/2, df=cc$df)
     }
     r <- (cc
         %>% transmute(conf.low=Estimate-mult*`Std. Error`,
                       conf.high=Estimate+mult*`Std. Error`)
     )
  }  else {
      r <- as.data.frame(confint(x, method=method, level=level, ...))
  }
  r <- (r
      %>% tibble()
      %>% setNames(c("conf.low", "conf.high"))
      )
  return(r)
}


## put specified columns (if they exist) as first columns in output, leave
##  other columns as is
reorder_frame <- function(x, first_cols = c("effect", "group", "term", "estimate")) {
  ## order of first arg to intersect() determines order of results ...
  first_cols <- intersect(first_cols, names(x))
  other_cols <- setdiff(names(x), first_cols)
  return(x[, c(first_cols, other_cols)])
}

## FIXME: store functions to run as a list of expressions,
##  allow user-specified 'skip' argument?
finish_glance <- function(ret = dplyr::tibble(), x) {
  stopifnot(length(ret) == 0 || nrow(ret) == 1)

  ## catch NULL, numeric(0), error responses

  tfun <- function(e) {
    tt <- tryCatch(eval(substitute(e)), error = function(e) NA)
    if (length(tt) == 0) tt <- NA
    return(tt)
  }

  newvals <- dplyr::tibble(
    nobs = tfun(stats::nobs(x)),
    sigma = tfun(stats::sigma(x)),
    logLik = tfun(as.numeric(stats::logLik(x))),
    AIC = tfun(stats::AIC(x)),
    BIC = tfun(stats::BIC(x)),
    deviance = suppressWarnings(tfun(stats::deviance(x))),
    df.residual = tfun(stats::df.residual(x))
  )
  ## drop NA values
  newvals <- newvals[!vapply(newvals, is.na, logical(1))]

  if (length(ret) == 0) {
    return(newvals)
  } else {
    return(dplyr::bind_cols(ret, newvals))
  }
}

######
## experimental finish_glance ...
f2 <- function(ret = data.frame(), x, skip_funs = character(0)) {
  tfun <- function(f) {
    tt <- tryCatch(f(x), error = function(e) NA)
    if (length(tt) == 0) tt <- NA
    return(tt)
  }

  stopifnot(length(ret) == 0 || nrow(ret) == 1)

  funs <- c("nobs", "logLik", "AIC", "BIC", "deviance", "df.residual")
  funs <- setdiff(funs, skip_funs)

  newvals <- lapply(funs, function(f) as.numeric(tfun(get(f, "package:stats"))))
  newvals <- as.data.frame(newvals)
  names(newvals) <- funs
  ## drop NA values
  newvals <- newvals[!vapply(newvals, is.na, logical(1))]
  if (length(ret) == 0) {
    return(unrowname(newvals))
  } else {
    return(unrowname(data.frame(ret, newvals)))
  }
}

## like process_lm, but without lm-specific confint stuff
## applied *downstream* (after CIs etc have already been added)
trans_coef <- function(ret, x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE,
                       trans = identity) {
  ## FIXME: should transform sds as well
  if (missing(trans)) {
    if (exponentiate) {
      if (is.null(x$family) || !grepl("log", x$family$link)) {
        warning(paste(
          "Exponentiating coefficients, ",
          "but model did not use ",
          "a (log, logit, cloglog) link function"
        ))
      }
      trans <- exp
    } else {
      trans <- identity
    }
  }
  ret <- (ret
  %>%
    mutate_at(intersect(c("term", "conf.low", "conf.high")), trans))
  return(ret)
}


## naming function
ran_pars_name <- function(x, ran_prefix) {
  v <- na.omit(unlist(x))
  if (length(v) == 0) v <- "Observation"
  p <- paste(v, collapse = ".")
  if (!identical(ran_prefix, NA)) {
    p <- paste(ran_prefix[length(v)], p,
      sep = getOption("broom.mixed.sep1")
    )
  }
  return(p)
}


## FIXME: 1. sds_..., sigma not properly translated
##        2. names of
## translate brms-style "terms" into standard broom.mixed
## term -> effect, group, term
trans_brms_params <- function(tidy_obj) {
  tt <- tidy_obj[["term"]]
  effcodes <- c("b", "sd", "cor", "s", "sigma", "sds", "r", "lp__")
  neweffcodes <- c(
    "fixed", "ran_pars", "ran_pars",
    "ran_vals", "ran_pars", "???", "ran_vals", "lp__"
  )
  effc2 <- effcodes
  effc2[4] <- "s(?!(igma))" ## negative lookahead ...
  effc2 <- paste0("^(", paste(effc2, collapse = "|"), ")")
  effects <- stringr::str_extract(tt, effc2)
  tt2 <- stringr::str_remove(tt, paste0(effc2, "_?"))
  ## keep r/s distinction a little longer
  ## https://stackoverflow.com/questions/42457189/greedy-regex-for-one-part-non-greedy-for-other?rq=1
  ## (.*?) go until FIRST occurence of next pattern
  ## (?= ...  ) lookahead -- don't include this stuff in the extracted string
  group <- stringr::str_extract(tt2, "(.*?)(?=(__|\\[))")
  grpvals <- effects %in% c("sd", "cor", "r")
  ## remove group__ for sd/cor
  tt2[grpvals] <- stringr::str_remove(tt2[grpvals], "(.*?)__")
  tt2[grpvals] <- stringr::str_remove(tt2[grpvals], "(.*?)(?=(\\[))")
  effects <- as.character(factor(effects,
    levels = effcodes,
    labels = neweffcodes
  ))
  ## replace 'term' (in place) with 'effect', 'group', 'term'
  term_col <- which(names(tidy_obj) == "term")
  prev_cols <- if (term_col > 1) seq(term_col - 1) else numeric(0)
  ## restore sd/cor to beginning of
  res <- bind_cols(tidy_obj[prev_cols],
    effect = effects,
    group = group, term = tt2,
    tidy_obj[(term_col + 1):ncol(tidy_obj)]
  )
  return(res)
}

## enforce consistent column order for *existing* columns
## should contain all possible column names
reorder_cols <- function(x) {
  all_cols <- c(
    "response","effect",
    "component", ## glmmTMB, brms
    "group", "level", "term", "index", "estimate",
    "std.error", "statistic",
    "df", "p.value",
    "conf.low", "conf.high", "rhat", "ess"
  )
  return(select(x, intersect(all_cols, names(x))))
}

rename_cols <- function(x,
                        from = c("Estimate", "Std. Error", "(z|Z|t) value", "Pr\\(>"),
                        to = c("estimate", "std.error", "statistic", "p.value")) {
  if (!is.data.frame(x)) x <- dplyr::as_tibble(x)
  for (i in seq_along(from)) {
    if (length(m <- grep(from[i], names(x))) > 0) {
      names(x)[m] <- to[i]
    }
  }
  return(x)
}

has_rownames <- function(df) {
  return (!tibble::is_tibble(df) &&
          any(rownames(df) != as.character(seq(nrow(df)))))
}

## previously from broom
## converts to tibble, adding non-trivial rownames and optionally renaming existing columns
fix_data_frame <- function(df, newnames=NULL, newcol="term") {
    rn <- rownames(df) ## grab rownames *before* df conversion
    ## must happen **AFTER** saving rownames
    df <- as_tibble(df, .name_repair="minimal")
    if (!is.null(newnames)) df <- setNames(df,newnames)
    ## add rownames as term **if necessary**
    if (!("term" %in% newnames) &&
        !("term" %in% names(df)) &&
        !is.null(rn)) {
        df <- tibble(rn,df)
        names(df)[1] <- newcol
    }
    return(df)
}

##' Retrieve all method/class combinations currently provided by the broom.mixed package
##' @examples print(get_methods(), n = Inf)
##' @importFrom tidyr complete pivot_wider separate
##' @export
get_methods <- function() {
  fun <- method <- provided <- NULL ## NSE code check
  ## TO DO: include associated package? not necessarily easy to find
  ## the package associated with a class ... (can look for print method
  ## with getAnywhere(), but would need package loaded ...)
  res <- (tibble(fun = ls(getNamespace("broom.mixed")))
    %>% filter(grepl("^(tidy|glance|augment)\\.", fun))
    %>% separate(fun, into = c("method", "class"), sep = "\\.", extra = "merge")
    %>% mutate(provided = TRUE)
    %>% complete(method, class, fill = list(provided = FALSE))
    %>% pivot_wider(names_from = method, values_from = provided)
    ## reorder
    %>% dplyr::select(class, tidy, glance, augment)
  )
  class(res) <- c("show_methods", class(res))
  return(res)
}

