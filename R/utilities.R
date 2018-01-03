## most of these are unexported (small) functions from broom;
## could be removed if these were exported

#' strip rownames from an object
#' 
#' @param x a data frame
unrowname <- function(x) {
    rownames(x) <- NULL
    return(x)
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
    indices[which(stats::complete.cases(original))] = seq_len(nrow(x))
    x[indices, ]
}

## list of regex matches for mixed-effect columns -> broom names
col_matches <- list(estimate="^(Estimate|Value)$",
                    std.error="Std\\. ?Error",
                    df = "df",
                    statistic = "(t|Z)[ -]value",
                    p.value = "(Pr\\(>|[tZ]\\)|p[ -]value)")

## like match(), but with a table of regexes
regex_match <- function(x,table) {
    r <- sapply(x,
           function(z) {
              m <- vapply(col_matches,grepl,x=z,ignore.case=TRUE,logical(1))
              if (any(m)) return(which(m)) else return(NA)
    })
    return(unname(r))
}

## rename columns according to regex matches
rename_regex_match <- function(x,table=col_matches) {
    names(x) <- names(table)[regex_match(names(x),table)]
    return(x)
}

## convert confint output to a data frame and relabel columns
cifun <- function(x,...) {
    r <- confint(x,...) %>%
        data.frame %>%
        setNames(c("conf.low","conf.high"))
    return(r)
}


## put specified columns (if they exist) as first columns in output, leave
##  other columns as is
reorder_frame <- function(x,first_cols=c("effect","group","term","estimate")) {
    ## order of first arg to intersect() determines order of results ...
    first_cols <- intersect(first_cols,names(x))
    other_cols <- setdiff(names(x),first_cols)
    return(x[,c(first_cols,other_cols)])
}

## FIXME: store functions to run as a list of expressions,
##  allow user-specified 'skip' argument?
finish_glance <- function(ret=data.frame(),x) {

    stopifnot(length(ret)==0 || nrow(ret)==1)
    
    ## catch NULL, numeric(0), error responses

    tfun <- function(e) {
        tt <- tryCatch(eval(substitute(e)),error=function(e) NA)
        if (length(tt)==0) tt <- NA
        return(tt)
    }
    
    newvals <- data.frame(sigma=tfun(sigma(x)),
                          logLik=tfun(as.numeric(stats::logLik(x))),
                          AIC=tfun(stats::AIC(x)),
                          BIC=tfun(stats::BIC(x)),
                          deviance=suppressWarnings(tfun(stats::deviance(x))),
                          df.residual=tfun(stats::df.residual(x)))
    ## drop NA values
    newvals <- newvals[!vapply(newvals,is.na,logical(1))]

    if (length(ret)==0) {
        return(unrowname(newvals))
    } else {
        return(unrowname(data.frame(ret,newvals)))
    }
}

######
## experimental finish_glance ...
f2 <- function(ret=data.frame(),x,skip_funs=character(0)) {

    tfun <- function(f) {
        tt <- tryCatch(f(x),error=function(e) NA)
        if (length(tt)==0) tt <- NA
        return(tt)
    }

    stopifnot(length(ret)==0 || nrow(ret)==1)

    funs <- c("logLik","AIC","BIC","deviance","df.residual")
    funs <- setdiff(funs,skip_funs)

    newvals <- lapply(funs, function(f) as.numeric(tfun(get(f,"package:stats"))))
    newvals <- as.data.frame(newvals)
    names(newvals) <- funs
    ## drop NA values
    newvals <- newvals[!vapply(newvals,is.na,logical(1))]
    if (length(ret)==0) {
        return(unrowname(newvals))
    } else {
        return(unrowname(data.frame(ret,newvals)))
    }
}
