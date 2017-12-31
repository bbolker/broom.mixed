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
          
