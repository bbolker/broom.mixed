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

