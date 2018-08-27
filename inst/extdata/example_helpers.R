save_file <- function(..., pkg, type = "rda") {
  f <- file.path("inst", "extdata", sprintf("%s_example.%s", pkg, type))
  if (type == "rda") {
    save(..., file = f)
  } else {
    saveRDS(..., file = f)
  }
  invisible(NULL)
}

pkg <- NULL
run_pkg <- function(pkg, e) {
  pkg <<- pkg
  if (require(pkg, character.only = TRUE)) {
    eval(e)
    return(TRUE)
  } else {
    cat(sprintf("%s examples not run\n", pkg))
    return(FALSE)
  }
}


hack_size <- function(x, ...) {
    UseMethod("hack_size")
}

hack_size.stanfit <- function(x) {
    x@stanmodel <- structure(numeric(0), class="stanmodel")
    x@.MISC <- new.env()
    return(x)
}

hack_size.brmsfit <- function(x) {
    x$fit <- hack_size(x$fit)
    return(x)
}

hack_size.stanreg <- function(x) {
    x$stanfit <- hack_size(x$stanfit)
    return(x)
}
