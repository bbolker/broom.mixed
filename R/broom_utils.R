## utilities stolen from broom/R/utilities.R (unexported)

#' @importFrom utils globalVariables
globalVariables(
    c(":=")
)


## a function that, given named arguments, will make a one-row
## tibble, switching out NULLs for the appropriate NA type.
as_glance_tibble <- function(..., na_types) {
  
  cols <- list(...)
  
  if (length(cols) != stringr::str_length(na_types)) {
    stop(
      "The number of columns provided does not match the number of ",
      "column types provided."
    )
  }
  
  na_types_long <- parse_na_types(na_types)
  
  entries <- purrr::map2(cols, 
                         na_types_long, 
                         function(.x, .y) {if (length(.x) == 0) .y else .x})
  
  tibble::as_tibble_row(entries)
  
}

#' Convert a data.frame or matrix to a tibble
#' 
#' This function is meant for use inside `tidy.*` methods.
#'
#' @param x a data.frame or matrix
#' @param new_names new column names, not including the rownames
#' @param new_col the name of the new rownames column
#'
#' @return a tibble with rownames moved into a column and new column
#' names assigned
#' @noRd
as_tidy_tibble <- function(x, new_names = NULL, new_column = "term") {
  
  if (!is.null(new_names) && length(new_names) != ncol(x)) {
    stop("newnames must be NULL or have length equal to number of columns")
  }

  ret <- x
  
  # matrices will take setNames by element, so rename conditionally
  if (!is.null(new_names)) {
    if (inherits(x, "data.frame")) {
      ret <- setNames(x, new_names)
    } else {
      colnames(ret) <- new_names
    }
  }
  
  if (all(rownames(x) == seq_len(nrow(x)))) {
    # don't need to move rownames into a new column
    tibble::as_tibble(ret)
  } else {
    # don't use tibble rownames to col because of name repairing
    dplyr::bind_cols(!!new_column := rownames(x), 
                     tibble::as_tibble(ret))
  }
}

na_types_dict <- list("r" = NA_real_,
                      "i" = rlang::na_int,
                      "c" = NA_character_,
                      "l" = rlang::na_lgl)

# a function that converts a string to a vector of NA types.
# e.g. "rri" -> c(NA_real_, NA_real_, rlang::na_int)
parse_na_types <- function(s) {
  
  positions <- purrr::map(
    stringr::str_split(s, pattern = ""), 
    match,
    table = names(na_types_dict)
  ) %>%
    unlist()
  
  na_types_dict[positions] %>%
    unlist() %>%
    unname()
}
