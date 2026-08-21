is_missing_value <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

has_code <- function(x, code) {
  any(trimws(as.character(x)) == as.character(code), na.rm = TRUE)
}

has_positive_count <- function(x) {
  values <- suppressWarnings(as.numeric(trimws(as.character(x))))
  any(values > 0, na.rm = TRUE)
}

has_count_at_least <- function(x, threshold) {
  values <- suppressWarnings(as.numeric(trimws(as.character(x))))
  any(values >= threshold, na.rm = TRUE)
}

first_non_missing <- function(x) {
  observed <- x[!is_missing_value(x)]
  if (length(observed) == 0) {
    return(x[NA_integer_][1])
  }
  observed[[1]]
}

all_existing_are <- function(x, codes) {
  observed <- trimws(as.character(x[!is_missing_value(x)]))
  length(observed) > 0 && all(observed %in% as.character(codes))
}

all_missing_or_code <- function(x, codes) {
  values <- trimws(as.character(x))
  all(is_missing_value(x) | values %in% as.character(codes))
}


apply_variable_spec <- function(data, spec) {
  for (source_name in names(spec)) {
    settings <- spec[[source_name]]
    new_name <- settings$new_name
    x <- data[[source_name]]

    x <- switch(
      settings$type,
      numeric = suppressWarnings(as.numeric(as.character(x))),
      character = as.character(x),
      factor = {
        if (!is.null(settings$levels)) {
          factor(
            x,
            levels = settings$levels,
            labels = settings$labels
          )
        } else {
          factor(x)
        }
      },
      stop("Unsupported variable type: ", settings$type)
    )

    data[[new_name]] <- x

    if (new_name != source_name) {
      data[[source_name]] <- NULL
    }
  }

  data
}