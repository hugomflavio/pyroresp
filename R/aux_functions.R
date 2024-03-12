#' check that the requested argument value is present in the data
#' 
#' @param arg the argument value
#' @param data A vector of the target data
#' @param name For messaging purposes; the type of data we're working with.
#' @param verbose Should a warning be sent if the argument does not fully
#' 	match the data?
#' 
#' @return the matching elements of arg
#' 
#' @keywords internal
#' 
check_arg_in_data <- function(arg, data, name, verbose = TRUE) {
	arg_check <- arg %in% data

	if (all(!arg_check)) {
		stop("Could not find requested ", name, " in the input", call. = FALSE)
	}
	if (verbose && any(!arg_check)) {
		if (sum(arg_check) < 5) {
			warning("Could not find ", name, " ",
					paste(arg[arg_check], collapse = ", "),
					" in the input. Discarding those.",
					immediate. = TRUE, call. = FALSE)
		} else {
			warning("Could not find ", sum(arg_check), 
					" of the requested ", name, " ",
					"in the input. Discarding those.",
					immediate. = TRUE, call. = FALSE)
		}
	}

	arg <- arg[arg_check]
	return(arg)
}

#' convert weight to volume in ml
#' 
#' @param w a units object in either grams or kilograms
#' @param d The density. defaults to 1 (1g = 1ml). Defaults to 1.
#'  Higher densities lead to lower volumes, e.g. if d = 2, then 1g = 0.5ml.
#' 
#' @return a units object in ml
#' 
conv_w_to_ml <- function(w, d = 1) {
  if (!inherits(w, "units")) {
    stop("w must be a units class object")
  }
  if (!(as.character(units(w)) %in% c("g", "kg"))) {
    stop("w must be either in grams or kilograms.")
  }

  x <- as.numeric(w) / d

  if (as.character(units(w)) == "kg") {
    x <- x/1000
  }

  units(x) <- "ml"

  return(x)
}
